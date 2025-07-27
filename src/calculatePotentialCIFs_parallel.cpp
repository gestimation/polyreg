#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <roptim.h>

using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;

// CIFTask remains unchanged
class CIFTask : public roptim::Functor {
public:
  CIFTask(double a1, vec b1, double a2, vec b2, const std::string& e1, const std::string& e2, double pb)
    : alpha_tmp_1(a1), beta_tmp_1(b1), alpha_tmp_2(a2), beta_tmp_2(b2),
      effect1(e1), effect2(e2), prob_bound(pb) {}

  double operator()(const arma::vec& log_p) override {
    auto clampLogP = [this](double x) {
      double ex = exp(x);
      if (ex < prob_bound) return log(prob_bound);
      if ((1 - ex) < prob_bound) return log(1 - prob_bound);
      return x;
    };

    vec clog_p(4), exp_lp(4);
    for (int i = 0; i < 4; ++i) {
      clog_p[i] = clampLogP(log_p[i]);
      exp_lp[i] = exp(clog_p[i]);
    }

    double lp0102;
    if ((1 - exp_lp[0] - exp_lp[1] < prob_bound) || (1 - exp_lp[2] - exp_lp[3] < prob_bound)) {
      lp0102 = log(prob_bound);
    } else {
      lp0102 = log(abs(1 - exp_lp[0] - exp_lp[1])) + log(abs(1 - exp_lp[2] - exp_lp[3]));
    }

    vec ret(4);
    if (effect1 == "RR") {
      ret[0] = alpha_tmp_1 - clog_p[0] - clog_p[2] + lp0102;
      ret[1] = beta_tmp_1[0] - clog_p[2] + clog_p[0];
    } else if (effect1 == "OR") {
      ret[0] = alpha_tmp_1 - clog_p[0] - clog_p[2] + lp0102;
      ret[1] = beta_tmp_1[0] - clog_p[2] + clog_p[0] +
        log(1 - exp_lp[2]) - log(1 - exp_lp[0]);
    } else if (effect1 == "SHR") {
      ret[0] = alpha_tmp_1 - clog_p[0] - clog_p[2] + lp0102;
      ret[1] = exp(beta_tmp_1[0]) - (log(1 - exp_lp[2]) / log(1 - exp_lp[0]));
    }

    if (effect2 == "RR") {
      ret[2] = alpha_tmp_2 - clog_p[1] - clog_p[3] + lp0102;
      ret[3] = beta_tmp_2[0] - clog_p[3] + clog_p[1];
    } else if (effect2 == "OR") {
      ret[2] = alpha_tmp_2 - clog_p[1] - clog_p[3] + lp0102;
      ret[3] = beta_tmp_2[0] - clog_p[3] + clog_p[1] +
        log(1 - exp_lp[3]) - log(1 - exp_lp[1]);
    } else if (effect2 == "SHR") {
      ret[2] = alpha_tmp_2 - clog_p[1] - clog_p[3] + lp0102;
      ret[3] = exp(beta_tmp_2[0]) - (log(1 - exp_lp[3]) / log(1 - exp_lp[1]));
    }

    return sum(square(ret));
  }

private:
  double alpha_tmp_1;
  vec beta_tmp_1;
  double alpha_tmp_2;
  vec beta_tmp_2;
  const std::string& effect1;
  const std::string& effect2;
  double prob_bound;
};

// Thread-safe Worker with RMatrix
struct CIFWorker : public RcppParallel::Worker {
  const arma::vec& alpha_tmp_1_vec;
  const arma::vec& alpha_tmp_2_vec;
  const arma::vec& beta_tmp_1;
  const arma::vec& beta_tmp_2;
  const std::string& effect1;
  const std::string& effect2;
  const double prob_bound;
  const std::string& method;
  const int maxit;
  const double reltol;
  const arma::vec log_p0;

  RMatrix<double> results;

  CIFWorker(const arma::vec& a1, const arma::vec& a2,
            const arma::vec& b1, const arma::vec& b2,
            const std::string& e1, const std::string& e2,
            double pb, const std::string& meth, int m, double tol,
            const arma::vec& lp0, Rcpp::NumericMatrix& res)
    : alpha_tmp_1_vec(a1), alpha_tmp_2_vec(a2),
      beta_tmp_1(b1), beta_tmp_2(b2),
      effect1(e1), effect2(e2), prob_bound(pb),
      method(meth), maxit(m), reltol(tol),
      log_p0(lp0), results(res) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      vec log_p = log_p0;
      CIFTask task(alpha_tmp_1_vec[i], beta_tmp_1,
                   alpha_tmp_2_vec[i], beta_tmp_2,
                   effect1, effect2, prob_bound);
      try {
        roptim::Roptim<CIFTask> opt(method);
        opt.control.maxit = maxit;
        opt.control.reltol = reltol;
        opt.control.trace = 0;
        opt.minimize(task, log_p);
        vec exp_log_p = exp(log_p);
        for (int j = 0; j < 4; ++j) {
          results(i, j) = exp_log_p[j];
        }
      } catch (...) {
        for (int j = 0; j < 4; ++j) {
          results(i, j) = NA_REAL;
        }
      }
    }
  }
};

// [[Rcpp::export]]
Rcpp::NumericMatrix calculatePotentialCIFs_parallel(
    arma::mat x_l,
    arma::vec log_p0,
    arma::vec alpha_tmp_1_vec,
    arma::vec alpha_tmp_2_vec,
    arma::vec beta_tmp_1,
    arma::vec beta_tmp_2,
    std::string effect1,
    std::string effect2,
    double prob_bound,
    std::string method = "BFGS",
    int maxit = 100,
    double reltol = 1e-8
) {
  int n = x_l.n_rows;
  Rcpp::NumericMatrix results(n, 4);
  CIFWorker worker(alpha_tmp_1_vec, alpha_tmp_2_vec,
                   beta_tmp_1, beta_tmp_2,
                   effect1, effect2, prob_bound,
                   method, maxit, reltol,
                   log_p0, results);
  parallelFor(0, n, worker);
  return results;
}
