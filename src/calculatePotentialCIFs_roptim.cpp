#include <RcppArmadillo.h>
#include <roptim.h>

using namespace Rcpp;
using namespace arma;

// Define CIFTask as Functor for roptim
class CIFTask : public roptim::Functor {
public:
  CIFTask(double a1, vec b1, double a2, vec b2, std::string e1, std::string e2, double pb)
  : alpha_tmp_1(a1), beta_tmp_1(b1), alpha_tmp_2(a2), beta_tmp_2(b2),
    effect1(e1), effect2(e2), prob_bound(pb) {}

  double operator()(const arma::vec& log_p) override {
    auto clampLogP = [this](double x) {
      double ex = exp(x);
      if (ex < prob_bound) return log(prob_bound);
      if (1 - ex < prob_bound) return log(1 - prob_bound);
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
  std::string effect1, effect2;
  double prob_bound;
};

// [[Rcpp::export]]
arma::vec calculatePotentialCIFs_roptim(
  arma::vec log_p0,
  double alpha_tmp_1,
  arma::vec beta_tmp_1,
  double alpha_tmp_2,
  arma::vec beta_tmp_2,
  std::string effect1,
  std::string effect2,
  double prob_bound,
  std::string method = "BFGS",
  int maxit = 100,
  double reltol = 1e-8
) {
  CIFTask task(alpha_tmp_1, beta_tmp_1, alpha_tmp_2, beta_tmp_2, effect1, effect2, prob_bound);
  roptim::Roptim<CIFTask> opt(method);
  opt.control.maxit = maxit;
  opt.control.reltol = reltol;
  opt.control.trace = 1;
  opt.minimize(task, log_p0);
  return log_p0;
}

// [[Rcpp::export]]
arma::mat calculatePotentialCIFs_roptim_all(
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
  arma::mat results(n, 4);

  for (int i = 0; i < n; ++i) {
    CIFTask task(alpha_tmp_1_vec[i], beta_tmp_1, alpha_tmp_2_vec[i], beta_tmp_2,
                 effect1, effect2, prob_bound);

    arma::vec log_p = log_p0;
    roptim::Roptim<CIFTask> opt(method);
    opt.control.maxit = maxit;
    opt.control.reltol = reltol;
    opt.control.trace = 0;

    opt.minimize(task, log_p);
    results.row(i) = log_p.t();
  }
  return arma::exp(results);
}
