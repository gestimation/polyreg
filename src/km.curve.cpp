#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List calculateKM_rcpp(Rcpp::NumericVector t, Rcpp::IntegerVector d,
                            Rcpp::NumericVector w = Rcpp::NumericVector::create(),
                            Rcpp::IntegerVector strata = Rcpp::IntegerVector::create(),
                            std::string error = "greenwood") {

                //  Rcpp::Rcout << "Method: " << error[0] << std::endl;

                Rcpp::List survfit_list;

                std::vector<double> combined_times;
                std::vector<double> combined_surv;
                std::vector<int> combined_n_risk;
                std::vector<int> combined_n_event;
                std::vector<int> combined_n_censor;
                std::vector<double> combined_std_err;
                std::vector<int> combined_n_stratum;
                std::vector<int> combined_u_stratum;

                if ((strata.size() == 0 || Rcpp::unique(strata).size() == 1) && (w.size() == 0 || (Rcpp::unique(w).size() == 1 && w[0] == 1))) {

                  Rcpp::NumericVector unique_times = Rcpp::unique(t);
                  std::sort(unique_times.begin(), unique_times.end());

                  int u = unique_times.size();
                  Rcpp::NumericVector km(u);
                  Rcpp::NumericVector km_i(u);
                  Rcpp::IntegerVector weighted_n_risk(u);
                  Rcpp::IntegerVector weighted_n_event(u);
                  Rcpp::IntegerVector weighted_n_censor(u);
                  Rcpp::NumericVector std_err(u);

                  int n_stratum = t.size();
                  int n = t.size();
                  for (int i = 0; i < u; ++i) {
                    double time = unique_times[i];
                    double weighted_n_i = 0;
                    for (int j = 0; j < n; ++j) {
                      if (t[j] >= time) {
                        weighted_n_i++;
                      }
                    }
                    weighted_n_risk[i] = weighted_n_i;

                    double weighted_events = 0;
                    for (int j = 0; j < n; ++j) {
                      if (t[j] == time && d[j] == 1) {
                        weighted_events++;
                      }
                    }

                    if (weighted_n_risk[i] > 0) {
                      km_i[i] = 1 - weighted_events / weighted_n_risk[i];
                      weighted_n_event[i] += weighted_events;
                    } else {
                      km_i[i] = 1;
                    }

                    if (i > 0) {
                      weighted_n_censor[i-1] = weighted_n_risk[i-1] - weighted_n_risk[i] - weighted_n_event[i-1];
                    }
                    if (i == u-1) {
                      weighted_n_censor[i] = weighted_n_risk[i] - weighted_n_event[i];
                    }

                    if (i == 0) {
                      km[i] = km_i[i];
                    } else {
                      km[i] = km[i - 1] * km_i[i];
                    }

                    double sum_se = 0;
                    for (int j = 0; j <= i; ++j) {
                      double n_i = weighted_n_risk[j];
                      double d_i = 0;
                      for (int k = 0; k < n; ++k) {
                        if (t[k] == unique_times[j] && d[k] == 1) {
                          d_i ++;
                        }
                      }
                      if (n_i > d_i) {
                        if (error == "tsiatis") {
                          sum_se += (d_i / (n_i * n_i));
                        } else if (error == "greenwood") {
                          sum_se += (d_i / (n_i * (n_i - d_i)));
                        }
                      } else {
                        sum_se = std::numeric_limits<double>::infinity();
                      }
                    }
                    std_err[i] = sqrt(sum_se);
                  }

                  combined_times.insert(combined_times.end(), unique_times.begin(), unique_times.end());
                  combined_surv.insert(combined_surv.end(), km.begin(), km.end());
                  combined_n_risk.insert(combined_n_risk.end(), weighted_n_risk.begin(), weighted_n_risk.end());
                  combined_n_event.insert(combined_n_event.end(), weighted_n_event.begin(), weighted_n_event.end());
                  combined_n_censor.insert(combined_n_censor.end(), weighted_n_censor.begin(), weighted_n_censor.end());
                  combined_std_err.insert(combined_std_err.end(), std_err.begin(), std_err.end());
                  combined_n_stratum.insert(combined_n_stratum.end(), n_stratum);

                } else if (strata.size() == 0 || Rcpp::unique(strata).size() == 1) {

                  Rcpp::NumericVector unique_times = Rcpp::unique(t);
                  std::sort(unique_times.begin(), unique_times.end());

                  int u = unique_times.size();
                  Rcpp::NumericVector km(u);
                  Rcpp::NumericVector km_i(u);
                  Rcpp::IntegerVector weighted_n_risk(u);
                  Rcpp::IntegerVector weighted_n_event(u);
                  Rcpp::IntegerVector weighted_n_censor(u);
                  Rcpp::NumericVector std_err(u);

                  int n_stratum = t.size();
                  int n = t.size();
                  for (int i = 0; i < u; ++i) {
                    double time = unique_times[i];
                    double weighted_n_i = 0;
                    for (int j = 0; j < n; ++j) {
                      if (t[j] >= time) {
                        weighted_n_i += w[j];
                      }
                    }
                    weighted_n_risk[i] = weighted_n_i;

                    double weighted_events = 0;
                    for (int j = 0; j < n; ++j) {
                      if (t[j] == time && d[j] == 1) {
                        weighted_events += w[j];
                      }
                    }

                    if (weighted_n_risk[i] > 0) {
                      km_i[i] = 1 - weighted_events / weighted_n_risk[i];
                      weighted_n_event[i] += weighted_events;
                    } else {
                      km_i[i] = 1;
                    }

                    if (i > 0) {
                      weighted_n_censor[i-1] = weighted_n_risk[i-1] - weighted_n_risk[i] - weighted_n_event[i-1];
                    }
                    if (i == u-1) {
                      weighted_n_censor[i] = weighted_n_risk[i] - weighted_n_event[i];
                    }

                    if (i == 0) {
                      km[i] = km_i[i];
                    } else {
                      km[i] = km[i - 1] * km_i[i];
                    }

                    double sum_se = 0;
                    for (int j = 0; j <= i; ++j) {
                      double n_i = weighted_n_risk[j];
                      double d_i = 0;
                      for (int k = 0; k < n; ++k) {
                        if (t[k] == unique_times[j] && d[k] == 1) {
                          d_i += w[k];
                        }
                      }
                      if (n_i > d_i) {
                        if (error == "tsiatis") {
                          sum_se += (d_i / (n_i * n_i));
                        } else if (error == "greenwood") {
                          sum_se += (d_i / (n_i * (n_i - d_i)));
                        }
                      } else {
                        sum_se = std::numeric_limits<double>::infinity();
                      }
                    }
                    std_err[i] = sqrt(sum_se);
                  }

                  combined_times.insert(combined_times.end(), unique_times.begin(), unique_times.end());
                  combined_surv.insert(combined_surv.end(), km.begin(), km.end());
                  combined_n_risk.insert(combined_n_risk.end(), weighted_n_risk.begin(), weighted_n_risk.end());
                  combined_n_event.insert(combined_n_event.end(), weighted_n_event.begin(), weighted_n_event.end());
                  combined_n_censor.insert(combined_n_censor.end(), weighted_n_censor.begin(), weighted_n_censor.end());
                  combined_std_err.insert(combined_std_err.end(), std_err.begin(), std_err.end());
                  combined_n_stratum.insert(combined_n_stratum.end(), n_stratum);

                } else if (w.size() == 0 || (Rcpp::unique(w).size() == 1 && w[0] == 1)) {

                  Rcpp::IntegerVector strata_vec = strata;

                  // Loop over each strata level and calculate the Kaplan-Meier estimate
                  for (int i = 0; i < Rcpp::max(strata_vec); ++i) {
                    // Select the data subset corresponding to the current strata level
                    Rcpp::LogicalVector strata_condition = (strata_vec == (i + 1));
                    Rcpp::NumericVector t_selected = t[strata_condition];
                    Rcpp::IntegerVector d_selected = d[strata_condition];

                    // Calculate weighted_n for this strata (sum of weights in the selected dataset)
                    int n_stratum = t_selected.size();
                    double weighted_n_stratum = 0;
                    int n_event_stratum = 0;
                    for (int j = 0; j < n_stratum; ++j) {
                      weighted_n_stratum += d_selected[j];
                      n_event_stratum += d_selected[j];
                    }

                    // Calculate Kaplan-Meier for the selected subset
                    Rcpp::NumericVector unique_times = Rcpp::unique(t_selected);
                    std::sort(unique_times.begin(), unique_times.end());
                    int u_stratum = unique_times.size();

                    Rcpp::NumericVector km(u_stratum);
                    Rcpp::NumericVector km_i(u_stratum);
                    Rcpp::IntegerVector weighted_n_risk(u_stratum);
                    Rcpp::IntegerVector weighted_n_event(u_stratum);
                    Rcpp::IntegerVector weighted_n_censor(u_stratum);
                    Rcpp::NumericVector std_err(u_stratum);

                    for (int j = 0; j < u_stratum; ++j) {
                      double time = unique_times[j];
                      double weighted_n_sub = 0;
                      for (int k = 0; k < t_selected.size(); ++k) {
                        if (t_selected[k] >= time) {
                          weighted_n_sub++;
                        }
                      }
                      weighted_n_risk[j] = weighted_n_sub;

                      double weighted_events = 0;
                      for (int k = 0; k < t_selected.size(); ++k) {
                        if (t_selected[k] == time && d_selected[k] == 1) {
                          weighted_events++;
                        }
                      }

                      if (weighted_n_risk[j] > 0) {
                        km_i[j] = 1 - weighted_events / weighted_n_risk[j];
                        weighted_n_event[j] += weighted_events;
                      } else {
                        km_i[j] = 1;
                      }

                      if (j > 0) {
                        weighted_n_censor[j-1] = weighted_n_risk[j-1] - weighted_n_risk[j] - weighted_n_event[j-1];
                      }
                      if (j == u_stratum-1) {
                        weighted_n_censor[j] = weighted_n_risk[j] - weighted_n_event[j];
                      }

                      km[j] = (j == 0) ? km_i[j] : km[j - 1] * km_i[j];

                      double sum_se = 0;
                      for (int k = 0; k <= j; ++k) {
                        double n_i = weighted_n_risk[k];
                        double d_i = 0;
                        for (int m = 0; m < t_selected.size(); ++m) {
                          if (t_selected[m] == unique_times[k] && d_selected[m] == 1) {
                            d_i += d_selected[m];
                          }
                        }
                        if (n_i > d_i) {
                          if (error == "tsiatis") {
                            sum_se += (d_i / (n_i * n_i));
                          } else if (error == "greenwood") {
                            sum_se += (d_i / (n_i * (n_i - d_i)));
                          }
                        } else {
                          sum_se = std::numeric_limits<double>::infinity();
                        }
                      }
                      std_err[j] = sqrt(sum_se);
                    }

                    combined_times.insert(combined_times.end(), unique_times.begin(), unique_times.end());
                    combined_surv.insert(combined_surv.end(), km.begin(), km.end());
                    combined_n_risk.insert(combined_n_risk.end(), weighted_n_risk.begin(), weighted_n_risk.end());
                    combined_n_event.insert(combined_n_event.end(), weighted_n_event.begin(), weighted_n_event.end());
                    combined_n_censor.insert(combined_n_censor.end(), weighted_n_censor.begin(), weighted_n_censor.end());
                    combined_std_err.insert(combined_std_err.end(), std_err.begin(), std_err.end());
                    combined_n_stratum.insert(combined_n_stratum.end(), n_stratum);
                    combined_u_stratum.insert(combined_u_stratum.end(), u_stratum);
                  }
                } else {

                  Rcpp::IntegerVector strata_vec = strata;

                  // Loop over each strata level and calculate the Kaplan-Meier estimate
                  for (int i = 0; i < Rcpp::max(strata_vec); ++i) {
                    // Select the data subset corresponding to the current strata level
                    Rcpp::LogicalVector strata_condition = (strata_vec == (i + 1));  // 1-based indexing for strata
                    Rcpp::NumericVector t_selected = t[strata_condition];
                    Rcpp::IntegerVector d_selected = d[strata_condition];
                    Rcpp::NumericVector w_selected = w[strata_condition];

                    // Calculate weighted_n for this strata (sum of weights in the selected dataset)
                    int n_stratum = t_selected.size();
                    double weighted_n_stratum = 0;
                    int n_event_stratum = 0;
                    for (int j = 0; j < n_stratum; ++j) {
                      weighted_n_stratum += w_selected[j];
                      n_event_stratum += d_selected[j];
                    }

                    // Calculate Kaplan-Meier for the selected subset
                    Rcpp::NumericVector unique_times = Rcpp::unique(t_selected);
                    std::sort(unique_times.begin(), unique_times.end());
                    int u_stratum = unique_times.size();

                    Rcpp::NumericVector km(u_stratum);
                    Rcpp::NumericVector km_i(u_stratum);
                    Rcpp::IntegerVector weighted_n_risk(u_stratum);
                    Rcpp::IntegerVector weighted_n_event(u_stratum);
                    Rcpp::IntegerVector weighted_n_censor(u_stratum);
                    Rcpp::NumericVector std_err(u_stratum);

                    for (int j = 0; j < u_stratum; ++j) {
                      double time = unique_times[j];

                      double weighted_n_sub = 0;
                      for (int k = 0; k < t_selected.size(); ++k) {
                        if (t_selected[k] >= time) {
                          weighted_n_sub += w_selected[k];
                        }
                      }
                      weighted_n_risk[j] = weighted_n_sub;

                      double weighted_events = 0;
                      for (int k = 0; k < t_selected.size(); ++k) {
                        if (t_selected[k] == time && d_selected[k] == 1) {
                          weighted_events += w_selected[k];
                        }
                      }

                      if (weighted_n_risk[j] > 0) {
                        km_i[j] = 1 - weighted_events / weighted_n_risk[j];
                        weighted_n_event[j] += weighted_events;
                      } else {
                        km_i[j] = 1;
                      }

                      if (j > 0) {
                        weighted_n_censor[j-1] = weighted_n_risk[j-1]-weighted_n_risk[j] - weighted_n_event[j-1];
                      }
                      if (j == u_stratum-1) {
                        weighted_n_censor[j] = weighted_n_risk[j] - weighted_n_event[j];
                      }

                      km[j] = (j == 0) ? km_i[j] : km[j - 1] * km_i[j];

                      double sum_se = 0;
                      for (int k = 0; k <= j; ++k) {
                        double n_i = weighted_n_risk[k];
                        double d_i = 0;
                        for (int m = 0; m < t_selected.size(); ++m) {
                          if (t_selected[m] == unique_times[k] && d_selected[m] == 1) {
                            d_i += w_selected[m];
                          }
                        }
                        if (n_i > d_i) {
                          if (error == "tsiatis") {
                            sum_se += (d_i / (n_i * n_i));
                          } else if (error == "greenwood") {
                            sum_se += (d_i / (n_i * (n_i - d_i)));
                          }
                        } else {
                          sum_se = std::numeric_limits<double>::infinity();
                        }
                      }
                      std_err[j] = sqrt(sum_se);
                    }

                    combined_times.insert(combined_times.end(), unique_times.begin(), unique_times.end());
                    combined_surv.insert(combined_surv.end(), km.begin(), km.end());
                    combined_n_risk.insert(combined_n_risk.end(), weighted_n_risk.begin(), weighted_n_risk.end());
                    combined_n_event.insert(combined_n_event.end(), weighted_n_event.begin(), weighted_n_event.end());
                    combined_n_censor.insert(combined_n_censor.end(), weighted_n_censor.begin(), weighted_n_censor.end());
                    combined_std_err.insert(combined_std_err.end(), std_err.begin(), std_err.end());
                    combined_n_stratum.insert(combined_n_stratum.end(), n_stratum);
                    combined_u_stratum.insert(combined_u_stratum.end(), u_stratum);
                  }
                }

                Rcpp::NumericVector all_times = Rcpp::wrap(combined_times);
                Rcpp::NumericVector all_surv = Rcpp::wrap(combined_surv);
                Rcpp::IntegerVector all_n_risk = Rcpp::wrap(combined_n_risk);
                Rcpp::IntegerVector all_n_event = Rcpp::wrap(combined_n_event);
                Rcpp::IntegerVector all_n_censor = Rcpp::wrap(combined_n_censor);
                Rcpp::NumericVector all_std_err = Rcpp::wrap(combined_std_err);
                Rcpp::IntegerVector all_n_stratum = Rcpp::wrap(combined_n_stratum);
                Rcpp::IntegerVector all_u_stratum = Rcpp::wrap(combined_u_stratum);

                // Create final combined output for all strata
                survfit_list = Rcpp::List::create(
                  Rcpp::_["time"] = all_times,
                  Rcpp::_["surv"] = all_surv,
                  Rcpp::_["n.risk"] = all_n_risk,
                  Rcpp::_["n"] = combined_n_stratum,
                  Rcpp::_["n.event"] = all_n_event,
                  Rcpp::_["n.censor"] = all_n_censor,
                  Rcpp::_["unweighted.n"] = all_n_stratum,
                  Rcpp::_["std.err"] = all_std_err,
                  Rcpp::_["high"] = R_NilValue,
                  Rcpp::_["low"] = R_NilValue,
                  Rcpp::_["conf.type"] = "log-log",
                  Rcpp::_["strata"] = all_u_stratum,
                  Rcpp::_["type"] = "kaplan-meier",
                  Rcpp::_["method"] = "Kaplan-Meier"
                );
                return survfit_list;
}

