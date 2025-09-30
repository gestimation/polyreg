test_that("km.curve yields the same outputs as survfit", {
  library(survival)
  library(ggsurvfit)
  library(Rcpp)
  df_test <- createTestData(200, 2, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE)
  e <- survfit(Surv(t, d)~strata, df_test, weight=w, conf.type = "log-log")
  t <- km.curve(Surv(t, d)~strata, df_test, weight="w", conf.type = "log-log", report.ggsurvfit = FALSE, report.survfit.std.err = TRUE)
  e$lower <- sapply(e$lower, function(x) ifelse(is.nan(x), NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(is.nan(x), NA, x))
  e$lower <- sapply(e$lower, function(x) ifelse(x>=1, NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(x>=1, NA, x))
  e$lower <- sapply(e$lower, function(x) ifelse(x<=0, NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(x<=0, NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(is.nan(x), NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(is.nan(x), NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(x>=1, NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(x>=1, NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(x<=0, NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(x<=0, NA, x))
  expected <- as.numeric(c(e$time, round(e$surv,digit=5), e$n, e$n.risk, e$n.event, e$n.censor, round(e$std.err,digit=5), round(e$lower,digit=5), round(e$upper,digit=5), e$strata))
  tested <- as.numeric(c(t$time, round(t$surv,digit=5), t$n, t$n.risk, t$n.event, t$n.censor, round(t$std.err,digit=5), round(t$lower,digit=5), round(t$upper,digit=5), t$strata))
  expect_equal(expected, tested)
})

test_that("km.curve yields the same outputs as survfit", {
  library(survival)
  library(ggsurvfit)
  library(Rcpp)
  df_test <- createTestData(200, 2, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE)
  e <- survfit(Surv(t, d)~strata, df_test, weight=w, conf.type = "log")
  t <- km.curve(Surv(t, d)~strata, df_test, weight="w", conf.type = "log", report.ggsurvfit = FALSE, report.survfit.std.err = TRUE)
  e$lower <- sapply(e$lower, function(x) ifelse(is.nan(x), NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(is.nan(x), NA, x))
  e$lower <- sapply(e$lower, function(x) ifelse(x>=1, NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(x>=1, NA, x))
  e$lower <- sapply(e$lower, function(x) ifelse(x<=0, NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(x<=0, NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(is.nan(x), NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(is.nan(x), NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(x>=1, NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(x>=1, NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(x<=0, NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(x<=0, NA, x))

  expected <- as.numeric(c(e$time, round(e$surv,digit=5), e$n, e$n.risk, e$n.event, e$n.censor, round(e$std.err,digit=5), round(e$lower,digit=5), round(e$upper,digit=5), e$strata))
  tested <- as.numeric(c(t$time, round(t$surv,digit=5), t$n, t$n.risk, t$n.event, t$n.censor, round(t$std.err,digit=5), round(t$lower,digit=5), round(t$upper,digit=5), t$strata))
  expect_equal(expected, tested)
})

test_that("km.curve yields the same outputs as survfit", {
  library(survival)
  library(ggsurvfit)
  library(Rcpp)
  df_test <- createTestData(200, 2, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE)
  e <- survfit(Surv(t, d)~strata, df_test, weight=w, conf.type = "a")
  t <- km.curve(Surv(t, d)~strata, df_test, weight="w", conf.type = "a", report.ggsurvfit = FALSE, report.survfit.std.err = TRUE)
  e$lower <- sapply(e$lower, function(x) ifelse(is.nan(x), NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(is.nan(x), NA, x))
  e$lower <- sapply(e$lower, function(x) ifelse(x>=1, NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(x>=1, NA, x))
  e$lower <- sapply(e$lower, function(x) ifelse(x<=0, NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(x<=0, NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(is.nan(x), NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(is.nan(x), NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(x>=1, NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(x>=1, NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(x<=0, NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(x<=0, NA, x))

  expected <- as.numeric(c(e$time, round(e$surv,digit=5), e$n, e$n.risk, e$n.event, e$n.censor, round(e$std.err,digit=5), round(e$lower,digit=5), round(e$upper,digit=5), e$strata))
  tested <- as.numeric(c(t$time, round(t$surv,digit=5), t$n, t$n.risk, t$n.event, t$n.censor, round(t$std.err,digit=5), round(t$lower,digit=5), round(t$upper,digit=5), t$strata))
  expect_equal(expected, tested)
})

test_that("km.curve yields the same outputs as survfit", {
  library(survival)
  library(ggsurvfit)
  library(Rcpp)
  df_test <- createTestData(200, 2, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE)
  e <- survfit(Surv(t, d)~strata, df_test, weight=w, conf.type = "plain")
  t <- km.curve(Surv(t, d)~strata, df_test, weight="w", conf.type = "plain", report.ggsurvfit = FALSE, report.survfit.std.err = TRUE)
  e$lower <- sapply(e$lower, function(x) ifelse(is.nan(x), NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(is.nan(x), NA, x))
  e$lower <- sapply(e$lower, function(x) ifelse(x>=1, NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(x>=1, NA, x))
  e$lower <- sapply(e$lower, function(x) ifelse(x<=0, NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(x<=0, NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(is.nan(x), NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(is.nan(x), NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(x>=1, NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(x>=1, NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(x<=0, NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(x<=0, NA, x))

  expected <- as.numeric(c(e$time, round(e$surv,digit=5), e$n, e$n.risk, e$n.event, e$n.censor, round(e$std.err,digit=5), round(e$lower,digit=5), round(e$upper,digit=5), e$strata))
  tested <- as.numeric(c(t$time, round(t$surv,digit=5), t$n, t$n.risk, t$n.event, t$n.censor, round(t$std.err,digit=5), round(t$lower,digit=5), round(t$upper,digit=5), t$strata))
  expect_equal(expected, tested)
})

test_that("empinf in boot package yields expected pseudo observations", {
  library(boot)
  fn1 <- function(data) {
    out_km <- calculateKM_rcpp(data$t, data$d, data$w, as.integer(data$strata), "none")
    return(out_km$surv[5])
  }
  fn2 <- function(data, i) {
    data <- data[i,]
    out_km <- calculateKM_rcpp(data$t, data$d, data$w, as.integer(data$strata), "none")
    #    out_km <- calculateKM_rcpp(data$t, data$d, w, as.integer(data$strata), "none")
    return(out_km$surv[5])
  }
  df_test <- createTestData(20, 1, first_zero=FALSE, last_zero=TRUE, subset_present=FALSE, logical_strata=TRUE, na_strata=FALSE)
  jk1 <- calculateJackKnife(df_test, fn1)
  e <- jk1$pseudo_observations - 0.2
  t <- empinf(data=df_test, statistic=fn2, type="jack", stype="i")
  expected <- round(e, digit=5)
  tested <- round(t, digit=5)
  expect_equal(expected, tested)
})



#test_that("calculateJackKnifeSE yields expected pseudo observations", {
#  fn1 <- function(data) {
#    out_km <- calculateKM_rcpp(data$t, data$d, data$w, as.integer(data$strata), "none")
#    return(out_km$surv[1])
#  }
#  fn2 <- function(data) {
#    out_km <- calculateKM_rcpp(data$t, data$d, data$w, as.integer(data$strata), "none")
#    return(out_km$surv[5])
#  }
#  fn3 <- function(data) {
#    out_km <- calculateKM_rcpp(data$t, data$d, data$w, as.integer(data$strata), "none")
#    return(out_km$surv[9])
#  }
#  out_readSurv <- createTestData(20, 2, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE)
#  jk1 <- calculateJackKnife(out_readSurv, fn1)
#  jk2 <- calculateJackKnife(out_readSurv, fn2)
#  jk3 <- calculateJackKnife(out_readSurv, fn3)
#  expected <- c(round(jk1$pseudo_observations, digit=5), round(jk2$pseudo_observations, digit=5), round(jk3$pseudo_observations, digit=5))

#  out_km <- calculateKM_rcpp(out_readSurv$t, out_readSurv$d, out_readSurv$w, as.integer(out_readSurv$strata), "none")
#  if (error=="jackknife") {
#    fn <- function(data) {
#      out_km <- calculateKM_rcpp(data$t, data$d, data$w, as.integer(data$strata), "none")
#      return(out_km)
#    }
#    out_jk <- calculateJackKnifeSE_wo_sort(out_readSurv, fn)
#    out_km$std.err <- out_jk$std.err
#  }
#  tested <- c(round(out_jk$pseudo.obs[,1], digit=5), round(out_jk$pseudo.obs[,5], digit=5), round(out_jk$pseudo.obs[,9], digit=5))
#  expect_equal(expected, tested)
#})
