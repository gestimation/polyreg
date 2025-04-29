#test_that("calculateKM_rcpp and get_surv produced expected KM in diabetes.complications", {
#  library(Rcpp)
#  data("diabetes.complications")
#  diabetes.complications$d <- as.numeric(diabetes.complications$epsilon==0)
#  resC <- suppressWarnings(phreg(Surv(t, d) ~ 1, data=diabetes.complications))
#  out_predict1 <- suppressWarnings(stats::predict(resC, newdata = diabetes.complications, type = "survival", times = diabetes.complications$t, individual.time = TRUE, se = FALSE, km = TRUE, tminus = TRUE))
#  expected <- out_predict1$surv
#  expected <- round(expected[1:10,],digit=3)
#  out_km <- calculateKM_rcpp(diabetes.complications$t, diabetes.complications$d)
#  tested <- as.matrix(get_surv(diabetes.complications$t, out_km$surv, out_km$time))
#  tested <- round(tested[1:10,], digit=3)
#  expect_equal(tested, expected)
#})

#test_that("calculateKM_rcpp and get_surv produced expected stratified KM in diabetes.complications", {
#  library(Rcpp)
#  data("diabetes.complications")
#  diabetes.complications$d <- as.numeric(diabetes.complications$epsilon==0)
#  resC <- suppressWarnings(phreg(Surv(t, d) ~ strata(strata), data=diabetes.complications))
#  out_predict1 <- suppressWarnings(stats::predict(resC, newdata = diabetes.complications, type = "survival", times = diabetes.complications$t, individual.time = TRUE, se = FALSE, km = TRUE, tminus = TRUE))
#  expected <- out_predict1$surv
#  expected <- round(expected[1:10,],digit=3)
#  out_km <- calculateKM_rcpp(t=diabetes.complications$t, d=diabetes.complications$d, strata=diabetes.complications$strata)
#  tested <- as.matrix(get_surv(diabetes.complications$t, out_km$surv, out_km$time, diabetes.complications$strata, out_km$strata))
#  tested <- round(tested[1:10,], digit=3)
#  expect_equal(tested, expected)
#})

