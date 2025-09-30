test_that("ci.curve yields the same outputs as cif of mets", {
  testdata <- createTestData(20, 1, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=TRUE, na_strata=FALSE)
  library(mets)
  library(Rcpp)
  Surv <- km.curve:::Surv
  e <- cif(Event(t,epsilon)~1, data=testdata, cause=1)
  e_surv <- 1-e$mu
  e_time <- e$times
  expected <- as.numeric(c(e_surv[c(1,2,4,6,8,10,12,14)], e_time[c(1,2,4,6,8,10,12,14)]))
  t <- ci.curve(Surv(t, epsilon)~1, data=testdata, code.event=c(1,2), report.ggsurvfit = FALSE, conf.type = "none")
  tested <- as.numeric(c(t$surv[c(2:9)], t$time[c(2:9)]))
  expect_equal(tested, expected)
})

test_that("ci.curve by strata yields the same outputs as subsetting", {
  testdata <- createTestData(20, 1, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=TRUE, na_strata=FALSE)
  testdata_f <- testdata[testdata$strata==FALSE,]
  library(Rcpp)
  Surv <- km.curve:::Surv
  e <- cif(Event(t,epsilon)~1, data=testdata, cause=1)
  e <- ci.curve(Surv(t, epsilon)~1, data=testdata_f, code.event=c(1,2), report.ggsurvfit = FALSE, conf.type = "none")
  expected <- c(e$surv, e$time)
  t <- ci.curve(Surv(t, epsilon)~strata, data=testdata, code.event=c(1,2), report.ggsurvfit = FALSE, conf.type = "none")
  tested <- c(t$surv[c(1:5)], t$time[c(1:5)])
  expect_equal(tested, expected)
})
