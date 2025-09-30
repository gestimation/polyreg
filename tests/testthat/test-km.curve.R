test_that("Surv yields the same outputs as Surv of survfit", {
  library(survival)
  testdata <- createTestData(100, 1, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE)
  expected <- survival::Surv(testdata$t, testdata$d)
  tested <- Surv(testdata$t, testdata$d)
  expect_equal(expected, tested)
})

test_that("Surv yields the same outputs as Surv of survfit with NA", {
  library(survival)
  testdata <- createTestData(100, 1, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE)
  testdata$t[1] <- NA
  testdata$d[1] <- NA
  expected <- survival::Surv(testdata$t, testdata$d)
  tested <- Surv(testdata$t, testdata$d)
  expect_equal(expected, tested)
})

test_that("Surv yields the same outputs as Surv of survfit with a factor", {
  library(survival)
  testdata <- createTestData(100, 1, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE)
  f <- as.factor(testdata$d)
  expected <- survival::Surv(testdata$t, f)
  tested <- Surv(testdata$t, f)
  expect_equal(expected[,2], tested[,2])
})

test_that("km.curve yields the same outputs as survfit", {
  library(survival)
  library(ggsurvfit)
  library(Rcpp)
  testdata <- createTestData(20, 2, first_zero=TRUE, last_zero=FALSE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE)
  e <- survfit(Surv(t, d)~1, testdata, weight=w, na.action=na.omit, conf.type = "none")
  t <- km.curve(Surv(t, d)~1, testdata, weight="w", na.action=na.omit, conf.type = "none", report.ggsurvfit = FALSE)
  expected <- as.numeric(c(e$time, round(e$surv,digit=5), e$n, e$n.risk, e$n.event, e$n.censor, round(e$std.err,digit=5), e$strata))
  tested <- as.numeric(c(t$time, round(t$surv,digit=5), t$n, t$n.risk, t$n.event, t$n.censor, round(t$std.err,digit=5), t$strata))
  expect_equal(expected, tested)
})

test_that("km.curve yields the same outputs as survfit when strata is logical", {
  library(survival)
  library(ggsurvfit)
  library(Rcpp)
  testdata <- createTestData(20, 1, first_zero=FALSE, last_zero=TRUE, subset_present=FALSE, logical_strata=TRUE, na_strata=FALSE)
  e <- survfit(Surv(t, d)~strata, testdata, weight=w, conf.type = "none")
  t <- km.curve(Surv(t, d)~strata, testdata, weight="w", conf.type = "none", report.ggsurvfit = FALSE)
  expected <- as.numeric(c(e$time, round(e$surv,digit=5), e$n, e$n.risk, e$n.event, e$n.censor, round(e$std.err,digit=5), e$strata))
  tested <- as.numeric(c(t$time, round(t$surv,digit=5), t$n, t$n.risk, t$n.event, t$n.censor, round(t$std.err,digit=5), t$strata))
  expect_equal(expected, tested)
})

test_that("Jack Knife standard error of km.curve yields the same outputs as separate jack knifing when strata is present", {
  library(ggsurvfit)
  library(Rcpp)
  testdata <- createTestData(20, 1, first_zero=TRUE, last_zero=FALSE, subset_present=FALSE, logical_strata=TRUE, na_strata=FALSE)
  testdata1 <- subset(testdata, strata==FALSE)
  testdata2 <- subset(testdata, strata==TRUE)
  t <- km.curve(Surv(t, d)~strata, testdata, weight="w", error = "jackknife", report.ggsurvfit = FALSE)
  e1 <- km.curve(Surv(t, d)~1, testdata1, weight="w", error = "jackknife", report.ggsurvfit = FALSE)
  e2 <- km.curve(Surv(t, d)~1, testdata2, weight="w", error = "jackknife", report.ggsurvfit = FALSE)
  expected <- c(e1$std.err, e2$std.err)
  tested <- t$std.err
  expect_equal(expected, tested)
})

