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

test_that("survival.curve yields the same outputs as survfit when strata is logical", {
  library(survival)
  library(ggsurvfit)
  library(Rcpp)
  testdata <- createTestData(20, 1, first_zero=FALSE, last_zero=TRUE, subset_present=FALSE, logical_strata=TRUE, na_strata=FALSE)
  e <- survfit(Surv(t, d)~strata, testdata, weight=w, conf.type = "none")
  t <- survival.curve(Surv(t, d) ~ strata, data = testdata, weight="w", report.ggsurvfit = FALSE, conf.type = "none", outcome.type = "SURVIVAL")
  #  expected <- as.numeric(c(e$time, round(e$surv,digit=5), e$n, e$n.risk, e$n.event, e$n.censor, round(e$std.err,digit=5), e$strata))
  #  tested <- as.numeric(c(t$time, round(t$surv,digit=5), t$n, t$n.risk, t$n.event, t$n.censor, round(t$std.err,digit=5), t$strata))
  expected <- as.numeric(c(e$time, round(e$surv,digit=5), e$n, e$n.risk, e$n.event, e$n.censor, e$lower, e$strata))
  tested <- as.numeric(c(t$time, round(t$surv,digit=5), t$n, t$n.risk, t$n.event, t$n.censor, t$lower, t$strata))
  expect_equal(expected, tested)
})

test_that("survival.curve matches ci.curve for competing risks output", {
  data(diabetes.complications)
  ci_fit <- ci.curve(Surv(t, epsilon) ~ strata, data = diabetes.complications, report.ggsurvfit = FALSE, conf.type = "none")
  surv_fit <- survival.curve(Surv(t, epsilon) ~ strata, data = diabetes.complications, report.ggsurvfit = FALSE, conf.type = "none", outcome.type = "COMPETING-RISK")
  expect_equal(surv_fit$surv, ci_fit$surv)
  expect_equal(surv_fit$std.err, ci_fit$std.err)
  expect_equal(surv_fit$type, "Aalen-Johansen")
})

test_that("survival.curve produced expected survfit under coding other than the default", {
  data(diabetes.complications)
  diabetes.complications$epsilon1 <- diabetes.complications$epsilon + 1
  expected <- survival.curve(Surv(t, epsilon) ~ strata, data = diabetes.complications, code.event1=1, code.event2=2, code.censoring=0, conf.type = "none", outcome.type="C")
  tested <- survival.curve(Surv(t, epsilon1) ~ strata, data = diabetes.complications, code.event1=2, code.event2=3, code.censoring=1, conf.type = "none", outcome.type="C")
  expect_equal(tested$surv, expected$surv)
  expect_equal(tested$std.err, expected$std.err)
})

test_that("survival.curve produced expected survfit under coding other than the default", {
  data(diabetes.complications)
  diabetes.complications$d <- as.numeric(diabetes.complications$epsilon>0)
  diabetes.complications$d1 <- as.numeric(diabetes.complications$epsilon>0) + 1
  expected <- survival.curve(Surv(t, d) ~ strata, data = diabetes.complications, code.censoring=0, conf.type = "none", outcome.type="S")
  tested <- survival.curve(Surv(t, d1) ~ strata, data = diabetes.complications, code.censoring=1, conf.type = "none", outcome.type="S")
  expect_equal(tested$surv, expected$surv)
  expect_equal(tested$std.err, expected$std.err)
})
