test_that("survival.curve matches km.curve for survival output", {
  testdata <- createTestData(20, 1, first_zero = TRUE, last_zero = TRUE, subset_present = FALSE, logical_strata = TRUE, na_strata = FALSE)
  km_fit <- km.curve(Surv(t, d) ~ strata, data = testdata, report.ggsurvfit = FALSE, conf.type = "none")
  surv_fit <- survival.curve(Surv(t, d) ~ strata, data = testdata, report.ggsurvfit = FALSE, conf.type = "none", outcome.type = "SURVIVAL")
  expect_equal(surv_fit$surv, km_fit$surv)
  expect_equal(surv_fit$time, km_fit$time)
  expect_equal(surv_fit$type, "kaplan-meier")
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
