test_that("survival.curve matches km.curve for survival output", {
  testdata <- createTestData(20, 1, first_zero = TRUE, last_zero = TRUE, subset_present = FALSE, logical_strata = TRUE, na_strata = FALSE)
  Surv <- km.curve:::Surv
  km_fit <- km.curve(Surv(t, d) ~ strata, data = testdata, report.ggsurvfit = FALSE, conf.type = "none")
  surv_fit <- survival.curve(Surv(t, d) ~ strata, data = testdata, report.ggsurvfit = FALSE, conf.type = "none", output.type = "SURVIVAL")
  expect_equal(surv_fit$surv, km_fit$surv)
  expect_equal(surv_fit$time, km_fit$time)
  expect_equal(surv_fit$type, "kaplan-meier")
})

test_that("survival.curve matches ci.curve for competing risks output", {
  testdata <- createTestData(20, 1, first_zero = TRUE, last_zero = TRUE, subset_present = FALSE, logical_strata = TRUE, na_strata = FALSE)
  Surv <- km.curve:::Surv
  ci_fit <- ci.curve(Surv(t, epsilon) ~ strata, data = testdata, code.event = c(1, 2), report.ggsurvfit = FALSE, conf.type = "none")
  surv_fit <- survival.curve(Surv(t, epsilon) ~ strata, data = testdata, code.event = c(1, 2), report.ggsurvfit = FALSE, conf.type = "none", output.type = "COMPETING-RISK")
  expect_equal(surv_fit$surv, ci_fit$surv)
  expect_equal(surv_fit$time, ci_fit$time)
  expect_equal(surv_fit$type, "Aalen-Johansen")
})

test_that("survival.curve returns numeric confidence bounds without CI", {
  library(survival)
  library(ggsurvfit)
  library(Rcpp)

  testdata <- createTestData(
    30,
    2,
    first_zero = TRUE,
    last_zero = FALSE,
    subset_present = FALSE,
    logical_strata = TRUE,
    na_strata = FALSE
  )

  expect_silent({
    fit <- survival.curve(
      Surv(t, d) ~ strata,
      data = testdata,
      weight = "w",
      conf.type = "none",
      report.ggsurvfit = FALSE
    )
  })

  expect_true(is.numeric(fit$upper))
  expect_true(is.numeric(fit$lower))
  expect_length(fit$upper, length(fit$surv))
  expect_length(fit$lower, length(fit$surv))
})
