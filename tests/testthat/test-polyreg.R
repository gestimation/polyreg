test_that("polyreg produced expected coefficients and variance covariance matrix from competing risks data in diabetes.complications", {
  data(diabetes.complications)
  output <- polyreg(nuisance.model = Event(t,epsilon)~+1, exposure = 'fruitq1', data = diabetes.complications, effect.measure1='RR', effect.measure2='RR', time.point=8, outcome.type='C')
  tested_coefficient <- round(output$coefficient,digit=3)
  tested_cov <- round(output$cov[1,],digit=3)
  tested <- as.vector(cbind(tested_coefficient,tested_cov))
#  expected <- c(-1.383, 0.300, -3.991, 0.076, 0.007, -0.005, -0.001, 0.005)
  expected <- c(-1.383, 0.300, -3.991, 0.076, 0.008, -0.005, 0.003, -0.002)
  expect_equal(tested, expected)
})

test_that("polyreg produced expected coefficients and variance covariance matrix for categorical exposure in diabetes.complications", {
  data(diabetes.complications)
  output <- polyreg(nuisance.model = Event(t,epsilon)~+1, exposure = 'strata', data = diabetes.complications, effect.measure1='RR', effect.measure2='RR', time.point=8, outcome.type='C', code.exposure.ref = 1)
  tested_coefficient <- round(output$coefficient,digit=3)
  tested_cov <- round(output$cov[1,],digit=3)
  tested <- as.vector(cbind(tested_coefficient,tested_cov))
  expected <- c(-2.297, 1.143, -0.197, 1.213, -4.817, 0.132, 1.506, 0.163, 0.040, -0.024, -0.031, -0.025, 0.009, -0.006, -0.006, -0.006)
  expect_equal(tested, expected)
})

test_that("polyreg produced expected coefficients and variance covariance matrix when stratified IPCW used", {
  data(diabetes.complications)
  output <- polyreg(nuisance.model = Event(t,epsilon)~+1, exposure = 'fruitq1', data = diabetes.complications, strata = 'strata', effect.measure1='SHR', effect.measure2='SHR', time.point=8, outcome.type='C')
  tested_coefficient <- round(output$coefficient,digit=3)
  tested_cov <- round(output$cov[1,],digit=3)
  tested <- as.vector(cbind(tested_coefficient,tested_cov))
  #  expected <- c(-1.383, 0.300, -3.991, 0.076, 0.007, -0.005, -0.001, 0.005)
  #  expected <- c(-1.383, 0.300, -3.988, 0.078, 0.007, -0.005, -0.001, 0.005)
  expected <- c(-1.383, 0.363, -3.988, 0.081, 0.009, -0.006, 0.004, -0.003)
  expect_equal(tested, expected)
})

test_that("polyreg produced expected coefficients and variance covariance matrix under coding other than the default", {
  data(diabetes.complications)
  diabetes.complications$epsilon1 <- diabetes.complications$epsilon + 1
  output <- polyreg(nuisance.model = Event(t,epsilon1)~+1, exposure = 'fruitq1', data = diabetes.complications,
                    code.event1=2, code.event2=3, code.censoring=1, code.exposure.ref = 1,
                    effect.measure1='RR', effect.measure2='RR', time.point=8, outcome.type='C')
  tested_coefficient <- round(output$coefficient,digit=3)
  tested_cov <- round(output$cov[1,],digit=3)
  tested <- as.vector(cbind(tested_coefficient,tested_cov))
#  expected <- c(-1.383, -0.300, -3.991, -0.076, 0.016, -0.012, 0.004, -0.004)
#  expected <- c(-1.393, -0.301, -4.006, -0.074, 0.016, -0.012, 0.004, -0.004)
  expected <- c(-1.393, -0.301, -4.006, -0.074, 0.017, -0.012, 0.010, -0.008)
  expect_equal(tested, expected)
})

test_that("polyreg produced expected coefficients and variance covariance matrix from survival data in diabetes.complications", {
  data(diabetes.complications)
  diabetes.complications$d <- as.numeric(diabetes.complications$epsilon>0)
  output <- polyreg(nuisance.model = Event(t,d)~+1, exposure = 'fruitq1', data = diabetes.complications, effect.measure1='RR', effect.measure2='RR', time.point=8, outcome.type='S')
  tested_coefficient <- round(output$coefficient,digit=3)
  tested_cov <- round(output$cov[1,],digit=3)
  tested <- as.vector(cbind(tested_coefficient,tested_cov))
  expected <- c(-0.901, 0.252, 0.003, -0.003)
  expect_equal(tested, expected)
})

test_that("polyreg produced expected coefficients and variance covariance matrix when binomial regression is applied to diabetes.complications", {
  data(diabetes.complications)
  diabetes.complications$d <- as.numeric(diabetes.complications$epsilon>0)
  output <- polyreg(nuisance.model = d~+1, exposure = 'fruitq1', data = diabetes.complications, effect.measure1='RR', effect.measure2='RR', time.point=8, outcome.type='B')
  tested_coefficient <- round(output$coefficient,digit=3)
  tested_cov <- round(output$cov[1,],digit=3)
  tested <- as.vector(cbind(tested_coefficient,tested_cov))
  expected <- c(-0.877, 0.249, 0.003, -0.003)
  expect_equal(tested, expected)
})

test_that("polyreg produced expected coefficients and variance covariance matrix from survival data with missing data", {
  data(diabetes.complications)
  diabetes.complications$d <- as.numeric(diabetes.complications$epsilon>0)

  expected_df <- diabetes.complications[-(1:10), ]
  expected_output <- polyreg(nuisance.model = Event(t,d)~sex, exposure = 'fruitq1', strata = 'strata', data = expected_df, effect.measure1='RR', effect.measure2='RR', time.point=8, outcome.type='S')
  expected <- round(expected_output$coefficient,digit=3)

  diabetes.complications$t[1:2] <- NA
  diabetes.complications$d[3:4] <- NA
  diabetes.complications$sex[5:6] <- NA
  diabetes.complications$fruitq1[7:8] <- NA
  diabetes.complications$strata[9:10] <- NA
  output <- polyreg(nuisance.model = Event(t,d)~sex, exposure = 'fruitq1', strata = 'strata', data = diabetes.complications, effect.measure1='RR', effect.measure2='RR', time.point=8, outcome.type='S')
  tested <- round(output$coefficient,digit=3)
  expect_equal(tested, expected)
})


test_that("polyreg produced expected coefficients and variance covariance matrix from survival data in men", {
  data(diabetes.complications)
  nuisance.model <- Event(t,d)~1
  diabetes.complications$d <- as.numeric(diabetes.complications$epsilon>0)
  subset_condition1="(diabetes.complications$sex == 1)"
  subset_condition2="(sex == 1)"
  subset_condition3="(data$sex == 1)"
  subset_condition4=NULL

  other.variables.analyzed <- c("fruitq1", "strata")
  all_vars <- c(all.vars(nuisance.model), other.variables.analyzed)
  expected <- diabetes.complications[, all_vars, drop = FALSE]
  expected_df <- subset(expected, eval(parse(text = subset_condition1)))
  #  print(nrow(expected_df))

  expected_output <- polyreg(nuisance.model = Event(t,d)~1, exposure = 'fruitq1', strata = 'strata', data = expected_df, effect.measure1='RR', effect.measure2='RR', time.point=8, outcome.type='S')
  expected <- round(expected_output$coefficient,digit=3)

  output <- polyreg(nuisance.model = Event(t,d)~1, exposure = 'fruitq1', strata = 'strata', data = diabetes.complications, subset.condition="(diabetes.complications$sex == 1)", effect.measure1='RR', effect.measure2='RR', time.point=8, outcome.type='S')
  tested <- round(output$coefficient,digit=3)
  expect_equal(tested, expected)
})

test_that("polyreg produced expected common effects at 1:5 in prostate", {
  library(dplyr)
  library(janitor)
  library(boot)
  data(prostate)
  prostate <- prostate %>% mutate(epsilon=ifelse(status=="alive",0,
                                                 ifelse(status=="dead - prostatic ca",1,
                                                        ifelse(status=="dead - other ca",1,
                                                               ifelse(status=="dead - heart or vascular",2,
                                                                      ifelse(status=="dead - cerebrovascular",2,2)
                                                               )))))
  prostate$epsilon <- as.numeric(prostate$epsilon)
  prostate$a <- as.numeric((prostate$rx=="placebo"))
  prostate$t <- prostate$dtime/12
  output <- polyreg(nuisance.model = Event(t,epsilon) ~ +1, exposure = 'a', strata='stage', data = prostate,
                    effect.measure1='RR', effect.measure2='RR', time.point=1:5, outcome.type='POLY-PROPORTIONAL', report.boot.conf=FALSE)
  tested <- round(output$coefficient,digit=3)
#  expected <- c(-2.246, -2.710, -1.420, -0.669, -0.109, -0.034, -3.702, -2.144, -0.882, -0.129,  0.506,  0.058)
  expected <- c(-2.289, -1.753, -1.461, -1.363, -1.254, -0.034, -2.056, -1.609, -1.245, -1.067, -0.971, 0.058)
  expect_equal(tested, expected)
})

test_that("polyreg produced expected common effects in prostate", {
  library(dplyr)
  library(janitor)
  library(boot)
  data(prostate)
  prostate <- prostate %>% mutate(d=ifelse(status=="alive",0,
                                           ifelse(status=="dead - prostatic ca",1,
                                                  ifelse(status=="dead - other ca",1,
                                                         ifelse(status=="dead - heart or vascular",1,
                                                                ifelse(status=="dead - cerebrovascular",1,1)
                                                         )))))
  prostate$d <- as.numeric(prostate$d)
  prostate$a <- as.numeric((prostate$rx=="placebo"))
  prostate$t <- prostate$dtime/12
  output <- polyreg(nuisance.model = Event(t,d) ~ +1, exposure = 'a', strata='stage', data = prostate,
                    effect.measure1='RR', outcome.type='PROPORTIONAL', report.boot.conf=FALSE)
  tested <- round(output$coefficient[1:4], digit=3)
  #expected <- c(-15.988, -15.550, -7.569, -7.560)
  #expected <- c(-6.835, -6.019, -5.016, -4.902)
  expected <- c(-3.442, -2.882, -2.631, -2.526)
  expect_equal(tested, expected)
})

#test_that("polyreg produced expected coefficients and variance covariance matrix from bmt dataset", {
#  library(mets)
#  library(nleqslv)
#  data(bmt)
#  output <- polyreg(nuisance.model = Event(time, cause)~age+tcell, exposure = 'platelet', cens.model = Event(time,cause==0)~+1, data = bmt, effect.measure1='RR', effect.measure2='RR', time.point=24, outcome.type='COMPETINGRISK')
#  tested_coefficient <- round(output$coefficient,digit=2)
#  tested_cov <- round(output$cov[1,],digit=2)
#  tested <- as.vector(cbind(tested_coefficient,tested_cov))
#  expected <- c(-0.45, 1.01, -1.32, -0.48, -1.68, 0.41, 0.63, 0.12, 0.02, 0.01, -0.02, -0.01, -0.00, 0.05, -0.02, 0.01)
#  expect_equal(tested, expected)
#})


#test_that("polyreg produced expected bootstrap confidence intervals in prostate", {
#  library(dplyr)
#  library(janitor)
#  library(boot)
#  data(prostate)
#  prostate <- prostate %>% mutate(epsilon=ifelse(status=="alive",0,
#                                                 ifelse(status=="dead - prostatic ca",1,
#                                                        ifelse(status=="dead - other ca",1,
#                                                               ifelse(status=="dead - heart or vascular",2,
#                                                                      ifelse(status=="dead - cerebrovascular",2,2)
#                                                               )))))
#  prostate$epsilon <- as.numeric(prostate$epsilon)
#  prostate$a <- as.numeric((prostate$rx=="placebo"))
#  prostate$t <- prostate$dtime/12
#  output <- polyreg(nuisance.model = Event(t,epsilon) ~ +1, exposure = 'a', strata='stage', data = prostate,
#                    effect.measure1='RR', effect.measure2='RR', time.point=1:5, outcome.type='PROPORTIONAL', boot.parameter1=1000)
#  tested_conf.low <- round(output$summary$event1$tidy$conf.low,digit=3)
#  tested_conf.high <- round(output$summary$event1$tidy$conf.high,digit=3)
#  tested <- as.vector(cbind(tested_conf.low, tested_conf.high))
#  expected <- c(-0.362, 0.335)
#  expected <- c(-0.380, 0.323)
#  expect_equal(tested, expected)
#})


test_that("calculateD produced the same D matrix for covariance as calculateD_old", {

  calculateD_old <- function(potential.CIFs, x_a, x_l, estimand, prob.bound) { #分類曝露未対応
    CIF1 <- ifelse(potential.CIFs[, 1] == 0, prob.bound, ifelse(potential.CIFs[, 1] == 1, 1 - prob.bound, potential.CIFs[, 1]))
    CIF2 <- ifelse(potential.CIFs[, 2] == 0, prob.bound, ifelse(potential.CIFs[, 2] == 1, 1 - prob.bound, potential.CIFs[, 2]))
    CIF3 <- ifelse(potential.CIFs[, 3] == 0, prob.bound, ifelse(potential.CIFs[, 3] == 1, 1 - prob.bound, potential.CIFs[, 3]))
    CIF4 <- ifelse(potential.CIFs[, 4] == 0, prob.bound, ifelse(potential.CIFs[, 4] == 1, 1 - prob.bound, potential.CIFs[, 4]))

    calculateA <- function(effect_measure, exposed.CIFs, unexposed.CIFs, a) {
      if (effect_measure == "RR") {
        return(a * (1 / exposed.CIFs) + (1 - a) * (1 / unexposed.CIFs))
      } else if (effect_measure == "OR") {
        return(a * (1 / exposed.CIFs + 1 / (1 - exposed.CIFs)) + (1 - a) * (1 / unexposed.CIFs + 1 / (1 - unexposed.CIFs)))
      } else if (effect_measure == "SHR") {
        tmp1_exposed <- -1 / (1 - exposed.CIFs)
        tmp1_unexposed <- -1 / (1 - unexposed.CIFs)
        tmp2_exposed <- log(1 - exposed.CIFs)
        tmp2_unexposed <- log(1 - unexposed.CIFs)
        return(a * (tmp1_exposed / tmp2_exposed) + (1 - a) * (tmp1_unexposed / tmp2_unexposed))
      } else {
        stop("Invalid effect_measure. Must be RR, OR or SHR.")
      }
    }

    a <- as.vector(x_a)
    a11 <- a12 <- a22 <- NULL
    a11 <- calculateA(estimand$effect.measure1, CIF2, CIF1, a)
    a22 <- calculateA(estimand$effect.measure2, CIF4, CIF3, a)
    #  a11 <- calculateA(estimand$effect.measure1, CIF3, CIF1, a)
    #  a22 <- calculateA(estimand$effect.measure2, CIF4, CIF2, a)
    a12 <- matrix(0, nrow = length(x_a), ncol = 1)

    d_ey_d_beta_11 <- a22 / (a11 * a22 - a12 * a12)
    d_ey_d_beta_12 <- -a12 / (a11 * a22 - a12 * a12)
    d_ey_d_beta_22 <- a11 / (a11 * a22 - a12 * a12)
    c12 <- a * (1 / (1 - CIF2 - CIF4)) + (1 - a) * (1 / (1 - CIF1 - CIF3))
    c11 <- a11 + c12
    c22 <- a22 + c12
    d_ey_d_alpha_11 <- c22 / (c11 * c22 - c12 * c12)
    d_ey_d_alpha_12 <- -c12 / (c11 * c22 - c12 * c12)
    d_ey_d_alpha_22 <- c11 / (c11 * c22 - c12 * c12)
    d_11 <- cbind((d_ey_d_alpha_11 * x_l), (d_ey_d_beta_11 * a))
    d_12 <- cbind((d_ey_d_alpha_12 * x_l), (d_ey_d_beta_12 * a))
    d_22 <- cbind((d_ey_d_alpha_22 * x_l), (d_ey_d_beta_22 * a))
    return(list(d_11 = d_11, d_12 = d_12, d_22 = d_22))
  }

set.seed(123)
n <- 10
x_a <- matrix(rbinom(n, 1, 0.5), ncol = 1)
x_l <- matrix(rnorm(n*2), ncol = 2)
estimand <- list(exposure.levels = 2, effect.measure1 = "SHR", effect.measure2 = "SHR")
prob.bound <- 1e-6

CIF1 <- runif(n, 0.05, 0.4)  # event1, level1
CIF2 <- runif(n, 0.05, 0.4)  # event1, level2
CIF3 <- runif(n, 0.05, 0.4)  # event2, level1
CIF4 <- runif(n, 0.05, 0.4)  # event2, level2
potential.CIFs <- cbind(CIF1, CIF2, CIF3, CIF4)

tested <- calculateD(potential.CIFs, x_a, x_l, estimand, prob.bound)
expected <- calculateD_old(potential.CIFs, x_a, x_l, estimand, prob.bound)

expect_equal(expected, tested)
})
