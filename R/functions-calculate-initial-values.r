getInitialValues <- function(formula, data, outcome.type, exposure, estimand, specific.time, data.initial.values, prob.bound) {
  cl <- match.call()
  mf <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "cluster", "offset")
  Terms <- terms(formula, special, data = data)
  mf$formula <- Terms
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Y <- model.extract(mf, "response")

  if (!inherits(Y, c("Event", "Surv"))) {
    if (outcome.type %in% c('COMPETINGRISK', 'SURVIVAL', 'PROPORTIONAL', 'POLY-PROPORTIONAL')) {
      stop("Expected a 'Surv' or 'Event'-object when outcome.type is COMPETINGRISK, SURVIVAL, PROPORTIONAL or POLY-PROPORTIONAL. ")
    } else {
      t <- rep(0, length(Y))
      epsilon <- Y
      if (!all(epsilon %in% c(estimand$code.event1, estimand$code.censoring))) {
        stop("Invalid event codes. Must be 0 or 1, if event codes are not specified. ")
      }
    }
  } else {
    t <- as.numeric(Y[, 1])
    if (any(t<0))
      stop("Invalid time variable. Expected non-negative values. ")
    if (any(is.na(t)))
      stop("Time variable contains NA values")
    if (any(is.na(Y[, 2]))) {
      stop("Event variable contains NA values")
    } else {
      epsilon <- Y[, 2]
    }
    if (!all(epsilon %in% c(estimand$code.event1, estimand$code.censoring)) & (outcome.type == 'SURVIVAL')) {
      stop("Invalid event codes. Must be 0 or 1, with 0 representing censoring, if event codes are not specified. ")
    } else if (!all(epsilon %in% c(estimand$code.event1, estimand$code.event2, estimand$code.censoring))) {
      stop("Invalid event codes. Must be 0, 1 or 2, with 0 representing censoring, if event codes are not specified. ")
    }
  }
  if (!is.null(offsetpos <- attributes(Terms)$specials$offset)) {
    ts <- survival::untangle.specials(Terms, "offset")
    if (length(ts$vars) > 0) {
      Terms <- Terms[-ts$terms]
      offset <- mf[[ts$vars]]
    } else {
      offset <- rep(0, nrow(mf))
    }
  } else {
    offset <- rep(0, nrow(mf))
  }

  a_ <- as.factor(data[[exposure]])
  if (estimand$code.exposure.ref==0) {
    x_a <- as.matrix(model.matrix(~ a_)[, -1])
  } else {
    x_a <- as.matrix(rep(1,length(t)) - model.matrix(~ a_)[, -1])
  }
  a <- as.vector(x_a)
  x_l <- model.matrix(Terms, mf)
  n_para_1 <- ncol(x_l)
  if (any(is.na(data[[exposure]])))
    stop("Exposure variable contains NA values")
  if (any(is.na(x_l)))
    stop("Covariates contain NA values")

  if (!is.null(data.initial.values)) {
    x_l <- model.matrix(Terms, mf)
    if (!(1+ncol(x_l))*2 == length(data.initial.values))
      stop("Invalid initial value dataset. Must contain the same number of initial values as parameters")
    return(data.initial.values)
  }

  binarizeIfContinuous <- function(x) {
    if (outcome.type == 'SURVIVAL' | outcome.type == 'BINOMIAL') {
      if (is.numeric(x) & length(unique(x)) > 2) {
        l <- as.numeric((x >= median(x)) == TRUE)
        out <- calculateInitialValuesSurvival(t, epsilon, a, l, estimand, specific.time, prob.bound)
        return(out)
      } else if (length(unique(x)) == 2) {
        l <- x
        out <- calculateInitialValuesSurvival(t, epsilon, a, l, estimand, specific.time, prob.bound)
        return(out)
      }
    } else {
      if (is.numeric(x) & length(unique(x)) > 2) {
        l <- as.numeric((x >= median(x)) == TRUE)
        out <- calculateInitialValuesCompetingRisk(t, epsilon, a, l, estimand, specific.time, prob.bound)
        return(out)
      } else if (length(unique(x)) == 2) {
        l <- x
        out <- calculateInitialValuesCompetingRisk(t, epsilon, a, l, estimand, specific.time, prob.bound)
        return(out)
      }
    }
  }

  if (all(epsilon %in% c(estimand$code.event1, estimand$code.censoring))) {
    if (n_para_1>1) {
      out_bic_1 <- t(binarizeIfContinuous(x_l[,2])[1,1:2])
      if (n_para_1>2) {
        for (i_para in 3:n_para_1) {
          out_bic_i <- binarizeIfContinuous(x_l[,i_para])
          out_bic_1 <- cbind(out_bic_1, t(out_bic_i[1,2]))
        }
        out_bic_1 <- cbind(out_bic_1, t(out_bic_i[1,3]))
      } else {
        out_bic_1 <- cbind(out_bic_1, t(binarizeIfContinuous(x_l[,2])[1,3]))
      }
      init_vals <- out_bic_1
    } else {
      l <- NULL
      init_vals <- calculateInitialValuesSurvival(t, epsilon, a, l, estimand, specific.time, prob.bound)
    }
  } else {
    if (n_para_1>1) {
      out_bic_1 <- t(binarizeIfContinuous(x_l[,2])[1,1:2])
      out_bic_2 <- t(binarizeIfContinuous(x_l[,2])[1,4:5])
      if (n_para_1>2) {
        for (i_para in 3:n_para_1) {
          out_bic_i <- binarizeIfContinuous(x_l[,i_para])
          out_bic_1 <- cbind(out_bic_1, t(out_bic_i[1,2]))
          out_bic_2 <- cbind(out_bic_2, t(out_bic_i[1,5]))
        }
        out_bic_1 <- cbind(out_bic_1, t(out_bic_i[1,3]))
        out_bic_2 <- cbind(out_bic_2, t(out_bic_i[1,6]))
      } else {
        out_bic_1 <- cbind(out_bic_1, t(binarizeIfContinuous(x_l[,2])[1,3]))
        out_bic_2 <- cbind(out_bic_2, t(binarizeIfContinuous(x_l[,2])[1,6]))
      }
      init_vals <- cbind(out_bic_1, out_bic_2)
    } else {
      l <- NULL
      init_vals <- calculateInitialValuesCompetingRisk(t, epsilon, a, l, estimand, specific.time, prob.bound)
    }
  }
  return(init_vals)
}

calculateInitialValuesCompetingRisk <- function(t, epsilon, a, l = NULL, estimand, specific.time, prob.bound) {
  epsilon0 <- epsilon[a == 0]
  epsilon1 <- epsilon[a == 1]
  t0 <- t[a == 0]
  t1 <- t[a == 1]

  p_10 <- (sum(epsilon0 == estimand$code.event1 & t0 <= specific.time) / length(epsilon0)) + prob.bound
  p_20 <- (sum(epsilon0 == estimand$code.event2 & t0 <= specific.time) / length(epsilon0)) + prob.bound
  p_00 <- 1 - p_10 - p_20
  p_11 <- (sum(epsilon1 == estimand$code.event1 & t1 <= specific.time) / length(epsilon1)) + prob.bound
  p_21 <- (sum(epsilon1 == estimand$code.event2 & t1 <= specific.time) / length(epsilon1)) + prob.bound
  p_01 <- 1 - p_11 - p_21

  alpha_1 <- log((p_10 * p_11) / (p_00 * p_01))
  alpha_2 <- log((p_20 * p_21) / (p_00 * p_01))
  if(!is.null(l)){
    epsilon00 <- epsilon[a == 0 & l == 0]
    epsilon10 <- epsilon[a == 1 & l == 0]
    epsilon01 <- epsilon[a == 0 & l == 1]
    epsilon11 <- epsilon[a == 1 & l == 1]
    t00 <- t[a == 0 & l == 0]
    t10 <- t[a == 1 & l == 0]
    t01 <- t[a == 0 & l == 1]
    t11 <- t[a == 1 & l == 1]

    p_100 <- (sum(epsilon00 == estimand$code.event1 & t00 <= specific.time) / length(epsilon00)) + prob.bound
    p_200 <- (sum(epsilon00 == estimand$code.event2 & t00 <= specific.time) / length(epsilon00)) + prob.bound
    p_000 <- 1 - p_100 - p_200
    p_110 <- (sum(epsilon10 == estimand$code.event1 & t10 <= specific.time) / length(epsilon10)) + prob.bound
    p_210 <- (sum(epsilon10 == estimand$code.event2 & t10 <= specific.time) / length(epsilon10)) + prob.bound
    p_010 <- 1 - p_110 - p_210
    p_101 <- (sum(epsilon01 == estimand$code.event1 & t01 <= specific.time) / length(epsilon01)) + prob.bound
    p_201 <- (sum(epsilon01 == estimand$code.event2 & t01 <= specific.time) / length(epsilon01)) + prob.bound
    p_001 <- 1 - p_101 - p_201
    p_111 <- (sum(epsilon11 == estimand$code.event1 & t11 <= specific.time) / length(epsilon11)) + prob.bound
    p_211 <- (sum(epsilon11 == estimand$code.event2 & t11 <= specific.time) / length(epsilon11)) + prob.bound
    p_011 <- 1 - p_111 - p_211

    alpha_10 <- log((p_100 * p_110) / (p_000 * p_010))
    alpha_20 <- log((p_200 * p_210) / (p_000 * p_010))
    alpha_11 <- log((p_101 * p_111) / (p_000 * p_011)) - alpha_10
    alpha_21 <- log((p_201 * p_211) / (p_000 * p_011)) - alpha_20
  }

  if (estimand$effect.measure1 == 'RR') {
    beta_1 <- log(p_11 / p_10)
  } else if (estimand$effect.measure1 == 'OR') {
    beta_1 <- log((p_11 / (1 - p_11)) / (p_10 / (1 - p_10)))
  } else if (estimand$effect.measure1 == 'SHR') {
    beta_1 <- log(log(1 - p_11) / log(1 - p_10))
  } else {
    stop("Invalid effect measure code. Must be RR, OR or SHR.")
  }

  if (estimand$effect.measure2 == 'RR') {
    beta_2 <- log(p_21 / p_20)
  } else if (estimand$effect.measure2 == 'OR') {
    beta_2 <- log((p_21 / (1 - p_21)) / (p_20 / (1 - p_20)))
  } else if (estimand$effect.measure2 == 'SHR') {
    beta_2 <- log(log(1 - p_21) / log(1 - p_20))
  } else {
    stop("Invalid effect measure code. Must be RR, OR or SHR.")
  }

  alpha_beta <- if (is.null(l)) {
    cbind(alpha_1, beta_1, alpha_2, beta_2)
  } else {
    cbind(alpha_10, alpha_11, beta_1, alpha_20, alpha_21, beta_2)
  }
  return(alpha_beta)
}

calculateInitialValuesSurvival <- function(t, epsilon, a, l = NULL, estimand, specific.time, prob.bound) {
  if (is.null(l)) {
    epsilon0 <- epsilon[a == 0]
    epsilon1 <- epsilon[a == 1]
    t0 <- t[a == 0]
    t1 <- t[a == 1]
    p_10 <- (sum(epsilon0 == estimand$code.event1 & t0 <= specific.time) / length(epsilon0)) + prob.bound
    p_00 <- 1 - p_10
    p_11 <- (sum(epsilon1 == estimand$code.event1 & t1 <= specific.time) / length(epsilon1)) + prob.bound
    p_01 <- 1 - p_11

    if (estimand$effect.measure1 == 'RR') {
      beta_1 <- log(p_11 / p_10)
    } else if (estimand$effect.measure1 == 'OR') {
      beta_1 <- log((p_11 / (1 - p_11)) / (p_10 / (1 - p_10)))
    } else if (estimand$effect.measure1 == 'SHR') {
      beta_1 <- log(log(1 - p_11) / log(1 - p_10))
    } else {
      stop("Invalid effect measure code. Must be RR, OR or SHR.")
    }

    alpha_1 <- log((p_10 * p_11) / (p_00 * p_01))
    return(cbind(alpha_1, beta_1))
  } else {
    epsilon00 <- epsilon[a == 0 & l == 0]
    epsilon10 <- epsilon[a == 1 & l == 0]
    epsilon01 <- epsilon[a == 0 & l == 1]
    epsilon11 <- epsilon[a == 1 & l == 1]
    t00 <- t[a == 0 & l == 0]
    t10 <- t[a == 1 & l == 0]
    t01 <- t[a == 0 & l == 1]
    t11 <- t[a == 1 & l == 1]

    p_100 <- (sum(epsilon00 == estimand$code.event1 & t00 <= specific.time) / length(epsilon00)) + prob.bound
    p_000 <- 1 - p_100
    p_110 <- (sum(epsilon10 == estimand$code.event1 & t10 <= specific.time) / length(epsilon10)) + prob.bound
    p_010 <- 1 - p_110
    p_101 <- (sum(epsilon01 == estimand$code.event1 & t01 <= specific.time) / length(epsilon01)) + prob.bound
    p_001 <- 1 - p_101
    p_111 <- (sum(epsilon11 == estimand$code.event1 & t11 <= specific.time) / length(epsilon11)) + prob.bound
    p_011 <- 1 - p_111

    if (any(c(p_100, p_000, p_110, p_010, p_101, p_001, p_111, p_011) == 0, na.rm = TRUE)) {
      stop("Complete separation detected in initial value search")
    }

    alpha_10 <- log((p_100 * p_110) / (p_000 * p_010))
    alpha_11 <- log((p_101 * p_111) / (p_000 * p_011)) - alpha_10

    epsilon0 <- epsilon[a == 0]
    epsilon1 <- epsilon[a == 1]
    t0 <- t[a == 0]
    t1 <- t[a == 1]
    p_10 <- (sum(epsilon0 == estimand$code.event1 & t0 <= specific.time) / length(epsilon0)) + prob.bound
    p_00 <- 1 - p_10
    p_11 <- (sum(epsilon1 == estimand$code.event1 & t1 <= specific.time) / length(epsilon1)) + prob.bound
    p_01 <- 1 - p_11

    if (estimand$effect.measure1 == 'RR') {
      beta_1 <- log(p_11 / p_10)
    } else if (estimand$effect.measure1 == 'OR') {
      beta_1 <- log((p_11 / (1 - p_11)) / (p_10 / (1 - p_10)))
    } else if (estimand$effect.measure1 == 'SHR') {
      beta_1 <- log(log(1 - p_11) / log(1 - p_10))
    } else {
      stop("Invalid effect measure code. Must be RR, OR or SHR.")
    }

    return(cbind(alpha_10, alpha_11, beta_1))
  }
}

getInitialValuesProportional <- function(formula, data, outcome.type, exposure, estimand, data.initial.values, prob.bound, out_normalizeCovariate) {
  cl <- match.call()
  mf <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "cluster", "offset")
  Terms <- terms(formula, special, data = data)
  mf$formula <- Terms
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Y <- model.extract(mf, "response")
  if (!inherits(Y, c("Event", "Surv")))
    stop("Expected a 'Surv' or 'Event'-object")
  if (ncol(Y) == 2) {
    t <- Y[, 1]
    epsilon <- Y[, 2]
  }

  n_para_1 <- out_normalizeCovariate$n_covariate + 1
  n_para_2 <- out_normalizeCovariate$n_covariate + 2
  n_para_3 <- out_normalizeCovariate$n_covariate + 3
  n_para_4 <- 2*out_normalizeCovariate$n_covariate + 3
  n_para_5 <- 2*out_normalizeCovariate$n_covariate + 4
  n_para_6 <- length(estimand$time.point) * (n_para_5 - 2) + 2

  if (outcome.type == "PROPORTIONAL") {
    alpha_beta_0 <- rep(NA_real_, n_para_6/2)
  } else {
    alpha_beta_0 <- rep(NA_real_, n_para_6)
  }

  sum1 <- 0
  sum2 <- 0
  i_time <- 0

  for (specific.time in estimand$time.point) {
    specific.time <- as.numeric(specific.time)
    out_getInitialValues <- getInitialValues(
      formula = formula,
      data = data,
      outcome.type = outcome.type,
      exposure = exposure,
      estimand = estimand,
      specific.time = specific.time,
      data.initial.values = data.initial.values,
      prob.bound = prob.bound
    )

    i_time <- i_time + 1
    i_para <- n_para_1*(i_time - 1) + 1

    tmp1 <- out_getInitialValues[seq.int(1, length.out = n_para_1)]
    if (outcome.type == "PROPORTIONAL") {
      idx1 <- seq.int(i_para, length.out = n_para_1)
      alpha_beta_0[idx1] <- tmp1
    } else {
      idx1 <- seq.int(i_para, length.out = n_para_1)
      alpha_beta_0[idx1] <- tmp1
    }
    sum1 <- sum1 + out_getInitialValues[n_para_2]

    if (outcome.type == "POLY-PROPORTIONAL") {
      tmp2 <- out_getInitialValues[seq.int(n_para_3, n_para_4)]
      idx2_start <- (n_para_6 %/% 2) + i_para
      idx2 <- seq.int(idx2_start, length.out = n_para_1)
      alpha_beta_0[idx2] <- tmp2
      sum2 <- sum2 + out_getInitialValues[n_para_5]
    }
  }
  if (outcome.type == "PROPORTIONAL") {
    alpha_beta_0[n_para_6 %/% 2] <- sum1 / length(estimand$time.point)
  } else {
    alpha_beta_0[n_para_6 %/% 2] <- sum1 / length(estimand$time.point)
    alpha_beta_0[n_para_6]       <- sum2 / length(estimand$time.point)
  }
  return(alpha_beta_0)
}
