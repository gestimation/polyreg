calculatePotentialCIFs <- function(
    alpha_beta_tmp,
    x_a,
    x_l,
    offset,
    epsilon,
    estimand,
    optim.method,
    prob.bound,
    initial.CIFs = NULL
) {
  i_parameter <- rep(NA, 7)
  i_parameter <- calculateIndexForParameter(i_parameter,x_l,x_a)
  alpha_1 <- alpha_beta_tmp[1:i_parameter[1]]
  alpha_tmp_1 <- x_l %*% as.matrix(alpha_1) + offset
  beta_tmp_1  <- alpha_beta_tmp[i_parameter[2]:i_parameter[3]]
  alpha_2 <- alpha_beta_tmp[i_parameter[4]:i_parameter[5]]
  alpha_tmp_2 <- x_l %*% as.matrix(alpha_2) + offset
  beta_tmp_2  <- alpha_beta_tmp[i_parameter[6]:i_parameter[7]]

  n <- length(epsilon)
  p0 <- c(
    sum(epsilon == estimand$code.event1) / n + prob.bound,
    sum(epsilon == estimand$code.event2) / n + prob.bound,
    sum(epsilon == estimand$code.event1) / n + prob.bound,
    sum(epsilon == estimand$code.event2) / n + prob.bound
  )
  log_p0 <- log(p0)
  list.CIFs <- vector("list", nrow(x_l))
  previous.CIFs <- NULL
  for (i_x in seq_len(nrow(x_l))) {
    # Skip the calculation if the current observation is the same as the previous one
    if (!is.null(previous.CIFs) && all(x_l[i_x, ] == x_l[i_x-1, ])) {
      list.CIFs[[i_x]] <- previous.CIFs
      next
    }

    # Use the previous prediction value of observation i_x if initial.CIFs is not NULL
    if (!is.null(initial.CIFs)) {
      log_p0 <- log(initial.CIFs[i_x, ])
    }

    if (optim.method$inner.optim.method == 'optim' | optim.method$inner.optim.method == 'BFGS') {
      eq_fn <- function(lp) {
        estimating_equation_CIFs(
          log_p          = lp,
          alpha_tmp_1    = alpha_tmp_1[i_x],
          beta_tmp_1     = beta_tmp_1,
          alpha_tmp_2    = alpha_tmp_2[i_x],
          beta_tmp_2     = beta_tmp_2,
          estimand       = estimand,
          optim.method   = optim.method,
          prob.bound     = prob.bound
        )
      }
      sol <- optim(par = log_p0, fn = eq_fn, method = "BFGS", control = list(maxit=optim.method$optim.parameter7, reltol=optim.method$optim.parameter6))
      list.CIFs[[i_x]] <- exp(sol$par)
      previous.CIFs <- list.CIFs[[i_x]]
    } else if (optim.method$inner.optim.method == 'SANN') {
      eq_fn <- function(lp) {
        estimating_equation_CIFs(
          log_p          = lp,
          alpha_tmp_1    = alpha_tmp_1[i_x],
          beta_tmp_1     = beta_tmp_1,
          alpha_tmp_2    = alpha_tmp_2[i_x],
          beta_tmp_2     = beta_tmp_2,
          estimand       = estimand,
          optim.method   = optim.method,
          prob.bound     = prob.bound
        )
      }
      sol <- optim(par = log_p0, fn = eq_fn, method = "SANN", control = list(maxit=optim.method$optim.parameter7, reltol=optim.method$optim.parameter6))
      list.CIFs[[i_x]] <- exp(sol$par)
      previous.CIFs <- list.CIFs[[i_x]]
    }
  }
  potential.CIFs <- do.call(rbind, list.CIFs)
  return(potential.CIFs)
}

estimating_equation_CIFs <- function(
    log_p,
    alpha_tmp_1,
    beta_tmp_1,
    alpha_tmp_2,
    beta_tmp_2,
    estimand,
    optim.method,
    prob.bound
) {
  objective_function <- function(log_p) {
    clampLogP <- function(x) {
      exp_x <- exp(x)
      if (exp_x < prob.bound) {
        return(log(prob.bound))
      } else if ((1 - exp_x) < prob.bound) {
        return(log(1 - prob.bound))
      } else {
        return(x)
      }
    }
    clog_p <- rep(NA,4)
    clog_p[1] <- clampLogP(log_p[1])
    clog_p[2] <- clampLogP(log_p[2])
    clog_p[3] <- clampLogP(log_p[3])
    clog_p[4] <- clampLogP(log_p[4])
    exp_lp <- exp(clog_p)

    if ((1 - exp_lp[1] - exp_lp[2] < prob.bound) |
        (1 - exp_lp[3] - exp_lp[4] < prob.bound)) {
      lp0102 <- log(prob.bound)
    } else {
      lp01   <- log(abs(1 - exp_lp[1] - exp_lp[2]))
      lp02   <- log(abs(1 - exp_lp[3] - exp_lp[4]))
      lp0102 <- lp01 + lp02
    }

    ret <- numeric(4)
    if (estimand$effect.measure1 == 'RR') {
      ret[1] <- alpha_tmp_1 - clog_p[1] - clog_p[3] + lp0102
      ret[2] <- beta_tmp_1  - clog_p[3] + clog_p[1]
    } else if (estimand$effect.measure1 == 'OR') {
      ret[1] <- alpha_tmp_1 - clog_p[1] - clog_p[3] + lp0102
      ret[2] <- beta_tmp_1  - clog_p[3] + clog_p[1] +
        log(1 - exp_lp[3]) - log(1 - exp_lp[1])
    } else if (estimand$effect.measure1 == 'SHR') {
      ret[1] <- alpha_tmp_1 - clog_p[1] - clog_p[3] + lp0102
      ret[2] <- exp(beta_tmp_1) - ( log(1 - exp_lp[3]) / log(1 - exp_lp[1]) )
    } else {
      stop("Invalid effect_measure. Must be RR, OR or SHR.")
    }

    if (estimand$effect.measure2 == 'RR') {
      ret[3] <- alpha_tmp_2 - clog_p[2] - clog_p[4] + lp0102
      ret[4] <- beta_tmp_2  - clog_p[4] + clog_p[2]
    } else if (estimand$effect.measure2 == 'OR') {
      ret[3] <- alpha_tmp_2 - clog_p[2] - clog_p[4] + lp0102
      ret[4] <- beta_tmp_2  - clog_p[4] + clog_p[2] +
        log(1 - exp_lp[4]) - log(1 - exp_lp[2])
    } else if (estimand$effect.measure2 == 'SHR') {
      ret[3] <- alpha_tmp_2 - clog_p[2] - clog_p[4] + lp0102
      ret[4] <- exp(beta_tmp_2) - ( log(1 - exp_lp[4]) / log(1 - exp_lp[2]) )
    } else {
      stop("Invalid effect_measure. Must be RR, OR or SHR.")
    }
    return(sum(ret^2))
  }
  return(objective_function(log_p))
}

calculatePotentialRisk <- function(alpha_beta, x_l, offset, estimand) {
  if (estimand$effect.measure1 == 'RR') {
    one <- rep(1, nrow(x_l))
    n_para_1 <- ncol(x_l)
    n_para_2 <- ncol(x_l) + 1
    alpha_beta_ <- as.matrix(as.vector(alpha_beta))
    if (ncol(alpha_beta_ == 1)) {
      alpha_beta_ <- t(alpha_beta_)
    }
    phi <- x_l %*% alpha_beta_[, 1:n_para_1] + offset
    theta <- one * alpha_beta_[, n_para_2]
    expphi <- exp(phi)
    exptheta <- exp(theta)
    if (all(phi == 0)) {
      p_10 <- one / (one + exptheta)
      p_11 <- exptheta * p_10
    } else {
      denomi_1 <- -(exptheta + one) * expphi
      denomi_2 <- sqrt(exp(2 * phi) * (exptheta + one) * (exptheta + one) + 4 * exp(theta + phi) * (one - expphi))
      denomi <- denomi_1 + denomi_2
      numera <- 2 * exptheta * (one - expphi)
      p_10 <- denomi / numera
      p_11 <- exptheta * p_10
    }
  } else if (estimand$effect.measure1 == 'OR') {
    one <- rep(1, nrow(x_l))
    n_para_1 <- ncol(x_l)
    n_para_2 <- ncol(x_l) + 1
    phi <- x_l %*% as.matrix(alpha_beta)[, 1:n_para_1] + offset
    theta <- one * as.matrix(alpha_beta)[, n_para_2]
    sqrt1 <- sqrt(exp(-theta - phi))
    sqrt2 <- sqrt(exp(theta - phi))
    if (all(phi == theta)) {
      p_10 <- 0.5 * one
      p_11 <- one/(one + sqrt1)
    } else {
      p_10 <- one/(one + sqrt2)
      p_11 <- one/(one + sqrt1)
    }
  } else {
    stop("Invalid effect_measure. Must be RR or OR.")
  }
  return(cbind(p_10, p_11))
}


