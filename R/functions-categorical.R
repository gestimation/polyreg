calculatePotentialCIFs_categorical <- function(
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

  ###########################################################
  alpha_1 <- alpha_beta_tmp[1:i_parameter[1]]
  beta_1  <- alpha_beta_tmp[i_parameter[2]:i_parameter[3]]
  alpha_2 <- alpha_beta_tmp[i_parameter[4]:i_parameter[5]]
  beta_2  <- alpha_beta_tmp[i_parameter[6]:i_parameter[7]]
  alpha_tmp_1 <- x_l %*% as.matrix(alpha_1) + offset
  alpha_tmp_2 <- x_l %*% as.matrix(alpha_2) + offset
  beta_tmp_1  <- x_a %*% as.matrix(beta_1)
  beta_tmp_2  <- x_a %*% as.matrix(beta_2)
  #  alpha_1 <- alpha_beta_tmp[1:i_parameter[1]]
  #  alpha_tmp_1 <- x_l %*% as.matrix(alpha_1) + offset
  #  beta_tmp_1  <- alpha_beta_tmp[i_parameter[2]:i_parameter[3]]
  #  alpha_2 <- alpha_beta_tmp[i_parameter[4]:i_parameter[5]]
  #  alpha_tmp_2 <- x_l %*% as.matrix(alpha_2) + offset
  #  beta_tmp_2  <- alpha_beta_tmp[i_parameter[6]:i_parameter[7]]
  ###########################################################

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
        estimating_equation_CIFs_categorical(
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
        estimating_equation_CIFs_categorical(
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
  ################################################
  #  potential.CIFs <- do.call(rbind, list.CIFs)
  #  return(potential.CIFs)

  observed.CIFs <- matrix(NA,nrow(x_l),2)
  for (i_x in seq_len(nrow(x_l))) {
    tmp1 <- list.CIFs[[i_x]]
    tmp2 <- cbind(1-sum(x_a[i_x,]),x_a[i_x,])
    observed.CIFs[i_x,1] <- tmp2 %*% tmp1[1:(length(tmp1)/2)]
    observed.CIFs[i_x,2] <- tmp2 %*% tmp1[(length(tmp1)/2+1):length(tmp1)]
  }
  return(observed.CIFs)
}


estimating_equation_CIFs_categorical <- function(
    log_p,      # c(log_p10, log_p20, log_p11, log_p21) > c(log_p10, log_p11, log_p20, log_p21)
    alpha_tmp_1, beta_tmp_1,
    alpha_tmp_2, beta_tmp_2,
    estimand, optim.method, prob.bound
) {
  objective_function <- function(log_p) {
    ret <- numeric(length(log_p))
    clog_p <- numeric(length(log_p))

    clamp_log_p <- function(lp) {
      val <- exp(lp)
      if (val < prob.bound) {
        return(log(prob.bound))
      } else if ((1 - val) < prob.bound) {
        return(log(1 - prob.bound))
      } else {
        return(lp)
      }
    }
    for (i in seq_len(length(log_p))) {
      clog_p[i] <- clamp_log_p(log_p[i])
    }

    log_rr1 <- log_p[2:(length(clog_p)/2)]/clog_p[1]
    log_rr2 <- log_p[(length(clog_p)/2+2):length(clog_p)]/clog_p[length(clog_p)/2+1]
    plor_1n <- sum(clog_p[1:(length(clog_p)/2)])
    plor_2n <- sum(clog_p[(length(clog_p)/2+1):length(clog_p)])
    p <- exp(clog_p)
    p01 <- 1-sum(p[1:(length(clog_p)/2)])
    p02 <- 1-sum(p[(length(clog_p)/2+1):length(clog_p)])
    log_p01 <- log(p01)
    log_p02 <- log(p02)
    plor_d <- log_p01 + log_p02

    if ((1 - sum(p[1:(length(clog_p)/2)]) < prob.bound) |
        (1 - sum(p[(length(clog_p)/2+1):length(clog_p)]) < prob.bound)) {
      plor_d <- log(prob.bound)
      #      lp0102 <- log(prob.bound)
    } else {
      p01 <- 1-sum(p[1:(length(clog_p)/2)])
      p02 <- 1-sum(p[(length(clog_p)/2+1):length(clog_p)])
      log_p01 <- log(p01)
      log_p02 <- log(p02)
      plor_d <- log_p01 + log_p02
    }
    plor_1 <- plor_1n - plor_d
    plor_2 <- plor_2n - plor_d

    log_p1 <- (log_p[1])
    log_p2 <- (log_p[2])
    log_p3 <- (log_p[3])
    log_p4 <- (log_p[4])

    if (estimand$effect.measure1 == 'RR') {
      #ret[1] <- alpha_tmp_1 - plor_1
      #ret[2] <- beta_tmp_1  - log_rr1
      ret[1] <- alpha_tmp_1 - log_p1 - log_p2 + plor_d
      ret[2] <- beta_tmp_1  - log_p2 + log_p1
      #      ret[1] <- alpha_tmp_1 - log_p1 - log_p3 + lp0102
      #      ret[2] <- beta_tmp_1  - log_p3 + log_p1
    } else if (estimand$effect.measure1 == 'OR') {
      ret[1] <- alpha_tmp_1 - plor_1
      ret[2] <- beta_tmp_1  - log_p3 + log_p1 + log(1 - exp_lp3) - log(1 - exp_lp1)
    } else if (estimand$effect.measure1 == 'SHR') {
      ret[1] <- alpha_tmp_1 - plor_1
      ret[2] <- exp(beta_tmp_1) - ( log(1 - exp_lp3) / log(1 - exp_lp1) )
    } else {
      stop("Invalid effect_measure. Must be RR, OR or SHR.")
    }

    if (estimand$effect.measure2 == 'RR') {
      #ret[3] <- alpha_tmp_2 - plor_2
      #ret[4] <- beta_tmp_2  - log_rr2
      ret[3] <- alpha_tmp_1 - log_p3 - log_p4 + plor_d
      ret[4] <- beta_tmp_1  - log_p4 + log_p3
    } else if (estimand$effect.measure2 == 'OR') {
      ret[3] <- alpha_tmp_2 - plor_2
      ret[4] <- beta_tmp_2  - log_p4 + log_p2 + log(1 - exp_lp4) - log(1 - exp_lp2)
    } else if (estimand$effect.measure2 == 'SHR') {
      ret[3] <- alpha_tmp_2 - plor_2
      ret[4] <- exp(beta_tmp_2) - ( log(1 - exp_lp4) / log(1 - exp_lp2) )
    } else {
      stop("Invalid effect_measure. Must be RR, OR or SHR.")
    }
    if (optim.method$inner.optim.method == 'multiroot'){
      return(ret)
    } else {
      return(sum(ret^2))
    }
  }
  return(objective_function(log_p))
}


estimating_equation_ipcw_categorical <- function(
    formula,
    data,
    exposure,
    ip.weight,
    alpha_beta,
    estimand,
    optim.method,
    prob.bound,
    initial.CIFs = NULL
) {
  cl <- match.call()
  mf <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "cluster", "offset")
  Terms <- terms(formula, special, data = data)
  mf$formula <- Terms
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Y <- model.extract(mf, "response")
  t <- Y[, 1]  # time variable
  epsilon <- Y[, 2]  # status variable
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

  y_0 <- ifelse(epsilon == estimand$code.censoring | t > estimand$time.point, 1, 0)
  y_1 <- ifelse(epsilon == estimand$code.event1 & t <= estimand$time.point, 1, 0)
  y_2 <- ifelse(epsilon == estimand$code.event2 & t <= estimand$time.point, 1, 0)
  y_0_ <- ifelse(epsilon == estimand$code.censoring, 1, 0)
  y_1_ <- ifelse(epsilon == estimand$code.event1, 1, 0)
  y_2_ <- ifelse(epsilon == estimand$code.event2, 1, 0)

  a_ <- as.factor(data[[exposure]])
  if (estimand$code.exposure.ref==0) {
    x_a <- as.matrix(model.matrix(~ a_)[, -1])
  } else {
    x_a <- as.matrix(rep(1,length(t)) - model.matrix(~ a_)[, -1])
  }
  x_l <- model.matrix(Terms, mf)
  x_la <- cbind(x_l, x_a)
  i_parameter <- rep(NA, 7)
  i_parameter <- calculateIndexForParameter(i_parameter,x_l,x_a)

  #################################################################
  #  potential.CIFs <- calculatePotentialCIFs(alpha_beta,x_a,x_l,offset,epsilon,estimand,optim.method,prob.bound,initial.CIFs)
  #  one <- rep(1, nrow(x_l))
  #  a <- as.vector(x_a)
  #  ey_1 <- potential.CIFs[,3]*a + potential.CIFs[,1]*(one - a)
  #  ey_2 <- potential.CIFs[,4]*a + potential.CIFs[,2]*(one - a)
  observed.CIFs <- calculatePotentialCIFs_categorical(alpha_beta,x_a,x_l,offset,epsilon,estimand,optim.method,prob.bound,initial.CIFs)
  ey_1 <- observed.CIFs[,1]
  ey_2 <- observed.CIFs[,2]
  #################################################################

  v11 <- ey_1 * (1 - ey_1)
  v12 <- -ey_1 * ey_2
  v22 <- ey_2 * (1 - ey_2)
  denom <- v11*v22 - v12*v12
  w11 <- v22 / denom
  w12 <- -v12 / denom
  w22 <- v11 / denom
  wy_1 <- ip.weight * y_1
  wy_2 <- ip.weight * y_2
  wy_1ey_1 <- w11*(wy_1 - ey_1) + w12*(wy_2 - ey_2)
  wy_2ey_2 <- w12*(wy_1 - ey_1) + w22*(wy_2 - ey_2)

  x_la <- cbind(x_l, x_a)
  zero <- matrix(0, nrow=nrow(x_la), ncol=ncol(x_la))
  tmp1 <- cbind(x_la, zero)
  tmp2 <- cbind(zero, x_la)
  d    <- rbind(tmp1, tmp2)
  residual <- c(wy_1ey_1, wy_2ey_2)
  ret <- as.vector(t(d) %*% residual / nrow(x_l))

  n_col_d <- ncol(d)
  score.matrix <- matrix(NA, nrow=nrow(d), ncol=n_col_d)
  for (j in seq_len(n_col_d)) {
    score.matrix[,j] <- d[,j]*residual
  }
  colnames(potential.CIFs) <- c("p10", "p20", "p11", "p21")
  out <- list(
    ret   = ret,
    score = score.matrix,
    ey_1  = ey_1,
    ey_2  = ey_2,
    w11   = w11,
    w12   = w12,
    w22   = w22,
    t     = t,
    y_0   = y_0,
    y_1   = y_1,
    y_2   = y_2,
    y_0_  = y_0_,
    y_1_  = y_1_,
    y_2_  = y_2_,
    x_a   = x_a,
    x_l = x_l,
    potential.CIFs = potential.CIFs,
    ip.weight = ip.weight,
    i_parameter = i_parameter
  )
  return(out)
}

