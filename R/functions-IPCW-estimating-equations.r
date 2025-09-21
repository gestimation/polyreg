estimating_equation_ipcw <- function(
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
  i_parameter <- calculateIndexForParameter(NA,x_l,x_a)

  if (optim.method$computation.order.method=="SEQUENTIAL") {
    potential.CIFs <- calculatePotentialCIFs_old(alpha_beta,x_a,x_l,offset,epsilon,estimand,optim.method,prob.bound,initial.CIFs)
  } else if (optim.method$computation.order.method=="PARALLEL") {
#    potential.CIFs <- calculatePotentialCIFs_parallel(alpha_beta,x_a,x_l,offset,epsilon,estimand,i_parameter,optim.method,prob.bound,initial.CIFs)
    potential.CIFs <- calculatePotentialCIFs_LM(alpha_beta,x_a,x_l,offset,epsilon,estimand,optim.method,prob.bound,initial.CIFs)
  } else {
    potential.CIFs <- calculatePotentialCIFs_tinyLM(alpha_beta,x_a,x_l,offset,epsilon,estimand,optim.method,prob.bound,initial.CIFs)
  }
  one <- rep(1, nrow(x_l))
  a <- as.vector(x_a)
  ey_1 <- potential.CIFs[,3]*a + potential.CIFs[,1]*(one - a)
  ey_2 <- potential.CIFs[,4]*a + potential.CIFs[,2]*(one - a)

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
  ret <- drop(crossprod(d, residual)) / nrow(x_l)
  #ret <- as.vector(t(d) %*% residual / nrow(x_l))

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

calculateCov <- function(objget_results, estimand, prob.bound)
{
  score <- objget_results$score
  ey_1 <- objget_results$ey_1
  ey_2 <- objget_results$ey_2
  w11 <- objget_results$w11
  w12 <- objget_results$w12
  w22 <- objget_results$w22
  t <- objget_results$t
  y_0_ <- objget_results$y_0_
  y_1 <- objget_results$y_1
  y_2 <- objget_results$y_2
  x_a <- objget_results$x_a
  x_l <- objget_results$x_l
  potential.CIFs <- objget_results$potential.CIFs
  i_parameter <- objget_results$i_parameter
  n <- length(t)

  censoring_dna	<- calculateNelsonAalen(t,y_0_)
  censoring_martingale <- diag(y_0_) - (outer(t, t, ">") * censoring_dna)
  censoring_km <- calculateKaplanMeier(t, y_0_)
  censoring_km[censoring_km == 0] <- 1e-5
  censoring_mkm <- censoring_martingale / censoring_km
  y_12 <- (y_1 + y_2 > 0)
  survival_km <- calculateKaplanMeier(t, y_12)
  wy_1 <- w11 * (y_1 - ey_1) + w12 * (y_2 - ey_2)
  wy_2 <- w12 * (y_1 - ey_1) + w22 * (y_2 - ey_2)
  x_la <- cbind(x_l, x_a)
  AB1 <- score[1:n, 1:i_parameter[2]]
  AB2 <- score[(n + 1):(2 * n), i_parameter[3]:i_parameter[5]]
  for (i_para in 1:i_parameter[2]) {
    tmp0 <- x_la[, i_para]
    use <- (t <= estimand$time.point)
    tmp1 <- colSums((use * tmp0) * (outer(t, t, ">=") * wy_1))
    tmp1 <- tmp1 / survival_km / n
    integrand1 <- tmp1 * censoring_mkm
    tmp2 <- colSums((use * tmp0) * (outer(t, t, ">=") * wy_2))
    tmp2 <- tmp2 / survival_km / n
    integrand2 <- tmp2 * censoring_mkm
    for (i_score in 1:n) {
      integral1 <- cumsum(integrand1[i_score, ])
      AB1[i_score, i_para] <- AB1[i_score, i_para] + integral1[n]
      integral2 <- cumsum(integrand2[i_score, ])
      AB2[i_score, i_para] <- AB2[i_score, i_para] + integral2[n]
    }
  }

  out_calculateD <- calculateD(potential.CIFs, x_a, x_l, estimand, prob.bound)
  hesse_d11 <- crossprod(x_la, w11 * out_calculateD$d_11) / n
  hesse_d12 <- crossprod(x_la, w12 * out_calculateD$d_12) / n
  hesse_d22 <- crossprod(x_la, w22 * out_calculateD$d_22) / n
  #hesse_d11 <- t(x_la) %*% (w11 * out_calculateD$d_11) / n
  #hesse_d12 <- t(x_la) %*% (w12 * out_calculateD$d_12) / n
  #hesse_d22 <- t(x_la) %*% (w22 * out_calculateD$d_22) / n

  hesse_d1 <- cbind(hesse_d11, hesse_d12)
  hesse_d2 <- cbind(hesse_d12, hesse_d22)
  hesse <- rbind(hesse_d1, hesse_d2)

  total_score <- cbind(AB1, AB2)
  influence.function <- t(solve(hesse, t(total_score)))
  cov_estimated <- crossprod(influence.function) / n^2
  #influence.function <- total_score %*% t(solve(hesse))
  #cov_estimated <- t(influence.function) %*% influence.function / n / n
  if (ncol(influence.function)==4) {
    colnames(influence.function) <- c("intercept", "exposure", "intercept", "exposure")
  } else if (ncol(influence.function)>4) {
    colnames(influence.function) <- c("intercept", paste("covariate", 1:(ncol(influence.function)/2 - 2), sep = ""), "exposure",
                                      "intercept", paste("covariate", 1:(ncol(influence.function)/2 - 2), sep = ""), "exposure")
  }
  return(list(cov_estimated = cov_estimated, score.function = total_score, influence.function = influence.function))
}

calculateD <- function(potential.CIFs, x_a, x_l, estimand, prob.bound) {
  CIF1 <- ifelse(potential.CIFs[, 1] == 0, prob.bound, ifelse(potential.CIFs[, 1] == 1, 1 - prob.bound, potential.CIFs[, 1]))
  CIF2 <- ifelse(potential.CIFs[, 2] == 0, prob.bound, ifelse(potential.CIFs[, 2] == 1, 1 - prob.bound, potential.CIFs[, 2]))
  CIF3 <- ifelse(potential.CIFs[, 3] == 0, prob.bound, ifelse(potential.CIFs[, 3] == 1, 1 - prob.bound, potential.CIFs[, 3]))
  CIF4 <- ifelse(potential.CIFs[, 4] == 0, prob.bound, ifelse(potential.CIFs[, 4] == 1, 1 - prob.bound, potential.CIFs[, 4]))

  calculateA <- function(effect_measure, exposed.CIFs, unexposed.CIFs, a) {
    if (effect_measure == 'RR') {
      return(a * (1 / exposed.CIFs) + (1 - a) * (1 / unexposed.CIFs))
    } else if (effect_measure == 'OR') {
      return(a * (1 / exposed.CIFs + 1 / (1 - exposed.CIFs)) + (1 - a) * (1 / unexposed.CIFs + 1 / (1 - unexposed.CIFs)))
    } else if (effect_measure == 'SHR') {
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
  n <- length(x_a)
  a11 <- a12 <- a22 <- NULL
  a11 <- calculateA(estimand$effect.measure1, CIF3, CIF1, a)
  a22 <- calculateA(estimand$effect.measure2, CIF4, CIF2, a)
  a12 <- matrix(0, nrow = n, ncol = 1)

  d_ey_d_beta_11 <- a22 / (a11 * a22 - a12 * a12)
  d_ey_d_beta_12 <- -a12 / (a11 * a22 - a12 * a12)
  d_ey_d_beta_22 <- a11 / (a11 * a22 - a12 * a12)
  c12 <- a * (1 / (1 - CIF3 - CIF4)) + (1 - a) * (1 / (1 - CIF1 - CIF2))
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

estimating_equation_survival <- function(
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

  if (!inherits(Y, c("Event", "Surv"))) {
    t <- rep(0, length(Y))
    epsilon <- Y
  } else {
    t <- Y[, 1]
    epsilon <- Y[, 2]
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

  y_0 <- ifelse(epsilon == estimand$code.censoring | t > estimand$time.point, 1, 0)
  y_1 <- ifelse(epsilon == estimand$code.event1 & t <= estimand$time.point, 1, 0)
  y_0_ <- ifelse(epsilon == estimand$code.censoring, 1, 0)
  y_1_ <- ifelse(epsilon == estimand$code.event1, 1, 0)

  a_ <- as.factor(data[[exposure]])
  if (estimand$code.exposure.ref==0) {
    x_a <- as.matrix(model.matrix(~ a_)[, -1])
  } else {
    x_a <- as.matrix(rep(1,length(t)) - model.matrix(~ a_)[, -1])
  }
  x_l <- model.matrix(Terms, mf)
  x_la <- cbind(x_l, x_a)
  i_parameter <- rep(NA, 7)
  i_parameter <- calculateIndexForParameter(NA,x_l,x_a)

  potential.CIFs <- calculatePotentialRisk(alpha_beta, x_a, x_l, offset, estimand)
  one <- rep(1, nrow(x_l))
  a <- as.vector(x_a)
  ey_1 <- potential.CIFs[,2]*a + potential.CIFs[,1]*(one - a)
  w11 <- 1 / (ey_1 * (1 - ey_1))
  wy_1 <- ip.weight * y_1
  wy_1ey_1 <- w11*(wy_1 - ey_1)
  d <- cbind(x_l, x_a)
  residual <- wy_1ey_1
  ret <- drop(crossprod(d, residual)) / nrow(x_l)
  #ret <- as.vector(t(d) %*% residual / nrow(x_l))

  n_col_d <- ncol(d)
  score.matrix <- matrix(NA, nrow=nrow(d), ncol=n_col_d)
  for (j in seq_len(n_col_d)) {
    score.matrix[,j] <- d[,j]*residual
  }
  out <- list(
    ret   = ret,
    score = score.matrix,
    ey_1  = ey_1,
    w11   = w11,
    t     = t,
    y_0   = y_0,
    y_1   = y_1,
    y_0_  = y_0_,
    y_1_  = y_1_,
    x_a   = x_a,
    x_l   = x_l,
    potential.CIFs = potential.CIFs,
    ip.weight = ip.weight,
    i_parameter = i_parameter
  )
  return(out)
}

calculateCovSurvival <- function(objget_results, estimand, prob.bound)
{
  score <- objget_results$score
  ey_1 <- objget_results$ey_1
  w11 <- objget_results$w11
  t <- objget_results$t
  y_0_ <- objget_results$y_0_
  y_1 <- objget_results$y_1
  x_a <- objget_results$x_a
  x_l <- objget_results$x_l
  potential.CIFs <- objget_results$potential.CIFs
  i_parameter <- objget_results$i_parameter
  n <- length(t)

  censoring_dna	<- calculateNelsonAalen(t,y_0_)
  censoring_martingale <- diag(y_0_) - (outer(t, t, ">") * censoring_dna)
  censoring_km <- calculateKaplanMeier(t, y_0_)
  censoring_km[censoring_km == 0] <- 1e-5
  censoring_mkm <- censoring_martingale / censoring_km
  y_12 <- (y_1 > 0)
  survival_km <- calculateKaplanMeier(t, y_12)
  wy_1 <- w11 * (y_1 - ey_1)
  x_la <- cbind(x_l, x_a)
  AB1 <- score[1:n, 1:i_parameter[2]]
  for (i_para in 1:i_parameter[2]) {
    tmp0 <- x_la[, i_para]
    use <- (t <= estimand$time.point)
    tmp1 <- colSums((use * tmp0) * (outer(t, t, ">=") * wy_1))
    tmp1 <- tmp1 / survival_km / n
    integrand1 <- tmp1 * censoring_mkm
    for (i_score in 1:n) {
      integral1 <- cumsum(integrand1[i_score, ])
      AB1[i_score, i_para] <- AB1[i_score, i_para] + integral1[n]
    }
  }
  out_calculateDSurvival <- calculateDSurvival(potential.CIFs, x_a, x_l, estimand, prob.bound)
#  hesse <- t(x_la) %*% (w11 * out_calculateDSurvival) / n
  hesse <- crossprod(x_la, w11 * out_calculateDSurvival) / n

  total_score <- AB1
  #influence.function <- total_score %*% t(solve(hesse))
  #cov_estimated <- t(influence.function) %*% influence.function / n / n
  influence.function <- t(solve(hesse, t(total_score)))
  cov_estimated <- crossprod(influence.function) / n^2
  if (ncol(influence.function)==2) {
    colnames(influence.function) <- c("intercept", "exposure")
  } else if (ncol(influence.function)>2) {
    colnames(influence.function) <- c("intercept", paste("covariate", 1:(ncol(influence.function) - 2), sep = ""), "exposure")
  }
  return(list(cov_estimated = cov_estimated, score.function = total_score, influence.function = influence.function))
}

calculateDSurvival <- function(potential.CIFs, x_a, x_l, estimand, prob.bound) {
  CIF1 <- ifelse(potential.CIFs[, 1] == 0, prob.bound, ifelse(potential.CIFs[, 1] == 1, 1 - prob.bound, potential.CIFs[, 1]))
  CIF2 <- ifelse(potential.CIFs[, 2] == 0, prob.bound, ifelse(potential.CIFs[, 2] == 1, 1 - prob.bound, potential.CIFs[, 2]))

  calculateA <- function(effect_measure, exposed.CIFs, unexposed.CIFs, a) {
    if (effect_measure == 'RR') {
      return(a * (1 / exposed.CIFs) + (1 - a) * (1 / unexposed.CIFs))
    } else if (effect_measure == 'OR') {
      return(a * (1 / exposed.CIFs + 1 / (1 - exposed.CIFs)) + (1 - a) * (1 / unexposed.CIFs + 1 / (1 - unexposed.CIFs)))
    } else if (effect_measure == 'SHR') {
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
  n <- length(x_a)
  a11 <-NULL
  a11 <- calculateA(estimand$effect.measure1, CIF2, CIF1, a)

  d_ey_d_beta_11 <- 1 / a11
  d_ey_d_alpha_11 <- 1 / a11
  d_11 <- cbind((d_ey_d_alpha_11 * x_l), (d_ey_d_beta_11 * x_a))
  return(d_11)
}

estimating_equation_proportional <- function(
    formula,
    data,
    exposure,
    ip.weight.matrix,
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
  if (!inherits(Y, c("Event", "Surv"))) {
    t <- rep(0, length(Y))
    epsilon <- Y
  } else {
    t <- Y[, 1]
    epsilon <- Y[, 2]
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
  if (is.null(estimand$time.point)) {
    time.point <- with(data, t[epsilon > 0])
    time.point <- unique(time.point)
  } else {
    time.point <- estimand$time.point
  }

  y_0_ <- ifelse(epsilon == estimand$code.censoring, 1, 0)
  y_1_ <- ifelse(epsilon == estimand$code.event1, 1, 0)

  a_ <- as.factor(data[[exposure]])
  if (estimand$code.exposure.ref==0) {
    x_a <- as.matrix(model.matrix(~ a_)[, -1])
  } else {
    x_a <- as.matrix(rep(1,length(t)) - model.matrix(~ a_)[, -1])
  }
  x_l <- model.matrix(Terms, mf)
  x_la <- cbind(x_l, x_a)
  i_parameter <- rep(NA, 7)
  i_parameter <- calculateIndexForParameter(NA,x_l,x_a,length(time.point))

  one <- rep(1, nrow(x_l))
  a <- as.vector(x_a)
  score_beta <- NULL
  score_alpha1 <- 0
  score_alpha2 <- 0
  i_time <- 0
  alpha_beta_i <- rep(NA, i_parameter[7])
  for (specific.time in time.point) {
    i_time <- i_time + 1
    i_para <- i_parameter[1]*(i_time-1)+1
    alpha_beta_i[1:i_parameter[1]]              <- alpha_beta[i_para:(i_para+i_parameter[1]-1)]
    alpha_beta_i[i_parameter[2]:i_parameter[3]] <- alpha_beta[i_parameter[8]/2]
    y_0 <- ifelse(epsilon == estimand$code.censoring | t > specific.time, 1, 0)
    y_1 <- ifelse(epsilon == estimand$code.event1 & t <= specific.time, 1, 0)

    if (optim.method$computation.order.method=="OLD") {
      potential.CIFs <- calculatePotentialCIFs_old(alpha_beta,x_a,x_l,offset,epsilon,estimand,optim.method,prob.bound,initial.CIFs)
    } else {
      potential.CIFs <- calculatePotentialCIFs_parallel(alpha_beta,x_a,x_l,offset,epsilon,estimand,optim.method,prob.bound,initial.CIFs)
    }
    one <- rep(1, nrow(x_l))
    a <- as.vector(x_a)
    ey_1 <- potential.CIFs[,2]*a + potential.CIFs[,1]*(one - a)
    w11 <- 1 / (ey_1 * (1 - ey_1))
#    wy_1 <- ip.weight * y_1
    wy_1 <- ip.weight.matrix[,i_time] * y_1
    wy_1ey_1 <- w11*(wy_1 - ey_1)
    d <- cbind(x_l, x_a)
    residual <- wy_1ey_1
#    ret <- as.vector(t(d) %*% residual / nrow(x_l))
#    n_col_d <- ncol(d)
#    score.matrix <- matrix(NA, nrow=nrow(d), ncol=n_col_d)
#    for (j in seq_len(n_col_d)) {
#      score.matrix[,j] <- d[,j]*residual
#    }

    subscore <- as.vector(t(d) %*% residual / nrow(x_l))
    tmp1 <- t(subscore[1:i_parameter[1]])
    score_beta <- cbind(score_beta, tmp1)
    score_alpha1 <- score_alpha1 + subscore[i_parameter[2]:i_parameter[3]]
#    tmp1 <- t(subscore[1:i_parameter[1],])
#    tmp2 <- t(subscore[i_parameter[4]:i_parameter[5],])
#    score_beta <- cbind(score_beta, tmp1, tmp2)
#    score_alpha1 <- score_alpha1 + subscore[i_parameter[2]:i_parameter[3],]
#    score_alpha2 <- score_alpha2 + subscore[i_parameter[6]:i_parameter[7],]
  }
#  score <- cbind(score_beta, score_alpha1, score_alpha2)
  score <- cbind(score_beta, score_alpha1)
  out <- list(
    ret   = score,
    t     = t,
    y_0   = y_0,
    y_1   = y_1,
    y_0_  = y_0_,
    y_1_  = y_1_,
    x_a   = x_a,
    x_l   = x_l,
    i_parameter = i_parameter
  )
  return(out)
}

estimating_equation_pproportional <- function(
  formula,
  data,
  exposure,
  ip.weight.matrix,
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
  t <- Y[, 1]
  epsilon <- Y[, 2]
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
  if (is.null(estimand$time.point)) {
    time.point <- with(data, t[epsilon > 0])
    time.point <- unique(time.point)
  } else {
    time.point <- estimand$time.point
  }

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
  i_parameter <- calculateIndexForParameter(NA,x_l,x_a,length(time.point))

  #n_para_1 <- ncol(x_l)
  #n_para_2 <- n_para_1 + 1
  #n_para_3 <- n_para_1 + 2
  #n_para_4 <- 2*n_para_1 + 1
  #n_para_5 <- 2*n_para_1 + 2
  #n_para_6 <- length(time.point)*(n_para_5-2) + 2

  one <- rep(1, nrow(x_l))
  a <- as.vector(x_a)
  score_beta <- NULL
  score_alpha1 <- 0
  score_alpha2 <- 0
  i_time <- 0
  alpha_beta_i <- rep(NA, i_parameter[7])
  #alpha_beta_i <- rep(NA, n_para_5)
  for (specific.time in time.point) {
    i_time <- i_time + 1
    i_para <- i_parameter[1]*(i_time-1)+1
#    i_para <- n_para_1*(i_time-1)+1
    alpha_beta_i[1:i_parameter[1]]              <- alpha_beta[i_para:(i_para+i_parameter[1]-1)]
    alpha_beta_i[i_parameter[2]:i_parameter[3]] <- alpha_beta[i_parameter[8]/2]
    alpha_beta_i[i_parameter[4]:i_parameter[5]] <- alpha_beta[(i_parameter[8]/2+i_para):(i_parameter[8]/2+i_para+i_parameter[1]-1)]
    alpha_beta_i[i_parameter[6]:i_parameter[7]] <- alpha_beta[i_parameter[8]]

    #alpha_beta_i[1:n_para_1]        <- alpha_beta[i_para:(i_para+n_para_1-1)]
    #alpha_beta_i[n_para_2]          <- alpha_beta[n_para_6/2]
    #alpha_beta_i[n_para_3:n_para_4] <- alpha_beta[(n_para_6/2+i_para):(n_para_6/2+i_para+n_para_1-1)]
    #alpha_beta_i[n_para_5]          <- alpha_beta[n_para_6]

    y_0 <- ifelse(epsilon == 0 | t > specific.time, 1, 0)
    y_1 <- ifelse(epsilon == 1 & t <= specific.time, 1, 0)
    y_2 <- ifelse(epsilon == 2 & t <= specific.time, 1, 0)

    #  potential.CIFs <- calculatePotentialCIFs(alpha_beta,x_a,x_l,offset,epsilon,estimand,optim.method,prob.bound,initial.CIFs)
    potential.CIFs <- calculatePotentialCIFs_parallel(alpha_beta,x_a,x_l,offset,epsilon,estimand,optim.method,prob.bound,initial.CIFs)
    ey_1 <- potential.CIFs[,3]*a + potential.CIFs[,1]*(one - a)
    ey_2 <- potential.CIFs[,4]*a + potential.CIFs[,2]*(one - a)

    v11 <- ey_1 * (1 - ey_1)
    v12 <- -ey_1 * ey_2
    v22 <- ey_2 * (1 - ey_2)
    denom <- v11*v22 - v12*v12

    w11 <- v22 / denom
    w12 <- -v12 / denom
    w22 <- v11 / denom

    wy_1 <- ip.weight.matrix[,i_time] * y_1
    wy_2 <- ip.weight.matrix[,i_time] * y_2

    wy_1ey_1 <- w11*(wy_1 - ey_1) + w12*(wy_2 - ey_2)
    wy_2ey_2 <- w12*(wy_1 - ey_1) + w22*(wy_2 - ey_2)

    x_la <- cbind(x_l, x_a)
    zero <- matrix(0, nrow=nrow(x_la), ncol=ncol(x_la))
    tmp1 <- cbind(x_la, zero)
    tmp2 <- cbind(zero, x_la)
    d    <- rbind(tmp1, tmp2)

    residual <- c(wy_1ey_1, wy_2ey_2)
    #subscore <- as.matrix(t(d) %*% residual / nrow(x_l))
    subscore <- crossprod(d, residual) / nrow(x_l)
    #tmp1 <- t(subscore[1:n_para_1,])
    #tmp2 <- t(subscore[n_para_3:n_para_4,])
    #score_beta <- cbind(score_beta, tmp1, tmp2)
    #score_alpha1 <- score_alpha1 + subscore[n_para_2,]
    #score_alpha2 <- score_alpha2 + subscore[n_para_5,]

    tmp1 <- t(subscore[1:i_parameter[1],])
    tmp2 <- t(subscore[i_parameter[4]:i_parameter[5],])
    score_beta <- cbind(score_beta, tmp1, tmp2)
    score_alpha1 <- score_alpha1 + subscore[i_parameter[2]:i_parameter[3],]
    score_alpha2 <- score_alpha2 + subscore[i_parameter[6]:i_parameter[7],]
  }
  score <- cbind(score_beta, score_alpha1, score_alpha2)
  out <- list(
    ret   = score,
    t     = t,
    y_0   = y_0,
    y_1   = y_1,
    y_2   = y_2,
    y_0_  = y_0_,
    y_1_  = y_1_,
    y_2_  = y_2_,
    x_a   = x_a,
    x_l   = x_l,
    i_parameter = i_parameter
  )
#  out <- list(ret = score)
  return(out)
}

createAtRiskMatrix <- function(t) {
  atrisk <- outer(t, t, "<=")
  return(atrisk)
}

calculateKaplanMeier <- function(t, d){
  n <- length(t)
  data <- data.frame(t = t, d = d, id = 1:n)
  sorted_data <- data[order(data$t), ]
  sorted_t <- sorted_data$t
  sorted_d <- sorted_data$d
  sorted_id <- sorted_data$id
  atrisk <- createAtRiskMatrix(sorted_t)
  n_atrisk <- rowSums(atrisk)
  s <- 1 - sorted_d / n_atrisk
  km <- cumprod(s)
  data <- data.frame(id=sorted_id, km=km)
  sorted_data <- data[order(data$id), ]
  km <- sorted_data$km
  return(km)
}

calculateNelsonAalen <- function(t, d) {
  atrisk <- createAtRiskMatrix(t)
  n_atrisk <- rowSums(atrisk)
  na <- d / n_atrisk
  return(na)
}

## ===== 0) ユーティリティ（必要なら流用） ==========================
`%||%` <- function(x, y) if (is.null(x)) y else x

clampP <- function(p, bound) {
  # 確率ベクトル p を [bound, 1-bound] にクランプしつつ、(1-p1-p2) も安全化
  p <- pmax(bound, pmin(1 - bound, p))
  p
}
clampLogP <- function(lp, bound = 0) {
  # log確率に軽い下駄（必要なら bound を利用）
  # 既存の clampLogP があるならそちらを使ってOK
  lp
}

## ===== 1) effect.measure の小カーネルを作るファクトリー ===========
## idx = (i,j) は (1,3) or (2,4) のペアを想定
make_effect_kernel <- function(measure, idx) {
  i <- idx[1]; j <- idx[2]
  if (measure == "RR") {
    res <- function(p, clogp, beta) {
      beta - clogp[j] + clogp[i]
    }
    jac <- function(p, clogp, beta) {
      g <- numeric(4); g[i] <-  1; g[j] <- -1; g
    }
  } else if (measure == "OR") {
    res <- function(p, clogp, beta) {
      beta - clogp[j] + clogp[i] + log1p(-p[j]) - log1p(-p[i])
    }
    jac <- function(p, clogp, beta) {
      g <- numeric(4)
      g[i] <-  1 +  p[i]/(1 - p[i])
      g[j] <- -1 + (-p[j]/(1 - p[j]))
      g
    }
  } else if (measure == "SHR" || identical(measure, "")) {
    res <- function(p, clogp, beta) {
      exp(beta) - (log1p(-p[j]) / log1p(-p[i]))
    }
    jac <- function(p, clogp, beta) {
      # ∂/∂log p ： d p / d log p = p
      g <- numeric(4)
      Bi <- log1p(-p[i])  # < 0
      Bj <- log1p(-p[j])
      # d/d log p_i of [-log(1-p_j)/log(1-p_i)] =  (Bj * (-p_i/(1-p_i))) / Bi^2
      g[i] <-  ( Bj * (-p[i]/(1 - p[i])) ) / (Bi^2)
      # d/d log p_j of [-log(1-p_j)/log(1-p_i)] =  p_j / ((1-p_j) * Bi)
      g[j] <-   p[j] / ((1 - p[j]) * Bi)
      g
    }
  } else {
    stop("Invalid measure: must be 'RR','OR','SHR' (or '' as SHR).")
  }
  list(res = res, jac = jac, idx = c(i, j))
}

## ===== 2) 残差・ヤコビアン（本体は分岐なし、カーネル差し替え） ====
residuals_CIFs_generic <- function(log_p, alpha1, beta1, alpha2, beta2,
                                   estimand, prob.bound, k1, k2) {
  clogp <- clampLogP(as.numeric(log_p))
  if (length(clogp) < 4L) clogp <- rep(clogp, 4L)
  p <- exp(clogp)

  rem12 <- max(1 - p[1] - p[2], prob.bound)
  rem34 <- max(1 - p[3] - p[4], prob.bound)
  lp0102 <- log(rem12) + log(rem34)

  r <- numeric(4)
  # 共通部（ret1, ret3）
  r[1] <- alpha1 - clogp[1] - clogp[3] + lp0102
  r[3] <- alpha2 - clogp[2] - clogp[4] + lp0102
  # 差分部（ret2, ret4）はカーネルで
  r[2] <- k1$res(p, clogp, beta1)
  r[4] <- k2$res(p, clogp, beta2)
  r
}

jacobian_CIFs_generic <- function(log_p, alpha1, beta1, alpha2, beta2,
                                  estimand, prob.bound, k1, k2) {
  clogp <- clampLogP(as.numeric(log_p))
  if (length(clogp) < 4L) clogp <- rep(clogp, 4L)
  p <- exp(clogp)

  rem12 <- max(1 - p[1] - p[2], prob.bound)
  rem34 <- max(1 - p[3] - p[4], prob.bound)

  dlp12 <- c(-p[1]/rem12, -p[2]/rem12)  # ∂/∂logp log(1-p1-p2)
  dlp34 <- c(-p[3]/rem34, -p[4]/rem34)

  J <- matrix(0.0, 4, 4)
  # 共通部（ret1, ret3 のヤコビアン）
  J[1,] <- c(-1 + dlp12[1], dlp12[2], -1 + dlp34[1], dlp34[2])
  J[3,] <- c(dlp12[1], -1 + dlp12[2], dlp34[1], -1 + dlp34[2])
  # 差分部（ret2, ret4）
  J[2,] <- k1$jac(p, clogp, beta1)
  J[4,] <- k2$jac(p, clogp, beta2)
  J
}

## ===== 3) 4×4専用・極小LM（base Rのみ） =============================
## ── 4×4 専用 Levenberg–Marquardt（堅牢版）──
## * nls.lm 相当の挙動（ρで λ を更新）
## * line search 付き
## * base R のみ（chol/backsolve/forwardsolve/crossprod）
tiny_lm4_robust <- function(start,
                            res_fun, jac_fun,
                            maxit = 100L,
                            ftol  = 1e-10,
                            ptol  = 1e-10,
                            lambda0    = 1e-3,
                            lambda_min = 1e-12,
                            lambda_max = 1e12,
                            eta_inc    = 2.0,
                            eta_dec    = 0.5,
                            do_linesearch = TRUE,
                            ls_c      = 1e-4,
                            ls_shrink = 0.5,
                            verbose   = FALSE) {
  ## 初期化
  lp   <- as.numeric(start)
  r    <- res_fun(lp)
  f2   <- drop(crossprod(r))            # ||r||^2
  J    <- jac_fun(lp)
  g    <- drop(crossprod(J, r))         # J^T r
  lam  <- lambda0
  prev_f2 <- Inf                        # 直前の目的関数値

  if (isTRUE(verbose)) {
    hist <- list(it=integer(), f2=double(),
                 grad_inf=double(), lambda=double(),
                 step_norm=double(), rho=double())
  }

  for (k in seq_len(maxit)) {
    ## 停止1: 勾配ノルム
    grad_inf <- max(abs(g))
    if (is.finite(grad_inf) && grad_inf < 1e-6) break

    ## (J^T J + λI) δ = -J^T r
    A <- crossprod(J)
    diag(A) <- diag(A) + lam
    R <- try(chol(A), silent = TRUE)
    if (inherits(R, "try-error")) {
      lam <- min(lam * 10, lambda_max)
      next
    }
    delta <- -backsolve(R, forwardsolve(t(R), g))
    step_norm <- max(abs(delta))

    ## 停止2: ステップが小さい
    if (is.finite(step_norm) &&
        step_norm <= ptol * (ptol + max(1, max(abs(lp))))) break

    ## 予測減少（LM モデル）
    pred <- 0.5 * sum(delta * (lam * delta - g))
    if (!is.finite(pred) || pred <= 0) pred <- .Machine$double.eps

    ## 候補点
    lp_try <- lp + delta
    r_try  <- res_fun(lp_try)
    f2_try <- drop(crossprod(r_try))

    ## ρ = 実際の減少 / 予測減少
    rho <- (f2 - f2_try) / pred
    accepted <- FALSE

    if (is.finite(rho) && rho > 0) {
      ## 受容 → λを減らす
      lp <- lp_try; r <- r_try; f2 <- f2_try
      lam <- max(lambda_min, lam * max(eta_dec, 1/(1 + rho)))
      accepted <- TRUE
    } else {
      ## 不受容 → λを増やす
      lam <- min(lambda_max, lam * (eta_inc * (1 + abs(ifelse(is.finite(rho), rho, 0)))))
      if (do_linesearch) {
        ## Armijo 型ラインサーチ
        t <- 1.0; gTd <- sum(g * delta)
        while (t > 1e-6) {
          lp_ls <- lp + t * delta
          r_ls  <- res_fun(lp_ls)
          f2_ls <- drop(crossprod(r_ls))
          if (is.finite(f2_ls) &&
              f2_ls <= f2 + ls_c * t * gTd) {
            lp <- lp_ls; r <- r_ls; f2 <- f2_ls
            accepted <- TRUE
            break
          }
          t <- t * ls_shrink
        }
      }
    }

    ## 次イテレーション用に更新
    J <- jac_fun(lp)
    g <- drop(crossprod(J, r))

    ## 停止3: 目的関数の変化が小さい
    if (k > 1L &&
        is.finite(prev_f2) && is.finite(f2) &&
        abs(prev_f2 - f2) <= ftol * (abs(f2) + ftol)) break

    ## ログ
    if (isTRUE(verbose)) {
      hist$it        <- c(hist$it, k)
      hist$f2        <- c(hist$f2, f2)
      hist$grad_inf  <- c(hist$grad_inf, grad_inf)
      hist$lambda    <- c(hist$lambda, lam)
      hist$step_norm <- c(hist$step_norm, step_norm)
      hist$rho       <- c(hist$rho, if (exists("rho")) rho else NA_real_)
    }

    prev_f2 <- f2
  }

  if (isTRUE(verbose)) attr(lp, "history") <- hist
  lp
}

## ===== 4) 観測1件用の小ラッパ ================================
# NULL 合体演算子（未定義ならこれも定義しておく）
`%||%` <- function(x, y) if (is.null(x)) y else x

# 方法A：明示指定で tiny_lm4_robust を呼ぶ安全版ラッパ
solve_one_CIFs_tinyLM <- function(start_lp,
                                  alpha1, beta1, alpha2, beta2,
                                  estimand, prob.bound,
                                  k1, k2,
                                  ctrl = list()) {
  # 残差・ヤコビアン（分岐なし、本体は共通）
  res_fun <- function(lp) residuals_CIFs_generic(lp, alpha1, beta1, alpha2, beta2,
                                                 estimand, prob.bound, k1, k2)
  jac_fun <- function(lp) jacobian_CIFs_generic(lp, alpha1, beta1, alpha2, beta2,
                                                estimand, prob.bound, k1, k2)

  # 必要なパラメータだけ明示的に拾って渡す（他は無視）
  lp_hat <- tiny_lm4_robust(
    start         = start_lp,
    res_fun       = res_fun,
    jac_fun       = jac_fun,
    maxit         = ctrl$maxit        %||% 100L,
    ftol          = ctrl$ftol         %||% 1e-10,
    ptol          = ctrl$ptol         %||% 1e-10,
    lambda0       = ctrl$lambda0      %||% 1e-3,
    lambda_min    = ctrl$lambda_min   %||% 1e-12,
    lambda_max    = ctrl$lambda_max   %||% 1e12,
    eta_inc       = ctrl$eta_inc      %||% 2.0,     # ρ<0 等でλを増やす倍率
    eta_dec       = ctrl$eta_dec      %||% 0.5,     # ρ>0 でλを減らす倍率
    do_linesearch = ctrl$do_linesearch%||% TRUE,
    ls_c          = ctrl$ls_c         %||% 1e-4,
    ls_shrink     = ctrl$ls_shrink    %||% 0.5,
    verbose       = ctrl$verbose      %||% FALSE
  )

  lp_hat
}

## ===== 5) drop-in 置換: calculatePotentialCIFs_* =====================
## 旧 calculatePotentialCIFs_old と同じ役割＆返り値
calculatePotentialCIFs_tinyLM <- function(
    alpha_beta_tmp,
    x_a, x_l, offset, epsilon,
    estimand, optim.method, prob.bound,
    initial.CIFs = NULL,
    lm_ctrl = list(maxit=30L, ftol=1e-10, ptol=1e-10,
                   lambda0=1e-3, lambda_up=10, lambda_down=0.1)
) {
  # 0) effect.measure を一度だけ固定（以降分岐なし）
  k1 <- make_effect_kernel(estimand$effect.measure1, c(1,3))
  k2 <- make_effect_kernel(estimand$effect.measure2 %||% "SHR", c(2,4))

  # 1) パラメータ分割
  i_parameter <- rep(NA_integer_, 7L)
  i_parameter <- calculateIndexForParameter(i_parameter, x_l, x_a)
  alpha_1    <- alpha_beta_tmp[seq_len(i_parameter[1])]
  beta_tmp_1 <- alpha_beta_tmp[seq.int(i_parameter[2], i_parameter[3])]
  alpha_2    <- alpha_beta_tmp[seq.int(i_parameter[4], i_parameter[5])]
  beta_tmp_2 <- alpha_beta_tmp[seq.int(i_parameter[6], i_parameter[7])]

  # 2) 線形予測子
  alpha_tmp_1 <- as.numeric(x_l %*% matrix(alpha_1, ncol = 1) + offset)
  alpha_tmp_2 <- as.numeric(x_l %*% matrix(alpha_2, ncol = 1) + offset)

  # 3) 初期確率（安全化）
  n  <- length(epsilon)
  p0 <- c(
    sum(epsilon == estimand$code.event1) / n + prob.bound,
    sum(epsilon == estimand$code.event2) / n + prob.bound,
    sum(epsilon == estimand$code.event1) / n + prob.bound,
    sum(epsilon == estimand$code.event2) / n + prob.bound
  )
  p0     <- clampP(p0, prob.bound)
  log_p0 <- log(p0)

  # 4) 全行ユニーク化キャッシュ（baseのみ）
  keys <- apply(x_l, 1, function(r) paste0(r, collapse = "\r"))
  uniq <- match(keys, unique(keys))
  cache_lp <- vector("list", length = max(uniq))

  out <- matrix(NA_real_, nrow = nrow(x_l), ncol = 4L)

  for (i in seq_len(nrow(x_l))) {
    k <- uniq[i]

    if (!is.null(cache_lp[[k]])) {
      lp_hat <- cache_lp[[k]]
    } else {
      ip <- if (!is.null(initial.CIFs)) as.numeric(initial.CIFs[i, 1:4, drop = FALSE]) else exp(log_p0)
      start_lp <- log(clampP(ip, prob.bound))

      lp_hat <- solve_one_CIFs_tinyLM(
        start_lp  = start_lp,
        alpha1    = alpha_tmp_1[i],
        beta1     = beta_tmp_1,
        alpha2    = alpha_tmp_2[i],
        beta2     = beta_tmp_2,
        estimand  = estimand,
        prob.bound = prob.bound,
        k1 = k1, k2 = k2,
        ctrl = lm_ctrl
      )
      cache_lp[[k]] <- lp_hat
    }

    out[i, ] <- clampP(exp(lp_hat), prob.bound)
  }

  colnames(out) <- c("p10", "p20", "p11", "p21")
  out
}
