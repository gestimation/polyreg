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
  i_parameter <- calculateIndexForParameter(i_parameter,x_l,x_a)

  potential.CIFs <- calculatePotentialCIFs(alpha_beta,x_a,x_l,offset,epsilon,estimand,optim.method,prob.bound,initial.CIFs)
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
  hesse_d11 <- t(x_la) %*% (w11 * out_calculateD$d_11) / n
  hesse_d12 <- t(x_la) %*% (w12 * out_calculateD$d_12) / n
  hesse_d22 <- t(x_la) %*% (w22 * out_calculateD$d_22) / n

  hesse_d1 <- cbind(hesse_d11, hesse_d12)
  hesse_d2 <- cbind(hesse_d12, hesse_d22)
  hesse <- rbind(hesse_d1, hesse_d2)

  total_score <- cbind(AB1, AB2)
  influence.function <- total_score %*% t(solve(hesse))
  cov_estimated <- t(influence.function) %*% influence.function / n / n
  if (ncol(influence.function)==4) {
    colnames(influence.function) <- c("intercept", "exposure", "intercept", "exposure")
  } else {
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
  i_parameter <- calculateIndexForParameter(i_parameter,x_l,x_a)

  potential.CIFs <- calculatePotentialRisk(alpha_beta, x_l, offset, estimand)
  one <- rep(1, nrow(x_l))
  a <- as.vector(x_a)
  ey_1 <- potential.CIFs[,2]*a + potential.CIFs[,1]*(one - a)
  w11 <- 1 / (ey_1 * (1 - ey_1))
  wy_1 <- ip.weight * y_1
  wy_1ey_1 <- w11*(wy_1 - ey_1)
  d <- cbind(x_l, x_a)
  residual <- wy_1ey_1
  ret <- as.vector(t(d) %*% residual / nrow(x_l))

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
  hesse <- t(x_la) %*% (w11 * out_calculateDSurvival) / n
  total_score <- AB1
  influence.function <- total_score %*% t(solve(hesse))
  cov_estimated <- t(influence.function) %*% influence.function / n / n
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
  time.point,
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
  i_parameter <- calculateIndexForParameter(i_parameter,x_l,x_a,length(time.point))

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

    potential.CIFs <- calculatePotentialCIFs(alpha_beta_i,x_a,x_l,offset,epsilon,estimand,optim.method,prob.bound,initial.CIFs)
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
    subscore <- as.matrix(t(d) %*% residual / nrow(x_l))
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

calculateIndexForParameter <- function(i_parameter,x_l,x_a,length.time.point=1) {
  i_parameter[1] <- ncol(x_l)
  i_parameter[2] <- i_parameter[1] + 1
  i_parameter[3] <- i_parameter[1] + ncol(x_a)
  i_parameter[4] <- i_parameter[1] + ncol(x_a) + 1
  i_parameter[5] <- 2 * i_parameter[1] + ncol(x_a)
  i_parameter[6] <- 2 * i_parameter[1] + ncol(x_a) + 1
  i_parameter[7] <- 2 * i_parameter[1] + 2 * ncol(x_a)
  i_parameter[8] <- length.time.point*(2 * i_parameter[1]) + 2 * ncol(x_a)
  return(i_parameter)
}

