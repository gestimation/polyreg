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
  out_terms <- terms(formula, special, data = data)
  mf$formula <- out_terms
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Y <- model.extract(mf, "response")
  t <- Y[, 1]
  epsilon <- Y[, 2]
  if (!is.null(offsetpos <- attributes(out_terms)$specials$offset)) {
    ts <- survival::untangle.specials(out_terms, "offset")
    if (length(ts$vars) > 0) {
      out_terms <- out_terms[-ts$terms]
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

  out_defineExposureDesign <- defineExposureDesign(data, exposure, estimand$code.exposure.ref)
  x_a <- out_defineExposureDesign$x_a
  x_l <- model.matrix(out_terms, mf)
  x_la <- cbind(x_l, x_a)

  potential.CIFs <- calculatePotentialCIFs(alpha_beta,x_a,x_l,offset,epsilon,estimand,optim.method,prob.bound,initial.CIFs)
  ey <- calculateEY(potential.CIFs, x_a)
  ey_1 <- ey$ey_1
  ey_2 <- ey$ey_2

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

  n_col_d <- ncol(d)
  score.matrix <- matrix(NA, nrow=nrow(d), ncol=n_col_d)
  for (j in seq_len(n_col_d)) {
    score.matrix[,j] <- d[,j]*residual
  }
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
    ip.weight = ip.weight
  )
  return(out)
}

calculateEY <- function(potential.CIFs, x_a) {
  n <- nrow(x_a)
  K <- ncol(x_a) + 1
  if (nrow(potential.CIFs) != n) stop("Row dimension of potential.CIFs is more than rows of x_a")
  if (ncol(potential.CIFs) != 2*K) stop("Column dimension of potential.CIFs is not 2*K")
  index <- max.col(cbind(0, x_a), ties.method = "first")
  CIF1 <- potential.CIFs[, seq_len(K), drop = FALSE]
  CIF2 <- potential.CIFs[, K + seq_len(K), drop = FALSE]
  ey_1 <- CIF1[cbind(seq_len(n), index)]
  ey_2 <- CIF2[cbind(seq_len(n), index)]
  return(list(ey_1 = ey_1, ey_2 = ey_2))
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
  index.vector <- estimand$index.vector
  n <- length(t)

  censoring_dna	<- calculateNelsonAalen(t,y_0_)
  censoring_martingale <- diag(y_0_) - (outer(t, t, ">") * censoring_dna)
  censoring_km <- calculateKaplanMeier(t, y_0_)
  censoring_km[censoring_km == 0] <- prob.bound
  censoring_mkm <- censoring_martingale / censoring_km
  y_12 <- (y_1 + y_2 > 0)
  survival_km <- calculateKaplanMeier(t, y_12)
  wy_1 <- w11 * (y_1 - ey_1) + w12 * (y_2 - ey_2)
  wy_2 <- w12 * (y_1 - ey_1) + w22 * (y_2 - ey_2)
  x_la <- cbind(x_l, x_a)
  #AB1 <- score[1:n, 1:index.vector[2]]
  #AB2 <- score[(n + 1):(2 * n), index.vector[3]:index.vector[5]]
  AB1 <- score[1:n, 1:index.vector[3]]
  AB2 <- score[(n + 1):(2 * n), index.vector[4]:index.vector[7]]
  for (i_para in 1:index.vector[2]) {
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

  hesse_d1 <- cbind(hesse_d11, hesse_d12)
  hesse_d2 <- cbind(hesse_d12, hesse_d22)
  hesse <- rbind(hesse_d1, hesse_d2)

  total_score <- cbind(AB1, AB2)
  influence.function <- t(solve(hesse, t(total_score)))
  cov_estimated <- crossprod(influence.function) / n^2
  if (ncol(influence.function)==4) {
    colnames(influence.function) <- c("intercept", "exposure", "intercept", "exposure")
  } else if (ncol(influence.function)>4) {
    colnames(influence.function) <- c("intercept", paste("covariate", 1:(ncol(influence.function)/2 - 2), sep = ""), "exposure",
                                      "intercept", paste("covariate", 1:(ncol(influence.function)/2 - 2), sep = ""), "exposure")
  }
  return(list(cov_estimated = cov_estimated, score.function = total_score, influence.function = influence.function))
}

calculateD <- function(potential.CIFs, x_a, x_l, estimand, prob.bound) {
  n <- nrow(x_l)
  K <- estimand$exposure.levels
  if (ncol(potential.CIFs) != 2 * K) stop("Column dimension of potential.CIFs is not 2 * K.")

  seq_n <- seq_len(n)
  idx <- max.col(cbind(0, x_a), ties.method = "first")

  CIF1_mat <- clampP(potential.CIFs[, seq_len(K), drop = FALSE], prob.bound)
  CIF2_mat <- clampP(potential.CIFs[, K + seq_len(K), drop = FALSE], prob.bound)
  survival_mat <- pmax(1 - CIF1_mat - CIF2_mat, prob.bound)

  CIF1_sel <- CIF1_mat[cbind(seq_n, idx)]
  CIF2_sel <- CIF2_mat[cbind(seq_n, idx)]
  survival_sel <- survival_mat[cbind(seq_n, idx)]

  calculateA <- function(effect.measure, CIFs_selected) {
    CIFs_selected <- clampP(CIFs_selected, prob.bound)
    if (effect.measure == "RR") {
      return(1 / CIFs_selected)
    } else if (effect.measure == "OR") {
      denom <- clampP(1 - CIFs_selected, prob.bound)
      return(1 / CIFs_selected + 1 / denom)
    } else if (effect.measure == "SHR") {
      denom <- clampP(1 - CIFs_selected, prob.bound)
      tmp1 <- -1 / denom
      tmp2 <- log(denom)
      return(tmp1 / tmp2)
    } else {
      stop("Invalid effect_measure. Must be RR, OR or SHR.")
    }
  }

  a11 <- calculateA(estimand$effect.measure1, CIF1_sel)
  a22 <- calculateA(estimand$effect.measure2, CIF2_sel)
  a12 <- rep(0, length.out = n)

  denom_beta <- a11 * a22 - a12 * a12
  d_ey_d_beta_11 <- a22 / denom_beta
  d_ey_d_beta_12 <- -a12 / denom_beta
  d_ey_d_beta_22 <- a11 / denom_beta

  c12 <- 1 / survival_sel
  c11 <- a11 + c12
  c22 <- a22 + c12

  denom_alpha <- c11 * c22 - c12 * c12
  d_ey_d_alpha_11 <- c22 / denom_alpha
  d_ey_d_alpha_12 <- -c12 / denom_alpha
  d_ey_d_alpha_22 <- c11 / denom_alpha

  scale_matrix <- function(mat, vec) {
    if (is.null(mat) || ncol(mat) == 0) {
      return(NULL)
    }
    sweep(mat, 1, vec, `*`)
  }

  beta_11 <- scale_matrix(x_a, d_ey_d_beta_11)
  beta_12 <- scale_matrix(x_a, d_ey_d_beta_12)
  beta_22 <- scale_matrix(x_a, d_ey_d_beta_22)

  d_11 <- cbind(sweep(x_l, 1, d_ey_d_alpha_11, `*`), beta_11)
  d_12 <- cbind(sweep(x_l, 1, d_ey_d_alpha_12, `*`), beta_12)
  d_22 <- cbind(sweep(x_l, 1, d_ey_d_alpha_22, `*`), beta_22)
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

  out_defineExposureDesign <- defineExposureDesign(data, exposure, estimand$code.exposure.ref)
  x_a <- out_defineExposureDesign$x_a
  x_l <- model.matrix(Terms, mf)
  x_la <- cbind(x_l, x_a)

  potential.CIFs <- calculatePotentialRisk(alpha_beta, x_a, x_l, offset, estimand)
  one <- rep(1, nrow(x_l))
  a <- as.vector(x_a)
  ey_1 <- potential.CIFs[,2]*a + potential.CIFs[,1]*(one - a) #曝露分類未対応
  w11 <- 1 / (ey_1 * (1 - ey_1))
  wy_1 <- ip.weight * y_1
  wy_1ey_1 <- w11*(wy_1 - ey_1)
  d <- cbind(x_l, x_a)
  residual <- wy_1ey_1
  ret <- drop(crossprod(d, residual)) / nrow(x_l)

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
    ip.weight = ip.weight
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
  index.vector <- estimand$index.vector
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
#  AB1 <- score[1:n, 1:index.vector[2]]
  AB1 <- score[1:n, 1:index.vector[3]]
  for (i_para in 1:index.vector[2]) {
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
  hesse <- crossprod(x_la, w11 * out_calculateDSurvival) / n

  total_score <- AB1
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
  n <- nrow(x_l)
  K <- estimand$exposure.levels
  if (ncol(potential.CIFs) < K) stop("Column dimension of potential.CIFs is insufficient for the number of exposure levels.")

  seq_n <- seq_len(n)
  idx <- max.col(cbind(0, x_a), ties.method = "first")

  CIF_mat <- clampP(potential.CIFs[, seq_len(K), drop = FALSE], prob.bound)
  CIF_sel <- CIF_mat[cbind(seq_n, idx)]

  calculateA <- function(effect_measure, CIFs_selected) {
    CIFs_selected <- clampP(CIFs_selected, prob.bound)
    if (effect_measure == "RR") {
      return(1 / CIFs_selected)
    } else if (effect_measure == "OR") {
      denom <- clampP(1 - CIFs_selected, prob.bound)
      return(1 / CIFs_selected + 1 / denom)
    } else if (effect_measure == "SHR") {
      denom <- clampP(1 - CIFs_selected, prob.bound)
      tmp1 <- -1 / denom
      tmp2 <- log(denom)
      return(tmp1 / tmp2)
    } else {
      stop("Invalid effect_measure. Must be RR, OR or SHR.")
    }
  }

  a11 <- calculateA(estimand$effect.measure1, CIF_sel)
  d_ey_d_beta_11 <- 1 / a11
  d_ey_d_alpha_11 <- 1 / a11

  beta_11 <- if (ncol(x_a) == 0) NULL else sweep(x_a, 1, d_ey_d_beta_11, `*`)
  d_11 <- cbind(sweep(x_l, 1, d_ey_d_alpha_11, `*`), beta_11)
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

  out_defineExposureDesign <- defineExposureDesign(data, exposure, estimand$code.exposure.ref)
  x_a <- out_defineExposureDesign$x_a
  x_l <- model.matrix(Terms, mf)
  x_la <- cbind(x_l, x_a)
  index.vector <- estimand$index.vector

  one <- rep(1, nrow(x_l))
  a <- as.vector(x_a)
  score_beta <- NULL
  score_alpha1 <- 0
  score_alpha2 <- 0
  i_time <- 0
  alpha_beta_i <- rep(NA, index.vector[7])
  for (specific.time in time.point) {
    i_time <- i_time + 1                    #時間のインデックス
    i_para <- index.vector[1]*(i_time-1)+1  #パラメータのうち時間依存性切片項のインデックス
    alpha_beta_i[seq_len(index.vector[1])]                  <- alpha_beta[seq.int(i_para, i_para+index.vector[1]-1)]  #時間依存性切片項から共変量回帰係数までのインデックス
    alpha_beta_i[seq.int(index.vector[2], index.vector[3])] <- alpha_beta[index.vector[8]/2] #パラメータのうち曝露のインデックス, 曝露分類未対応

    y_0 <- ifelse(epsilon == estimand$code.censoring | t > specific.time, 1, 0)
    y_1 <- ifelse(epsilon == estimand$code.event1 & t <= specific.time, 1, 0)

    potential.CIFs <- calculatePotentialRisk(alpha_beta, x_a, x_l, offset, estimand)
    ey_1 <- potential.CIFs[,2]*a + potential.CIFs[,1]*(one - a) #曝露分類未対応
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
    tmp1 <- t(subscore[1:index.vector[1]])
    score_beta <- cbind(score_beta, tmp1)
    score_alpha1 <- score_alpha1 + subscore[index.vector[2]:index.vector[3]]
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
    x_l   = x_l
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

  out_defineExposureDesign <- defineExposureDesign(data, exposure, estimand$code.exposure.ref)
  x_a <- out_defineExposureDesign$x_a
  x_l <- model.matrix(Terms, mf)
  x_la <- cbind(x_l, x_a)
  index.vector <- estimand$index.vector

  #n_para_1 <- ncol(x_l)
  #n_para_2 <- n_para_1 + 1
  #n_para_3 <- n_para_1 + 2
  #n_para_4 <- 2*n_para_1 + 1
  #n_para_5 <- 2*n_para_1 + 2
  #n_para_6 <- length(time.point)*(n_para_5-2) + 2

  one <- rep(1, nrow(x_l))
#  a <- as.vector(x_a)
  score_beta <- NULL
  score_alpha1 <- 0
  score_alpha2 <- 0
  i_time <- 0
  alpha_beta_i <- rep(NA, index.vector[7])
  for (specific.time in time.point) {
    i_time <- i_time + 1
    i_para <- index.vector[1]*(i_time-1)+1
    alpha_beta_i[seq_len(index.vector[1])]                 <- alpha_beta[seq.int(i_para, i_para+index.vector[1]-1)]
    alpha_beta_i[seq.int(index.vector[2], index.vector[3])] <- alpha_beta[index.vector[8]/2]
    alpha_beta_i[seq.int(index.vector[4], index.vector[5])] <- alpha_beta[seq.int((index.vector[8]/2+i_para),(index.vector[8]/2+i_para+index.vector[1]-1))]
    alpha_beta_i[seq.int(index.vector[6], index.vector[7])] <- alpha_beta[index.vector[8]]
    #alpha_beta_i[1:n_para_1]        <- alpha_beta[i_para:(i_para+n_para_1-1)]
    #alpha_beta_i[n_para_2]          <- alpha_beta[n_para_6/2]
    #alpha_beta_i[n_para_3:n_para_4] <- alpha_beta[(n_para_6/2+i_para):(n_para_6/2+i_para+n_para_1-1)]
    #alpha_beta_i[n_para_5]          <- alpha_beta[n_para_6]

    y_0 <- ifelse(epsilon == 0 | t > specific.time, 1, 0)
    y_1 <- ifelse(epsilon == 1 & t <= specific.time, 1, 0)
    y_2 <- ifelse(epsilon == 2 & t <= specific.time, 1, 0)

    potential.CIFs <- calculatePotentialCIFs(alpha_beta,x_a,x_l,offset,epsilon,estimand,optim.method,prob.bound,initial.CIFs)
    ey <- calculateEY(potential.CIFs, x_a)
    ey_1 <- ey$ey_1
    ey_2 <- ey$ey_2

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

    tmp1 <- t(subscore[1:index.vector[1],])
    tmp2 <- t(subscore[index.vector[4]:index.vector[5],])
    score_beta <- cbind(score_beta, tmp1, tmp2)
    score_alpha1 <- score_alpha1 + subscore[index.vector[2]:index.vector[3],]
    score_alpha2 <- score_alpha2 + subscore[index.vector[6]:index.vector[7],]
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
    x_l   = x_l
  )
#  out <- list(ret = score)
  return(out)
}

solveEstimatingEquation <- function(
    nuisance.model,
    exposure,
    strata,
    normalized_data,
    outcome.type,
    estimand,
    optim.method,
    out_normalizeCovariate,
    prob.bound,
    alpha_beta_0
) {
  if (outcome.type == "COMPETING-RISK" | outcome.type == "SURVIVAL") {
    ip.weight <- calculateIPCW(nuisance.model, normalized_data, estimand$code.censoring, strata, estimand$time.point)
  } else if (outcome.type == "BINOMIAL") {
    ip.weight <- rep(1,nrow(normalized_data))
  } else if (outcome.type == "PROPORTIONAL" | outcome.type == "POLY-PROPORTIONAL") {
    ip.weight.matrix <- calculateIPCWMatrix(nuisance.model, normalized_data, estimand$code.censoring, strata, estimand, out_normalizeCovariate)
  }

  makeObjectiveFunction <- function() {
    out_ipcw <- list()
    initial.CIFs <- NULL

    call_and_capture <- function(fun, ...) {
      out_ipcw <<- do.call(fun, list(...))
      out_ipcw$ret
    }

    estimating_equation_i <- function(p) call_and_capture(
      estimating_equation_ipcw,
      formula = nuisance.model, data = normalized_data, exposure = exposure,
      ip.weight = ip.weight, alpha_beta = p, estimand = estimand,
      optim.method = NULL, prob.bound = prob.bound,
      initial.CIFs = initial.CIFs
    )

    estimating_equation_s <- function(p) call_and_capture(
      estimating_equation_survival,
      formula = nuisance.model, data = normalized_data, exposure = exposure,
      ip.weight = ip.weight, alpha_beta = p, estimand = estimand,
      prob.bound = prob.bound, initial.CIFs = initial.CIFs
    )

    estimating_equation_p <- function(p) call_and_capture(
      estimating_equation_proportional,
      formula = nuisance.model, data = normalized_data, exposure = exposure,
      ip.weight.matrix = ip.weight.matrix, alpha_beta = p, estimand = estimand,
      optim.method = NULL, prob.bound = prob.bound,
      initial.CIFs = initial.CIFs
    )

    estimating_equation_pp <- function(p) call_and_capture(
      estimating_equation_pproportional,
      formula = nuisance.model, data = normalized_data, exposure = exposure,
      ip.weight.matrix = ip.weight.matrix, alpha_beta = p, estimand = estimand,
      optim.method = NULL, prob.bound = prob.bound,
      initial.CIFs = initial.CIFs
    )

    setInitialCIFs <- function(new.CIFs) initial.CIFs <<- new.CIFs
    getResults     <- function() out_ipcw

    list(
      estimating_equation_i  = estimating_equation_i,
      estimating_equation_s  = estimating_equation_s,
      estimating_equation_p  = estimating_equation_p,
      estimating_equation_pp = estimating_equation_pp,
      setInitialCIFs = setInitialCIFs,
      getResults     = getResults
    )
  }

  assessConvergence <- function(new_params, current_params, current_obj_value, optim.parameter1, optim.parameter2, optim.parameter3) {
    if (any(abs(new_params) > optim.parameter3)) {
      stop("Estimates are either too large or too small, and convergence might not be achieved.")
    }
    param_diff <- abs(new_params - current_params)
    max.absolute.difference <- max(param_diff)
    relative.difference <- assessRelativeDifference(new_params, current_params)
    obj_value <- get_obj_value(new_params)
    converged <- (relative.difference <= optim.parameter1) || (obj_value <= optim.parameter2) || is_stalled(c(current_obj_value, obj_value))

    criteria1 <- (relative.difference <= optim.parameter1)
    criteria2 <- (obj_value <= optim.parameter2)
    criteria3 <- is_stalled(c(current_obj_value, obj_value))
    converged  <- (criteria1 || criteria2 || criteria3)
    converged.by <- if (!converged) NA_character_
    else if (criteria1)   "Converged in relative difference"
    else if (criteria2) "Converged in objective function"
    else "Stalled"

    list(converged = converged, converged.by=converged.by, relative.difference = relative.difference, max.absolute.difference = max.absolute.difference, obj_value = obj_value)
  }

  assessRelativeDifference <- function(new, old) {
    max(abs(new - old) / pmax(1, abs(old)))
  }

  is_stalled <- function(x, stall_patience=3, stall_eps=1e-3) {
    if (length(x) < stall_patience) return(FALSE)
    recent <- tail(x, stall_patience)
    (diff(range(recent)) / max(1e-12, mean(recent))) <= stall_eps
  }

  choose_estimating_equation <- function(outcome.type, obj) {
    if (outcome.type == "COMPETING-RISK") {
      obj$estimating_equation_i
    } else if (outcome.type == "SURVIVAL" | outcome.type == "BINOMIAL") {
      obj$estimating_equation_s
    } else if (outcome.type == "PROPORTIONAL") {
      obj$estimating_equation_p
    } else if (outcome.type == "POLY-PROPORTIONAL") {
      obj$estimating_equation_pp
    } else {
      stop("Unknown outcome.type: ", outcome.type)
    }
  }

  choose_nleqslv_method <- function(nleqslv.method) {
    if (nleqslv.method %in% c("nleqslv", "Broyden")) {
      "Broyden"
    } else if (nleqslv.method == "Newton") {
      "Newton"
    } else {
      stop("Unsupported nleqslv.method without optim(): ", nleqslv.method)
    }
  }

  get_obj_value <- switch(outcome.type,
                          "COMPETING-RISK"   = function(p) drop(crossprod(obj$estimating_equation_i(p))),
                          "SURVIVAL"         = function(p) drop(crossprod(obj$estimating_equation_s(p))),
                          "BINOMIAL"         = function(p) drop(crossprod(obj$estimating_equation_s(p))),
                          "PROPORTIONAL"     = function(p) drop(crossprod(obj$estimating_equation_p(p))),
                          "POLY-PROPORTIONAL"= function(p) drop(crossprod(obj$estimating_equation_pp(p))),
                          stop("Unknown outcome.type: ", outcome.type)
  )

  obj <- makeObjectiveFunction()
  estimating_fun <- choose_estimating_equation(outcome.type, obj)
  nleqslv_method <- choose_nleqslv_method(optim.method$nleqslv.method)

  iteration <- 0L
  max.absolute.difference <- Inf
  out_nleqslv <- NULL
  current_params <- alpha_beta_0
  current_obj_value <- numeric(0)
  trace_df  <- NULL

  computation.time0 <- proc.time()

  while ((iteration < optim.method$optim.parameter4) & (max.absolute.difference > optim.method$optim.parameter1)) {
    iteration <- iteration + 1L
    prev_params <- current_params

    out_nleqslv <- nleqslv(
      prev_params,
      estimating_fun,
      method  = optim.method$nleqslv_method,
      control = list(maxit = optim.method$optim.parameter5, allowSingular = FALSE)
    )
    new_params <- out_nleqslv$x

    current_obj_value <- get_obj_value(new_params)
    obj$setInitialCIFs(obj$getResults()$potential.CIFs)

    ac <- assessConvergence(
      new_params, prev_params, current_obj_value,
      optim.method$optim.parameter1, optim.method$optim.parameter2, optim.method$optim.parameter3
    )
    current_params <- new_params
    max.absolute.difference <- ac$max.absolute.difference
    if (ac$converged) break
  }
  return(current_params)
}


