## K-level residuals (length = 2K)
## log_p: [p1_ref, p1_L2..p1_LK, p2_ref, p2_L2..p2_LK]
## alpha1, alpha2: numeric length K-1  (ℓ = 2..K に対応)
## beta1, beta2: numeric scalar
## effect.measure*: "RR" | "OR" | "SHR" | ""(=SHR)
residuals_CIFs_K <- function(log_p, alpha1, beta1, alpha2, beta2, estimand, prob.bound, c1,c2) {
  effect.measure1 <- estimand$effect.measure1
  effect.measure2 <- estimand$effect.measure2
  clogp <- clampLogP(as.numeric(log_p))
  m <- length(clogp); if (m %% 2L != 0L) stop("log_p length must be even (2K).")
  K <- m %/% 2L
  if (length(alpha1) != K-1 || length(alpha2) != K-1)
    stop("alpha1/alpha2 must have length K-1.")
  if (is.null(effect.measure1) || identical(effect.measure1, "")) effect.measure1 <- "SHR"
  if (is.null(effect.measure2) || identical(effect.measure2, "")) effect.measure2 <- "SHR"

  i1 <- 1:K                      # indices for event1 block
  i2 <- (K+1):(2*K)              # indices for event2 block
  p1 <- exp(clogp[i1])
  p2 <- exp(clogp[i2])

  ## 残り確率は「各レベルごと」に計算（ref と各ℓのペア）
  rem <- pmax(1 - p1 - p2, prob.bound)
  lrem <- log(rem)

  r <- numeric(2*K)

  ## ---- α constraints: 各 ℓ=2..K で2本（event1/event2）
  for (ell in 2:K) {
    ## event1: α1_ℓ-1 - log p1_ref - log p1_ℓ + log rem_ref + log rem_ℓ = 0
    r[ell - 1]        <- alpha1[ell - 1] - clogp[i1[1]] - clogp[i1[ell]] + lrem[1] + lrem[ell]
    ## event2: α2_ℓ-1 - log p2_ref - log p2_ℓ + log rem_ref + log rem_ℓ = 0
    r[(K - 1) + (ell - 1)] <- alpha2[ell - 1] - clogp[i2[1]] - clogp[i2[ell]] + lrem[1] + lrem[ell]
  }

  ## ---- β constraints（各イベント1本に集約）
  ## RR:      beta - log p_L + log p_ref
  ## OR:      beta - [logit(p_L) - logit(p_ref)]
  ## SHR:     exp(beta) - [ log(1-p_L) / log(1-p_ref) ]
  agg_res <- function(meas, pref, pL, clogpref, clogpL, beta) {
    if (meas == "RR") {
      sum(beta - clogpL + clogpref)
    } else if (meas == "OR") {
      sum(beta - (log(pL) - log1p(-pL)) + (log(pref) - log1p(-pref)))
    } else if (meas == "SHR") {
      Bi <- log1p(-pref)        # < 0
      Bj <- log1p(-pL)
      sum(exp(beta) - (Bj / Bi))
    } else stop("effect.measure must be RR/OR/SHR.")
  }
  r[2*K - 1] <- agg_res(effect.measure1, p1[1], p1[2:K], clogp[i1[1]], clogp[i1[2:K]], beta1)
  r[2*K]     <- agg_res(effect.measure2, p2[1], p2[2:K], clogp[i2[1]], clogp[i2[2:K]], beta2)

  r
}

jacobian_CIFs_K <- function(log_p, alpha1, beta1, alpha2, beta2, estimand, prob.bound, c1,c2) {
  effect.measure1 <- estimand$effect.measure1
  effect.measure2 <- estimand$effect.measure2
  clogp <- clampLogP(as.numeric(log_p))
  m <- length(clogp); if (m %% 2L != 0L) stop("log_p length must be even (2K).")
  K <- m %/% 2L
  i1 <- 1:K; i2 <- (K+1):(2*K)
  p1 <- exp(clogp[i1]); p2 <- exp(clogp[i2])
  rem <- pmax(1 - p1 - p2, prob.bound)

  dlogrem_dlp1 <- (-p1) / rem
  dlogrem_dlp2 <- (-p2) / rem

  J <- matrix(0.0, nrow = 2*K, ncol = 2*K)

  ## α 行
  for (ell in 2:K) {
    r1 <- ell - 1
    J[r1,   i1[1]]   <- -1 + dlogrem_dlp1[1]
    J[r1,   i2[1]]   <-       dlogrem_dlp2[1]
    J[r1,   i1[ell]] <- -1 + dlogrem_dlp1[ell]
    J[r1,   i2[ell]] <-       dlogrem_dlp2[ell]

    r2 <- (K - 1) + (ell - 1)
    J[r2,   i2[1]]   <- -1 + dlogrem_dlp2[1]
    J[r2,   i1[1]]   <-       dlogrem_dlp1[1]
    J[r2,   i2[ell]] <- -1 + dlogrem_dlp2[ell]
    J[r2,   i1[ell]] <-       dlogrem_dlp1[ell]
  }

  add_effect_grad <- function(meas, pref, pL, row, idx_ref, idx_L) {
    if (meas == "RR") {
      ## d/d log p_ref: +1 を (K-1) 回分加算、d/d log p_L: -1
      J[row, idx_ref] <- J[row, idx_ref] + (K - 1)
      for (k in idx_L) J[row, k] <- J[row, k] - 1
    } else if (meas == "OR") {
      ## d/d log p = d/dp * dp/dlogp = g'(p) * p
      ## g(p)=logit(p)=log p - log(1-p) → g'(p)=1/p + 1/(1-p)
      ## ref: +[1 + p_ref/(1-p_ref)] を (K-1) 回
      J[row, idx_ref] <- J[row, idx_ref] + (1 + pref/(1 - pref)) * (K - 1)
      for (k in seq_along(idx_L)) {
        pl <- pL[k]
        J[row, idx_L[k]] <- J[row, idx_L[k]] + (-1 - pl/(1 - pl))
      }
    } else if (meas == "SHR") {
      Bi <- log1p(-pref)   # <0
      ## ref への寄与: sum_ℓ (Bj * (-pref/(1-pref))) / Bi^2
      coef_ref <- sum(log1p(-pL)) * ( -pref/(1 - pref) ) / (Bi^2)
      J[row, idx_ref] <- J[row, idx_ref] + coef_ref
      ## 各 L: pL / ((1 - pL) * Bi)
      for (k in seq_along(idx_L)) {
        pl <- pL[k]
        J[row, idx_L[k]] <- J[row, idx_L[k]] + pl / ((1 - pl) * Bi)
      }
    } else stop("effect.measure must be RR/OR/SHR.")
  }

  row_e1 <- 2*K - 1
  row_e2 <- 2*K
  add_effect_grad(effect.measure1, pref = p1[1], pL = p1[2:K], row = row_e1, idx_ref = i1[1], idx_L = i1[2:K])
  add_effect_grad(effect.measure2, pref = p2[1], pL = p2[2:K], row = row_e2, idx_ref = i2[1], idx_L = i2[2:K])

  J
}

chooseEffectMeasureLevenbergMarquardt <- function(effect.measure, index) {
  i <- index[1]; j <- index[2]
  if (effect.measure == "RR") {
    res <- function(p, clogp, beta) {
      beta - clogp[j] + clogp[i]
    }
    jac <- function(p, clogp, beta) {
      g <- numeric(4); g[i] <-  1; g[j] <- -1; g
    }
  } else if (effect.measure == "OR") {
    res <- function(p, clogp, beta) {
      beta - clogp[j] + clogp[i] + log1p(-p[j]) - log1p(-p[i])
    }
    jac <- function(p, clogp, beta) {
      g <- numeric(4)
      g[i] <-  1 +  p[i]/(1 - p[i])
      g[j] <- -1 + (-p[j]/(1 - p[j]))
      g
    }
  } else if (effect.measure == "SHR" || identical(effect.measure, "")) {
    res <- function(p, clogp, beta) {
      exp(beta) - (log1p(-p[j]) / log1p(-p[i]))
    }
    jac <- function(p, clogp, beta) {
      g <- numeric(4)
      Bi <- log1p(-p[i])  # < 0
      Bj <- log1p(-p[j])
      g[i] <-  ( Bj * (-p[i]/(1 - p[i])) ) / (Bi^2)
      g[j] <-   p[j] / ((1 - p[j]) * Bi)
      g
    }
  } else {
    stop("Invalid measure: must be 'RR','OR','SHR' (or '' as SHR).")
  }
  return(list(res = res, jac = jac, index = c(i, j)))
}

residuals_CIFs_generic <- function(log_p, alpha1, beta1, alpha2, beta2, estimand, prob.bound, cemlm1, cemlm2) {
  clogp <- clampLogP(as.numeric(log_p))
  if (length(clogp) < 4L) clogp <- rep(clogp, 4L)
  p <- exp(clogp)
  rem12 <- max(1 - p[1] - p[2], prob.bound)
  rem34 <- max(1 - p[3] - p[4], prob.bound)
  lp0102 <- log(rem12) + log(rem34)
  r <- numeric(4)
  r[1] <- alpha1 - clogp[1] - clogp[3] + lp0102
  r[3] <- alpha2 - clogp[2] - clogp[4] + lp0102
  r[2] <- cemlm1$res(p, clogp, beta1)
  r[4] <- cemlm2$res(p, clogp, beta2)
  return(r)
}

jacobian_CIFs_generic <- function(log_p, alpha1, beta1, alpha2, beta2, estimand, prob.bound, cemlm1, cemlm2) {
  clogp <- clampLogP(as.numeric(log_p))
  if (length(clogp) < 4L) clogp <- rep(clogp, 4L)
  p <- exp(clogp)
  rem12 <- max(1 - p[1] - p[2], prob.bound)
  rem34 <- max(1 - p[3] - p[4], prob.bound)
  dlp12 <- c(-p[1]/rem12, -p[2]/rem12)
  dlp34 <- c(-p[3]/rem34, -p[4]/rem34)
  J <- matrix(0.0, 4, 4)
  J[1,] <- c(-1 + dlp12[1], dlp12[2], -1 + dlp34[1], dlp34[2])
  J[3,] <- c(dlp12[1], -1 + dlp12[2], dlp34[1], -1 + dlp34[2])
  J[2,] <- cemlm1$jac(p, clogp, beta1)
  J[4,] <- cemlm2$jac(p, clogp, beta2)
  return(J)
}

LevenbergMarquardt <- function(start,
                               res_fun, jac_fun,
                               maxit = 30,
                               ftol  = 1e-10,
                               ptol  = 1e-10,
                               lambda0    = 1e-3,
                               lambda_max = 10,
                               lambda_min = 0.1,
                               lambda_increment    = 2.0,
                               lambda_decrement    = 0.5,
                               do_linesearch = TRUE,
                               ls_c      = 1e-4,
                               ls_shrink = 0.5,
                               verbose   = FALSE) {
  lp   <- as.numeric(start)
  r    <- res_fun(lp)
  f2   <- drop(crossprod(r))
  J    <- jac_fun(lp)
  g    <- drop(crossprod(J, r))
  lam  <- lambda0
  prev_f2 <- Inf

  if (isTRUE(verbose)) {
    hist <- list(it=integer(), f2=double(),
                 grad_inf=double(), lambda=double(),
                 step_norm=double(), rho=double())
  }

  for (k in seq_len(maxit)) {
    grad_inf <- max(abs(g))
    if (is.finite(grad_inf) && grad_inf < 1e-6) break
    A <- crossprod(J)
    diag(A) <- diag(A) + lam
    R <- try(chol(A), silent = TRUE)
    if (inherits(R, "try-error")) {
      lam <- min(lam * 10, lambda_max)
      next
    }
    delta <- -backsolve(R, forwardsolve(t(R), g))
    step_norm <- max(abs(delta))

    if (is.finite(step_norm) &&
        step_norm <= ptol * (ptol + max(1, max(abs(lp))))) break

    pred <- 0.5 * sum(delta * (lam * delta - g))
    if (!is.finite(pred) || pred <= 0) pred <- .Machine$double.eps

    lp_try <- lp + delta
    r_try  <- res_fun(lp_try)
    f2_try <- drop(crossprod(r_try))

    rho <- (f2 - f2_try) / pred
    accepted <- FALSE

    if (is.finite(rho) && rho > 0) {
      lp <- lp_try; r <- r_try; f2 <- f2_try
      lam <- max(lambda_min, lam * max(lambda_decrement, 1/(1 + rho)))
      accepted <- TRUE
    } else {
      lam <- min(lambda_max, lam * (lambda_increment * (1 + abs(ifelse(is.finite(rho), rho, 0)))))
      if (do_linesearch) {
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

    J <- jac_fun(lp)
    g <- drop(crossprod(J, r))

    if (k > 1L &&
        is.finite(prev_f2) && is.finite(f2) &&
        abs(prev_f2 - f2) <= ftol * (abs(f2) + ftol)) break

    if (isTRUE(verbose)) { # ログの記録
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
  return(lp)
}

callLevenbergMarquardt <- function(log_CIFs0, alpha1, beta1, alpha2, beta2, optim.method, estimand, prob.bound, cemlm1, cemlm2)
{
  #  res_fun <- function(lp) residuals_CIFs_K(lp, alpha1, beta1, alpha2, beta2, estimand, prob.bound, cemlm1, cemlm2)
  #  jac_fun <- function(lp) jacobian_CIFs_K(lp, alpha1, beta1, alpha2, beta2, estimand, prob.bound, cemlm1, cemlm2)
  res_fun <- function(lp) residuals_CIFs_generic(lp, alpha1, beta1, alpha2, beta2, estimand, prob.bound, cemlm1, cemlm2)
  jac_fun <- function(lp) jacobian_CIFs_generic(lp, alpha1, beta1, alpha2, beta2, estimand, prob.bound, cemlm1, cemlm2)
  out_LevenbergMarquardt <- LevenbergMarquardt(
    start         = log_CIFs0,
    res_fun       = res_fun,
    jac_fun       = jac_fun,
    maxit         = optim.method$optim.parameter6   %||% 30,
    ftol          = optim.method$optim.parameter7   %||% 1e-10,
    ptol          = optim.method$optim.parameter8   %||% 1e-10,
    lambda0       = optim.method$optim.parameter9   %||% 1e-3,
    lambda_max    = optim.method$optim.parameter10  %||% 10,
    lambda_min    = optim.method$optim.parameter11  %||% 0.1,
    lambda_increment  = optim.method$optim.parameter12  %||% 2.0,
    lambda_decrement  = optim.method$optim.parameter13  %||% 0.5,
    do_linesearch = TRUE,
    ls_c          = 1e-4,
    ls_shrink     = 0.5,
    verbose       = FALSE
  )
  return(out_LevenbergMarquardt)
}

calculatePotentialCIFs <- function(alpha_beta_tmp, x_a, x_l, offset, epsilon, estimand, optim.method, prob.bound, initial.CIFs = NULL) {

  i_parameter <- rep(NA_integer_, 7L)
  i_parameter <- calculateIndexForParameter(i_parameter, x_l, x_a)
  alpha_1    <- alpha_beta_tmp[seq_len(i_parameter[1])]
  beta_tmp_1 <- alpha_beta_tmp[seq.int(i_parameter[2], i_parameter[3])]
  alpha_2    <- alpha_beta_tmp[seq.int(i_parameter[4], i_parameter[5])]
  beta_tmp_2 <- alpha_beta_tmp[seq.int(i_parameter[6], i_parameter[7])]

  alpha_tmp_1 <- as.numeric(x_l %*% matrix(alpha_1, ncol = 1) + offset)
  alpha_tmp_2 <- as.numeric(x_l %*% matrix(alpha_2, ncol = 1) + offset)

  n  <- length(epsilon)
  p0 <- c(
    sum(epsilon == estimand$code.event1) / n + prob.bound,
    sum(epsilon == estimand$code.event2) / n + prob.bound,
    sum(epsilon == estimand$code.event1) / n + prob.bound,
    sum(epsilon == estimand$code.event2) / n + prob.bound
  )
  p0     <- clampP(p0, prob.bound)
  log_p0 <- log(p0)

  keys <- apply(x_l, 1, function(r) paste0(r, collapse = "\r"))
  uniq <- match(keys, unique(keys))
  cache_log_CIFs <- vector("list", length = max(uniq))

  cemlm1 <- chooseEffectMeasureLevenbergMarquardt(estimand$effect.measure1, c(1,3))
  cemlm2 <- chooseEffectMeasureLevenbergMarquardt(estimand$effect.measure2 %||% "SHR", c(2,4))
  potential.CIFs <- matrix(NA_real_, nrow = nrow(x_l), ncol = 4L)

  for (i in seq_len(nrow(x_l))) {
    k <- uniq[i]

    if (!is.null(cache_log_CIFs[[k]])) {
      log_CIFs <- cache_log_CIFs[[k]]
    } else {
      CIFs0 <- if (!is.null(initial.CIFs)) as.numeric(initial.CIFs[i, 1:4, drop = FALSE]) else exp(log_p0)
      log_CIFs0 <- log(clampP(CIFs0, prob.bound))

      log_CIFs <- callLevenbergMarquardt(
        log_CIFs0  = log_CIFs0,
        alpha1    = alpha_tmp_1[i],
        beta1     = beta_tmp_1,
        alpha2    = alpha_tmp_2[i],
        beta2     = beta_tmp_2,
        estimand  = estimand,
        optim.method = optim.method,
        prob.bound = prob.bound,
        cemlm1 = cemlm1,
        cemlm2 = cemlm2
      )
      cache_log_CIFs[[k]] <- log_CIFs
    }
    potential.CIFs[i, ] <- clampP(exp(log_CIFs), prob.bound)
  }
  colnames(potential.CIFs) <- c("p10", "p20", "p11", "p21")
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
    clog_p <- clampLogP(as.numeric(log_p))
    if (length(clog_p) < 4L) clog_p <- rep(clog_p, length.out = 4L)
    exp_lp <- exp(clog_p)

    ## 対の残り確率（log(1 - p1 - p2)）も安全に
    rem12 <- 1 - exp_lp[1] - exp_lp[2]
    rem34 <- 1 - exp_lp[3] - exp_lp[4]
    if (rem12 < prob.bound || rem34 < prob.bound) {
      lp0102 <- log(prob.bound)
    } else {
      lp01   <- log1p(-(exp_lp[1] + exp_lp[2]))
      lp02   <- log1p(-(exp_lp[3] + exp_lp[4]))
      lp0102 <- lp01 + lp02
    }

    ret <- numeric(4L)
    if (estimand$effect.measure1 == "RR") {
      ret[1] <- alpha_tmp_1 - clog_p[1] - clog_p[3] + lp0102
      ret[2] <- beta_tmp_1  - clog_p[3] + clog_p[1]
    } else if (estimand$effect.measure1 == "OR") {
      ret[1] <- alpha_tmp_1 - clog_p[1] - clog_p[3] + lp0102
      ret[2] <- beta_tmp_1  - clog_p[3] + clog_p[1] +
        log1p(-exp_lp[3]) - log1p(-exp_lp[1])
    } else if (estimand$effect.measure1 == "SHR") {
      ret[1] <- alpha_tmp_1 - clog_p[1] - clog_p[3] + lp0102
      ret[2] <- exp(beta_tmp_1) - ( log1p(-exp_lp[3]) / log1p(-exp_lp[1]) )
    } else {
      stop("Invalid effect_measure1. Must be RR, OR or SHR.")
    }

    if (estimand$effect.measure2 == "RR") {
      ret[3] <- alpha_tmp_2 - clog_p[2] - clog_p[4] + lp0102
      ret[4] <- beta_tmp_2  - clog_p[4] + clog_p[2]
    } else if (estimand$effect.measure2 == "OR") {
      ret[3] <- alpha_tmp_2 - clog_p[2] - clog_p[4] + lp0102
      ret[4] <- beta_tmp_2  - clog_p[4] + clog_p[2] +
        log1p(-exp_lp[4]) - log1p(-exp_lp[2])
    } else if (estimand$effect.measure2 == '') {
      ret[3] <- alpha_tmp_2 - clog_p[2] - clog_p[4] + lp0102
      ret[4] <- exp(beta_tmp_2) - ( log1p(-exp_lp[4]) / log1p(-exp_lp[2]) )
    } else {
      stop("Invalid effect_measure2. Must be RR, OR or SHR.")
    }
    sum(ret^2)
  }
  objective_function(log_p)
}

calculatePotentialRisk <- function(alpha_beta, x_a, x_l, offset, estimand) {
  i_parameter <- rep(NA_integer_, 7L)
  i_parameter <- calculateIndexForParameter(i_parameter, x_l, x_a)

  one <- rep(1, nrow(x_l))
  alpha_beta_ <- as.matrix(as.vector(alpha_beta))
  if (ncol(alpha_beta_) == 1L) alpha_beta_ <- t(alpha_beta_)
  phi <- x_l %*% alpha_beta_[, seq_len(i_parameter[1])] + offset
  theta <- one * alpha_beta_[, i_parameter[2]]

  if (estimand$effect.measure1 == "RR") {
    expphi <- exp(phi)
    exptheta <- exp(theta)
    tol <- 1e-8
    if (all(abs(phi) < tol)) {
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
  } else if (estimand$effect.measure1 == "OR") {
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

# 残差ベクトル ret を返す（sum(ret^2) ではなく ret を直接返す）
residuals_CIFs <- function(
    log_p, alpha_tmp_1, beta_tmp_1, alpha_tmp_2, beta_tmp_2, estimand, prob.bound
) {
  clog_p <- clampLogP(as.numeric(log_p))
  if (length(clog_p) < 4L) clog_p <- rep(clog_p, length.out = 4L)
  p1 <- exp(clog_p[1]); p2 <- exp(clog_p[2]); p3 <- exp(clog_p[3]); p4 <- exp(clog_p[4])

  # 安全な残り確率
  rem12 <- max(1 - p1 - p2, prob.bound)
  rem34 <- max(1 - p3 - p4, prob.bound)
  lp0102 <- log(rem12) + log(rem34)

  ret <- numeric(4L)

  # effect1
  if (estimand$effect.measure1 == "RR") {
    ret[1] <- alpha_tmp_1 - clog_p[1] - clog_p[3] + lp0102
    ret[2] <- beta_tmp_1  - clog_p[3] + clog_p[1]
  } else if (estimand$effect.measure1 == "OR") {
    ret[1] <- alpha_tmp_1 - clog_p[1] - clog_p[3] + lp0102
    ret[2] <- beta_tmp_1  - clog_p[3] + clog_p[1] +
      log1p(-p3) - log1p(-p1)
  } else if (estimand$effect.measure1 == "SHR") {
    ret[1] <- alpha_tmp_1 - clog_p[1] - clog_p[3] + lp0102
    ret[2] <- exp(beta_tmp_1) - ( log1p(-p3) / log1p(-p1) )
  } else stop("Invalid effect_measure1.")

  # effect2
  if (estimand$effect.measure2 == "RR") {
    ret[3] <- alpha_tmp_2 - clog_p[2] - clog_p[4] + lp0102
    ret[4] <- beta_tmp_2  - clog_p[4] + clog_p[2]
  } else if (estimand$effect.measure2 == "OR") {
    ret[3] <- alpha_tmp_2 - clog_p[2] - clog_p[4] + lp0102
    ret[4] <- beta_tmp_2  - clog_p[4] + clog_p[2] +
      log1p(-p4) - log1p(-p2)
  } else if (estimand$effect.measure2 == "SHR" || estimand$effect.measure2 == "") {
    ret[3] <- alpha_tmp_2 - clog_p[2] - clog_p[4] + lp0102
    ret[4] <- exp(beta_tmp_2) - ( log1p(-p4) / log1p(-p2) )
  } else stop("Invalid effect_measure2.")

  ret
}
jacobian_CIFs <- function(
    log_p, alpha_tmp_1, beta_tmp_1, alpha_tmp_2, beta_tmp_2, estimand, prob.bound
) {
  clog_p <- clampLogP(as.numeric(log_p))
  if (length(clog_p) < 4L) clog_p <- rep(clog_p, length.out = 4L)
  p1 <- exp(clog_p[1]); p2 <- exp(clog_p[2]); p3 <- exp(clog_p[3]); p4 <- exp(clog_p[4])
  rem12 <- max(1 - p1 - p2, prob.bound)
  rem34 <- max(1 - p3 - p4, prob.bound)

  # dl/dlogp for logs of remaining mass
  dlp12_dlp1 <- -p1 / rem12
  dlp12_dlp2 <- -p2 / rem12
  dlp34_dlp3 <- -p3 / rem34
  dlp34_dlp4 <- -p4 / rem34

  J <- matrix(0.0, nrow = 4L, ncol = 4L)
  # ret1 = alpha1 - logp1 - logp3 + log(rem12)+log(rem34)
  J[1,1] <- -1 + dlp12_dlp1               # ∂ret1/∂logp1
  J[1,2] <-      dlp12_dlp2               # ∂ret1/∂logp2
  J[1,3] <- -1 + dlp34_dlp3               # ∂ret1/∂logp3
  J[1,4] <-      dlp34_dlp4               # ∂ret1/∂logp4

  if (estimand$effect.measure1 == "RR") {
    # ret2 = beta1 - logp3 + logp1
    J[2,1] <-  1
    J[2,3] <- -1
  } else if (estimand$effect.measure1 == "OR") {
    # ret2 = beta1 - logp3 + logp1 + log(1-p3) - log(1-p1)
    J[2,1] <-  1 +  p1/(1-p1)
    J[2,3] <- -1 + (-p3/(1-p3))
  } else if (estimand$effect.measure1 == "SHR") {
    # ret2 = exp(beta1) - log(1-p3)/log(1-p1)
    A <- log1p(-p3); B <- log1p(-p1)              # B<0
    J[2,1] <-  -( -A * ( -p1/(1-p1) ) ) / (B^2)   # =  A * (-p1/(1-p1)) / B^2
    J[2,3] <-  - ( (-p3/(1-p3)) / B )             # =  p3/((1-p3)*B)
  }

  # ret3 = alpha2 - logp2 - logp4 + log(rem12)+log(rem34)
  J[3,1] <-      dlp12_dlp1
  J[3,2] <- -1 + dlp12_dlp2
  J[3,3] <-      dlp34_dlp3
  J[3,4] <- -1 + dlp34_dlp4

  if (estimand$effect.measure2 == "RR") {
    # ret4 = beta2 - logp4 + logp2
    J[4,2] <-  1
    J[4,4] <- -1
  } else if (estimand$effect.measure2 == "OR") {
    # ret4 = beta2 - logp4 + logp2 + log(1-p4) - log(1-p2)
    J[4,2] <-  1 +  p2/(1-p2)
    J[4,4] <- -1 + (-p4/(1-p4))
  } else if (estimand$effect.measure2 == "SHR" || estimand$effect.measure2 == "") {
    # ret4 = exp(beta2) - log(1-p4)/log(1-p2)
    A <- log1p(-p4); B <- log1p(-p2)
    J[4,2] <-  -( -A * ( -p2/(1-p2) ) ) / (B^2)   # =  A * (-p2/(1-p2)) / B^2
    J[4,4] <-  - ( (-p4/(1-p4)) / B )             # =  p4/((1-p4)*B)
  }

  J
}
