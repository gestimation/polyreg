calculatePotentialCIFs_parallel <- function(
    alpha_beta_tmp,
    x_a,
    x_l,
    offset,
    epsilon,
    estimand,
    i_parameter,
    optim.method,
    prob.bound,
    initial.CIFs = NULL
) {

  safe_exp <- function(x) {
    x <- ifelse(is.finite(x), x, 700) # exp(700) ~ 1e304
    exp(x)
  }

  clampCIFs <- function(p, eps = 1e-5) {
    p <- pmin(pmax(p, eps), 1 - eps)
    s0 <- p[1] + p[2]
    s1 <- p[3] + p[4]
    max_sum <- 1 - eps
    if (s0 > max_sum) {
      p[1:2] <- p[1:2] * (max_sum / s0)
    }
    if (s1 > max_sum) {
      p[3:4] <- p[3:4] * (max_sum / s1)
    }
    return(p)
  }

  index1 <- seq_len(i_parameter[1])
  index23 <- seq.int(i_parameter[2], i_parameter[3])
  index45 <- seq.int(i_parameter[4], i_parameter[5])
  index67 <- seq.int(i_parameter[6], i_parameter[7])
  alpha_1 <- alpha_beta_tmp[index1]
  beta_tmp_1  <- alpha_beta_tmp[index23]
  alpha_2 <- alpha_beta_tmp[index45]
  beta_tmp_2  <- alpha_beta_tmp[index67]
  alpha_tmp_1 <- as.numeric(x_l %*% matrix(alpha_1, ncol = 1) + offset)                                  # ### INDEX-SAFE
  alpha_tmp_2 <- as.numeric(x_l %*% matrix(alpha_2, ncol = 1) + offset)                                  # ### INDEX-SAFE
  n <- length(epsilon)
  freq1 <- sum(epsilon == estimand$code.event1) / n + prob.bound
  freq2 <- sum(epsilon == estimand$code.event2) / n + prob.bound
  log_p0 <- log(c(freq1, freq2, freq1, freq2))

  # ---- 重複検出（キー）----
  x_l_key <- apply(x_l, 1, paste0, collapse = "_")
  unique_keys <- !duplicated(x_l_key)
  key_to_index <- match(x_l_key, x_l_key[unique_keys])

  # ---- バッチ構成 ----
  unique_idx <- which(unique_keys)
  batch.size <- if (!is.null(optim.method$computation.order.batch.size)) {
    optim.method$computation.order.batch.size
  } else {
    cores <- parallel::detectCores()
    max(10, floor(length(unique_idx) / (cores * 2))) * 2
  }
  batch_indices <- split(unique_idx, ceiling(seq_along(unique_idx) / batch.size))

  # 目的関数ラッパ（log_p 入力、内部で安全化）
  obj_fun <- function(lp, a1, b1, a2, b2, est, bound) {
    # lp が非有限なら log_p0 に置換
    if (any(!is.finite(lp))) lp <- log_p0
    estimating_equation_CIFs(
      log_p       = lp,
      alpha_tmp_1 = a1,
      beta_tmp_1  = b1,
      alpha_tmp_2 = a2,
      beta_tmp_2  = b2,
      estimand    = est,
      prob.bound  = bound
    )
  }

  # ---- 1点 i_x を解く関数（堅牢化版）----
  solve_one <- function(i_x) {
    # 初期点
    lp_init <- if (is.null(initial.CIFs)) log_p0 else log(initial.CIFs[i_x, ])

    try_paths <- list(
      list(method = "BFGS", par = lp_init,        maxit = optim.method$optim.parameter7, reltol = optim.method$optim.parameter6),
      list(method = "BFGS", par = log_p0,         maxit = optim.method$optim.parameter7 * 2, reltol = optim.method$optim.parameter6 * 0.1),
      list(method = "Nelder-Mead", par = lp_init, maxit = optim.method$optim.parameter7, reltol = optim.method$optim.parameter6),
      list(method = "Nelder-Mead", par = log_p0,  maxit = optim.method$optim.parameter7 * 2, reltol = optim.method$optim.parameter6)
    )

    if (optim.method$optim.parameter9 > 0) {
      for (k in seq_len(optim.method$optim.parameter9)) {
        try_paths <- c(try_paths, list(
          list(method = "BFGS",
               par = lp_init + rnorm(4, sd = optim.method$optim.parameter8),
               maxit = optim.method$optim.parameter7, reltol = optim.method$optim.parameter6)
        ))
      }
      for (k in seq_len(optim.method$optim.parameter9)) {
        try_paths <- c(try_paths, list(
          list(method = "Nelder-Mead",
               par = log_p0 + rnorm(4, sd = optim.method$optim.parameter8),
               maxit = optim.method$optim.parameter7, reltol = optim.method$optim.parameter6)
        ))
      }
    }

    best_val <- Inf
    best_par <- lp_init
    succeeded <- FALSE

    for (cfg in try_paths) {
      sol <- tryCatch(
        optim(
          par     = cfg$par,
          fn      = function(lp) obj_fun(lp,
                                         a1 = alpha_tmp_1[i_x],
                                         b1 = beta_tmp_1,
                                         a2 = alpha_tmp_2[i_x],
                                         b2 = beta_tmp_2,
                                         est = estimand,
                                         bound = prob.bound),
          method  = cfg$method,
          control = list(maxit = cfg$maxit, reltol = cfg$reltol)
        ),
        error = function(e) NULL
      )

      if (is.null(sol)) next

      # sol$par 非有限ならスキップ
      if (!all(is.finite(sol$par))) next

      # 最良値の更新（convergence==0 を最優先）
      if (sol$convergence == 0) {
        best_par <- sol$par
        best_val <- sol$value
        succeeded <- TRUE
        break
      } else {
        # 失敗でも値が改善してたら暫定採用
        if (is.finite(sol$value) && sol$value < best_val) {
          best_par <- sol$par
          best_val <- sol$value
        }
      }
    }

    p <- safe_exp(best_par)
    if (!all(is.finite(p))) {
      p <- safe_exp(log_p0)
    }
    p <- clampCIFs(p, prob.bound)
    p <- pmin(pmax(p, prob.bound), 1 - prob.bound)
    p
  }

  solve_CIF_batch <- function(idx) {
    mat <- matrix(
      data = unlist(lapply(idx, function(i_x) {
        out <- tryCatch(solve_one(i_x), error = function(e) NULL)
        if (is.null(out)) {
          out <- clampCIFs(safe_exp(log_p0), prob.bound)
        }
        out
      })),
      ncol = 4,
      byrow = TRUE
    )
    mat
  }

  result_list <- future.apply::future_lapply(batch_indices, solve_CIF_batch, future.seed = TRUE)
  unique_CIFs <- do.call(rbind, result_list)

  # 万一の NA/Inf を最終保険で埋める
  bad <- !is.finite(unique_CIFs)
  if (any(bad)) {
    fallback <- matrix(rep(clampCIFs(safe_exp(log_p0), prob.bound),
                           length.out = length(unique_idx) * 4),
                       ncol = 4, byrow = TRUE)
    unique_CIFs[bad] <- fallback[bad]
  }

  CIFs_all <- unique_CIFs[key_to_index, , drop = FALSE]
  CIFs_all
}

calculatePotentialCIFs_old <- function(
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
  i_parameter <- rep(NA_integer_, 7L)
  i_parameter <- calculateIndexForParameter(i_parameter, x_l, x_a)
  index1 <- seq_len(i_parameter[1])
  index23 <- seq.int(i_parameter[2], i_parameter[3])
  index45 <- seq.int(i_parameter[4], i_parameter[5])
  index67 <- seq.int(i_parameter[6], i_parameter[7])
  alpha_1 <- alpha_beta_tmp[index1]
  beta_tmp_1  <- alpha_beta_tmp[index23]
  alpha_2 <- alpha_beta_tmp[index45]
  beta_tmp_2  <- alpha_beta_tmp[index67]

  alpha_tmp_1 <- as.numeric(x_l %*% matrix(alpha_1, ncol = 1) + offset)                                  # ### INDEX-SAFE
  alpha_tmp_2 <- as.numeric(x_l %*% matrix(alpha_2, ncol = 1) + offset)                                  # ### INDEX-SAFE

  n <- length(epsilon)
  p0 <- c(
    sum(epsilon == estimand$code.event1) / n + prob.bound,
    sum(epsilon == estimand$code.event2) / n + prob.bound,
    sum(epsilon == estimand$code.event1) / n + prob.bound,
    sum(epsilon == estimand$code.event2) / n + prob.bound
  )
  p0     <- clampP(p0,prob.bound)
  log_p0 <- log(p0)

  list.CIFs     <- vector("list", nrow(x_l))
  previous.CIFs <- NULL

  for (i_x in seq_len(nrow(x_l))) {
    if (i_x > 1) {                                                                                       # ### INDEX-SAFE
      same_row <- isTRUE(all.equal(
        as.numeric(x_l[i_x, , drop = FALSE]),                                                            # ### INDEX-SAFE
        as.numeric(x_l[i_x - 1L, , drop = FALSE])
      ))
      if (!is.null(previous.CIFs) && same_row) {
        list.CIFs[[i_x]] <- previous.CIFs
        next
      }
    }

    if (!is.null(initial.CIFs)) {
      ip <- as.numeric(initial.CIFs[i_x, , drop = FALSE])
      if (length(ip) >= 4) ip <- ip[seq_len(4)]                                                          # ### INDEX-SAFE
      log_p0 <- log(clampP(ip, prob.bound))
    }

    eq_fn <- function(lp) {
      estimating_equation_CIFs(
        log_p        = lp,
        alpha_tmp_1  = alpha_tmp_1[i_x],
        beta_tmp_1   = beta_tmp_1,
        alpha_tmp_2  = alpha_tmp_2[i_x],
        beta_tmp_2   = beta_tmp_2,
        estimand     = estimand,
        optim.method = optim.method,
        prob.bound   = prob.bound
      )
    }

    if (optim.method$inner.optim.method %in% c("optim", "BFGS")) {
      sol <- optim(par = log_p0, fn = eq_fn, method = "BFGS",
                   control = list(maxit = optim.method$optim.parameter7,
                                  reltol = optim.method$optim.parameter6))
    } else if (optim.method$inner.optim.method == "SANN") {
      sol <- optim(par = log_p0, fn = eq_fn, method = "SANN",
                   control = list(maxit = optim.method$optim.parameter7,
                                  reltol = optim.method$optim.parameter6))
    } else {
      stop("Unsupported inner.optim.method")
    }

    probs <- clampP(exp(sol$par), prob.bound)
    list.CIFs[[i_x]] <- probs
    previous.CIFs     <- probs
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
calculatePotentialCIFs_LM <- function(
    alpha_beta_tmp, x_a, x_l, offset, epsilon, estimand, optim.method, prob.bound, initial.CIFs = NULL
) {
  stopifnot(requireNamespace("minpack.lm", quietly = TRUE))

  i_parameter <- rep(NA_integer_, 7L)
  i_parameter <- calculateIndexForParameter(i_parameter, x_l, x_a)
  alpha_1   <- alpha_beta_tmp[seq_len(i_parameter[1])]
  beta_tmp_1 <- alpha_beta_tmp[seq.int(i_parameter[2], i_parameter[3])]
  alpha_2   <- alpha_beta_tmp[seq.int(i_parameter[4], i_parameter[5])]
  beta_tmp_2 <- alpha_beta_tmp[seq.int(i_parameter[6], i_parameter[7])]

  alpha_tmp_1 <- as.numeric(x_l %*% matrix(alpha_1, ncol = 1) + offset)
  alpha_tmp_2 <- as.numeric(x_l %*% matrix(alpha_2, ncol = 1) + offset)

  n <- length(epsilon)
  p0 <- c(
    sum(epsilon == estimand$code.event1) / n + prob.bound,
    sum(epsilon == estimand$code.event2) / n + prob.bound,
    sum(epsilon == estimand$code.event1) / n + prob.bound,
    sum(epsilon == estimand$code.event2) / n + prob.bound
  )
  p0     <- clampP(p0, prob.bound)
  log_p0 <- log(p0)

  # 重複行スキップ（同じ x_l 行の連続重複をキャッシュ）
  out_mat <- matrix(NA_real_, nrow = nrow(x_l), ncol = 4L)
  previous_key <- NULL
  previous_sol <- NULL

  ctrl <- minpack.lm::nls.lm.control(
    maxiter = optim.method$optim.parameter7 %||% 100L,
    ftol    = optim.method$optim.parameter6 %||% 1e-10,
    ptol    = 1e-10
  )

  for (i_x in seq_len(nrow(x_l))) {
    key <- paste0(x_l[i_x, ], collapse = "\r")

    if (i_x > 1L && !is.null(previous_sol) && isTRUE(key == previous_key)) {
      out_mat[i_x, ] <- previous_sol
      next
    }

    ip <- if (!is.null(initial.CIFs)) as.numeric(initial.CIFs[i_x, 1:4, drop = FALSE]) else exp(log_p0)
    start <- log(clampP(ip, prob.bound))

    # 観測 i の残差とJac
    res_fun <- function(lp) residuals_CIFs(
      log_p       = lp,
      alpha_tmp_1 = alpha_tmp_1[i_x],
      beta_tmp_1  = beta_tmp_1,
      alpha_tmp_2 = alpha_tmp_2[i_x],
      beta_tmp_2  = beta_tmp_2,
      estimand    = estimand,
      prob.bound  = prob.bound
    )
    jac_fun <- function(lp) jacobian_CIFs(
      log_p       = lp,
      alpha_tmp_1 = alpha_tmp_1[i_x],
      beta_tmp_1  = beta_tmp_1,
      alpha_tmp_2 = alpha_tmp_2[i_x],
      beta_tmp_2  = beta_tmp_2,
      estimand    = estimand,
      prob.bound  = prob.bound
    )

    fit <- minpack.lm::nls.lm(
      par = start,
      fn  = res_fun,
      jac = jac_fun,
      control = ctrl
    )

    sol <- clampP(exp(fit$par), prob.bound)
    out_mat[i_x, ] <- sol
    previous_sol <- sol
    previous_key <- key
  }

  colnames(out_mat) <- c("p10", "p20", "p11", "p21")
  out_mat
}

