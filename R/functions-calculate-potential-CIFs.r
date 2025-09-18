calculatePotentialCIFs_parallel <- function(
    alpha_beta_tmp,
    x_a,
    x_l,
    offset,
    epsilon,
    estimand,
    optim.method,
    prob.bound,
    initial.CIFs = NULL,
    use.parallel = TRUE
) {
  # パラメータ
  i_parameter <- calculateIndexForParameter(NA, x_l, x_a)
  alpha_1     <- alpha_beta_tmp[1:i_parameter[1]]
  beta_tmp_1  <- alpha_beta_tmp[i_parameter[2]:i_parameter[3]]
  alpha_2     <- alpha_beta_tmp[i_parameter[4]:i_parameter[5]]
  beta_tmp_2  <- alpha_beta_tmp[i_parameter[6]:i_parameter[7]]

  # 線形予測子
  alpha_tmp_1 <- x_l %*% alpha_1 + offset
  alpha_tmp_2 <- x_l %*% alpha_2 + offset

  # 初期log(p)
  n <- length(epsilon)
  freq1 <- sum(epsilon == estimand$code.event1) / n + prob.bound
  freq2 <- sum(epsilon == estimand$code.event2) / n + prob.bound
  log_p0 <- log(c(freq1, freq2, freq1, freq2))

  # 重複検出（キー作成）
  x_l_key <- apply(x_l, 1, paste0, collapse = "_")
  unique_keys <- !duplicated(x_l_key)
  key_to_index <- match(x_l_key, x_l_key[unique_keys])  # optionally: fastmatch::fmatch()

  # バッチ構成
  unique_idx <- which(unique_keys)
  batch.size <- if (!is.null(optim.method$computation.order.batch.size)) {
    optim.method$computation.order.batch.size
  } else {
    cores <- parallel::detectCores()
    max(10, floor(length(unique_idx) / (cores * 2))) * 2
  }
  batch_indices <- split(unique_idx, ceiling(seq_along(unique_idx) / batch.size))

  # CIF計算関数
  solve_CIF_batch <- function(idx) {
    matrix(
      data = unlist(lapply(idx, function(i_x) {
        log_p_i <- if (is.null(initial.CIFs)) log_p0 else log(initial.CIFs[i_x, ])
        sol <- optim(
          par     = log_p_i,
          fn      = function(lp) estimating_equation_CIFs(
            log_p       = lp,
            alpha_tmp_1 = alpha_tmp_1[i_x],
            beta_tmp_1  = beta_tmp_1,
            alpha_tmp_2 = alpha_tmp_2[i_x],
            beta_tmp_2  = beta_tmp_2,
            estimand    = estimand,
            prob.bound  = prob.bound
          ),
          method  = "BFGS",
          control = list(
            maxit  = optim.method$optim.parameter7,
            reltol = optim.method$optim.parameter6
          )
        )
        exp(sol$par)
      })),
      ncol = 4,
      byrow = TRUE
    )
  }

  # 並列 or 逐次実行
  result_list <- if (use.parallel) {
    future.apply::future_lapply(batch_indices, solve_CIF_batch, future.seed = TRUE)
  } else {
    lapply(batch_indices, solve_CIF_batch)
  }

  # 結合・復元
  unique_CIFs <- do.call(rbind, result_list)
  CIFs_all <- unique_CIFs[key_to_index, , drop = FALSE]
  return(CIFs_all)
}

calculatePotentialCIFs_sequential <- function(
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

    if (optim.method$inner.optim.method %in% c('optim', 'BFGS')) {
      sol <- optim(par = log_p0, fn = eq_fn, method = "BFGS",
                   control = list(maxit = optim.method$optim.parameter7,
                                  reltol = optim.method$optim.parameter6))
    } else if (optim.method$inner.optim.method == 'SANN') {
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
    if (estimand$effect.measure1 == 'RR') {
      ret[1] <- alpha_tmp_1 - clog_p[1] - clog_p[3] + lp0102
      ret[2] <- beta_tmp_1  - clog_p[3] + clog_p[1]
    } else if (estimand$effect.measure1 == 'OR') {
      ret[1] <- alpha_tmp_1 - clog_p[1] - clog_p[3] + lp0102
      ret[2] <- beta_tmp_1  - clog_p[3] + clog_p[1] +
        log1p(-exp_lp[3]) - log1p(-exp_lp[1])
    } else if (estimand$effect.measure1 == 'SHR') {
      ret[1] <- alpha_tmp_1 - clog_p[1] - clog_p[3] + lp0102
      ret[2] <- exp(beta_tmp_1) - ( log1p(-exp_lp[3]) / log1p(-exp_lp[1]) )
    } else {
      stop("Invalid effect_measure1. Must be RR, OR or SHR.")
    }

    if (estimand$effect.measure2 == 'RR') {
      ret[3] <- alpha_tmp_2 - clog_p[2] - clog_p[4] + lp0102
      ret[4] <- beta_tmp_2  - clog_p[4] + clog_p[2]
    } else if (estimand$effect.measure2 == 'OR') {
      ret[3] <- alpha_tmp_2 - clog_p[2] - clog_p[4] + lp0102
      ret[4] <- beta_tmp_2  - clog_p[4] + clog_p[2] +
        log1p(-exp_lp[4]) - log1p(-exp_lp[2])
    } else if (estimand$effect.measure2 == 'SHR') {
      ret[3] <- alpha_tmp_2 - clog_p[2] - clog_p[4] + lp0102
      ret[4] <- exp(beta_tmp_2) - ( log1p(-exp_lp[4]) / log1p(-exp_lp[2]) )
    } else {
      stop("Invalid effect_measure2. Must be RR, OR or SHR.")
    }
    sum(ret^2)
  }
  objective_function(log_p)
}


estimating_equation_CIFs2 <- function(
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
    clog_p <- rep(NA,4)
    clog_p[1] <- clampLogP(log_p[1],prob.bound)
    clog_p[2] <- clampLogP(log_p[2],prob.bound)
    clog_p[3] <- clampLogP(log_p[3],prob.bound)
    clog_p[4] <- clampLogP(log_p[4],prob.bound)
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

calculatePotentialRisk <- function(alpha_beta, x_a, x_l, offset, estimand) {
  i_parameter <- rep(NA_integer_, 7L)
  i_parameter <- calculateIndexForParameter(i_parameter, x_l, x_a)

  one <- rep(1, nrow(x_l))
  alpha_beta_ <- as.matrix(as.vector(alpha_beta))
  if (ncol(alpha_beta_) == 1L) alpha_beta_ <- t(alpha_beta_)
  phi <- x_l %*% alpha_beta_[, seq_len(i_parameter[1])] + offset
  theta <- one * alpha_beta_[, i_parameter[2]]

  if (estimand$effect.measure1 == 'RR') {
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
  } else if (estimand$effect.measure1 == 'OR') {
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


