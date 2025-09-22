calculatePotentialCIFs_parallel20250921 <- function(
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

  alpha_1    <- alpha_beta_tmp[seq_len(i_parameter[1])]
  beta_tmp_1 <- alpha_beta_tmp[seq.int(i_parameter[2], i_parameter[3])]
  alpha_2    <- alpha_beta_tmp[seq.int(i_parameter[4], i_parameter[5])]
  beta_tmp_2 <- alpha_beta_tmp[seq.int(i_parameter[6], i_parameter[7])]

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

calculateIPCW_20250401 <- function(formula, data, code.censoring, strata_name, specific.time) {
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
    t <- Y[, 1]  # time variable
    epsilon <- Y[, 2] # status variable
    if (any(t<0))
      stop("Expected non-negative time variable")
  }

  censoring.model <- createCensoringFormula(formula=formula, code.censoring.updated=code.censoring, strata_name=strata_name)
  resC <- phreg(censoring.model, data)
  if (resC$p > 0) kmt <- FALSE
  kmt <- TRUE
  out_predict1 <- suppressWarnings(stats::predict(resC, newdata = data, type = "survival", times = t, individual.time = TRUE, se = FALSE, km = kmt, tminus = TRUE))
  out_predict2 <- suppressWarnings(stats::predict(resC, newdata = data, type = "survival", times = specific.time, individual.time = TRUE, se = FALSE, km = kmt, tminus = TRUE))
  km1 <- out_predict1$surv
  km2 <- out_predict2$surv[1]
  tmp1 <- ifelse(km1 > 0, 1 / km1, 0)
  tmp2 <- ifelse(km2 > 0, 1 / km2, 0)
  censoring_status <- as.numeric(epsilon==code.censoring)
  tmp3 <- (t <= specific.time) * (censoring_status==0) * tmp1
  tmp4 <- (t > specific.time) * tmp2
  ip.weight <- tmp3 + tmp4
  if (any(is.na(ip.weight)))
    stop("Inverse probability weights contain NA values")
  return(ip.weight)
}

calculatePotentialCIFs20250918 <- function(
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
      #      sol <- optim(par = log_p0, fn = eq_fn, method = "BFGS", control = list(maxit = optim.method$optim.parameter7, reltol = optim.method$optim.parameter6))
      sol <- optim(par = log_p0, fn = eq_fn, method = "BFGS", control = list(maxit = 200, reltol = 1e-8))
    } else if (optim.method$inner.optim.method == "SANN") {
      #      sol <- optim(par = log_p0, fn = eq_fn, method = "SANN", control = list(maxit = optim.method$optim.parameter7, reltol = optim.method$optim.parameter6))
      sol <- optim(par = log_p0, fn = eq_fn, method = "SANN", control = list(maxit = 200, reltol = 1e-8))
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

sortByCovariate <- function(formula, data, optim.method, n_covariate) {
  if (optim.method$computation.order.method == "SEQUENTIAL" & n_covariate>0) {
    terms_obj <- terms(formula)
    covariate_names <- attr(terms_obj, "term.labels")
    missing_vars <- setdiff(covariate_names, names(data))
    if (length(missing_vars) > 0) {
      stop("The following covariates are missing: ", paste(missing_vars, collapse = ", "))
    }
    sorted_data <- data[do.call(order, data[covariate_names]), , drop = FALSE]
    return(sorted_data)
  } else {
    return(data)
  }
}

checkDependentPackages <- function(computation.order.method = c("SEQUENTIAL", "PARALLEL")) {
  computation.order.method <- match.arg(computation.order.method)

  # 1) 依存確認（存在だけチェック）
  required_pkgs <- c("ggsurvfit", "Rcpp", "nleqslv", "boot")
  parallel_pkgs <- c("future", "future.apply")
  pkgs <- if (computation.order.method=="PARALLEL") {
    c(required_pkgs, parallel_pkgs)
  } else {
    required_pkgs
  }

  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Required packages not installed: ", paste(missing_pkgs, collapse = ", "))
  }

  # 2) 並列プラン（必要なときだけ設定）
  if (computation.order.method=="PARALLEL") {
    # すでに multisession でなければ設定
    current_plan <- future::plan("list")[[1]]  # 現行ストラテジ関数
    if (!inherits(current_plan, "multisession")) {
      future::plan(future::multisession)
    }
    # RNG 再現性が必要なら（任意）
    # RNGkind("L'Ecuyer-CMRG")  # 呼び出し側で行うなら削除してOK
  }

  invisible(TRUE)
}

checkDependentPackages_old <- function() {
  if (requireNamespace("ggsurvfit", quietly = TRUE) & requireNamespace("Rcpp", quietly = TRUE)) {
    suppressWarnings(library(ggsurvfit))
    suppressWarnings(library(Rcpp))
  } else {
    stop("Required packages 'ggsurvfit' and/or 'Rcpp' are not installed.")
  }
}

