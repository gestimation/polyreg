getInitialValues <- function(formula, data, outcome.type, exposure, estimand,
                             specific.time, data.initial.values = NULL, prob.bound = 1e-8) {

  ## ---------- 1) モデル枠の構築と基本チェック ----------
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
      stop("Expected a 'Surv' or 'Event' object for the specified outcome.type.")
    } else {
      t <- rep(0, length(Y))
      epsilon <- Y
      if (!all(epsilon %in% c(estimand$code.event1, estimand$code.censoring))) {
        stop("Invalid event codes (need 0/1).")
      }
    }
  } else {
    t <- as.numeric(Y[, 1])
    if (any(t < 0)) stop("Invalid time variable (negative values).")
    if (any(is.na(t))) stop("Time variable contains NA.")
    if (any(is.na(Y[, 2]))) stop("Event variable contains NA.")
    epsilon <- Y[, 2]
    if (outcome.type == 'SURVIVAL') {
      if (!all(epsilon %in% c(estimand$code.event1, estimand$code.censoring)))
        stop("SURVIVAL requires event codes {censoring,event1}.")
    } else {
      if (!all(epsilon %in% c(estimand$code.event1, estimand$code.event2, estimand$code.censoring)))
        stop("COMPETINGRISK requires event codes {censoring,event1,event2}.")
    }
  }

  ## ---------- 2) offset と x_l（共変量デザイン） ----------
  if (!is.null(attributes(Terms)$specials$offset)) {
    ts <- survival::untangle.specials(Terms, "offset")
    if (length(ts$vars) > 0) {
      Terms2 <- Terms[-ts$terms]             # offset を式から除去
      offset <- mf[[ts$vars]]
      x_l <- model.matrix(Terms2, mf)        # 切片列は含まれる（デフォルト）
    } else {
      offset <- rep(0, nrow(mf))
      x_l <- model.matrix(Terms, mf)
    }
  } else {
    offset <- rep(0, nrow(mf))
    x_l <- model.matrix(Terms, mf)
  }
  if (any(is.na(x_l))) stop("Covariate design (x_l) contains NA.")

  ## ---------- 3) x_a（曝露デザインの多値対応） ----------
  if (any(is.na(data[[exposure]]))) stop("Exposure variable contains NA.")
  a_raw <- data[[exposure]]

  build_x_a <- function(a_vec, ref_code) {
    if (is.factor(a_vec) || is.character(a_vec)) {
      a_fac <- as.factor(a_vec)
      M <- model.matrix(~ a_fac)           # 切片 + ダミー
      Xa <- M[, -1, drop = FALSE]          # 基準水準を除外（K-1 列）
      if (!is.null(ref_code) && ref_code != 0) {
        # ref の扱いを簡易に切り替えたい場合の軽実装
        Xa <- 1 - Xa
      }
      colnames(Xa) <- sub("^a_fac", exposure, colnames(Xa))
      Xa
    } else if (is.numeric(a_vec)) {
      # 数値はそのまま 1 列（連続曝露も想定）
      Xa <- as.matrix(a_vec)
      colnames(Xa) <- exposure
      Xa
    } else {
      stop("Unsupported exposure type. Use factor/character/numeric.")
    }
  }

  x_a <- build_x_a(a_raw, estimand$code.exposure.ref)
  if (ncol(x_a) == 0) stop("Exposure design x_a has zero columns.")

  p_l <- ncol(x_l)
  p_a <- ncol(x_a)

  ## ---------- 4) 既定の初期値（手入力）があればそのまま返す ----------
  if (!is.null(data.initial.values)) {
    expected_len <- if (outcome.type == 'SURVIVAL' ||
                        all(epsilon %in% c(estimand$code.event1, estimand$code.censoring))) {
      p_l + p_a
    } else {
      2 * (p_l + p_a)
    }
    if (length(data.initial.values) != expected_len) {
      stop(sprintf("Invalid initial values length: expected %d, got %d.",
                   expected_len, length(data.initial.values)))
    }
    return(data.initial.values)
  }

  ## ---------- 5) 効果指標→リンク関数 ----------
  link_from_measure <- function(measure) {
    switch(toupper(measure),
           "RR"  = "log",
           "OR"  = "logit",
           "SHR" = "cloglog",
           stop("effect.measure must be RR, OR, or SHR"))
  }

  ## ---------- 6) GLM で初期値を作る小関数 ----------
  fit_init <- function(y01, x_l, x_a, link, offset = NULL, pb = 1e-8) {
    # 分離対策：y を [pb, 1-pb] にスムージング
    y_tilde <- as.numeric(y01)
    y_tilde <- y_tilde * (1 - 2 * pb) + pb

    X <- cbind(x_l, x_a)
    df <- data.frame(y = y_tilde, X)
    fam <- binomial(link = link)

    fit <- suppressWarnings(
      try(glm(y ~ . - 1, data = df,
              family = fam,
              weights = rep(1, nrow(df)),
              offset = offset),
          silent = TRUE)
    )

    if (inherits(fit, "try-error")) {
      co <- rep(0, ncol(X))                 # フォールバック：全部0
      names(co) <- colnames(X)
    } else {
      co <- coef(fit)
      co[is.na(co)] <- 0
    }
    list(alpha = unname(co[seq_len(p_l)]),
         beta  = unname(co[p_l + seq_len(p_a)]))
  }

  ## ---------- 7) 目的変数（時点まで発生）と推定 ----------
  y1 <- as.integer(epsilon == estimand$code.event1 & t <= specific.time)

  if (outcome.type == 'SURVIVAL' ||
      all(epsilon %in% c(estimand$code.event1, estimand$code.censoring))) {
    link1 <- link_from_measure(estimand$effect.measure1)
    est1  <- fit_init(y1, x_l, x_a, link1, offset = offset, pb = prob.bound)
    return(c(est1$alpha, est1$beta))
  } else {
    # 競合リスク：イベント2も別途フィット
    y2    <- as.integer(epsilon == estimand$code.event2 & t <= specific.time)
    link1 <- link_from_measure(estimand$effect.measure1)
    link2 <- link_from_measure(estimand$effect.measure2)

    est1  <- fit_init(y1, x_l, x_a, link1, offset = offset, pb = prob.bound)
    est2  <- fit_init(y2, x_l, x_a, link2, offset = offset, pb = prob.bound)

    return(c(est1$alpha, est1$beta, est2$alpha, est2$beta))
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
