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

