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

Surv <- function(time, event) {
  if (missing(time))
    stop("Must have a time argument")
  if (!is.numeric(time))
    stop("Time variable is not numeric")
  if (any(is.na(time)))
    warning("Invalid time variable. NA values included")
  #  if (any(time<0))
  #    warning("Invalid time variable. Non-negative values included")
  if (length(event) != length(time))
    stop("Time and event variables are different lengths")
  if (missing(event))
    stop("Must have an event argument")
  if (is.numeric(event)) {
    if (any(is.na(event)))
      warning("Invalid event variable. NA values included")
    status <- event
  } else if (is.logical(event)) {
    if (any(is.na(event)))
      warning("Invalid event variable. NA values included")
    status <- as.numeric(event)
    warning("Event variable is logical, converted to numearic")
  } else if (is.factor(event)) {
    status <- as.numeric(as.factor(event)) - 1
    warning("Event variable is a factor, converted to numearic")
  } else stop("Invalid status value, must be logical or numeric")
  if (nlevels(as.factor(event)) > 3)
    warning("Event variable should not have more than three levels")
  ss <- cbind(time = time, status = status)
  type <- "right"

  inputAttributes <- list()
  if (!is.null(attributes(time)))
    inputAttributes$time <- attributes(time)
  cname <- dimnames(ss)[[2]]
  if (length(cname) == 0) {
    cname <- c("time", "status")
  }
  dimnames(ss) <- list(NULL, cname)
  attr(ss, "type") <- type
  if (length(inputAttributes) > 0)
    attr(ss, "inputAttributes") <- inputAttributes
  class(ss) <- "Surv"
  ss
}

Event <- function(time, event) {
  if (missing(time))
    stop("A time argument is required")
  if (!is.numeric(time))
    stop("Time variable is not numeric")
  if (any(is.na(time)))
    warning("Invalid time variable. NA values included")
  #  if (any(time<0))
  #    warning("Invalid time variable. Non-negative values included")
  if (length(event) != length(time))
    stop("Time and event variables are different lengths")
  if (missing(event))
    stop("An event argument is required")
  if (is.numeric(event)) {
    if (any(is.na(event)))
      warning("Invalid event variable. NA values included")
    status <- event
  } else if (is.is.logical(event)) {
    if (any(is.na(event)))
      warning("Invalid event variable. NA values included")
    status <- as.numeric(event)
    warning("Event variable is logical, converted to numearic")
  } else if (is.factor(event)) {
    status <- as.numeric(as.factor(event)) - 1
    warning("Event variable is a factor, converted to numearic")
  } else stop("Invalid status value, must be logical or numeric")
  if (nlevels(as.factor(event)) > 3)
    warning("Event variable should not have more than three levels")
  ss <- cbind(time = time, status = status)
  type <- "right"

  inputAttributes <- list()
  if (!is.null(attributes(time)))
    inputAttributes$time <- attributes(time)
  cname <- dimnames(ss)[[2]]
  if (length(cname) == 0) {
    cname <- c("time", "status")
  }
  dimnames(ss) <- list(NULL, cname)
  attr(ss, "type") <- type
  if (length(inputAttributes) > 0)
    attr(ss, "inputAttributes") <- inputAttributes
  class(ss) <- "Event"
  ss
}

read_time.point <- function(formula, data, outcome.type, code.censoring, time.point) {
  if (outcome.type %in% c("COMPETING-RISK","SURVIVAL")) {
    if (is.null(time.point) || !length(time.point)) stop("time.point is required when outcome.type is COMPETING-RISK or SURVIVAL.")
    tp <- suppressWarnings(max(time.point, na.rm = TRUE))
    if (!is.finite(tp) || tp < 0) stop("time.point must be non-negative and finite when outcome.type is COMPETING-RISK or SURVIVAL.")
    return(tp)
  } else if (outcome.type == "BINOMIAL") {
    tp <- Inf
    return(tp)
  } else if (outcome.type %in% c("PROPORTIONAL","POLY-PROPORTIONAL") & is.null(time.point)) {
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
    tp <- data$t[!is.na(data$epsilon) & data$epsilon != code.censoring]
    tp <- sort(unique(tp[is.finite(tp)]))
    return(tp)
  } else {
    return(time.point)
  }
}

createAnalysisDataset <- function(formula, data, other.variables.analyzed=NULL, subset.condition=NULL, na.action=na.pass) {
  if (!is.null(subset.condition)) {
    analysis_dataset <- subset(data, eval(parse(text = subset.condition)))
  } else {
    analysis_dataset <- data
  }
  all_vars <- c(all.vars(formula), other.variables.analyzed)
  analysis_dataset <- analysis_dataset[, all_vars, drop = FALSE]
  return(na.action(analysis_dataset))
}

clampP <- function(p, eps = 1e-5) pmin(pmax(p, eps), 1 - eps)
clampLogP <- function(x, eps = 1e-5) {
  ex <- exp(x)
  x[ex < eps]         <- log(eps)
  x[(1 - ex) < eps]   <- log(1 - eps)
  x
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

create_rr_text <- function(coefficient, cov, index, omit.conf.int=TRUE, conf.int=0.95) {
  alpha <- 1 - conf.int
  critical_value <- qnorm(1 - alpha / 2)
  coef <- coefficient[index]
  coef_se <- sqrt(diag(cov)[index])
  conf_low <- coef - critical_value * coef_se
  conf_high <- coef + critical_value * coef_se
  p_value <- floor(2 * (1 - pnorm(abs(coef) / coef_se)))
  if (omit.conf.int==TRUE) {
    if (p_value<0.01) text <- paste0("RR=", round(exp(coef), digit=2), ", p<0.01")
    else text <- paste0("RR=", round(exp(coef), digit=2), ", p=", p_value)
  } else {
    if (p_value<0.01) text <- paste0("RR=", round(exp(coef), digit=2), " (", round(exp(conf_low), digit=2), " to ", round(exp(conf_high), digit=2), ", p<0.01", ")")
    else text <- paste0("RR=", round(exp(coef), digit=2), " (", round(exp(conf_low), digit=2), " to ", round(exp(conf_high), digit=2), ", p=", p_value, ")")
  }
  return(text)
}

createTestData <- function(n, w, first_zero=FALSE, last_zero=FALSE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE) {
  one <- rep(1, n)
  t <- c(1:(n/2), 1:(n/2))
  epsilon <- rep(1, n)
  epsilon[2] <- 2
  epsilon[3] <- 2
  if (first_zero==TRUE) {
    epsilon[1] <- 0
    epsilon[n/2+1] <- 0
  }
  if (last_zero==TRUE) {
    epsilon[n/2] <- 0
    epsilon[n] <- 0
  }
  w <- rep(w, n)
  if (logical_strata==TRUE) {
    strata <- (t %% 2 == 0)
  } else {
    strata <- as.factor((t %% 2 == 0))
  }
  if (na_strata==TRUE) {
    strata[1] <- NA
  }
  subset <- rep(1, n)
  if (subset_present==TRUE) {
    subset[1] <- 0
  }
  d <- as.numeric(epsilon>0)
  return(data.frame(id = 1:n, t = t, epsilon = epsilon, d = d, w = w, strata = strata, subset=subset))
}


normalizeCovariate <- function(formula, data, should.normalize.covariate, outcome.type) {
  mf <- model.frame(formula, data)
  Y <- model.extract(mf, "response")
  response_term <- formula[[2]]
  if (inherits(mf[[1]], "Surv") | inherits(mf[[1]], "Event")) {
    response_vars <- all.vars(response_term)
    covariate_cols <- setdiff(all.vars(formula), response_vars)  # Remove time and event
  } else {
    covariate_cols <- all.vars(formula)[-1]  # Exclude the response variable
  }
  normalized_data <- data
  range_vector <- 1
  if (should.normalize.covariate == TRUE & length(covariate_cols)>0) {
    for (col in covariate_cols) {
      x <- normalized_data[[col]]
      range <- max(x)-min(x)
      normalized_data[[col]] <- x/range
      range_vector <- cbind(range_vector,range)
    }
    if (outcome.type == 'PROPORTIONAL') {
      range_vector <- cbind(range_vector)
    } else if (outcome.type == 'SURVIVAL' | outcome.type == 'BINOMIAL') {
      range_vector <- cbind(range_vector,1)
    } else {
      range_vector <- cbind(range_vector,1,range_vector,1)
    }
  } else {
    if (outcome.type == "PROPORTIONAL") {
      range_vector <- rep(1, (length(covariate_cols)+1))
    } else if (outcome.type == "SURVIVAL" | outcome.type == "BINOMIAL") {
      range_vector <- rep(1, (length(covariate_cols)+2))
    } else {
      range_vector <- rep(1, (2*length(covariate_cols)+4))
    }
  }
  n_covariate <- length(covariate_cols)
  out <- list(normalized_data=normalized_data, range=range_vector, n_covariate=n_covariate, n=nrow(data))
  return(out)
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

sortByCovariate_old <- function(formula, data, should.sort.data, n_covariate) {
  if (should.sort.data == TRUE & n_covariate>0) {
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

checkSpell_new <- function(outcome.type, effect.measure1, effect.measure2) {
  # 必須パッケージの存在チェック（attachしない）
  required_packages <- c("nleqslv", "boot")
  miss   <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Required packages not installed: ", paste(miss, collapse = ", "))

  # 正規化：大文字化＋空白/ハイフン等の除去
  canon <- function(x) {
    s <-
      if (is.character(x)) x[1L]
    else if (is.factor(x)) as.character(x)[1L]
    else paste(deparse(x), collapse = " ")  # 関数・式・数値なども全部OK
    toupper(gsub("[^A-Z]", "", trimws(s)))
  }

  ot_map <- c(
    "C"="COMPETINGRISK","CR"="COMPETINGRISK","COMPETINGRISK"="COMPETINGRISK","COMPETINGRISKS"="COMPETINGRISK",
    "S"="SURVIVAL","SURVIVAL"="SURVIVAL",
    "B"="BINOMIAL","BINOMIAL"="BINOMIAL",
    "P"="PROPORTIONAL","PROPORTIONAL"="PROPORTIONAL",
    "PP"="POLY-PROPORTIONAL","POLYPROPORTIONAL"="POLY-PROPORTIONAL"
  )
  ot <- ot_map[[canon(outcome.type)]]
  if (is.null(ot)) stop("Invalid outcome.type. Use COMPETINGRISK, SURVIVAL, BINOMIAL, PROPORTIONAL, POLY-PROPORTIONAL.")

  em_map <- c(
    "RR"="RR","RISKRATIO"="RR",
    "OR"="OR","ODDSRATIO"="OR",
    "SHR"="SHR","HR"="SHR",
    "SUBDISTRIBUTIONHAZARDRATIO"="SHR"
  )
  e1 <- em_map[[canon(effect.measure1)]]
  e2 <- em_map[[canon(effect.measure2)]]
  if (is.null(e1)) stop("Invalid effect.measure1. Use: RR, OR, SHR.")
  if (is.null(e2)) stop("Invalid effect.measure2. Use: RR, OR, SHR.")
  list(outcome.type = ot, effect.measure1 = e1, effect.measure2 = e2)
}

checkInput_new <- function(outcome.type, conf.level, report.boot.conf, outer.optim.method, inner.optim.method) {
  if (!is.numeric(conf.level) || length(conf.level) != 1L || conf.level <= 0 || conf.level >= 1)
    stop("conf.level must be a single number in (0, 1).")

  outer_choices <- c("nleqslv","Newton","Broyden","optim","BFGS","SANN","multiroot")
  outer.optim.method <- match.arg(outer.optim.method, choices = outer_choices)

  inner_choices <- c("optim","BFGS","SANN","multiroot")
  inner.optim.method <- match.arg(inner.optim.method, choices = inner_choices)

  if (is.null(report.boot.conf)) {
    report.boot.conf <- outcome.type %in% c("PROPORTIONAL","POLY-PROPORTIONAL")
  } else {
    report.boot.conf <- isTRUE(report.boot.conf)
  }

  list(
    conf.level = conf.level,
    report.boot.conf = report.boot.conf,
    outer.optim.method = outer.optim.method,
    inner.optim.method = inner.optim.method
  )
}

checkInput_old <- function(outcome.type, time.point, conf.level, report.boot.conf, outer.optim.method, inner.optim.method) {
  if (outcome.type == "COMPETING-RISK" | outcome.type == "SURVIVAL") {
    if (length(time.point)>1) {
      time.point <<- max(time.point)
    } else if (is.null(time.point)) {
      stop("Invalid input for time.point when outcome.type is COMPETING-RISK or SURVIVAL.")
    } else if (min(time.point)<0) {
      stop("time.point should be positive when outcome.type is COMPETING-RISK or SURVIVAL.")
    }
  } else if (outcome.type == "BINOMIAL") {
    time.point.corrected <<- Inf
  } else {
    time.point.corrected <<- time.point
  }
  if (conf.level <= 0 | conf.level >= 1)
    stop("Confidence level must be between 0 and 1")
  if (outcome.type == "COMPETINGRISK" & !outer.optim.method %in% c("nleqslv","Newton","Broyden","optim","BFGS","SANN","multiroot","partial")) {
    stop("Invalid input for 'optimization'. Choose 'nleqslv', 'Newton', 'Broyden', 'optim', 'BFGS', 'SANN', 'multiroot' or 'partial'.")
  }
  if (outcome.type == "SURVIVAL" & outer.optim.method == "partial") {
    stop("Invalid input for 'optimization'. Choose 'nleqslv', 'Newton', 'Broyden', 'optim', 'BFGS', 'SANN' or 'multiroot'.")
  }
  if (!inner.optim.method %in% c("optim","BFGS","SANN","multiroot")) {
    stop("Invalid input for 'optimization'. Choose 'optim', 'BFGS', 'SANN' or 'multiroot'.")
  }
  if (is.null(report.boot.conf) & (outcome.type == "PROPORTIONAL" | outcome.type == "POLY-PROPORTIONAL")) {
    report.boot.conf.corrected <<- TRUE
  } else if (is.null(report.boot.conf)) {
    report.boot.conf.corrected <<- FALSE
  } else {
    report.boot.conf.corrected <<- report.boot.conf
  }
}

checkInput2 <- function(data, formula, code.event1, code.event2, code.censoring, outcome.type, conf.level, report.boot.conf, computation.order.method, outer.optim.method, inner.optim.method) {
  cl <- match.call()
  if (missing(formula)) stop("A formula argument is required")
  mf <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "offset", "cluster")
  out_terms <- terms(formula, special, data = data)
  if (!is.null(attr(out_terms, "specials")$strata))
    stop("strata() cannot appear in formula")
  if (!is.null(attr(out_terms, "specials")$offset))
    stop("offset() cannot appear in formula")
  if (!is.null(attr(out_terms, "specials")$cluster))
    stop("cluster() cannot appear in formula")
  if (outcome.type %in% c("COMPETING-RISK","SURVIVAL","POLY-PROPORTIONAL","POLY-PROPORTIONAL")) {
    mf$formula <- out_terms
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    Y <- model.extract(mf, "response")
    if (!inherits(Y, c("Event", "Surv"))) {
      stop("Surv- or Event-object is expected")
    } else {
      t <- Y[, 1]
      if (any(t<0)) stop("Invalid time variable. Expected non-negative values. ")
      if (!all(Y[, 2] %in% c(code.event1, code.event2, code.censoring))) stop("Invalid event codes. Must be 0 or 1 for survival and 0, 1 or 2 for competing risks, with 0 representing censoring, if event codes are not specified. ")
    }
  }
  if (!is.numeric(conf.level) || length(conf.level) != 1L || conf.level <= 0 || conf.level >= 1)
    stop("conf.level must be a single number between 0 and 1")
  if (is.null(report.boot.conf) & (outcome.type == "PROPORTIONAL" | outcome.type == "POLY-PROPORTIONAL")) {
    report.boot.conf.corrected <<- TRUE
  } else if (is.null(report.boot.conf)) {
    report.boot.conf.corrected <<- FALSE
  } else {
    report.boot.conf.corrected <<- report.boot.conf
  }
  order_choices <- c("PARALLEL","SEQUENTIAL")
  match.arg(computation.order.method, choices = order_choices)
  outer_choices <- c("nleqslv","Newton","Broyden","optim","BFGS","SANN","multiroot")
  match.arg(outer.optim.method, choices = outer_choices)
  inner_choices <- c("optim","BFGS","SANN","multiroot")
  match.arg(inner.optim.method, choices = inner_choices)

  #  if (outcome.type == "COMPETING-RISK" & !outer.optim.method %in% c("nleqslv","Newton","Broyden","optim","BFGS","SANN","multiroot","partial")) {
  #   stop("Invalid input for 'optimization'. Choose 'nleqslv', 'Newton', 'Broyden', 'optim', 'BFGS', 'SANN', 'multiroot' or 'partial'.")
  #  }
  #  if (outcome.type == "SURVIVAL" & outer.optim.method == "partial") {
  #    stop("Invalid input for 'optimization'. Choose 'nleqslv', 'Newton', 'Broyden', 'optim', 'BFGS', 'SANN' or 'multiroot'.")
  #  }
  #  if (!inner.optim.method %in% c("optim","BFGS","SANN","multiroot")) {
  #    stop("Invalid input for 'optimization'. Choose 'optim', 'BFGS', 'SANN' or 'multiroot'.")
  #  }
}

checkInput <- function(data, formula, code.event1, code.event2, code.censoring, outcome.type, conf.level, report.boot.conf, computation.order.method, outer.optim.method, inner.optim.method) {
  cl <- match.call()
  if (missing(formula)) stop("A formula argument is required")
  mf <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "offset", "cluster")
  out_terms <- terms(formula, special, data = data)
  if (!is.null(attr(out_terms, "specials")$strata))
    stop("strata() cannot appear in formula")
  if (!is.null(attr(out_terms, "specials")$offset))
    stop("offset() cannot appear in formula")
  if (!is.null(attr(out_terms, "specials")$cluster))
    stop("cluster() cannot appear in formula")
  if (outcome.type %in% c("COMPETING-RISK","SURVIVAL","POLY-PROPORTIONAL","POLY-PROPORTIONAL")) {
    mf$formula <- out_terms
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    Y <- model.extract(mf, "response")
    if (!inherits(Y, c("Event", "Surv"))) {
      stop("Surv- or Event-object is expected")
    } else {
      t <- Y[, 1]
      if (any(t<0)) stop("Invalid time variable. Expected non-negative values. ")
      if (!all(Y[, 2] %in% c(code.event1, code.event2, code.censoring))) stop("Invalid event codes. Must be 0 or 1 for survival and 0, 1 or 2 for competing risks, with 0 representing censoring, if event codes are not specified. ")
    }
  }
  if (!is.numeric(conf.level) || length(conf.level) != 1L || conf.level <= 0 || conf.level >= 1)
    stop("conf.level must be a single number between 0 and 1")
  if (is.null(report.boot.conf) & (outcome.type == 'PROPORTIONAL' | outcome.type == 'POLY-PROPORTIONAL')) {
    report.boot.conf.corrected <- TRUE
  } else if (is.null(report.boot.conf)) {
    report.boot.conf.corrected <- FALSE
  } else {
    report.boot.conf.corrected <- report.boot.conf
  }
  order_choices <- c("PARALLEL","SEQUENTIAL")
  computation.order.method <- match.arg(computation.order.method, choices = order_choices)
  outer_choices <- c("nleqslv","Newton","Broyden","optim","BFGS","SANN","multiroot")
  outer.optim.method <- match.arg(outer.optim.method, choices = outer_choices)
  inner_choices <- c("optim","BFGS","SANN","multiroot")
  inner.optim.method <- match.arg(inner.optim.method, choices = inner_choices)
  return(list(report.boot.conf = report.boot.conf.corrected))
  #  if (outcome.type == "COMPETING-RISK" & !outer.optim.method %in% c("nleqslv","Newton","Broyden","optim","BFGS","SANN","multiroot","partial")) {
  #   stop("Invalid input for 'optimization'. Choose 'nleqslv', 'Newton', 'Broyden', 'optim', 'BFGS', 'SANN', 'multiroot' or 'partial'.")
  #  }
  #  if (outcome.type == "SURVIVAL" & outer.optim.method == "partial") {
  #    stop("Invalid input for 'optimization'. Choose 'nleqslv', 'Newton', 'Broyden', 'optim', 'BFGS', 'SANN' or 'multiroot'.")
  #  }
  #  if (!inner.optim.method %in% c("optim","BFGS","SANN","multiroot")) {
  #    stop("Invalid input for 'optimization'. Choose 'optim', 'BFGS', 'SANN' or 'multiroot'.")
  #  }
}

checkSpell <- function(outcome.type, effect.measure1, effect.measure2) {
  if (outcome.type %in% c("COMPETING-RISK", "COMPETINGRISK", "C", "CR", "COMPETING RISK", "COMPETING-RISKS", "COMPETINGRISKS", "COMPETING RISKS", "Competingrisk", "Competing-risk", "Competing risk", "Competingrisks", "Competing-risks", "Competing risks", "competing-risk", "competingrisk", "competing risk", "competing-risks", "competingrisks", "competing risks")) {
    outcome.type.corrected <- "COMPETING-RISK"
  } else if (outcome.type %in% c("SURVIVAL", "S", "Survival", "Survival")) {
    outcome.type.corrected <- "SURVIVAL"
  } else if (outcome.type %in% c("POLY-PROPORTIONAL", "PP", "Poly-proportional", "poly-proportional")) {
    outcome.type.corrected <- "POLY-PROPORTIONAL"
  } else if (outcome.type %in% c("PROPORTIONAL", "P", "Proportional", "proportional")) {
    outcome.type.corrected <- "PROPORTIONAL"
  } else if (outcome.type %in% c("BINOMIAL", "B", "Binomial", "binomial")) {
    outcome.type.corrected <- "BINOMIAL"
  } else {
    stop("Invalid input for outcome.type, Choose 'COMPETING-RISK', 'SURVIVAL', 'BINOMIAL', 'PROPORTIONAL', or 'POLY-PROPORTIONAL'.")
  }
  if (effect.measure1 %in% c("RR", "rr", "RISK RATIO", "Risk ratio", "risk ratio")) {
    effect.measure1.corrected <- "RR"
  } else if (effect.measure1 %in% c("OR", "or", "ODDS RATIO", "Odds ratio", "odds ratio")) {
    effect.measure1.corrected <- "OR"
  } else if (effect.measure1 %in% c("SHR", "shr", "HR", "hr", "SUBDISTRIBUTION HAZARD RATIO",
                                    "Subdistibution hazard ratio", "subdistibution hazard ratio")) {
    effect.measure1.corrected <- "SHR"
  } else {
    stop("Invalid input for effect.measure1, Choose 'RR', 'OR', or 'SHR'.")
  }
  if (effect.measure2 %in% c("RR", "rr", "RISK RATIO", "Risk ratio", "risk ratio")) {
    effect.measure2.corrected <- "RR"
  } else if (effect.measure2 %in% c("OR", "or", "ODDS RATIO", "Odds ratio", "odds ratio")) {
    effect.measure2.corrected <- "OR"
  } else if (effect.measure2 %in% c("SHR", "shr", "HR", "hr", "SUBDISTRIBUTION HAZARD RATIO",
                                    "Subdistibution hazard ratio", "subdistibution hazard ratio")) {
    effect.measure2.corrected <- "SHR"
  } else {
    stop("Invalid input for effect.measure2, Choose 'RR', 'OR', or 'SHR'.")
  }
  return(list(outcome.type = outcome.type.corrected, effect.measure1 = effect.measure1.corrected, effect.measure2 = effect.measure2.corrected))
}

checkSpell_ <- function(outcome.type, effect.measure1, effect.measure2) {
  if (requireNamespace("nleqslv", quietly = TRUE)
      & requireNamespace("boot", quietly = TRUE)) {
    suppressWarnings(library(nleqslv))
    suppressWarnings(library(boot))
  } else {
    stop("Required packages 'nleqslv' and/or 'boot' are not installed.")
  }
  if (outcome.type %in% c("COMPETING-RISK", "COMPETINGRISK", "C", "CR", "COMPETING RISK", "COMPETING-RISKS", "COMPETINGRISKS", "COMPETING RISKS", "Competingrisk", "Competing-risk", "Competing risk", "Competingrisks", "Competing-risks", "Competing risks", "competing-risk", "competingrisk", "competing risk", "competing-risks", "competingrisks", "competing risks")) {
    outcome.type.corrected <<- "COMPETING-RISK"
  } else if (outcome.type %in% c("SURVIVAL", "S", "Survival", "Survival")) {
    outcome.type.corrected <<- "SURVIVAL"
  } else if (outcome.type %in% c("POLY-PROPORTIONAL", "PP", "Poly-proportional", "poly-proportional")) {
    outcome.type.corrected <<- "POLY-PROPORTIONAL"
  } else if (outcome.type %in% c("PROPORTIONAL", "P", "Proportional", "proportional")) {
    outcome.type.corrected <<- "PROPORTIONAL"
  } else if (outcome.type %in% c("BINOMIAL", "B", "Binomial", "binomial")) {
    outcome.type.corrected <<- "BINOMIAL"
  } else {
    stop("Invalid input for outcome.type, Choose 'COMPETING-RISK', 'SURVIVAL', 'BINOMIAL', 'PROPORTIONAL', or 'POLY-PROPORTIONAL'.")
  }
  if (effect.measure1 %in% c("RR", "rr", "RISK RATIO", "Risk ratio", "risk ratio")) {
    effect.measure1.corrected <<- "RR"
  } else if (effect.measure1 %in% c("OR", "or", "ODDS RATIO", "Odds ratio", "odds ratio")) {
    effect.measure1.corrected <<- "OR"
  } else if (effect.measure1 %in% c("SHR", "shr", "HR", "hr", "SUBDISTRIBUTION HAZARD RATIO",
                                    "Subdistibution hazard ratio", "subdistibution hazard ratio")) {
    effect.measure1.corrected <<- "SHR"
  } else {
    stop("Invalid input for effect.measure1, Choose 'RR', 'OR', or 'SHR'.")
  }
  if (effect.measure2 %in% c("RR", "rr", "RISK RATIO", "Risk ratio", "risk ratio")) {
    effect.measure2.corrected <<- "RR"
  } else if (effect.measure2 %in% c("OR", "or", "ODDS RATIO", "Odds ratio", "odds ratio")) {
    effect.measure2.corrected <<- "OR"
  } else if (effect.measure2 %in% c("SHR", "shr", "HR", "hr", "SUBDISTRIBUTION HAZARD RATIO",
                                    "Subdistibution hazard ratio", "subdistibution hazard ratio")) {
    effect.measure2.corrected <<- "SHR"
  } else {
    stop("Invalid input for effect.measure2, Choose 'RR', 'OR', or 'SHR'.")
  }
}
