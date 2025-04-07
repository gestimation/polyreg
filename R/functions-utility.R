checkDependentPackages <- function() {
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
  if (!is.na(any(time)))
    warning("Invalid time variable. NA values included")
  #  if (any(time<0))
  #    warning("Invalid time variable. Non-negative values included")
  if (length(event) != length(time))
    stop("Time and event variables are different lengths")
  if (missing(event))
    stop("Must have an event argument")
  if (is.numeric(event)) {
    if (!is.na(any(event)))
      warning("Invalid event variable. NA values included")
    status <- event
  } else if (is.logical(event)) {
    if (!is.na(any(event)))
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
  if (!is.na(any(time)))
    warning("Invalid time variable. NA values included")
  #  if (any(time<0))
  #    warning("Invalid time variable. Non-negative values included")
  if (length(event) != length(time))
    stop("Time and event variables are different lengths")
  if (missing(event))
    stop("An event argument is required")
  if (is.numeric(event)) {
    if (!is.na(any(event)))
      warning("Invalid event variable. NA values included")
    status <- event
  } else if (is.logical(event)) {
    if (!is.na(any(event)))
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

readSurv <- function(formula, data, weights, code.event, code.censoring, subset.condition, na.action) {
  data <- createAnalysisDataset(formula, data, weights, subset.condition, na.action)
  cl <- match.call()
  if (missing(formula))
    stop("A formula argument is required")
  mf <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "offset", "cluster")
  out_terms <- terms(formula, special, data = data)
  if (!is.null(attr(out_terms, "specials")$strata))
    stop("strata() cannot appear in formula")
  if (!is.null(attr(out_terms, "specials")$offset))
    stop("offset() cannot appear in formula")
  if (!is.null(attr(out_terms, "specials")$cluster))
    stop("cluster() cannot appear in formula")
  mf$formula <- out_terms
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Y <- model.extract(mf, "response")
  if (!inherits(Y, c("Event", "Surv"))) {
    stop("A 'Surv' or 'Event' object is expected")
  } else {
    t <- Y[, 1]
    if (any(t<0)) {
      stop("Invalid time variable. Expected non-negative values. ")
    }
    if (!all(Y[, 2] %in% c(code.event, code.censoring))) {
      stop("Invalid event codes. Must be 0 or 1 for survival and 0, 1 or 2 for competing risks, with 0 representing censoring, if event codes are not specified. ")
    } else {
      epsilon <- Y[, 2]
      d <- ifelse(Y[, 2] == code.censoring, 0, 1)
      d0 <- ifelse(Y[, 2] == code.censoring, 1, 0)
    }
  }
  if (is.na(all.vars(out_terms)[3])) {
    strata <- rep(1, nrow(data))
    strata_name <- NULL
  } else {
    strata_name <- all.vars(out_terms)[3]
    strata <- as.factor(data[[strata_name]])
  }
  if (is.null(weights)) {
    w <- rep(1, nrow(data))
  } else {
    w <- data[[weights]]
    if (!is.numeric(w))
      stop("Weights must be numeric")
    if (any(!is.finite(w)))
      stop("Weights must be finite")
    if (any(w < 0))
      stop("Weights must be non-negative")
    if (any(is.na(w)))
      stop("Weights contain NA values")
  }
  return(list(t = t, epsilon = epsilon, d = d, d0 = d0, strata = strata, strata_name = strata_name, w=w))
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







