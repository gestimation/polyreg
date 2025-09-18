calculateIPCW <- function(formula, data, code.censoring, strata_name, specific.time) {
  cl <- match.call()
  mf <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "cluster", "offset")
  Terms <- terms(formula, special, data = data)
  mf$formula <- Terms
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Y <- model.extract(mf, "response")
  if (ncol(Y) == 2) {
    t <- Y[, 1]
    epsilon <- Y[, 2]
  }

  #  censoring.model <- createCensoringFormula(formula=formula, code.censoring.updated=code.censoring, strata_name=strata_name)
  #  resC <- phreg(censoring.model, data)
  #  if (resC$p > 0) kmt <- FALSE
  #  kmt <- TRUE
  #  out_predict1 <- suppressWarnings(stats::predict(resC, newdata = data, type = "survival", times = t, individual.time = TRUE, se = FALSE, km = kmt, tminus = TRUE))
  #  out_predict2 <- suppressWarnings(stats::predict(resC, newdata = data, type = "survival", times = specific.time, individual.time = TRUE, se = FALSE, km = kmt, tminus = TRUE))
  #  km1 <- out_predict1$surv
  #  km2 <- out_predict2$surv[1]

  d <- as.integer(epsilon==0)
  if (is.null(strata_name)) {
    strata <- rep(1,nrow(data))
  } else {
    strata <- data[[strata_name]]
    strata <- as.integer(strata)
  }
  out_km <- calculateKM_rcpp(t=t, d=d, strata=strata, error="none")
  s <- rep(specific.time, length(t))
  km1 <- get_surv(t, out_km$surv, out_km$time, strata, out_km$strata)
  km2 <- get_surv(s, out_km$surv, out_km$time, strata, out_km$strata)

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

calculateIPCWMatrix <- function(formula, data, code.censoring, strata_name, estimand, out_normalizeCovariate) {
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
  i_time <- 0
  ip.weight.matrix <- matrix(NA, nrow=out_normalizeCovariate$n, ncol=length(estimand$time.point))
  for (specific.time in estimand$time.point) {
    i_time <- i_time + 1
    ip.weight.matrix[,i_time] <- calculateIPCW(formula, data, code.censoring, strata_name, specific.time)
  }
  return(ip.weight.matrix)
}

get_surv <- function(predicted.time, estimated.surv, estimated.time, predicted.strata = NULL, estimated.strata = NULL) {
  predicted.surv <- numeric(length(predicted.time))

  if (!is.null(predicted.strata) && length(unique(predicted.strata)) == 1) {
    predicted.strata <- NULL
  }
  if (!is.null(estimated.strata) && length(unique(estimated.strata)) == 1) {
    estimated.strata <- NULL
  }
  if (!is.null(estimated.strata)) {
    strata_start <- c(1, head(cumsum(estimated.strata), -1) + 1)
    strata_end <- cumsum(estimated.strata)
  }

  for (i in seq_along(predicted.time)) {
    t <- predicted.time[i]

    # strataがNULL、または predicted.strata[i] がNAの場合は共通の方法でサバイバル推定
    if (is.null(estimated.strata) || is.null(predicted.strata) || is.na(predicted.strata[i])) {
      time_until_t <- estimated.time[estimated.time < t]
      if (length(time_until_t) > 0) {
        time_index <- which.max(time_until_t)
        predicted.surv[i] <- estimated.surv[time_index]
      } else {
        predicted.surv[i] <- 1
      }
    } else {
      strata_index <- predicted.strata[i]

      # strata_index が範囲外でないことをチェック
      if (!is.na(strata_index) && strata_index > 0 && strata_index <= length(estimated.strata)) {
        strata_size <- estimated.strata[strata_index]

        if (!is.na(strata_size) && strata_size > 0) {
          strata_indices <- strata_start[strata_index]:strata_end[strata_index]
          strata_time <- estimated.time[strata_indices]
          strata_surv <- estimated.surv[strata_indices]

          time_until_t <- strata_time[strata_time < t]

          if (length(time_until_t) > 0) {
            time_index <- which.max(time_until_t)
            predicted.surv[i] <- strata_surv[time_index]
          } else {
            predicted.surv[i] <- 1
          }
        } else {
          predicted.surv[i] <- NA
        }
      } else {
        predicted.surv[i] <- NA
      }
    }
  }

  return(predicted.surv)
}

calculateIPCW_old <- function(formula, data, code.censoring, strata_name, specific.time) {
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

  censoring.model <- createCensoringFormula(formula=formula, code.censoring.updated=code.censoring, strata_name=strata_name)
  resC <- phreg(censoring.model, data)
  if (resC$p > 0) kmt <- FALSE
  kmt <- TRUE
  s <- rep(specific.time, length(t))
  out_predict1 <- suppressWarnings(stats::predict(resC, newdata = data, type = "survival", times = t, individual.time = TRUE, se = FALSE, km = kmt, tminus = TRUE))
  out_predict2 <- suppressWarnings(stats::predict(resC, newdata = data, type = "survival", times = s, individual.time = TRUE, se = FALSE, km = kmt, tminus = TRUE))
  km1 <- out_predict1$surv
  km2 <- out_predict2$surv

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

createCensoringFormula <- function(formula, code.censoring.updated, strata_name = NULL) {
  lhs <- formula[[2]]
  if (is.call(lhs) && length(lhs) >= 3) {
    lhs[[3]] <- call("==", lhs[[3]], code.censoring.updated)
  }
  if (is.null(strata_name)) {
    rhs <- quote(+1)  # Intercept-only model
  } else {
    if (!is.character(strata_name)) {
      strata_name <- deparse(substitute(strata_name))  # Convert non-character input to character
    }
    rhs <- as.call(list(as.symbol("strata"), as.symbol(strata_name)))
  }
  censoring.model <- as.formula(call("~", lhs, rhs))
  return(censoring.model)
}
