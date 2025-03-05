calculateIPW <- function(formula, data, specific.time) {
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
    epsilon <- Y[, 2]  # status variable
    if (any(t<0))
      stop("Expected non-negative time variable")
    if (!all(epsilon %in% c(0, 1, 2)))
      stop("Expected only 0, 1 or 2, with 0 representing censoring")
  } else {
    stop("Expected only right censored data")
  }
  if (!nrow(data) == nrow(mf))
    stop("Variables contain NA values")

  resC <- phreg(formula, data)
  cens.strata <- resC$strata[order(resC$ord)]
  cens.nstrata <- resC$nstrata
  if (resC$p > 0) kmt <- FALSE
  kmt <- TRUE
  out_predict1 <- suppressWarnings(predict(resC, newdata = data, type = "survival", times = t, individual.time = TRUE, se = FALSE, km = kmt, tminus = TRUE))
  out_predict2 <- suppressWarnings(predict(resC, newdata = data, type = "survival", times = specific.time, individual.time = TRUE, se = FALSE, km = kmt, tminus = TRUE))
  km1 <- out_predict1$surv
  km2 <- out_predict2$surv[1]
  tmp1 <- ifelse(km1 > 0, 1 / km1, 0)
  tmp2 <- ifelse(km2 > 0, 1 / km2, 0)
  tmp3 <- (t <= specific.time) * (epsilon==0) * tmp1
  tmp4 <- (t > specific.time) * tmp2
  ip.weight <- tmp3 + tmp4
  if (any(is.na(ip.weight)))
    stop("Inverse probability weights contain NA values")
  return(ip.weight)
}

