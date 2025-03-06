calculateIPCW <- function(formula, data, code.censoring, strata_name, specific.time) {
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
