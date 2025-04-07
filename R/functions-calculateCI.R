
calculateCI <- function(survfit_object, conf.int, conf.type, conf.lower) {
  if (conf.int <= 0 | conf.int >= 1)
    stop("Confidence level must be between 0 and 1")
  alpha <- 1 - conf.int
  critical_value <- qnorm(1 - alpha / 2)
  if (is.null(conf.type) | conf.type == "none" | conf.type == "n") {
    lower <- NULL
    upper <- NULL
  } else if (conf.type == "arcsine-square root" | conf.type == "arcsin" | conf.type == "a") {
    se <- survfit_object$surv*survfit_object$std.err/2/sqrt(survfit_object$surv * (1 - survfit_object$surv))
    lower <- sin(pmax(asin(sqrt(survfit_object$surv)) - critical_value*se, 0))^2
    upper <- sin(pmin(asin(sqrt(survfit_object$surv)) + critical_value*se, pi/2))^2
  } else if (conf.type == "plain" | conf.type == "p" | conf.type == "linear") {
    lower <- pmax(survfit_object$surv - critical_value*survfit_object$surv*survfit_object$std.err, 0)
    upper <- pmin(survfit_object$surv + critical_value*survfit_object$surv*survfit_object$std.err, 1)
  } else if (conf.type == "log") {
    se <- survfit_object$std.err
    lower <- survfit_object$surv * exp(-critical_value*se)
    upper <- pmin(survfit_object$surv * exp(critical_value*se), 1)
  } else if (conf.type == "log-log") {
    se <- survfit_object$std.err / log(survfit_object$surv)
    lower <- survfit_object$surv^exp(-critical_value*se)
    upper <- survfit_object$surv^exp(critical_value*se)
  } else if (conf.type == "logit") {
    se <- survfit_object$std.err/(1 - survfit_object$surv)
    lower <- survfit_object$surv / (survfit_object$surv + (1 - survfit_object$surv)*exp(critical_value*se))
    upper <- survfit_object$surv / (survfit_object$surv + (1 - survfit_object$surv)*exp(-critical_value*se))
  }
  lower <- sapply(lower, function(x) ifelse(is.nan(x), NA, x))
  upper <- sapply(upper, function(x) ifelse(is.nan(x), NA, x))
  lower <- sapply(lower, function(x) ifelse(x>=1, 1, x))
  upper <- sapply(upper, function(x) ifelse(x>=1, 1, x))
  lower <- sapply(lower, function(x) ifelse(x<=0, 0, x))
  upper <- sapply(upper, function(x) ifelse(x<=0, 0, x))
  return(list(upper=upper, lower=lower))
}
