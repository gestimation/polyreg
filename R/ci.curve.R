#' #' Title Estimate and plot cumulative incidence curves
#'
#' @param formula formula Model formula representing outcome and strata
#' @param data data.frame Input dataset containing survival data.
#' @param weights character Column name representing the weights. The weights must be nonnegative and it is strongly recommended that they be strictly positive, since zero weights are ambiguous, compared to use of the subset argument.
#' @param subset character Specifies a condition for subsetting the data. Defaults to NULL.
#' @param code.event integer Specifies the code of event. Defaults to 1.
#' @param code.censoring integer Specifies the code of censoring. Defaults to 0.
#' @param na.action character Specifies a missing-data filter function, applied to the model frame, after any subset argument has been used. Defaults to na.pass.
#' @param conf.int numeric The level for a two-sided confidence interval on the survival probabilities. Defaults to 0.95.
#' @param error character Specifies standard error calculation. "aalen" for the Aalen formula, "delta" for the delta method or "jackknife" for the jack knife method. Defaults to "delta".
#' @param conf.type character Specifies transformation used to construct the confidence interval on the probabilities. Defaults to "arcsine-square root".
#' @param report.survfit.std.err logical Report standard error of log of survival probabilities. If this is not specified, the SE of survival probabilities is stored in std.err, unlike the original survfitThis is  according to the original survfit. Defaults to FALSE.
#' @param report.ggsurvfit logical Draw a survival plot using ggsurvfit. Defaults to TRUE.
#' @param label.strata character Labels of strata. Defaults to NULL.
#' @param label.x character Labels of x axis. Defaults to "Cumulative incidence probability".
#' @param label.y character Labels of y axis. Defaults to "Time".
#' @param lims.x vector Range of x axis. Defaults to NULL.
#' @param lims.y vector Range of y axis. Defaults to c(0, 1).
#' @param font.family character Specifies font family of the plot. Defaults to "sans".
#' @param font.size numeric Specifies font size of the plot. Defaults to 14.
#' @param legend.position character Specifies position of the legend of curves. Defaults to "top".
#' @importFrom ggsurvfit ggsurvfit
#' @importFrom Rcpp sourceCpp
#' @useDynLib km.curve, .registration = TRUE
#' @returns An object consists of Aalen-Johansen estimator and related statistics. This object is formatted to conform to the survfit class, so note that surv contains estimates corresponds to 1-cumulative incidence. Some methods for the class (e.g. residuals.survfit) are not supported.
#' @export km.curve
#' @export ci.curve
#' @export Surv
#' @export Event
#' @export createTestData
#' @export calculateKM_rcpp
#'
#' @examples
#' library(km.curve)
#' library(gtsummary)
#' library(dplyr)
#' library(labelled)
#' data(prostate)
#' prostate <- prostate %>% mutate(epsilon=ifelse(status=="alive",0,
#'                                        ifelse(status=="dead - prostatic ca",1,
#'                                        ifelse(status=="dead - other ca",1,
#'                                        ifelse(status=="dead - heart or vascular",2,
#'                                        ifelse(status=="dead - cerebrovascular",2,2))))))
#' prostate <- prostate %>% mutate(a=ifelse(rx=="placebo","Placebo","Experimental"))
#' prostate$t <- prostate$dtime/12
#' attr(prostate$a, "label") <- "Treatment"
#' survfit_by_group <- ci.curve(Event(t, epsilon)~a, data=prostate, code.event = c(1,2), label.y = "Cumulative incidence of cancer-specific death", label.x = "Years from randomization")
#' ci.curve(Event(t, epsilon)~a, data=prostate, code.event = c(2,1), label.y = "Cumulative incidence of non-cancer-specific death", label.x = "Years from randomization", conf.type = "n")
#'
#' survfit_overall <- ci.curve(Event(t, epsilon)~1, data=prostate, code.event = c(1,2), report.ggsurvfit=FALSE)
#' survfit_list <- list(survfit_overall, survfit_by_group)
#' table_from_survfit <- tbl_survfit(survfit_list, times = c(2, 4, 6), label_header = "**{time} years**", type="risk") |>
#'   modify_spanning_header(all_stat_cols() ~ "**Cumulative incidence of cancer-specific death**")
#' table_from_survfit

ci.curve <- function(formula,
                     data,
                     weights = NULL,
                     subset = NULL,
                     code.event = c(1, 2),
                     code.censoring = 0,
                     na.action = na.pass,
                     conf.int = 0.95,
                     error = "delta",
                     conf.type = "arcsine-square root",
                     report.survfit.std.err = FALSE,
                     report.ggsurvfit = TRUE,
                     label.x = "Time",
                     label.y = "Cumulative incidence probability",
                     label.strata = NULL,
                     lims.x = NULL,
                     lims.y = c(0, 1),
                     font.family = "sans",
                     font.size = 14,
                     legend.position = "top"
) {
  checkDependentPackages()
  out_readSurv <- readSurv(formula, data, weights, code.event, code.censoring, subset, na.action)
  out_aj <- calculateAJ(out_readSurv)
  out_aj <- readStrata(out_readSurv, out_aj, label.strata)
  if (any(as.integer(out_readSurv$strata) != 1)) {
    n <- table(as.integer(out_readSurv$strata))
    rep_list <- mapply(rep, n, out_aj$strata1, SIMPLIFY = FALSE)
    n.risk <- do.call(c, rep_list) - out_aj$n.cum.censor - out_aj$n.cum.event1 - out_aj$n.cum.event2
  } else {
    n <- length(out_readSurv$strata)
    n.risk <- n - out_aj$n.cum.censor - out_aj$n.cum.event1 - out_aj$n.cum.event2
  }
  out_aj$std.err <- calculateAalenDeltaSE(out_aj$time1, out_aj$aj1, out_aj$n.event1, out_aj$n.event2, n.risk, out_aj$time0, out_aj$km0, error)
  out_aj$surv <- 1-out_aj$aj1
  out_ci <- calculateCI(out_aj, conf.int, conf.type, conf.lower)
  if (report.survfit.std.err == TRUE) {
    out_aj$std.err <- out_aj$std.err/out_aj$surv
  }

  survfit_object <- list(
    time = out_aj$time1,
    surv = out_aj$surv,
    n = n,
    n.risk = n.risk,
    n.event = out_aj$n.event1,
    n.censor = out_aj$n.censor,
    std.err = out_aj$std.err,
    upper = out_ci$upper,
    lower = out_ci$lower,
    conf.type = conf.type,
    call = match.call(),
    type = "Aalen-Johansen",
    method = "Aalen-Johansen"
  )
  if (any(as.integer(out_readSurv$strata) != 1)) {
    survfit_object$strata <- out_aj$strata1
  }
  class(survfit_object) <- c("survfit")
  if (report.ggsurvfit) {
    if (is.null(lims.x)) lims.x <- c(0, max(out_readSurv$t))
    if (conf.type == "none" | conf.type == "n" | length(survfit_object$strata)>2) {
      survfit_object_ <- survfit_object
      survfit_object_$std.err <- rep(1, length(survfit_object_$surv))
      survfit_object_$lower <- survfit_object_$surv
      survfit_object_$upper <- survfit_object_$surv
      out_ggsurvfit <- ggsurvfit(survfit_object_, type = "risk") +
        theme_classic()+
        theme(legend.position = legend.position,
              axis.title = element_text(size = (font.size+2), family = font.family),
              axis.text = element_text(size = font.size, family = font.family),
              legend.text = element_text(size = font.size, family = font.family))+
        labs(x = label.x,
             y = label.y) +
        lims(x = lims.x,
             y = lims.y) +
        theme(legend.position = "top")+
        add_risktable(risktable_stats = c("n.risk"))
    } else {
      if (length(survfit_object$time) != length(survfit_object$lower)) stop("time and upper/lower required for ggsurvfit are different lengths")
      if (length(survfit_object$time) != length(survfit_object$n.risk)) stop("time and n.risk used required ggsurvfit are different lengths")
      out_ggsurvfit <- ggsurvfit(survfit_object, type = "risk") +
        theme_classic()+
        theme(legend.position = legend.position,
              axis.title = element_text(size = (font.size+2), family = font.family),
              axis.text = element_text(size = font.size, family = font.family),
              legend.text = element_text(size = font.size, family = font.family))+
        labs(x = label.x,
             y = label.y) +
        lims(x = lims.x,
             y = lims.y) +
        add_confidence_interval() +
        add_risktable(risktable_stats = c("n.risk"))
    }
    print(out_ggsurvfit)
  }
  return(survfit_object)
}

calculateAJ <- function(data) {
  out_km0 <- calculateKM_rcpp(data$t, data$d0, data$w, as.integer(data$strata), "none")
  km0 <- get_surv(data$t, out_km0$surv, out_km0$time, as.integer(data$strata), out_km0$strata)
  ip.weight <- (data$d0==0) * ifelse(km0 > 0, 1 / km0, 0)
  d1_ipw <- as.matrix(data$w*data$d1*ip.weight)

  aj1 <- time1 <- n.cum.event1 <- n.cum.event2 <- n.cum.censor <- n.event1 <- n.event2 <- n.censor <-strata1 <- NULL
  for (level in sort(unique(as.integer(data$strata)))) {
    sub_d1_ipw <- d1_ipw[as.integer(data$strata) == level, ]
    sub_t <- data$t[as.integer(data$strata) == level]
    sub_d0 <- data$d0[as.integer(data$strata) == level]
    sub_d1 <- data$d1[as.integer(data$strata) == level]
    sub_d2 <- data$d2[as.integer(data$strata) == level]
    not_atrisk <- outer(sub_t, sub_t, ">=")
    #   not_atrisk <- matrix(as.integer(not_atrisk), nrow = nrow(not_atrisk), ncol = ncol(not_atrisk))
    sub_aj1 <- not_atrisk %*% sub_d1_ipw / length(sub_t)
    sub_n.censor <- not_atrisk %*% as.matrix(sub_d0)
    sub_n.event1 <- not_atrisk %*% as.matrix(sub_d1)
    sub_n.event2 <- not_atrisk %*% as.matrix(sub_d2)
    #    unique_aj1 <- unique(cbind(sub_t, sub_d1, sub_aj1))
    #    selected_aj1 <- unique_aj1[unique_aj1[, 2] == 1, ]
    #    sorted_aj1 <- selected_aj1[order(selected_aj1[, 1]), ]
    #    aj1 <- c(aj1, sorted_aj1[,3])
    #    time1 <- c(time1, sorted_aj1[,1])
    #    strata1 <- c(strata1, length(sorted_aj1[,2]))
    is_unique <- !duplicated(sub_t)
    unique_t <- sub_t[is_unique]
    unique_aj1 <- sub_aj1[is_unique]
    unique_n.censor <- sub_n.censor[is_unique]
    unique_n.event1 <- sub_n.event1[is_unique]
    unique_n.event2 <- sub_n.event2[is_unique]
    sorted_t <- sort(unique_t)
    sorted_aj1 <- unique_aj1[order(unique_t)]
    sorted_n.censor <- unique_n.censor[order(unique_t)]
    sorted_n.event1 <- unique_n.event1[order(unique_t)]
    sorted_n.event2 <- unique_n.event2[order(unique_t)]
    time1 <- c(time1, sorted_t)
    aj1 <- c(aj1, sorted_aj1)
    n.cum.censor <- c(n.cum.censor, 0, sorted_n.censor[-length(sorted_n.censor)])
    n.cum.event1 <- c(n.cum.event1, 0, sorted_n.event1[-length(sorted_n.event1)])
    n.cum.event2 <- c(n.cum.event2, 0, sorted_n.event2[-length(sorted_n.event2)])
    strata1 <- c(strata1, length(sorted_t))
  }
  n.censor[1] <- n.cum.censor[1]
  n.event1[1] <- n.cum.event1[1]
  n.event2[1] <- n.cum.event2[1]
  for (i in 2:length(n.cum.event1)) {
    n.censor[i] <- n.cum.censor[i]-n.cum.censor[i-1]
    n.event1[i] <- n.cum.event1[i]-n.cum.event1[i-1]
    n.event2[i] <- n.cum.event2[i]-n.cum.event2[i-1]
  }
  return(list(time1=time1, aj1=aj1, n.event1=n.event1, n.event2=n.event2, n.censor=n.censor, n.cum.event1=n.cum.event1, n.cum.event2=n.cum.event2, n.cum.censor=n.cum.censor, strata1=strata1, time0=out_km0$time, km0=out_km0$surv))
}

readStrata <- function(out_readSurv, out_aj, label.strata=NULL) {
  if (!all(as.integer(out_readSurv$strata) == 1) & (is.null(label.strata))) {
    names(out_aj$strata1) <- levels(as.factor(out_readSurv$strata))
  } else if (!all(as.integer(out_readSurv$strata) == 1)) {
    names(out_aj$strata1) <- label.strata
  }
  return(out_aj)
}
