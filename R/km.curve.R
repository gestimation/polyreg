#' #' Title Estimate and plot Kaplan-Meier curves
#'
#' @param formula formula Model formula representing outcome and strata
#' @param data data.frame Input dataset containing survival data.
#' @param weights character Column name representing the weights. The weights must be nonnegative and it is strongly recommended that they be strictly positive, since zero weights are ambiguous, compared to use of the subset argument.
#' @param subset character Specifies a condition for subsetting the data. Defaults to NULL.
#' @param code.event integer Specifies the code of event. Defaults to 1.
#' @param code.censoring integer Specifies the code of censoring. Defaults to 0.
#' @param na.action character Specifies a missing-data filter function, applied to the model frame, after any subset argument has been used. Defaults to na.pass.
#' @param conf.int numeric The level for a two-sided confidence interval on the survival probabilities. Defaults to 0.95.
#' @param error character Specifies standard error calculation. "greenwood" for the Greenwood formula, "tsiatis" for the Tsiatis formula or "jackknife" for the jack knife method. Defaults to "greenwood".
#' @param conf.type character Specifies transformation used to construct the confidence interval on the probabilities. Defaults to "arcsine-square root".
#' @param report.survfit.std.err logical Report standard error of log of survival probabilities. If this is not specified, the SE of survival probabilities is stored in std.err, unlike the original survfitThis is  according to the original survfit. Defaults to FALSE.
#' @param report.ggsurvfit logical Draw a survival plot using ggsurvfit. Defaults to TRUE.
#' @param label.strata character Labels of strata. Defaults to NULL.
#' @param label.x character Labels of x axis. Defaults to "Survival probability".
#' @param label.y character Labels of y axis. Defaults to "Time".
#' @param lims.x vector Range of x axis. Defaults to NULL.
#' @param lims.y vector Range of y axis. Defaults to c(0, 1).
#' @param font.family character Specifies font family of the plot. Defaults to "sans".
#' @param font.size numeric Specifies font size of the plot. Defaults to 14.
#' @param legend.position character Specifies position of the legend of curves. Defaults to "top".
#' @importFrom ggsurvfit ggsurvfit
#' @importFrom Rcpp sourceCpp
#' @useDynLib km.curve, .registration = TRUE
#' @returns An object consists of Kaplan-Meier estimator and related statistics. This object is formatted to conform to the survfit class. Some methods for the class (e.g. residuals.survfit) are not supported.
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
#' prostate <- prostate %>% mutate(d=ifelse(status=="alive",0,1))
#' prostate <- prostate %>% mutate(a=ifelse(rx=="placebo","Placebo","Experimental"))
#' prostate$t <- prostate$dtime/12
#' attr(prostate$a, "label") <- "Treatment"
#' survfit_by_group <- km.curve(Event(t, d)~a, data=prostate, label.x = "Years from randomization")
#' quantile(survfit_by_group)
#' print(survfit_by_group, rmean=6)
#'
#' Surv <- km.curve::Surv
#' survfit_overall <- km.curve(Surv(t, d)~1, data=prostate, report.ggsurvfit=FALSE)
#'
#' survfit_list <- list(survfit_overall, survfit_by_group)
#' table_from_survfit <- tbl_survfit(survfit_list, times = c(2, 4, 6), label_header = "**{time} years**") |>
#'   modify_spanning_header(all_stat_cols() ~ "**Overall survival**")
#' table_from_survfit
km.curve <- function(formula,
                     data,
                     weights = NULL,
                     subset = NULL,
                     code.event = c(1, 2),
                     code.censoring = 0,
                     na.action = na.pass,
                     conf.int = 0.95,
                     error = "greenwood",
                     conf.type = "arcsine-square root",
                     report.survfit.std.err = FALSE,
                     report.ggsurvfit = TRUE,
                     label.x = "Time",
                     label.y = "Survival probability",
                     label.strata = NULL,
                     lims.x = NULL,
                     lims.y = c(0, 1),
                     font.family = "sans",
                     font.size = 14,
                     legend.position = "top"
) {
  checkDependentPackages()
  out_readSurv <- readSurv(formula, data, weights, code.event, code.censoring, subset, na.action)
  out_km <- calculateKM_rcpp(out_readSurv$t, out_readSurv$d, out_readSurv$w, as.integer(out_readSurv$strata), error)
  if (error=="jackknife") {
    fn <- function(data) {
      out_km <- calculateKM_rcpp(data$t, data$d, data$w, as.integer(data$strata), "none")
      return(out_km)
    }
    out_km$std.err <- calculateJackKnifeSE(out_readSurv, fn)
  } else {
    out_km$std.err <- out_km$surv*out_km$std.err
  }
  out_ci <- calculateCI(out_km, conf.int, conf.type, conf.lower)
  if (!all(as.integer(out_readSurv$strata) == 1) & (is.null(label.strata))) {
    names(out_km$strata) <- levels(as.factor(out_readSurv$strata))
  } else if (!all(as.integer(out_readSurv$strata) == 1)) {
    names(out_km$strata) <- label.strata
  }
  if (report.survfit.std.err == TRUE) {
    out_km$std.err <- out_km$std.err/out_km$surv
  }

  survfit_object <- list(
    time = out_km$time,
    surv = out_km$surv,
    n = out_km$n,
    n.risk = out_km$n.risk,
    n.event = out_km$n.event,
    n.censor = out_km$n.censor,
    std.err = out_km$std.err,
    upper = out_ci$upper,
    lower = out_ci$lower,
    conf.type = conf.type,
    call = match.call(),
    type = "kaplan-meier",
    method = "Kaplan-Meier"
  )
  if (any(as.integer(out_readSurv$strata) != 1)) {
    survfit_object$strata <- out_km$strata
  }

  class(survfit_object) <- c("survfit")
  if (report.ggsurvfit) {
    if (is.null(lims.x)) lims.x <- c(0, max(out_readSurv$t))
    if (conf.type == "none" | conf.type == "n" | length(survfit_object$strata)>2) {
      survfit_object_ <- survfit_object
      survfit_object_$lower <- survfit_object_$surv
      survfit_object_$upper <- survfit_object_$surv
      out_ggsurvfit <- ggsurvfit(survfit_object_) +
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
        add_risktable(risktable_stats = c("n.risk")) +
        add_censor_mark()
    } else {
      if (length(survfit_object$time) != length(survfit_object$lower)) stop("time and upper/lower required for ggsurvfit are different lengths")
      if (length(survfit_object$time) != length(survfit_object$n.risk)) stop("time and n.risk used required ggsurvfit are different lengths")
      out_ggsurvfit <- ggsurvfit(survfit_object) +
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
        add_risktable(risktable_stats = c("n.risk")) +
        add_censor_mark()
    }
    print(out_ggsurvfit)
  }
  return(survfit_object)
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
      d1 <- ifelse(Y[, 2] == code.event[1], 1, 0)
      d2 <- ifelse(Y[, 2] == code.event[2], 1, 0)
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
  return(list(t = t, epsilon = epsilon, d = d, d0 = d0, d1 = d1, d2 = d2, strata = strata, strata_name = strata_name, w=w))
}

