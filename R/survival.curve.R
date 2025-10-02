#' @title Survival curves and cumulative incidence curves
#' @description Estimate and plot survival curves based on Kaplan-Meier estimator or cumulative incidence curves based on Aalen-Johansen estimators.
#'
#' @param formula Model formula representing outcome and strata.
#' @param data data A data frame containing the outcome and strata variable if necessary.
#' @param weights Column name representing weights for estimation. The weights must be nonnegative and it is strongly recommended that they be strictly positive, since zero weights are ambiguous, compared to use of he subset argument.
#' @param subset.condition Optional expression (as a character string) defining a subset of \code{data} to analyse. Defaults to \code{NULL}.
#' @param na.action A function specifying the action to take on missing values. The default is \code{\link[stats]{na.omit}}.
#' @param code.event1 Integer code corresponding to the first event of interest. Defaults to \code{1}.
#' @param code.event2 Integer code corresponding to the competing event. Defaults to \code{2}.
#' @param code.censoring Integer code representing censoring. Defaults to \code{0}.
#' @param outcome.type Character string determining whether the Kaplan-Meier (\code{"SURVIVAL"}) or Aalen-Johansen (\code{"COMPETING-RISK"}) estimator is returned.
#' @param conf.int The level for a two-sided confidence interval on the survival probabilities. Defaults to \code{0.95}.
#' @param error Specifies standard error calculation. Supported options depend on `outcome.type`.
#' @param conf.type Specifies transformation used to construct the confidence interval on the probabilities. Defaults to \code{"arcsine-square root"}.
#' @param report.survfit.std.err Report standard error of log of survival probabilities. If this is not specified, the SE of survival probabilities is stored in std.err, unlike the original survfit. Defauts to FALSE.
#' @param report.ggsurvfit Draw a survival plot using ggsurvfit. Defaults to \code{TRUE}.
#' @param label.strata Labels of strata. Defaults to \code{NULL}.
#' @param label.x Labels of x axis. Defaults to \code{"Time"}.
#' @param label.y Labels of y axis. Defaults to \code{"Survival probability"}.
#' @param lims.x Range of x axis. Defaults to \code{NULL}.
#' @param lims.y Range of y axis. Defaults to \code{c(0, 1)}.
#' @param font.family Specifies font family of the plot. Defaults to \code{"sans"}.
#' @param font.size Specifies font size of the plot. Defaults to \code{14}.
#' @param legend.position Specifies position of the legend of curves. Defaults to \code{"top"}.
#'
#' @importFrom ggsurvfit ggsurvfit
#' @returns A survfit object containing the Kaplan-Meier estimator or Aalen-Johansen estimator and related statistics. This object is formatted to conform to the survfit class, so note that surv contains estimates corresponds to survival probability or 1-cumulative incidence. Some methods for the class (e.g. residuals.survfit) are not supported.
#' @export
survival.curve <- function(formula,
                           data,
                           weights = NULL,
                           subset.condition = NULL,
                           na.action = na.omit,
                           code.event1 = 1,
                           code.event2 = 2,
                           code.censoring = 0,
                           outcome.type = "SURVIVAL",
                           conf.int = 0.95,
                           error = "greenwood",
                           conf.type = "arcsine-square root",
                           report.survfit.std.err = FALSE,
                           report.ggsurvfit = TRUE,
                           ggsurvfit.type = NULL,
                           AddCensoringMark = FALSE,
                           AddConfidenceInterval = TRUE,
                           AddRiskTable = TRUE,
                           label.x = "Time",
                           label.y = "Survival probability",
                           label.strata = NULL,
                           lims.x = NULL,
                           lims.y = c(0, 1),
                           font.family = "sans",
                           font.size = 14,
                           legend.position = "top") {
  checkDependentPackages()
  outcome.type <- check_outcome.type(outcome.type)
  out_readSurv <- readSurv(formula, data, weights, code.event1, code.event2, code.censoring, subset.condition, na.action)
  error <- check_error(error, outcome.type)

  if (outcome.type == "SURVIVAL") {
    out_km <- calculateKM_rcpp(out_readSurv$t, out_readSurv$d, out_readSurv$w, as.integer(out_readSurv$strata), error)
    if (error == "jackknife") {
      fn <- function(data) {
        calculateKM_rcpp(data$t, data$d, data$w, as.integer(data$strata), "none")
      }
      out_km$std.err <- calculateJackKnifeSE(out_readSurv, fn)
    } else {
      out_km$std.err <- out_km$surv * out_km$std.err
    }
    out_ci <- calculateCI(out_km, conf.int, conf.type, conf.lower=NULL)

    if (!all(as.integer(out_readSurv$strata) == 1) & is.null(label.strata)) {
      names(out_km$strata) <- levels(as.factor(out_readSurv$strata))
    } else if (!all(as.integer(out_readSurv$strata) == 1)) {
      names(out_km$strata) <- label.strata
    }
    if (report.survfit.std.err) {
      out_km$std.err <- out_km$std.err / out_km$surv
    }

    survfit_object <- list(
      time = out_km$time,
      surv = out_km$surv,
      n = out_km$n,
      n.risk = out_km$n.risk,
      n.event = out_km$n.event,
      n.censor = out_km$n.censor,
      std.err = out_km$std.err,
      upper = if (is.null(conf.type) | conf.type == "none" | conf.type == "n") NULL else out_ci$upper,
      lower = if (is.null(conf.type) | conf.type == "none" | conf.type == "n") NULL else out_ci$lower,
      conf.type = conf.type,
      call = match.call(),
      type = "kaplan-meier",
      method = "Kaplan-Meier"
    )
    if (any(as.integer(out_readSurv$strata) != 1)) {
      survfit_object$strata <- out_km$strata
    }
    class(survfit_object) <- c("survfit")
  } else {
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
    out_aj$std.err <- calculateAalenDeltaSE(out_aj$time1, out_aj$aj1, out_aj$n.event1, out_aj$n.event2, n.risk, out_aj$time0, out_aj$km0, out_aj$strata1, error)
    out_aj$surv <- 1 - out_aj$aj1
    out_ci <- calculateCI(out_aj, conf.int, conf.type, conf.lower=NULL)
    if (report.survfit.std.err) {
      out_aj$std.err <- out_aj$std.err / out_aj$surv
    }

    survfit_object <- list(
      time = out_aj$time1,
      surv = out_aj$surv,
      n = n,
      n.risk = n.risk,
      n.event = out_aj$n.event1,
      n.censor = out_aj$n.censor,
      std.err = out_aj$std.err,
      upper = if (is.null(conf.type) | conf.type == "none" | conf.type == "n") NULL else out_ci$upper,
      lower = if (is.null(conf.type) | conf.type == "none" | conf.type == "n") NULL else out_ci$lower,
      conf.type = conf.type,
      call = match.call(),
      type = "Aalen-Johansen",
      method = "Aalen-Johansen"
    )
    if (any(as.integer(out_readSurv$strata) != 1)) {
      survfit_object$strata <- out_aj$strata1
    }
    class(survfit_object) <- c("survfit")
  }

  if (report.ggsurvfit) {
    if (is.null(lims.x)) lims.x <- c(0, max(out_readSurv$t))

    call_ggsurvfit <- function(sfobj) {
      if (identical(ggsurvfit.type, "risk")) {
        ggsurvfit(sfobj, type = "risk")
      } else {
        ggsurvfit(sfobj)
      }
    }

    if (conf.type %in% c("none", "n") || length(survfit_object$strata) > 2) {
      survfit_object_ <- survfit_object
      survfit_object_$lower <- survfit_object_$surv
      survfit_object_$upper <- survfit_object_$surv

      out_ggsurvfit <- call_ggsurvfit(survfit_object_) +
        theme_classic() +
        theme(
          legend.position = legend.position,
          axis.title = element_text(size = (font.size + 2), family = font.family),
          axis.text  = element_text(size = font.size, family = font.family),
          legend.text = element_text(size = font.size, family = font.family)
        ) +
        labs(x = label.x, y = label.y) +
        lims(x = lims.x, y = lims.y) +
        theme(legend.position = "top")

      if (isTRUE(AddConfidenceInterval)) {
        out_ggsurvfit <- out_ggsurvfit + add_confidence_interval()
      }
      if (isTRUE(AddRiskTable)) {
        out_ggsurvfit <- out_ggsurvfit + add_risktable(risktable_stats = c("n.risk"))
      }
      if (isTRUE(AddCensoringMark)) {
        out_ggsurvfit <- out_ggsurvfit + add_censor_mark()
      }

    } else {
      if (length(survfit_object$time) != length(survfit_object$lower))
        stop("time and upper/lower required for ggsurvfit are different lengths")
      if (length(survfit_object$time) != length(survfit_object$n.risk))
        stop("time and n.risk used required ggsurvfit are different lengths")

      out_ggsurvfit <- call_ggsurvfit(survfit_object) +
        theme_classic() +
        theme(
          legend.position = legend.position,
          axis.title = element_text(size = (font.size + 2), family = font.family),
          axis.text  = element_text(size = font.size, family = font.family),
          legend.text = element_text(size = font.size, family = font.family)
        ) +
        labs(x = label.x, y = label.y) +
        lims(x = lims.x, y = lims.y)

      if (isTRUE(AddConfidenceInterval)) {
        out_ggsurvfit <- out_ggsurvfit + add_confidence_interval()
      }
      if (isTRUE(AddRiskTable)) {
        out_ggsurvfit <- out_ggsurvfit + add_risktable(risktable_stats = c("n.risk"))
      }
      if (isTRUE(AddCensoringMark)) {
        out_ggsurvfit <- out_ggsurvfit + add_censor_mark()
      }
    }
    print(out_ggsurvfit)
  }

  return(survfit_object)
}
