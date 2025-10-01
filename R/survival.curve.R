#' Estimate and plot survival or cumulative incidence curves
#'
#' @param formula Model formula representing outcome and strata.
#' @param data data.frame Input dataset containing survival data.
#' @param weights Column name representing the weights. The weights must be nonnegative and it is strongly recommended that they be strictly positive, since zero weights are ambiguous, compared to use of he subset argument.
#' @param subset Specifies a condition for subsetting the data. Defaults to NULL.
#' @param code.event Specifies the code of event. Defaults to 1. For competing risks, provide two event codes where the first is used for the event of interest.
#' @param code.censoring Specifies the code of censoring. Defaults to 0.
#' @param na.action Specifies a missing-data filter function, applied to the model frame, after any subset argument has been used. Defaults to na.pass.
#' @param conf.int The level for a two-sided confidence interval on the survival probabilities. Defaults to 0.95.
#' @param error Specifies standard error calculation. Supported options depend on `output.type`.
#' @param conf.type Specifies transformation used to construct the confidence interval on the probabilities. Defaults to "arcsine-square root".
#' @param report.survfit.std.err Report standard error of log of survival probabilities. If this is not specified, the SE of survival probabilities is stored in std.err, unlike the original survfit. Defauts to FALSE.
#' @param report.ggsurvfit Draw a survival plot using ggsurvfit. Defaults to TRUE.
#' @param label.strata Labels of strata. Defaults to NULL.
#' @param label.x Labels of x axis. Defaults to "Time".
#' @param label.y Labels of y axis. Defaults depend on `output.type`.
#' @param lims.x Range of x axis. Defaults to NULL.
#' @param lims.y Range of y axis. Defaults to c(0, 1).
#' @param font.family Specifies font family of the plot. Defaults to "sans".
#' @param font.size Specifies font size of the plot. Defaults to 14.
#' @param legend.position Specifies position of the legend of curves. Defaults to "top".
#' @param output.type Character string determining whether the Kaplan-Meier ("SURVIVAL") or Aalen-Johansen ("COMPETING-RISK") estimator is returned.
#'
#' @importFrom ggsurvfit ggsurvfit
#' @returns A survfit object containing the Kaplan-Meier estimator or Aalen-Johansen estimator and related statistics. This object is formatted to conform to the survfit class, so note that surv contains estimates corresponds to survival probability or 1-cumulative incidence. Some methods for the class (e.g. residuals.survfit) are not supported.
#' @export
survival.curve <- function(formula,
                           data,
                           weights = NULL,
                           subset = NULL,
                           code.event = c(1, 2),
                           code.censoring = 0,
                           na.action = na.pass,
                           conf.int = 0.95,
                           error = NULL,
                           conf.type = "arcsine-square root",
                           report.survfit.std.err = FALSE,
                           report.ggsurvfit = TRUE,
                           ggsurvfit.type = "risk",
                           AddCensoringMark = FALSE,
                           AddConfidenceInterval = TRUE,
                           AddRiskTable = TRUE,
                           label.x = "Time",
                           label.y = NULL,
                           label.strata = NULL,
                           lims.x = NULL,
                           lims.y = c(0, 1),
                           font.family = "sans",
                           font.size = 14,
                           legend.position = "top",
                           output.type = c("SURVIVAL", "COMPETING-RISK")) {
  checkDependentPackages()
  output.type <- match.arg(output.type)
  conf.lower <- NULL

  out_readSurv <- readSurv(formula, data, weights, code.event, code.censoring, subset, na.action)

  if (output.type == "SURVIVAL") {
    if (is.null(error)) error <- "greenwood"
    valid_errors <- c("greenwood", "tsiatis", "jackknife")
    if (!error %in% valid_errors) {
      stop("Invalid error method for SURVIVAL output. Supported options are 'greenwood', 'tsiatis', and 'jackknife'.")
    }
    if (is.null(label.y)) label.y <- "Survival probability"

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
  } else {
    # COMPETING-RISK branch
    if (is.null(error)) error <- "delta"
    valid_errors <- c("aalen", "delta", "jackknife")
    if (!error %in% valid_errors) {
      stop("Invalid error method for COMPETING-RISK output. Supported options are 'aalen', 'delta', and 'jackknife'.")
    }
    if (is.null(label.y)) label.y <- "Cumulative incidence probability"

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
    out_ci <- calculateCI(out_aj, conf.int, conf.type, conf.lower)
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
  }

  if (report.ggsurvfit) {
    if (is.null(lims.x)) lims.x <- c(0, max(out_readSurv$t))

    # ヘルパ：ggsurvfit呼び分け（riskタイプ対応）
    call_ggsurvfit <- function(sfobj) {
      if (identical(ggsurvfit.type, "risk")) {
        ggsurvfit(sfobj, type = "risk")
      } else {
        ggsurvfit(sfobj)
      }
    }

    if (conf.type %in% c("none", "n") || length(survfit_object$strata) > 2) {
      # CIなし（もしくは層が多い）ケース：上下限=生存確率を与えておく
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

      # --- オプションの追加 ---
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
      # CIを使用する通常ケース
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

      # --- オプションの追加 ---
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
