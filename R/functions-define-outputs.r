getCoef <- function(
    index,
    alpha_beta_estimated,
    cov_estimated,
    report.boot.conf,
    out_bootstrap,
    conf.level) {
  alpha <- 1 - conf.level
  critical_value <- qnorm(1 - alpha / 2)
  if (report.boot.conf == TRUE) {
    coef <- out_bootstrap$boot.coef[index]
    coef_se <- out_bootstrap$boot.coef_se[index]
    conf_low <- out_bootstrap$boot.conf_low[index]
    conf_high <- out_bootstrap$boot.conf_high[index]
    p_value <- out_bootstrap$boot.p_value[index]
  } else {
    coef <- alpha_beta_estimated[index]
    coef_se <- sqrt(diag(cov_estimated)[index])
    conf_low <- coef - critical_value * coef_se
    conf_high <- coef + critical_value * coef_se
    p_value <- 2 * (1 - pnorm(abs(coef) / coef_se))
  }
  list(coef = coef, coef_se = coef_se, conf_low = conf_low, conf_high = conf_high, p_value = p_value)
}

getMessage <- function(outer.optim.method, sol) {
  if (outer.optim.method %in% c("nleqslv", "Newton", "Broyden")) return(sol$message)
  if (outer.optim.method %in% c("optim", "BFGS", "SANN")) return(ifelse(!is.null(sol$message), sol$message, "None"))
  return("-")
}

describeEvents <- function(events, observations, exposed) {
  paste(events, "events in", observations, if (exposed) "exposed observations" else "unexposed observations")
}

tidy_df <- function(coef, text) {
  data.frame(
    term = text,
    estimate = coef$coef,
    std.error = coef$coef_se,
    p.value = coef$p_value,
    conf.low = coef$conf_low,
    conf.high = coef$conf_high
  )
}

glance_df <- function(
    report.optim.convergence, effect_text, events, exposed_events, unexposed_events, median_followup, param_len, iteration, max.estimate.difference, max.function.value, n.solver.iteration, message) {
  if (report.optim.convergence == FALSE) {
    data.frame(
      effect.measure = effect_text,
      n.events = events,
      n.events.exposed = exposed_events,
      n.events.unexposed = unexposed_events,
      median.follow.up = median_followup,
      n.parameters = param_len,
      optimization.message = message
    )
  } else {
    data.frame(
      effect.measure = effect_text,
      n.events = events,
      n.events.exposed = exposed_events,
      n.events.unexposed = unexposed_events,
      median.follow.up = median_followup,
      n.parameters = param_len,
      n.loop.iteration = iteration,
      max.estimate.difference = max.estimate.difference,
      max.function.value = max.function.value,
      n.solver.iteration = n.solver.iteration,
      optimization.message = message
    )
  }
}

reportEffects <- function(outcome.type,
                          report.nuisance.parameter,
                          report.optim.convergence,
                          report.boot.conf,
                          nuisance.model,
                          exposure,
                          estimand,
                          alpha_beta_estimated,
                          cov_estimated,
                          out_bootstrap,
                          out_getResults,
                          iteration,
                          param_diff,
                          sol,
                          conf.level,
                          outer.optim.method) {

  i_parameter <- out_getResults$i_parameter

  x_a_mat <- as.matrix(out_getResults$x_a)
  lev <- attr(out_getResults$x_a, "levels")
  ref <- attr(out_getResults$x_a, "reference")
  if (is.null(lev) || is.null(ref)) {
    nonref_lev <- if (ncol(x_a_mat) > 0L) paste0("lvl", seq_len(ncol(x_a_mat))) else character(0L)
    ref <- estimand$code.exposure.ref %||% "ref"
  } else {
    nonref_lev <- setdiff(lev, ref)
  }

  exposure_terms <- if (length(nonref_lev)) {
    paste0(exposure, " [", nonref_lev, " vs ", ref, "]")
  } else character(0L)

  nl_terms <- attr(stats::terms(nuisance.model), "term.labels")
  nuisance_terms <- c("Intercept", nl_terms)

  idx_beta1 <- if (length(nonref_lev)) seq.int(i_parameter[2], i_parameter[3]) else 0L
  idx_beta2 <- if (length(nonref_lev) && outcome.type == "COMPETING-RISK") seq.int(i_parameter[6], i_parameter[7]) else 0L

  if (isTRUE(report.nuisance.parameter)) {
    idx1 <- c(seq_len(i_parameter[1]), idx_beta1)
    terms_text1 <- c(nuisance_terms, exposure_terms)
  } else {
    idx1 <- idx_beta1
    terms_text1 <- exposure_terms
  }

  coef1 <- if (length(idx1)) {
    getCoef(idx1, alpha_beta_estimated, cov_estimated, report.boot.conf, out_bootstrap, conf.level)
  } else {
    list(coef = numeric(0), coef_se = numeric(0), conf_low = numeric(0), conf_high = numeric(0), p_value = numeric(0))
  }

  exposed_any <- as.integer(rowSums(x_a_mat) > 0L)
  y1 <- out_getResults$y_1
  y2 <- out_getResults$y_2

  exposed_events1   <- sum(y1 * exposed_any)
  unexposed_events1 <- sum(y1 * (1 - exposed_any))

  effect1_text <- paste(estimand$effect.measure1, "of", exposure, "at", estimand$time.point,
                        "( ref =", ref, ")")

  if (outcome.type == "COMPETING-RISK") {
    coef2 <- if (length(idx_beta2) || isTRUE(report.nuisance.parameter)) {
      idx2 <- if (isTRUE(report.nuisance.parameter)) c(seq.int(i_parameter[4], i_parameter[5]), idx_beta2) else idx_beta2
      terms_text2 <- if (isTRUE(report.nuisance.parameter)) c(nuisance_terms, exposure_terms) else exposure_terms
      getCoef(idx2, alpha_beta_estimated, cov_estimated, report.boot.conf, out_bootstrap, conf.level)
    } else {
      list(coef = numeric(0), coef_se = numeric(0), conf_low = numeric(0), conf_high = numeric(0), p_value = numeric(0))
    }

    exposed_events2   <- sum(y2 * exposed_any)
    unexposed_events2 <- sum(y2 * (1 - exposed_any))

    text_values <- list(
      terms_text1 = terms_text1,
      terms_text2 = if (isTRUE(report.nuisance.parameter)) c(nuisance_terms, exposure_terms) else exposure_terms,
      effect1_text = effect1_text,
      event1_summary = sum(y1),
      event1_summary_exposed = describeEvents(exposed_events1, sum(exposed_any), TRUE),
      event1_summary_unexposed = describeEvents(unexposed_events1, length(y1) - sum(exposed_any), FALSE),
      effect2_text = paste(estimand$effect.measure2, "of", exposure, "at", estimand$time.point,
                           "( ref =", ref, ")"),
      event2_summary = sum(y2),
      event2_summary_exposed = describeEvents(exposed_events2, sum(exposed_any), TRUE),
      event2_summary_unexposed = describeEvents(unexposed_events2, length(y2) - sum(exposed_any), FALSE),
      median_followup = paste(round(median(out_getResults$t), 2),
                              "[", round(min(out_getResults$t), 2), ",", round(max(out_getResults$t), 2), "]")
    )

    param_len1 <- length(idx1)
    param_len2 <- if (isTRUE(report.nuisance.parameter)) (i_parameter[7] - i_parameter[4] + 1L) else length(idx_beta2)

    tg <- list(
      event1 = list(
        tidy = tidy_df(coef1, text_values$terms_text1),
        glance = glance_df(report.optim.convergence,
                           text_values$effect1_text,
                           text_values$event1_summary,
                           text_values$event1_summary_exposed,
                           text_values$event1_summary_unexposed,
                           text_values$median_followup,
                           param_len1,
                           iteration,
                           param_diff,
                           max(sol$fvec),
                           sol$iter,
                           getMessage(outer.optim.method, sol))
      ),
      event2 = list(
        tidy = tidy_df(coef2, text_values$terms_text2),
        glance = glance_df(report.optim.convergence,
                           text_values$effect2_text,
                           text_values$event2_summary,
                           text_values$event2_summary_exposed,
                           text_values$event2_summary_unexposed,
                           text_values$median_followup,
                           param_len2,
                           iteration,
                           param_diff,
                           max(sol$fvec),
                           sol$iter,
                           "-")
      )
    )
    class(tg$event1) <- class(tg$event2) <- "modelsummary_list"
    return(tg)

  } else {
    text_values <- list(
      terms_text1 = terms_text1,
      effect1_text = effect1_text,
      event1_summary = sum(y1),
      event1_summary_exposed = describeEvents(exposed_events1, sum(exposed_any), TRUE),
      event1_summary_unexposed = describeEvents(unexposed_events1, length(y1) - sum(exposed_any), FALSE),
      median_followup = paste(round(median(out_getResults$t), 2),
                              "[", round(min(out_getResults$t), 2), ",", round(max(out_getResults$t), 2), "]")
    )

    param_len1 <- length(idx1)

    tg1 <- list(
      tidy = tidy_df(coef1, text_values$terms_text1),
      glance = glance_df(report.optim.convergence,
                         text_values$effect1_text,
                         text_values$event1_summary,
                         text_values$event1_summary_exposed,
                         text_values$event1_summary_unexposed,
                         text_values$median_followup,
                         param_len1,
                         iteration,
                         param_diff,
                         max(sol$fvec),
                         sol$iter,
                         getMessage(outer.optim.method, sol))
    )
    class(tg1) <- "modelsummary_list"
    return(list("event 1 (no competing risk)" = tg1))
  }
}

reportEffects_ <- function(outcome.type,
                          report.nuisance.parameter,
                          report.optim.convergence,
                          report.boot.conf,
                          nuisance.model,
                          exposure,
                          estimand,
                          alpha_beta_estimated,
                          cov_estimated,
                          out_bootstrap,
                          out_getResults,
                          iteration,
                          param_diff,
                          sol,
                          conf.level,
                          outer.optim.method) {

  i_parameter <- out_getResults$i_parameter
  if (report.nuisance.parameter==TRUE) {
    coef1 <- getCoef(1:i_parameter[3], alpha_beta_estimated, cov_estimated, report.boot.conf, out_bootstrap, conf.level)
    terms_text <- c("Intercept", attr(terms(nuisance.model), "term.labels"), paste(exposure, "( ref =", estimand$code.exposure.ref, ")"))
  } else {
    coef1 <- getCoef(i_parameter[3], alpha_beta_estimated, cov_estimated, report.boot.conf, out_bootstrap, conf.level)
    terms_text <- paste(exposure, "( ref =", estimand$code.exposure.ref, ")")
  }
  exposed_events1 <- sum(out_getResults$y_1 * out_getResults$x_a)
  unexposed_events1 <- sum(out_getResults$y_1 * (1 - out_getResults$x_a))

  if (outcome.type == "COMPETING-RISK") {
    if (report.nuisance.parameter==TRUE) {
      coef2 <- getCoef(i_parameter[4]:i_parameter[7], alpha_beta_estimated, cov_estimated, report.boot.conf, out_bootstrap, conf.level)
    } else {
      coef2 <- getCoef(i_parameter[7], alpha_beta_estimated, cov_estimated, report.boot.conf, out_bootstrap, conf.level)
    }
    exposed_events2 <- sum(out_getResults$y_2 * out_getResults$x_a)
    unexposed_events2 <- sum(out_getResults$y_2 * (1 - out_getResults$x_a))
    text_values <- list(
      terms_text               = terms_text,
      effect1_text             = paste(estimand$effect.measure1, "of", exposure, "at", estimand$time.point, "( ref =", estimand$code.exposure.ref, ")"),
      event1_summary           = sum(out_getResults$y_1),
      event1_summary_exposed   = describeEvents(exposed_events1, sum(out_getResults$x_a), TRUE),
      event1_summary_unexposed = describeEvents(unexposed_events1, length(out_getResults$y_1) - sum(out_getResults$x_a), FALSE),
      effect2_text             = paste(estimand$effect.measure2, "of", exposure, "at", estimand$time.point, "( ref =", estimand$code.exposure.ref, ")"),
      event2_summary           = sum(out_getResults$y_2),
      event2_summary_exposed   = describeEvents(exposed_events2, sum(out_getResults$x_a), TRUE),
      event2_summary_unexposed = describeEvents(unexposed_events2, length(out_getResults$y_2) - sum(out_getResults$x_a), FALSE),
      median_followup          = paste(round(median(out_getResults$t),digit=2), "[", round(min(out_getResults$t),digit=2), ",", round(max(out_getResults$t),digit=2), "]")
    )
  } else if (outcome.type == "SURVIVAL") {
    text_values <- list(
      terms_text               = terms_text,
      effect1_text             = paste(estimand$effect.measure1, "of", exposure, "at", estimand$time.point, "( ref =", estimand$code.exposure.ref, ")"),
      event1_summary           = sum(out_getResults$y_1),
      event1_summary_exposed   = describeEvents(exposed_events1, sum(out_getResults$x_a), TRUE),
      event1_summary_unexposed = describeEvents(unexposed_events1, length(out_getResults$y_1) - sum(out_getResults$x_a), FALSE),
      median_followup          = paste(round(median(out_getResults$t),digit=2), "[", round(min(out_getResults$t),digit=2), ",", round(max(out_getResults$t),digit=2), "]")
    )
  } else if (outcome.type == "BINOMIAL") {
    text_values <- list(
      terms_text               = terms_text,
      effect1_text             = paste(estimand$effect.measure1, "of", exposure, "( ref =", estimand$code.exposure.ref, ")"),
      event1_summary           = sum(out_getResults$y_1),
      event1_summary_exposed   = describeEvents(exposed_events1, sum(out_getResults$x_a), TRUE),
      event1_summary_unexposed = describeEvents(unexposed_events1, length(out_getResults$y_1) - sum(out_getResults$x_a), FALSE),
      median_followup          = paste(round(median(out_getResults$t),digit=2), "[", round(min(out_getResults$t),digit=2), ",", round(max(out_getResults$t),digit=2), "]")
    )
  }

  if (outcome.type == "COMPETING-RISK") {
    tg <- list(
      event1 = list(
        tidy = tidy_df(coef1, text_values$terms_text),
        glance = glance_df(report.optim.convergence,
                           text_values$effect1_text,
                           text_values$event1_summary,
                           text_values$event1_summary_exposed,
                           text_values$event1_summary_unexposed,
                           text_values$median_followup,
                           iteration,
                           param_diff,
                           max(sol$fvec),
                           sol$iter,
                           i_parameter[3],
                           getMessage(outer.optim.method, sol))
      ),
      event2 = list(
        tidy = tidy_df(coef2, text_values$terms_text),
        glance = glance_df(report.optim.convergence,
                           text_values$effect2_text,
                           text_values$event2_summary,
                           text_values$event2_summary_exposed,
                           text_values$event2_summary_unexposed,
                           text_values$median_followup,
                           iteration,
                           param_diff,
                           max(sol$fvec),
                           sol$iter,
                           i_parameter[7] - i_parameter[4] + 1,
                           "-")
      )
    )
    class(tg$event1) <- class(tg$event2) <- "modelsummary_list"
  } else {
    tg1 <- list(
      tidy = tidy_df(coef1, text_values$terms_text),
      glance = glance_df(report.optim.convergence,
                         text_values$effect1_text,
                         text_values$event1_summary,
                         text_values$event1_summary_exposed,
                         text_values$event1_summary_unexposed,
                         text_values$median_followup,
                         iteration,
                         param_diff,
                         max(sol$fvec),
                         sol$iter,
                         i_parameter[3],
                         getMessage(outer.optim.method, sol))
    )
    class(tg1) <- "modelsummary_list"
    tg <- list("event 1 (no competing risk)" = tg1)
  }
  return(tg)
}

reportConstantEffects <- function(nuisance.model,
                                  exposure,
                                  estimand,
                                  out_bootstrap,
                                  out_getResults,
                                  iteration,
                                  param_diff,
                                  sol,
                                  outer.optim.method) {
  if (is.null(out_bootstrap)) return(NULL)

  i_parameter <- out_getResults$i_parameter

  x_a_mat <- as.matrix(out_getResults$x_a)
  lev <- attr(out_getResults$x_a, "levels")
  ref <- attr(out_getResults$x_a, "reference")
  if (is.null(lev) || is.null(ref)) {
    nonref_lev <- if (ncol(x_a_mat) > 0L) paste0("lvl", seq_len(ncol(x_a_mat))) else character(0L)
    ref <- estimand$code.exposure.ref %||% "ref"
  } else {
    nonref_lev <- setdiff(lev, ref)
  }

  # 非参照ごとの項目名
  exposure_terms <- if (length(nonref_lev)) {
    paste0(exposure, " [", nonref_lev, " vs ", ref, "]")
  } else character(0L)

  # βの位置（定数効果ブートストラップは係数配列の同じ並びを仮定）
  idx_beta1 <- if (length(nonref_lev)) seq.int(i_parameter[2], i_parameter[3]) else integer(0L)
  idx_beta2 <- if (length(nonref_lev)) seq.int(i_parameter[6], i_parameter[7]) else integer(0L)

  # ベクトルで抽出
  pick_boot <- function(idx) {
    if (!length(idx)) return(data.frame(term = character(0), estimate = numeric(0),
                                        std.error = numeric(0), p.value = numeric(0),
                                        conf.low = numeric(0), conf.high = numeric(0)))
    data.frame(
      term      = exposure_terms,
      estimate  = out_bootstrap$boot.coef[idx],
      std.error = out_bootstrap$boot.coef_se[idx],
      p.value   = out_bootstrap$boot.p_value[idx],
      conf.low  = out_bootstrap$boot.conf_low[idx],
      conf.high = out_bootstrap$boot.conf_high[idx]
    )
  }

  exposed_any <- as.integer(rowSums(x_a_mat) > 0L)
  y1 <- out_getResults$y_1
  y2 <- out_getResults$y_2

  effect1_text <- paste(estimand$effect.measure1, "of", exposure, "over time (ref =", ref, ")")
  effect2_text <- paste(estimand$effect.measure2, "of", exposure, "over time (ref =", ref, ")")

  glance_df_local <- function(effect_text, events, exposed_events, unexposed_events, param_len, message) {
    data.frame(
      effect.measure = effect_text,
      n.events = events,
      n.events.exposed = exposed_events,
      n.events.unexposed = unexposed_events,
      median.follow.up = paste(round(median(out_getResults$t), 2),
                               "[", round(min(out_getResults$t), 2), ",", round(max(out_getResults$t), 2), "]"),
      n.parameters = param_len,
      optimization.message = message
    )
  }

  tg <- list(
    event1 = list(
      tidy = pick_boot(idx_beta1),
      glance = glance_df_local(effect1_text,
                               sum(y1),
                               describeEvents(sum(y1 * exposed_any), sum(exposed_any), TRUE),
                               describeEvents(sum(y1 * (1 - exposed_any)), length(y1) - sum(exposed_any), FALSE),
                               length(idx_beta1),
                               getMessage(outer.optim.method, sol))
    ),
    event2 = list(
      tidy = pick_boot(idx_beta2),
      glance = glance_df_local(effect2_text,
                               sum(y2),
                               describeEvents(sum(y2 * exposed_any), sum(exposed_any), TRUE),
                               describeEvents(sum(y2 * (1 - exposed_any)), length(y2) - sum(exposed_any), FALSE),
                               length(idx_beta2),
                               "-")
    )
  )
  class(tg$event1) <- class(tg$event2) <- "modelsummary_list"
  return(tg)
}


reportConstantEffects_ <- function(nuisance.model,
                               exposure,
                               estimand,
                               out_bootstrap,
                               out_getResults,
                               iteration,
                               param_diff,
                               sol,
                               outer.optim.method) {
  if (is.null(out_bootstrap)) return(NULL)

  i_parameter <- out_getResults$i_parameter
  coef1 <- list(coef = out_bootstrap$boot.coef[1], coef_se = out_bootstrap$boot.coef_se[1], conf_low = out_bootstrap$boot.conf_low[1], conf_high = out_bootstrap$boot.conf_high[1], p_value = out_bootstrap$boot.p_value[1])
  coef2 <- list(coef = out_bootstrap$boot.coef[2], coef_se = out_bootstrap$boot.coef_se[2], conf_low = out_bootstrap$boot.conf_low[2], conf_high = out_bootstrap$boot.conf_high[2], p_value = out_bootstrap$boot.p_value[2])

  text_values <- list(
    exposure_text   = paste(exposure, "( ref =", estimand$code.exposure.ref, ")"),
    effect1_text    = paste(estimand$effect.measure1, "of", exposure, "over time"),
    effect2_text    = paste(estimand$effect.measure2, "of", exposure, "over time"),
    median_followup = paste(round(median(out_getResults$t),digit=2), "[", round(min(out_getResults$t),digit=3), ",", round(max(out_getResults$t),digit=2), "]")
  )

  tidy_df <- function(coef, text) {
    data.frame(
      term = text,
      estimate = coef$coef,
      std.error = coef$coef_se,
      p.value = coef$p_value,
      conf.low = coef$conf_low,
      conf.high = coef$conf_high
    )
  }

  glance_df <- function(effect_text, events, exposed_events, unexposed_events, param_len, message) {
    data.frame(
      effect.measure = effect_text,
      n.events = events,
      n.events.exposed = exposed_events,
      n.events.unexposed = unexposed_events,
      median.follow.up = text_values$median_followup,
      n.parameters = param_len,
      optimization.message = message
    )
  }

  tg <- list(
    event1 = list(
      tidy = tidy_df(coef1, text_values$exposure_text),
      glance = glance_df(text_values$effect1_text,
                         sum(out_getResults$y_1),
                         describeEvents(sum(out_getResults$y_1 * out_getResults$x_a), sum(out_getResults$x_a), TRUE),
                         describeEvents(sum(out_getResults$y_1 * (1 - out_getResults$x_a)), length(out_getResults$y_1) - sum(out_getResults$x_a), FALSE),
                         i_parameter[3],
                         getMessage(outer.optim.method, sol))
    ),
    event2 = list(
      tidy = tidy_df(coef2, text_values$exposure_text),
      glance = glance_df(text_values$effect2_text,
                         sum(out_getResults$y_2), describeEvents(sum(out_getResults$y_2 * out_getResults$x_a), sum(out_getResults$x_a), TRUE),
                         describeEvents(sum(out_getResults$y_2 * (1 - out_getResults$x_a)), length(out_getResults$y_2) - sum(out_getResults$x_a), FALSE),
                         i_parameter[7] - i_parameter[4] + 1,
                         "-")
    )
  )
  class(tg$event1) <- class(tg$event2) <- "modelsummary_list"
  return(tg)
}
