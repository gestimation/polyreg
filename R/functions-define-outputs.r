reportCompetingRisk <- function(nuisance.model, exposure, estimand, alpha_beta_estimated, cov_estimated,
                                out_getResults, iteration, param_diff, sol, conf.level, outer.optim.method) {
  alpha <- 1 - conf.level
  critical_value <- qnorm(1 - alpha / 2)
  i_parameter <- out_getResults$i_parameter

  getCoefDetails <- function(index) {
    coef <- alpha_beta_estimated[index]
    coef_se <- sqrt(diag(cov_estimated)[index])
    conf_low <- coef - critical_value * coef_se
    conf_high <- coef + critical_value * coef_se
    p_value <- 2 * (1 - pnorm(abs(coef) / coef_se))
    list(coef = coef, coef_se = coef_se, conf_low = conf_low, conf_high = conf_high, p_value = p_value)
  }

  coef1 <- getCoefDetails(i_parameter[3])
  coef2 <- getCoefDetails(i_parameter[7])

  summary_text <- function(events, observations, exposed) {
    paste(events, "events in", observations, if (exposed) "exposed observations" else "unexposed observations")
  }

  text_values <- list(
    exposure_text   = paste(exposure, "( ref =", estimand$code.exposure.ref, ")"),
    effect1_text    = paste(estimand$effect.measure1, "of", exposure, "at", estimand$time.point),
    effect2_text    = paste(estimand$effect.measure2, "of", exposure, "at", estimand$time.point),
    median_followup = paste(median(out_getResults$t), "[", min(out_getResults$t), ",", max(out_getResults$t), "]")
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

  getMessage <- function() {
    if (outer.optim.method %in% c("nleqslv", "Newton", "Broyden")) return(sol$message)
    if (outer.optim.method %in% c("optim", "BFGS", "SANN")) return(ifelse(!is.null(sol$message), sol$message, "None"))
    return("-")
  }

  tg <- list(
    event1 = list(
      tidy = tidy_df(coef1, text_values$exposure_text),
      glance = glance_df(text_values$effect1_text,
                         sum(out_getResults$y_1),
                         summary_text(sum(out_getResults$y_1 * out_getResults$x_a), sum(out_getResults$x_a), TRUE),
                         summary_text(sum(out_getResults$y_1 * (1 - out_getResults$x_a)), length(out_getResults$y_1) - sum(out_getResults$x_a), FALSE),
                         i_parameter[3],
                         getMessage())
    ),
    event2 = list(
      tidy = tidy_df(coef2, text_values$exposure_text),
      glance = glance_df(text_values$effect2_text,
                         sum(out_getResults$y_2), summary_text(sum(out_getResults$y_2 * out_getResults$x_a), sum(out_getResults$x_a), TRUE),
                         summary_text(sum(out_getResults$y_2 * (1 - out_getResults$x_a)), length(out_getResults$y_2) - sum(out_getResults$x_a), FALSE),
                         i_parameter[7] - i_parameter[4] + 1,
                         "-")
    )
  )

  class(tg$event1) <- class(tg$event2) <- "modelsummary_list"
  return(tg)
}

reportSurvival <- function(nuisance.model, exposure, estimand, alpha_beta_estimated, cov_estimated,
                           out_getResults, iteration, param_diff, sol, conf.level, outer.optim.method) {
  alpha <- 1 - conf.level
  critical_value <- qnorm(1 - alpha / 2)
  i_parameter <- out_getResults$i_parameter

  getCoefDetails <- function(index) {
    coef <- alpha_beta_estimated[index]
    coef_se <- sqrt(diag(cov_estimated)[index])
    conf_low <- coef - critical_value * coef_se
    conf_high <- coef + critical_value * coef_se
    p_value <- 2 * (1 - pnorm(abs(coef) / coef_se))
    list(coef = coef, coef_se = coef_se, conf_low = conf_low, conf_high = conf_high, p_value = p_value)
  }

  summary_text <- function(events, observations, exposed) {
    paste(events, "events in", observations, if (exposed) "exposed observations" else "unexposed observations")
  }

  coef1 <- getCoefDetails(i_parameter[3])

  exposed_events <- sum(out_getResults$y_1 * out_getResults$x_a)
  unexposed_events <- sum(out_getResults$y_1 * (1 - out_getResults$x_a))

  text_values <- list(
    exposure_text           = paste(exposure, "( ref =", estimand$code.exposure.ref, ")"),
    effect_text             = paste(estimand$effect.measure1, "of", exposure, "at", estimand$time.point, "( ref =", estimand$code.exposure.ref, ")"),
    event_summary           = sum(out_getResults$y_1),
    event_summary_exposed   = summary_text(exposed_events, sum(out_getResults$x_a), TRUE),
    event_summary_unexposed = summary_text(unexposed_events, length(out_getResults$y_1) - sum(out_getResults$x_a), FALSE),
    median_followup         = paste(median(out_getResults$t), "[", min(out_getResults$t), ",", max(out_getResults$t), "]")
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

  getMessage <- function() {
    if (outer.optim.method %in% c("nleqslv", "Newton", "Broyden", "partial")) return(sol$message)
    if (outer.optim.method %in% c("optim", "BFGS", "SANN")) return(ifelse(!is.null(sol$message), sol$message, "None"))
    return("-")
  }

  glance_df <- data.frame(
    effect.measure = text_values$effect_text,
    n.events = text_values$event_summary,
    n.events.exposed = text_values$event_summary_exposed,
    n.events.unexposed = text_values$event_summary_unexposed,
    median.follow.up = text_values$median_followup,
    n.parameters = length(alpha_beta_estimated[1:i_parameter[3]]),
    optimization.message = getMessage()
  )

  tg1 <- list(tidy = tidy_df(coef1, text_values$exposure_text), glance = glance_df)
  class(tg1) <- "modelsummary_list"
  tg <- list("event 1 (no competing risk)" = tg1)

  return(tg)
}


reportCompetingRiskFull <- function(nuisance.model, exposure, estimand, alpha_beta_estimated, cov_estimated,
                                     out_getResults, iteration, param_diff, sol, conf.level, outer.optim.method) {
  alpha <- 1 - conf.level
  critical_value <- qnorm(1 - alpha / 2)
  i_parameter <- out_getResults$i_parameter

  coef_1 <- alpha_beta_estimated[1:i_parameter[3]]
  coef_2 <- alpha_beta_estimated[i_parameter[4]:i_parameter[7]]
  coef_se_1 <- sqrt(diag(cov_estimated)[1:i_parameter[3]])
  coef_se_2 <- sqrt(diag(cov_estimated)[i_parameter[4]:i_parameter[7]])
  coef_1_l <- alpha_beta_estimated[1:i_parameter[3]] - critical_value * coef_se_1
  coef_1_u <- alpha_beta_estimated[1:i_parameter[3]] + critical_value * coef_se_1
  coef_2_l <- alpha_beta_estimated[i_parameter[4]:i_parameter[7]] - critical_value * coef_se_2
  coef_2_u <- alpha_beta_estimated[i_parameter[4]:i_parameter[7]] + critical_value * coef_se_2
  coef_z_1 <- abs(alpha_beta_estimated[1:i_parameter[3]]) / coef_se_1
  coef_z_2 <- abs(alpha_beta_estimated[i_parameter[4]:i_parameter[7]]) / coef_se_2
  coef_p_1 <- 2 * (1 - pnorm(coef_z_1))
  coef_p_2 <- 2 * (1 - pnorm(coef_z_2))

  t <- out_getResults$t
  y_0_ <- out_getResults$y_0_
  y_1_ <- out_getResults$y_1_
  y_2_ <- out_getResults$y_2_
  y_0 <- out_getResults$y_0
  y_1 <- out_getResults$y_1
  y_2 <- out_getResults$y_2

  effect.modifier <- NULL
  if (!is.null(effect.modifier)) {
    text1 <- c("Intercept", attr(terms(nuisance.model), "term.labels"), exposure, "Interaction")
  } else {
    text1 <- c("Intercept", attr(terms(nuisance.model), "term.labels"), exposure)
  }
  text2 <- paste(estimand$effect.measure1, "of", exposure, "at", estimand$time.point, "( ref =", estimand$code.exposure.ref, ")")
  text3 <- paste(estimand$effect.measure2, "of", exposure, "at", estimand$time.point, "( ref =", estimand$code.exposure.ref, ")")
  text4 <- paste(sum(y_1), "events in", length(y_1), "observations")
  text5 <- paste(sum(y_2), "events in", length(y_2), "observations")
  text6 <- paste(median(t), "[", min(t), ",", max(t), "]")

  ti_1 <- data.frame(
    term = text1,
    estimate = coef_1,
    std.error = coef_se_1,
    p.value = coef_p_1,
    conf.low = coef_1_l,
    conf.high = coef_1_u
  )
  ti_2 <- data.frame(
    term = text1,
    estimate = coef_2,
    std.error = coef_se_2,
    p.value = coef_p_2,
    conf.low = coef_2_l,
    conf.high = coef_2_u
  )
  if (outer.optim.method == 'nleqslv' || outer.optim.method == 'Newton' || outer.optim.method == 'Broyden' || outer.optim.method == 'partial'){
    gl_1 <- data.frame(
      effect.measure = text2,
      n.events = text4,
      median.follow.up = text6,
      #    n.parameters = length(alpha_beta_estimated[1:n_para_2]),
      #    min.parameters = min(alpha_beta_estimated[1:n_para_2]),
      #    max.parameter = max(alpha_beta_estimated[1:n_para_2]),
      n.loop.iteration = iteration,
      #    max.estimate.difference = param_diff,
      n.optimization.iteration = sol$iter,
      max.function.value = max(sol$fvec),
      optimization.message = sol$message
    )
    gl_2 <- data.frame(
      effect.measure = text3,
      n.events = text5,
      median.follow.up = text6,
      #    n.parameters = length(alpha_beta_estimated[n_para_3:n_para_5]),
      #    min.parameters = min(alpha_beta_estimated[n_para_3:n_para_5]),
      #    max.parameter = max(alpha_beta_estimated[n_para_3:n_para_5]),
      n.loop.iteration = "-",
      #    max.estimate.difference = param_diff,
      n.optimization.iteration = "-",
      max.function.value = "-",
      optimization.message = "-"
    )
  } else if (outer.optim.method == 'optim' || outer.optim.method == 'BFGS' || outer.optim.method == 'SANN') {
    gl_1 <- data.frame(
      effect.measure = text2,
      n.events = text4,
      median.follow.up = text6,
      #n.parameters = length(alpha_beta_estimated[1:n_para_2]),
      #min.parameters = min(alpha_beta_estimated[1:n_para_2]),
      #max.parameter = max(alpha_beta_estimated[1:n_para_2]),
      n.loop.iteration = iteration,
      #max.estimate.difference = param_diff,
      n.optimization.iteration = as.vector(sol$counts)[2],
      max.function.value = max(sol$value),
      optimization.message = if (!is.null(sol$message)) sol$message else "None"
    )
    gl_2 <- data.frame(
      effect.measure = text3,
      n.events = text5,
      median.follow.up = text6,
      #n.parameters = length(alpha_beta_estimated[n_para_3:n_para_5]),
      #min.parameters = min(alpha_beta_estimated[n_para_3:n_para_5]),
      #max.parameter = max(alpha_beta_estimated[n_para_3:n_para_5]),
      n.loop.iteration = "-",
      #max.estimate.difference = param_diff,
      n.optimization.iteration = "-",
      max.function.value = "-",
      optimization.message = "-"
    )
  } else {
    gl_1 <- data.frame(
      effect.measure = text2,
      n.events = text4,
      median.follow.up = text6,
      #n.parameters = length(alpha_beta_estimated[1:n_para_2]),
      #min.parameters = min(alpha_beta_estimated[1:n_para_2]),
      #max.parameter = max(alpha_beta_estimated[1:n_para_2]),
      n.loop.iteration = iteration,
      #max.estimate.difference = param_diff,
      n.optimization.iteration = sol$iter,
      max.function.value = max(sol$f.root),
      optimization.message = "-"
    )
    gl_2 <- data.frame(
      effect.measure = text3,
      n.events = text5,
      median.follow.up = text6,
      #    n.parameters = length(alpha_beta_estimated[n_para_3:n_para_5]),
      #    min.parameters = min(alpha_beta_estimated[n_para_3:n_para_5]),
      #    max.parameter = max(alpha_beta_estimated[n_para_3:n_para_5]),
      n.loop.iteration = "-",
      #    max.estimate.difference = param_diff,
      n.optimization.iteration = "-",
      max.function.value = "-",
      optimization.message = "-"
    )
  }
  tg1 <- list(tidy = ti_1,glance = gl_1)
  tg2 <- list(tidy = ti_2,glance = gl_2)
  class(tg1) <- "modelsummary_list"
  class(tg2) <- "modelsummary_list"
  tg <- list()
  tg[[1]] <- tg1
  tg[[2]] <- tg2
  names(tg) <- c("event 1", "event 2")
  return(tg)
}

reportSurvivalFull <- function(nuisance.model, exposure, estimand, alpha_beta_estimated, cov_estimated,
                               out_getResults, iteration, param_diff, sol, conf.level, outer.optim.method) {
  alpha <- 1 - conf.level
  critical_value <- qnorm(1 - alpha / 2)
  i_parameter <- out_getResults$i_parameter

  coef_1 <- alpha_beta_estimated[1:i_parameter[3]]
  coef_se_1 <- sqrt(diag(cov_estimated)[1:i_parameter[3]])
  coef_1_l <- alpha_beta_estimated[1:i_parameter[3]] - critical_value * coef_se_1
  coef_1_u <- alpha_beta_estimated[1:i_parameter[3]] + critical_value * coef_se_1
  coef_z_1 <- abs(alpha_beta_estimated[1:i_parameter[3]]) / coef_se_1
  coef_p_1 <- 2 * (1 - pnorm(coef_z_1))

  t <- out_getResults$t
  y_0_ <- out_getResults$y_0_
  y_1_ <- out_getResults$y_1_
  y_0 <- out_getResults$y_0
  y_1 <- out_getResults$y_1
  effect.modifier <- NULL
  if (!is.null(effect.modifier)) {
    text1 <- c("Intercept", attr(terms(nuisance.model), "term.labels"), exposure, "Interaction")
  } else {
    text1 <- c("Intercept", attr(terms(nuisance.model), "term.labels"), exposure)
  }
  text2 <- paste(estimand$effect.measure1, "of", exposure, "at", estimand$time.point, "( ref =", estimand$code.exposure.ref, ")")
  text4 <- paste(sum(y_1), "events in", length(y_1), "observations")
  text6 <- paste(median(t), "[", min(t), ",", max(t), "]")

  ti_1 <- data.frame(
    term = text1,
    estimate = coef_1,
    std.error = coef_se_1,
    p.value = coef_p_1,
    conf.low = coef_1_l,
    conf.high = coef_1_u
  )

  if (outer.optim.method == 'nleqslv' || outer.optim.method == 'Newton' || outer.optim.method == 'Broyden' || outer.optim.method == 'partial'){
    gl_1 <- data.frame(
      effect.measure = text2,
      n.events = text4,
      median.follow.up = text6,
      #    n.parameters = length(alpha_beta_estimated[1:n_para_2]),
      #    min.parameters = min(alpha_beta_estimated[1:n_para_2]),
      #    max.parameter = max(alpha_beta_estimated[1:n_para_2]),
      n.loop.iteration = iteration,
      #    max.estimate.difference = param_diff,
      #    max.function.value = max(sol$fvec),
      n.optimization.iteration = sol$iter,
      optimization.message = sol$message
    )
  } else if (outer.optim.method == 'optim' || outer.optim.method == 'BFGS' || outer.optim.method == 'SANN') {
    gl_1 <- data.frame(
      effect.measure = text2,
      n.events = text4,
      median.follow.up = text6,
      #n.parameters = length(alpha_beta_estimated[1:n_para_2]),
      #min.parameters = min(alpha_beta_estimated[1:n_para_2]),
      #max.parameter = max(alpha_beta_estimated[1:n_para_2]),
      n.loop.iteration = iteration,
      #max.estimate.difference = param_diff,
      #max.function.value = max(sol$value),
      n.solver.iteration = as.vector(sol$counts)[2],
      optimization.message = if (!is.null(sol$message)) sol$message else "None"
    )
  } else {
    gl_1 <- data.frame(
      effect.measure = text2,
      n.events = text4,
      median.follow.up = text6,
      #n.parameters = length(alpha_beta_estimated[1:n_para_2]),
      #min.parameters = min(alpha_beta_estimated[1:n_para_2]),
      #max.parameter = max(alpha_beta_estimated[1:n_para_2]),
      n.loop.iteration = iteration,
      #max.estimate.difference = param_diff,
      #max.function.value = max(sol$f.root),
      n.optimization.iteration = sol$iter,
      optimization.message = "-"
    )
  }

  tg1 <- list(tidy = ti_1,glance = gl_1)
  class(tg1) <- "modelsummary_list"
  tg <- list()
  tg[[1]] <- tg1
  names(tg) <- c("event 1 (no competing risk)")
  return(tg)
}

reportPrediction <- function(formula,
                                 data,
                                 exposure,
                                 alpha_beta,
                                 cov = NULL,
                                 outcome.type,
                                 estimand,
                                 optim.method = NULL,
                                 prob.bound = 1e-5,
                                 initial_pred = NULL)
{
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
    if (any(t<0))
      stop("Expected non-negative time variable")
    epsilon <- Y[, 2]  # status variable
    if (!all(epsilon %in% c(0, 1, 2)))
      stop("Expected only 0, 1 or 2, with 0 representing censoring")
    time.name <- all.vars(formula)[1]
    status.name <- all.vars(formula)[2]
  } else {
    stop("Expected only right censored data")
  }
  if (!is.null(offsetpos <- attributes(Terms)$specials$offset)) {
    ts <- survival::untangle.specials(Terms, "offset")
    Terms <- Terms[-ts$terms]
    offset <- m[[ts$vars]]
  } else {
    offset <- rep(0, length(t))
  }

  a_ <- as.factor(data[[exposure]])
  if (estimand$code.exposure.ref==0) {
    x_a <- as.matrix(model.matrix(~ a_)[, -1])
  } else {
    x_a <- as.matrix(rep(1,length(t)) - model.matrix(~ a_)[, -1])
  }
  x_l <- model.matrix(Terms, mf)
  x_la <- cbind(x_l, x_a)

  if (outcome.type == 'SURVIVAL') {
    pred <- calculatePotentialRisk(alpha_beta,x_l,offset,estimand)
  } else {
    pred <- calculatePotentialCIFs(alpha_beta,x_a,x_l,offset,epsilon,estimand,optim.method,prob.bound,initial_pred)
  }
  return(pred)
}
