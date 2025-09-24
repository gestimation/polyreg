reportEffects_new <- function(outcome.type,
                          report.nuisance.parameter,
                          report.optim.convergence,
                          report.boot.conf,
                          nuisance.model,
                          exposure,
                          estimand,               # estimand$code.exposure.ref を参照
                          alpha_beta_estimated,
                          cov_estimated,
                          out_bootstrap,
                          out_getResults,         # ここに x_a (K-1), y_1, y_2, t を想定
                          iteration,
                          param_diff,
                          sol,
                          conf.level,
                          outer.optim.method) {

  index.vector <- estimand$index.vector
  ref_lab <- as.character(estimand$code.exposure.ref)

  # ----- 係数の取り出しと term ラベルの構築 -----

  # 1) 交絡（nuisance）パラメータ: 既存仕様を踏襲
  if (report.nuisance.parameter == TRUE) {
    coef1 <- getCoef(1:index.vector[3], alpha_beta_estimated, cov_estimated, report.boot.conf, out_bootstrap, conf.level)
    nuisance_terms <- c("Intercept", attr(terms(nuisance.model), "term.labels"))
  } else {
    # nuisance を報告しない場合は、coef1 は「曝露効果ブロック」のみを取得
    coef1 <- getCoef(index.vector[3], alpha_beta_estimated, cov_estimated, report.boot.conf, out_bootstrap, conf.level)
    nuisance_terms <- character(0)
  }

  # 2) 曝露効果の「係数ラベル」を作る（reference vs each level）
  #    - x_a の列名は treatment coding で “非参照レベル” のみ
  #    - モデルの係数順と列名の対応が取れている前提
  #    - terms_text_exposure は <exposure>: <level> vs <ref>（K-1 行）
  x_a_mat <- out_getResults$x_a
  if (is.null(dim(x_a_mat))) x_a_mat <- matrix(x_a_mat, ncol = 1)  # 念のため
  nonref_levels <- colnames(x_a_mat)
  terms_text_exposure <- if (length(nonref_levels)) {
    paste0(exposure, ": ", sub("^.*_", "", nonref_levels), " vs ", ref_lab)
  } else {
    character(0)
  }

  # 3) terms_text（nuisance + exposure）
  terms_text <- c(nuisance_terms, terms_text_exposure)

  # ----- by-level 集計の作成（K列へ拡張して集計） -----

  # K列ワンホットへ復元
  Xk <- asLevelMatrix(x_a_mat, ref_lab)                         # n x K
  level_names <- colnames(Xk)

  # 観測数 by level
  obs_by_level <- colSums(Xk)

  # 追跡期間メッセージ
  median_followup <- paste(
    round(median(out_getResults$t), 2), "[",
    round(min(out_getResults$t), 2), ",",
    round(max(out_getResults$t), 2), "]"
  )

  # イベント1の集計
  exposed_events_by_level_1 <- as.numeric(t(out_getResults$y_1) %*% Xk) # 長さK
  names(exposed_events_by_level_1) <- level_names
  event1_total <- sum(out_getResults$y_1)

  # COMPETING-RISK の場合はイベント2も
  if (outcome.type == 'COMPETING-RISK') {
    exposed_events_by_level_2 <- as.numeric(t(out_getResults$y_2) %*% Xk)
    names(exposed_events_by_level_2) <- level_names
    event2_total <- sum(out_getResults$y_2)
  }

  # ----- 出力（modelsummary_list） -----

  # event1: tidy（係数は coef1、ラベルは terms_text）
  tg_event1_tidy <- tidy_df(coef1, terms_text)

  tg_event1_glance <- glance_df(
    report.optim.convergence,
    effect_text          = paste(estimand$effect.measure1, "of", exposure, "at", estimand$time.point,
                                 "(ref =", ref_lab, ")"),
    events_total         = event1_total,
    events_by_level      = exposed_events_by_level_1,
    obs_by_level         = obs_by_level,
    median_followup      = median_followup,
    param_len            = if (report.nuisance.parameter) index.vector[3] else length(coef1$coef),
    iteration            = iteration,
    max.estimate.difference = param_diff,
    max.function.value   = max(sol$fvec),
    n.solver.iteration   = if (report.nuisance.parameter) index.vector[3] else length(coef1$coef),
    message              = getMessage(outer.optim.method, sol)
  )

  if (outcome.type == 'COMPETING-RISK') {
    # event2 系の係数
    if (report.nuisance.parameter == TRUE) {
      coef2 <- getCoef(index.vector[4]:index.vector[7], alpha_beta_estimated, cov_estimated, report.boot.conf, out_bootstrap, conf.level)
      # ラベルは nuisance_terms + exposure（同じでOK）
      tg_event2_tidy <- tidy_df(coef2, terms_text)
      npar2 <- index.vector[7] - index.vector[4] + 1
    } else {
      coef2 <- getCoef(index.vector[7], alpha_beta_estimated, cov_estimated, report.boot.conf, out_bootstrap, conf.level)
      tg_event2_tidy <- tidy_df(coef2, terms_text)  # exposureブロックのみになる想定
      npar2 <- length(coef2$coef)
    }

    tg_event2_glance <- glance_df(
      report.optim.convergence,
      effect_text          = paste(estimand$effect.measure2, "of", exposure, "at", estimand$time.point,
                                   "(ref =", ref_lab, ")"),
      events_total         = event2_total,
      events_by_level      = exposed_events_by_level_2,
      obs_by_level         = obs_by_level,
      median_followup      = median_followup,
      param_len            = npar2,
      iteration            = iteration,
      max.estimate.difference = param_diff,
      max.function.value   = max(sol$fvec),
      n.solver.iteration   = npar2,
      message              = "-"  # 既存仕様を踏襲
    )

    tg <- list(
      event1 = list(tidy = tg_event1_tidy, glance = tg_event1_glance),
      event2 = list(tidy = tg_event2_tidy, glance = tg_event2_glance)
    )
    class(tg$event1) <- class(tg$event2) <- "modelsummary_list"
    return(tg)
  }

  # SURVIVAL / BINOMIAL（event1のみ）
  tg1 <- list(
    tidy   = tg_event1_tidy,
    glance = tg_event1_glance
  )
  class(tg1) <- "modelsummary_list"
  tg <- list("event 1 (no competing risk)" = tg1)
  return(tg)
}


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

asLevelMatrix <- function(x_a, ref_lab) {
  if (is.null(dim(x_a))) {
    # ベクトルで来た場合でも行列に変換
    x_a <- matrix(x_a, ncol = 1)
  }
  # 参照カテゴリ列を復元
  ref_col <- 1 - rowSums(x_a)
  Xk <- cbind(ref_col, x_a)
  colnames(Xk)[1] <- ref_lab
  Xk
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

  index.vector <- estimand$index.vector
  if (report.nuisance.parameter==TRUE) {
    coef1 <- getCoef(1:index.vector[3], alpha_beta_estimated, cov_estimated, report.boot.conf, out_bootstrap, conf.level)
    terms_text <- c("Intercept", attr(terms(nuisance.model), "term.labels"), paste(exposure, "( ref =", estimand$code.exposure.ref, ")"))
  } else {
    coef1 <- getCoef(index.vector[3], alpha_beta_estimated, cov_estimated, report.boot.conf, out_bootstrap, conf.level)
    terms_text <- paste(exposure, "( ref =", estimand$code.exposure.ref, ")")
  }
  exposed_events1 <- sum(out_getResults$y_1 * out_getResults$x_a)
  unexposed_events1 <- sum(out_getResults$y_1 * (1 - out_getResults$x_a))

  if (outcome.type == 'COMPETING-RISK') {
    if (report.nuisance.parameter==TRUE) {
      coef2 <- getCoef(index.vector[4]:index.vector[7], alpha_beta_estimated, cov_estimated, report.boot.conf, out_bootstrap, conf.level)
    } else {
      coef2 <- getCoef(index.vector[7], alpha_beta_estimated, cov_estimated, report.boot.conf, out_bootstrap, conf.level)
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
  } else if (outcome.type == 'SURVIVAL') {
    text_values <- list(
      terms_text               = terms_text,
      effect1_text             = paste(estimand$effect.measure1, "of", exposure, "at", estimand$time.point, "( ref =", estimand$code.exposure.ref, ")"),
      event1_summary           = sum(out_getResults$y_1),
      event1_summary_exposed   = describeEvents(exposed_events1, sum(out_getResults$x_a), TRUE),
      event1_summary_unexposed = describeEvents(unexposed_events1, length(out_getResults$y_1) - sum(out_getResults$x_a), FALSE),
      median_followup          = paste(round(median(out_getResults$t),digit=2), "[", round(min(out_getResults$t),digit=2), ",", round(max(out_getResults$t),digit=2), "]")
    )
  } else if (outcome.type == 'BINOMIAL') {
    text_values <- list(
      terms_text               = terms_text,
      effect1_text             = paste(estimand$effect.measure1, "of", exposure, "( ref =", estimand$code.exposure.ref, ")"),
      event1_summary           = sum(out_getResults$y_1),
      event1_summary_exposed   = describeEvents(exposed_events1, sum(out_getResults$x_a), TRUE),
      event1_summary_unexposed = describeEvents(unexposed_events1, length(out_getResults$y_1) - sum(out_getResults$x_a), FALSE),
      median_followup          = paste(round(median(out_getResults$t),digit=2), "[", round(min(out_getResults$t),digit=2), ",", round(max(out_getResults$t),digit=2), "]")
    )
  }

  if (outcome.type == 'COMPETING-RISK') {
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
                           index.vector[3],
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
                           index.vector[7] - index.vector[4] + 1,
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
                         index.vector[3],
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

  index.vector <- estimand$index.vector
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
                         index.vector[3],
                         getMessage(outer.optim.method, sol))
    ),
    event2 = list(
      tidy = tidy_df(coef2, text_values$exposure_text),
      glance = glance_df(text_values$effect2_text,
                         sum(out_getResults$y_2), describeEvents(sum(out_getResults$y_2 * out_getResults$x_a), sum(out_getResults$x_a), TRUE),
                         describeEvents(sum(out_getResults$y_2 * (1 - out_getResults$x_a)), length(out_getResults$y_2) - sum(out_getResults$x_a), FALSE),
                         index.vector[7] - index.vector[4] + 1,
                         "-")
    )
  )
  class(tg$event1) <- class(tg$event2) <- "modelsummary_list"
  return(tg)
}
