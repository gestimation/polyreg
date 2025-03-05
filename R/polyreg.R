#' @title Direct polynomial regression for survival and competing risks analysis
#'
#' @param nuisance.model formula Model formula representing outcome, exposure and covariates
#' @param exposure character Column name representing the exposure (1 = exposed, 0 = not exposed).
#' @param cens.model formula Model formula representing censoring and covariates.
#' @param data data.frame Input dataset containing survival data.
#' @param effect.measure1 character Specifies the effect measure for event (RR, OR, SHR).
#' @param effect.measure2 character Specifies the effect measure for competing risk (RR, OR, SHR).
#' @param exposure.reference integer Specifies the code of the reference category of effects. Defaults to 1.
#' @param time.point numeric The time point for exposure effects to be estimated.
#' @param outcome.type character Specifies the type of outcome (COMPETINGRISK or SURVIVAL).
#' @param conf.level numeric The level for confidence intervals.
#' @param outer.optim.method character Specifies the method of optimization (nleqslv, optim, multiroot).
#' @param inner.optim.method character Specifies the method of optimization (nleqslv, optim, multiroot).
#' @param optim.parameter1 numeric Convergence threshold for outer loop. Defaults to 1e-5.
#' @param optim.parameter2 integer Maximum number of iterations. Defaults to 20.
#' @param optim.parameter3 numeric Constraint range for parameters. Defaults to 100.
#' @param optim.parameter4 numeric Convergence threshold for optim in outer loop. Defaults to 1e-5.
#' @param optim.parameter5 integer Maximum number of iterations for nleqslv/optim in outer loop. Defaults to 200.
#' @param optim.parameter6 numeric Convergence threshold for optim in inner loop. Defaults to 1e-10.
#' @param optim.parameter7 integer Maximum number of iterations for optim in inner loop. Defaults to 200.
#' @param data.initlal.values data.frame Optional dataset containing initial values. Defaults to NULL.
#' @param should.normalize.covariate logical Indicates whether covariates are normalized (TRUE = normalize, FALSE = otherwise). Defaults to TRUE.
#' @param should.sort.data logical Indicates whether data are initially sorted to reduce computation steps (TRUE = sort, FALSE = otherwise). Defaults to TRUE.
#' @param prob.bound numeric Small threshold for clamping probabilities. Defaults to 1e-5.
#' @importFrom mets
#' @importFrom nleqslv
#'
#' @return A list of results from direct polynomial regression. coefficient and cov are estimated regression coefficients of exposure and covariates and their variance covariance matrix. summary and summary.full meets requirement of msummary function.
#' @export
#'
#' @examples
#' #' data(bmt)
#' result <- polyreg(nuisance.model = Event(time, cause)~age+tcell, exposure = 'platelet',
#' cens.model = Event(time,cause==0)~+1, data = bmt, effect.measure1='RR', effect.measure2='RR', time.point=24, outcome.type='COMPETINGRISK')
#' msummary(result$out_summary, statistic = c("conf.int"), exponentiate = TRUE)
polyreg <- function(
    nuisance.model,
    exposure,
    cens.model,
    data,
    effect.measure1 = 'RR',
    effect.measure2 = 'RR',
    exposure.reference = 0,
    time.point,
    outcome.type = 'COMPETINGRISK',
    conf.level = 0.95,
    outer.optim.method = 'nleqslv',
    inner.optim.method = 'optim',
    optim.parameter1 = 1e-5,
    optim.parameter2 = 20,
    optim.parameter3 = 100,
    optim.parameter4 = 1e-5,
    optim.parameter5 = 200,
    optim.parameter6 = 1e-10,
    optim.parameter7 = 200,
    data.initial.values = NULL,
    should.normalize.covariate = TRUE,
    should.sort.data = TRUE,
    prob.bound = 1e-5
) {

  #######################################################################################################
  # 1. Pre-processing (function: checkSpell, checkInput1, normalizeCovariate, sortByCovariate)
  #######################################################################################################
  checkSpell(outcome.type, effect.measure1, effect.measure2)
  checkInput1(outcome.type, time.point, conf.level, outer.optim.method, inner.optim.method)
  outcome.type <- outcome.type.corrected
  effect.measure1 <- effect.measure1.corrected
  effect.measure2 <- effect.measure2.corrected
  time.point <- time.point.corrected

  estimand <- list(effect.measure1=effect.measure1, effect.measure2=effect.measure2, time.point=time.point, exposure.reference=exposure.reference)
  optim.method <- list(
    outer.optim.method = outer.optim.method,
    inner.optim.method = inner.optim.method,
    optim.parameter1 = optim.parameter1,
    optim.parameter2 = optim.parameter2,
    optim.parameter3 = optim.parameter3,
    optim.parameter4 = optim.parameter4,
    optim.parameter5 = optim.parameter5,
    optim.parameter6 = optim.parameter6,
    optim.parameter7 = optim.parameter7
  )

  out_normalizeCovariate <- normalizeCovariate(nuisance.model, data, should.normalize.covariate, outcome.type)
  normalized_data <- out_normalizeCovariate$normalized_data
  sorted_data <- sortByCovariate(nuisance.model, normalized_data, should.sort.data, out_normalizeCovariate$n_covariate)

  #######################################################################################################
  # 2. Pre-processing and Calculating initial values alpha_beta_0 (function: calculateInitialValues)
  #######################################################################################################
  if (outcome.type == 'PROPORTIONAL') {
    n_para_1 <- out_normalizeCovariate$n_covariate+1
    n_para_2 <- out_normalizeCovariate$n_covariate+2
    n_para_3 <- out_normalizeCovariate$n_covariate+3
    n_para_4 <- 2*out_normalizeCovariate$n_covariate+3
    n_para_5 <- 2*out_normalizeCovariate$n_covariate+4
    n_para_6 <- length(time.point)*(n_para_5-2) + 2
    alpha_beta_0 <- rep(NA, n_para_6) # initial parameters for event 1 and 2 over time points
    i_time <- 0
    sum1 <- 0
    sum2 <- 0
    for (specific.time in time.point) {
      specific.time <- as.numeric(specific.time)
      out_calculateInitialValues <- calculateInitialValues(
        formula = nuisance.model,
        data = sorted_data,
        exposure = exposure,
        data.initial.values = data.initial.values,
        estimand = estimand,
        specific.time = specific.time,
        outcome.type = outcome.type,
        prob.bound = prob.bound
      )
      i_time <- i_time + 1
      i_para <- n_para_1*(i_time-1)+1
      tmp1 <- out_calculateInitialValues[1:n_para_1]
      tmp2 <- out_calculateInitialValues[n_para_3:n_para_4]
      alpha_beta_0[i_para:(i_para+n_para_1-1)]                          <- tmp1
      alpha_beta_0[(n_para_6/2+i_para):(n_para_6/2+i_para+n_para_1-1)]  <- tmp2
      sum1 <- sum1 + out_calculateInitialValues[n_para_2]
      sum2 <- sum2 + out_calculateInitialValues[n_para_5]
    }
    alpha_beta_0[n_para_6/2]  <- sum1/length(time.point)
    alpha_beta_0[n_para_6]    <- sum2/length(time.point)
  } else {
    out_calculateInitialValues <- calculateInitialValues(
      formula = nuisance.model,
      data = sorted_data,
      exposure = exposure,
      data.initial.values = data.initial.values,
      estimand = estimand,
      specific.time = estimand$time.point,
      outcome.type = outcome.type,
      prob.bound = prob.bound
    )
    alpha_beta_0 <- out_calculateInitialValues
  }

  #######################################################################################################
  # 3. Calculating IPCW (function: calculateIPW)
  #######################################################################################################
  if (outcome.type == 'PROPORTIONAL') {
    i_time <- 0
    ip.weight.matrix <- matrix(NA, nrow=out_normalizeCovariate$n, ncol=length(time.point))
    for (specific.time in time.point) {
      i_time <- i_time + 1
      ip.weight.matrix[,i_time] <- calculateIPW(formula = cens.model, data = sorted_data, specific.time = specific.time)
    }
  } else {
    ip.weight <- calculateIPW(formula = cens.model, data = sorted_data, specific.time = estimand$time.point)
  }

  #######################################################################################################
  # 4. Parameter estimation (functions: estimating_equation_ipcw, _survival, _proportional)
  #######################################################################################################
  makeObjectiveFunction <- function() {
    out_ipcw <- list()
    initial.CIFs <- NULL
    estimating_equation_i <- function(p) {
      out_ipcw <- estimating_equation_ipcw(
        formula = nuisance.model,
        data = sorted_data,
        exposure = exposure,
        ip.weight = ip.weight,
        alpha_beta = p,
        estimand = estimand,
        optim.method = optim.method,
        prob.bound = prob.bound,
        initial.CIFs = initial.CIFs)
      out_ipcw <<- out_ipcw
      return(out_ipcw$ret)
    }
    estimating_equation_s <- function(p) {
      out_ipcw <- estimating_equation_survival(
        formula = nuisance.model,
        data = sorted_data,
        exposure = exposure,
        ip.weight = ip.weight,
        alpha_beta = p,
        estimand = estimand,
        prob.bound = prob.bound,
        initial.CIFs = initial.CIFs)
      out_ipcw <<- out_ipcw
      return(out_ipcw$ret)
    }
    estimating_equation_p <- function(p) {
      out_ipcw <- estimating_equation_proportional(
        formula = nuisance.model,
        data = sorted_data,
        exposure = exposure,
        ip.weight.matrix = ip.weight.matrix,
        alpha_beta = p,
        estimand = estimand,
        time.point = time.point,
        optim.method = optim.method,
        prob.bound = prob.bound,
        initial.CIFs = initial.CIFs)
      out_ipcw <<- out_ipcw
      return(out_ipcw$ret)
    }
    setInitialCIFs <- function(new.CIFs) {
      initial.CIFs <<- new.CIFs
    }
    get_results <- function() {
      out_ipcw
    }
    list(
      estimating_equation_i = estimating_equation_i,
      estimating_equation_s = estimating_equation_s,
      estimating_equation_p = estimating_equation_p,
      setInitialCIFs = setInitialCIFs,
      get_results = get_results
    )
  }

  iteration <- 0
  max_param_diff  <- Inf
  obj <- makeObjectiveFunction()
  sol_list <- list()
  diff_list <- list()

  if (outcome.type == 'COMPETINGRISK') {
    current_params <- alpha_beta_0
    while ((iteration < optim.parameter2) & (max_param_diff > optim.parameter1)) {
      iteration <- iteration + 1
      if (outer.optim.method == "nleqslv" | outer.optim.method == "Broyden"){
        sol <- nleqslv(current_params, obj$estimating_equation_i, method="Broyden", control=list(maxit=optim.parameter5, allowSingular=FALSE))
        new_params <- sol$x
      } else if (outer.optim.method == "Newton"){
        sol <- nleqslv(current_params, obj$estimating_equation_i, method="Newton", control=list(maxit=optim.parameter5, allowSingular=FALSE))
        new_params <- sol$x
      } else if (outer.optim.method == "multiroot") {
        sol <- multiroot(obj$estimating_equation_i, start = current_params, maxiter=optim.parameter5, rtol = optim.parameter4)
        new_params <- sol$root
      } else if (outer.optim.method == "optim" | outer.optim.method == "SANN"){
        sol <- optim(par = current_params,
                     fn = function(params) {
                       sum(obj$estimating_equation_i(params)^2)
                     },
                     method = "SANN",  control = list(maxit=optim.parameter5, reltol=optim.parameter4)
        )
        new_params <- sol$par
      } else if (outer.optim.method == "BFGS"){
        sol <- optim(par = current_params,
                     fn = function(params) {
                       sum(obj$estimating_equation_i(params)^2)
                     },
                     method = "BFGS",  control = list(maxit=optim.parameter5, reltol=optim.parameter4)
        )
        new_params <- sol$par
      }
      if (any(abs(new_params) > optim.parameter3)) {
        stop("Estimates are either too large or too small, and convergence might not be achieved.")
      }
      param_diff <- abs(new_params - current_params)
      max_param_diff   <- max(param_diff)
      current_params <- new_params

      obj$setInitialCIFs(obj$get_results()$potential.CIFs)
      sol_list[[iteration]] <- sol
      diff_list[[iteration]] <- max_param_diff
    }
  } else if (outcome.type == 'SURVIVAL') {
    current_params <- alpha_beta_0
    while ((iteration < optim.parameter2) & (max_param_diff > optim.parameter1)) {
      iteration <- iteration + 1
      if (outer.optim.method == "nleqslv" | outer.optim.method == "Broyden"){
        sol <- nleqslv(current_params, obj$estimating_equation_s, method="Broyden", control=list(maxit=optim.parameter5, allowSingular=FALSE))
        new_params <- sol$x
      } else if (outer.optim.method == "Newton"){
        sol <- nleqslv(current_params, obj$estimating_equation_s, method="Newton", control=list(maxit=optim.parameter5, allowSingular=FALSE))
        new_params <- sol$x
      } else if (outer.optim.method == "multiroot") {
        sol <- multiroot(obj$estimating_equation_s, start = current_params, maxiter=optim.parameter5, rtol = optim.parameter4)
        new_params <- sol$root
      } else if (outer.optim.method == "optim" | outer.optim.method == "SANN"){
        sol <- optim(par = current_params,
                     fn = function(params) {
                       sum(obj$estimating_equation_s(params)^2)
                     },
                     method = "SANN",  control = list(maxit=optim.parameter5, reltol=optim.parameter4)
        )
        new_params <- sol$par
      } else if (outer.optim.method == "BFGS"){
        sol <- optim(par = current_params,
                     fn = function(params) {
                       sum(obj$estimating_equation_s(params)^2)
                     },
                     method = "BFGS",  control = list(maxit=optim.parameter5, reltol=optim.parameter4)
        )
        new_params <- sol$par
      }
      if (any(abs(new_params) > optim.parameter3)) {
        stop("Estimates are either too large or too small, and convergence might not be achieved.")
      }
      param_diff <- abs(new_params - current_params)
      max_param_diff   <- max(param_diff)
      current_params <- new_params

      obj$setInitialCIFs(obj$get_results()$potential.CIFs)
      sol_list[[iteration]] <- sol
      diff_list[[iteration]] <- max_param_diff
    }
  } else if (outcome.type == 'PROPORTIONAL') {
    current_params <- alpha_beta_0
    while ((iteration < optim.parameter2) & (max_param_diff > optim.parameter1)) {
      iteration <- iteration + 1
      if (outer.optim.method == "nleqslv" | outer.optim.method == "Broyden"){
        sol <- nleqslv(current_params, obj$estimating_equation_p, method="Broyden", control=list(maxit=optim.parameter5, allowSingular=FALSE))
        new_params <- sol$x
      } else if (outer.optim.method == "Newton"){
        sol <- nleqslv(current_params, obj$estimating_equation_p, method="Newton", control=list(maxit=optim.parameter5, allowSingular=FALSE))
        new_params <- sol$x
      } else if (outer.optim.method == "multiroot") {
        sol <- multiroot(obj$estimating_equation_p, start = current_params, maxiter=optim.parameter5, rtol=optim.parameter4)
        new_params <- sol$root
      } else if (outer.optim.method == "optim" | outer.optim.method == "SANN"){
        sol <- optim(par = current_params,
                     fn = function(params) {
                       sum(obj$estimating_equation_p(params)^2)
                     },
                     method = "SANN",  control = list(maxit=optim.parameter5, reltol=optim.parameter4)
        )
        new_params <- sol$par
      } else if (outer.optim.method == "BFGS"){
        sol <- optim(par = current_params,
                     fn = function(params) {
                       sum(obj$estimating_equation_p(params)^2)
                     },
                     method = "BFGS",  control = list(maxit=optim.parameter5, reltol=optim.parameter4)
        )
        new_params <- sol$par
      }
      if (any(abs(new_params) > optim.parameter3)) {
        stop("Estimates are either too large or too small, and convergence might not be achieved.")
      }
      param_diff <- abs(new_params - current_params)
      max_param_diff   <- max(param_diff)
      current_params <- new_params

      obj$setInitialCIFs(obj$get_results()$potential.CIFs)
      sol_list[[iteration]] <- sol
      diff_list[[iteration]] <- max_param_diff
    }
  }

  #######################################################################################################
  # 5. Calculating variance (functions: calculateCov, calculateCovSurvival)
  #######################################################################################################
  objget_results <- obj$get_results()
  normalizeEstimate <- function(out_calculateCov, should.normalize.covariate, current_params, out_normalizeCovariate) {
    if (should.normalize.covariate == FALSE & !is.null(out_calculateCov$cov_estimated)) {
      alpha_beta_estimated <- current_params
      cov_estimated <- out_calculateCov$cov_estimated
    } else if (should.normalize.covariate == FALSE & is.null(out_calculateCov$cov_estimated)) {
      alpha_beta_estimated <- current_params
      cov_estimated <- NULL
    } else if (should.normalize.covariate == TRUE & is.null(out_calculateCov$cov_estimated)) {
      adj <- 1 / as.vector(out_normalizeCovariate$range)
      alpha_beta_estimated <- adj * current_params
      cov_estimated <- NULL
    } else {
      adj <- 1 / as.vector(out_normalizeCovariate$range)
      alpha_beta_estimated <- adj * current_params
      adj_matrix <- diag(adj, length(adj), length(adj))
      cov_estimated <- adj_matrix %*% out_calculateCov$cov_estimated %*% adj_matrix
    }
    return(list(alpha_beta_estimated = alpha_beta_estimated, cov_estimated = cov_estimated))
  }

  if (outcome.type == 'COMPETINGRISK') {
    out_calculateCov <- calculateCov(objget_results, estimand, prob.bound)
    out_normalizeEstimate <- normalizeEstimate(out_calculateCov, should.normalize.covariate, current_params, out_normalizeCovariate)
  } else if (outcome.type == 'SURVIVAL') {
    out_calculateCov <- calculateCovSurvival(objget_results, estimand, prob.bound)
    out_normalizeEstimate <- normalizeEstimate(out_calculateCov, should.normalize.covariate, current_params, out_normalizeCovariate)
  } else if (outcome.type == 'PROPORTIONAL') {
    out_calculateCov <- NULL
    out_normalizeEstimate <- normalizeEstimate(out_calculateCov, should.normalize.covariate, current_params, out_normalizeCovariate)
  }
  alpha_beta_estimated <- out_normalizeEstimate$alpha_beta_estimated
  cov_estimated <- out_normalizeEstimate$cov_estimated

  #######################################################################################################
  # 6. Output (functions: reportSurvival, reportCompetingRisk, reportPrediction)
  #######################################################################################################
  reportSummary <- list(
    SURVIVAL = reportSurvival,
    COMPETINGRISK = reportCompetingRisk
  )
  if (outcome.type %in% names(reportSummary)) {
    out_summary <- reportSummary[[outcome.type]](
      nuisance.model, exposure, estimand, alpha_beta_estimated,
      cov_estimated, objget_results, iteration, max_param_diff, sol,
      conf.level, optim.method$outer.optim.method
    )
    #    out_prediction <- reportPrediction(nuisance.model,data,exposure,alpha_beta_estimated,cov_estimated,outcome.type,estimand,optim.method,prob.bound)
  }
  reportSummaryFull <- list(
    SURVIVAL = reportSurvivalFull,
    COMPETINGRISK = reportCompetingRiskFull
  )
  if (outcome.type %in% names(reportSummaryFull)) {
    out_summary_full <- reportSummaryFull[[outcome.type]](
      nuisance.model, exposure, estimand, alpha_beta_estimated,
      cov_estimated, objget_results, iteration, max_param_diff, sol,
      conf.level, optim.method$outer.optim.method
    )
    #    out_prediction <- reportPrediction(nuisance.model,data,exposure,alpha_beta_estimated,cov_estimated,outcome.type,estimand,optim.method,prob.bound)
  }
  #  if (outcome.type == 'PROPORTIONAL') {
  #    out_summary <- NULL
  #    out_prediction <- NULL
  #    sorted_data$ip.weight <- objget_results$ip.weight
  #  } else {
  #    sorted_data$influence.function <- out_calculateCov$influence.function
  #    sorted_data$ip.weight <- objget_results$ip.weight
  #  }
  #  out <- list(summary = out_summary, prediction = out_prediction, coefficient=alpha_beta_estimated, coefficient.before.normalization=current_params, cov=cov_estimated, ipw.influence=sorted_data)
  out <- list(summary = out_summary, summary.full = out_summary_full, coefficient=alpha_beta_estimated, cov=cov_estimated)
  return(out)
}
