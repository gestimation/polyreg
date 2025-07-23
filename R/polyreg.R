#' @title Direct polynomial regression for competing risks, survival and binomial analysis
#'
#' @param nuisance.model formula Model formula representing outcome and covariates other than exposure
#' @param exposure character Column name representing the binary exposure variable.
#' @param strata character Column name representing the stratification variable for adjustment for dependent censoring. Defaults to NULL.
#' @param data data.frame Input dataset containing the outcome, the exposure and covariates.
#' @param subset.condition character Specifies a condition for subsetting the data. Defaults to NULL.
#' @param na.action character Specifies a missing-data filter function, applied to the model frame, after any subset argument has been used. Defaults to na.omit.
#' @param code.event1 integer Specifies the code of event 1. Defaults to 1.
#' @param code.event2 integer Specifies the code of event 2. Defaults to 2.
#' @param code.censoring integer Specifies the code of censoring. Defaults to 0.
#' @param code.exposure.ref integer Specifies the code of the reference category of exposure. Defaults to 0.
#' @param effect.measure1 character Specifies the effect measure for event (RR, OR, SHR).
#' @param effect.measure2 character Specifies the effect measure for competing risk (RR, OR, SHR).
#' @param time.point numeric The time point for exposure effects to be estimated.
#' @param outcome.type character Specifies the type of outcome (COMPETINGRISK, SURVIVAL, BINOMIAL, PROPORTIONAL and POLY-PROPORTIONAL).
#' @param conf.level numeric The level for confidence intervals.
#' @param report.nuisance.parameter logical Specifies contents of return. (TRUE = report estimates of nuisance parameters, FALSE = otherwise). Defaults to FALSE.
#' @param report.optim.convergence logical Specifies contents of return. (TRUE = report indicators of convergence of parameter estimation, FALSE = otherwise). Defaults to FALSE.
#' @param report.boot.conf logical Specifies contents of return. (TRUE = report bootstrap confidence intervals, FALSE = otherwise). Defaults to FALSE.
#' @param boot.bca logical Specifies the method of bootstrap confidence intervals (TRUE = BCA method, FALSE = normal approximation).
#' @param boot.parameter1 integer Number of replications for bootstrap confidence intervals.
#' @param boot.parameter2 numeric Seed used for bootstrap confidence intervals.
#' @param outer.optim.method character Specifies the method of optimization (nleqslv, Broyden, Newton, optim, BFGS, SANN).
#' @param inner.optim.method character Specifies the method of optimization (optim, BFGS, SANN, roptim).
#' @param optim.parameter1 numeric A threshold for initial value search in outer loop. Defaults to 1e-5.
#' @param optim.parameter2 integer Maximum number of iterations. Defaults to 20.
#' @param optim.parameter3 numeric Constraint range for parameters. Defaults to 100.
#' @param optim.parameter4 numeric A threshold for determining convergence in outer loop. Defaults to 1e-5.
#' @param optim.parameter5 integer Maximum number of iterations for nleqslv/optim in outer loop. Defaults to 200.
#' @param optim.parameter6 numeric A threshold for determining convergence in inner loop. Defaults to 1e-10.
#' @param optim.parameter7 integer Maximum number of iterations for optim in inner loop. Defaults to 200.
#' @param data.initlal.values data.frame A dataset containing initial values. Defaults to NULL.
#' @param should.normalize.covariate logical Indicates whether covariates are normalized (TRUE = normalize, FALSE = otherwise). Defaults to TRUE.
#' @param should.sort.data logical Indicates whether data are initially sorted to reduce computation steps (TRUE = sort, FALSE = otherwise). Defaults to TRUE.
#' @param prob.bound numeric A threshold for clamping probabilities. Defaults to 1e-5.
#' @importFrom nleqslv nleqslv
#' @importFrom boot boot boot.ci
#' @importFrom Rcpp sourceCpp
#' @useDynLib polyreg, .registration = TRUE
#'
#' @return A list of results from direct polynomial regression. coefficient and cov are estimated regression coefficients of exposure and covariates and their variance covariance matrix. summary meets requirement of msummary function. diagnosis.statistics in includes inverse probability weights, influence functions, and predicted potential outcomes
#' @export polyreg
#'
#' @examples
#' data(diabetes.complications)
#' output <- polyreg(nuisance.model = Event(t,epsilon)~+1, exposure = 'fruitq1', data = diabetes.complications, effect.measure1='RR', effect.measure2='RR', time.point=8, outcome.type='C')
#' library(modelsummary)
#' msummary(output$summary, statistic = c("conf.int"), exponentiate = TRUE)
polyreg <- function(
    nuisance.model,
    exposure,
    strata = NULL,
    data,
    subset.condition = NULL,
    na.action = na.omit,
    code.event1 = 1,
    code.event2 = 2,
    code.censoring = 0,
    code.exposure.ref = 0,
    effect.measure1 = 'RR',
    effect.measure2 = 'RR',
    time.point = NULL,
    outcome.type = 'COMPETINGRISK',
    conf.level = 0.95,
    report.nuisance.parameter = FALSE,
    report.optim.convergence = FALSE,
    report.boot.conf = FALSE,
    boot.bca = TRUE,
    boot.parameter1 = 200,
    boot.parameter2 = 46,
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
  # 1. Pre-processing (function: checkSpell, checkInput, normalizeCovariate, sortByCovariate)
  #######################################################################################################
  checkSpell(outcome.type, effect.measure1, effect.measure2)
  checkInput(outcome.type, time.point, conf.level, outer.optim.method, inner.optim.method)
  outcome.type <- outcome.type.corrected
  effect.measure1 <- effect.measure1.corrected
  effect.measure2 <- effect.measure2.corrected
  time.point <- time.point.corrected

  estimand <- list(
    effect.measure1=effect.measure1,
    effect.measure2=effect.measure2,
    time.point=time.point,
    code.event1=code.event1,
    code.event2=code.event2,
    code.censoring=code.censoring,
    code.exposure.ref=code.exposure.ref
  )
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

  data <- createAnalysisDataset(formula=nuisance.model, data=data, other.variables.analyzed=c(exposure, strata), subset.condition=subset.condition, na.action=na.action)
#  if (!is.null(subset)) {
#    subset_condition <- subset
#    data <- subset(data, eval(parse(text = subset_condition)))
#  }
#  data <- na.action(data)
  out_normalizeCovariate <- normalizeCovariate(nuisance.model, data, should.normalize.covariate, outcome.type)
  normalized_data <- out_normalizeCovariate$normalized_data
  sorted_data <- sortByCovariate(nuisance.model, normalized_data, should.sort.data, out_normalizeCovariate$n_covariate)

  #######################################################################################################
  # 2. Pre-processing and Calculating initial values alpha_beta_0 (function: calculateInitialValues)
  #######################################################################################################
  if (outcome.type == 'COMPETINGRISK' | outcome.type == 'SURVIVAL' | outcome.type == 'BINOMIAL') {
    alpha_beta_0 <- getInitialValues(
      formula = nuisance.model,
      data = sorted_data,
      exposure = exposure,
      data.initial.values = data.initial.values,
      estimand = estimand,
      specific.time = estimand$time.point,
      outcome.type = outcome.type,
      prob.bound = prob.bound
    )
  } else if (outcome.type == 'PROPORTIONAL') {
    n_para_1 <- out_normalizeCovariate$n_covariate+1
    n_para_2 <- out_normalizeCovariate$n_covariate+2
    n_para_3 <- out_normalizeCovariate$n_covariate+3
    n_para_4 <- 2*out_normalizeCovariate$n_covariate+3
    n_para_5 <- 2*out_normalizeCovariate$n_covariate+4
    n_para_6 <- length(time.point)*(n_para_5-2) + 2
#    alpha_beta_0 <- rep(NA, n_para_6)
    alpha_beta_0 <- rep(NA, n_para_6/2) # initial parameters for event 1 over time points
    i_time <- sum1 <- sum2 <- 0
    for (specific.time in time.point) {
      specific.time <- as.numeric(specific.time)
      out_getInitialValues <- getInitialValues(
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
      tmp1 <- out_getInitialValues[1:n_para_1]
#      tmp2 <- out_getInitialValues[n_para_3:n_para_4]
      alpha_beta_0[i_para:(i_para+n_para_1-1)]                          <- tmp1
#      alpha_beta_0[(n_para_6/2+i_para):(n_para_6/2+i_para+n_para_1-1)]  <- tmp2
      sum1 <- sum1 + out_getInitialValues[n_para_2]
#      sum2 <- sum2 + out_getInitialValues[n_para_5]
    }
    alpha_beta_0[n_para_6/2]  <- sum1/length(time.point)
#    alpha_beta_0[n_para_6]    <- sum2/length(time.point)
  } else if (outcome.type == 'POLY-PROPORTIONAL') {
    n_para_1 <- out_normalizeCovariate$n_covariate+1
    n_para_2 <- out_normalizeCovariate$n_covariate+2
    n_para_3 <- out_normalizeCovariate$n_covariate+3
    n_para_4 <- 2*out_normalizeCovariate$n_covariate+3
    n_para_5 <- 2*out_normalizeCovariate$n_covariate+4
    n_para_6 <- length(time.point)*(n_para_5-2) + 2
    alpha_beta_0 <- rep(NA, n_para_6) # initial parameters for event 1 and 2 over time points
    i_time <- sum1 <- sum2 <- 0
    for (specific.time in time.point) {
      specific.time <- as.numeric(specific.time)
      out_getInitialValues <- getInitialValues(
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
      tmp1 <- out_getInitialValues[1:n_para_1]
      tmp2 <- out_getInitialValues[n_para_3:n_para_4]
      alpha_beta_0[i_para:(i_para+n_para_1-1)]                          <- tmp1
      alpha_beta_0[(n_para_6/2+i_para):(n_para_6/2+i_para+n_para_1-1)]  <- tmp2
      sum1 <- sum1 + out_getInitialValues[n_para_2]
      sum2 <- sum2 + out_getInitialValues[n_para_5]
    }
    alpha_beta_0[n_para_6/2]  <- sum1/length(time.point)
    alpha_beta_0[n_para_6]    <- sum2/length(time.point)
  }

  #######################################################################################################
  # 3. Calculating IPCW (function: calculateIPCW)
  #######################################################################################################
  if (outcome.type == 'COMPETINGRISK' | outcome.type == 'SURVIVAL') {
    ip.weight <- calculateIPCW(formula = nuisance.model, data = sorted_data, code.censoring=code.censoring, strata=strata, specific.time = estimand$time.point)
  } else if (outcome.type == 'PROPORTIONAL') {
    i_time <- 0
    ip.weight.matrix <- matrix(NA, nrow=out_normalizeCovariate$n, ncol=length(time.point))
    for (specific.time in time.point) {
      i_time <- i_time + 1
      ip.weight.matrix[,i_time] <- calculateIPCW(formula = nuisance.model, data = sorted_data, code.censoring=code.censoring, strata=strata, specific.time = specific.time)
    }
  } else if (outcome.type == 'BINOMIAL') {
    ip.weight <- rep(1,nrow(sorted_data))
  }

  #######################################################################################################
  # 4. Parameter estimation (functions: estimating_equation_ipcw, _survival, _proportional)
  #######################################################################################################
  makeObjectiveFunction <- function() {
    out_ipcw <- list()
    initial.CIFs <- NULL
    estimating_equation_i <- function(p) {
    if (optim.method$inner.optim.method == "roptim") {
      out_ipcw <- estimating_equation_ipcw_roptim(
        formula = nuisance.model,
        data = sorted_data,
        exposure = exposure,
        ip.weight = ip.weight,
        alpha_beta = p,
        estimand = estimand,
        optim.method = optim.method,
        prob.bound = prob.bound,
        initial.CIFs = initial.CIFs)
    } else {
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
    }
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
    estimating_equation_pp <- function(p) {
      out_ipcw <- estimating_equation_pproportional(
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
    getResults <- function() {
      out_ipcw
    }
    list(
      estimating_equation_i = estimating_equation_i,
      estimating_equation_s = estimating_equation_s,
      estimating_equation_p = estimating_equation_p,
      estimating_equation_pp = estimating_equation_pp,
      setInitialCIFs = setInitialCIFs,
      getResults = getResults
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

      obj$setInitialCIFs(obj$getResults()$potential.CIFs)
      sol_list[[iteration]] <- sol
      diff_list[[iteration]] <- max_param_diff
    }
  } else if (outcome.type == 'SURVIVAL' | outcome.type == 'BINOMIAL') {
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
      param_diff     <- abs(new_params - current_params)
      max_param_diff <- max(param_diff)
      current_params <- new_params

      obj$setInitialCIFs(obj$getResults()$potential.CIFs)
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

      obj$setInitialCIFs(obj$getResults()$potential.CIFs)
      sol_list[[iteration]] <- sol
      diff_list[[iteration]] <- max_param_diff
    }
  } else if (outcome.type == 'POLY-PROPORTIONAL') {
    current_params <- alpha_beta_0
    while ((iteration < optim.parameter2) & (max_param_diff > optim.parameter1)) {
      iteration <- iteration + 1
      if (outer.optim.method == "nleqslv" | outer.optim.method == "Broyden"){
        sol <- nleqslv(current_params, obj$estimating_equation_pp, method="Broyden", control=list(maxit=optim.parameter5, allowSingular=FALSE))
        new_params <- sol$x
      } else if (outer.optim.method == "Newton"){
        sol <- nleqslv(current_params, obj$estimating_equation_pp, method="Newton", control=list(maxit=optim.parameter5, allowSingular=FALSE))
        new_params <- sol$x
      } else if (outer.optim.method == "multiroot") {
        sol <- multiroot(obj$estimating_equation_pp, start = current_params, maxiter=optim.parameter5, rtol=optim.parameter4)
        new_params <- sol$root
      } else if (outer.optim.method == "optim" | outer.optim.method == "SANN"){
        sol <- optim(par = current_params,
                     fn = function(params) {
                       sum(obj$estimating_equation_pp(params)^2)
                     },
                     method = "SANN",  control = list(maxit=optim.parameter5, reltol=optim.parameter4)
        )
        new_params <- sol$par
      } else if (outer.optim.method == "BFGS"){
        sol <- optim(par = current_params,
                     fn = function(params) {
                       sum(obj$estimating_equation_pp(params)^2)
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

      obj$setInitialCIFs(obj$getResults()$potential.CIFs)
      sol_list[[iteration]] <- sol
      diff_list[[iteration]] <- max_param_diff
    }
  }
  out_getResults <- obj$getResults()

  #######################################################################################################
  # 5. Calculating variance (functions: calculateCov, calculateCovSurvival)
  #######################################################################################################
  normalizeEstimate <- function(out_calculateCov, should.normalize.covariate, report.boot.conf, current_params, out_normalizeCovariate) {
    if (should.normalize.covariate == TRUE & report.boot.conf == FALSE & !is.null(out_calculateCov$cov_estimated)) {
        adj <- 1 / as.vector(out_normalizeCovariate$range)
      alpha_beta_estimated <- adj * current_params
      adj_matrix <- diag(adj, length(adj), length(adj))
      cov_estimated <- adj_matrix %*% out_calculateCov$cov_estimated %*% adj_matrix
    } else if (should.normalize.covariate == FALSE & report.boot.conf == FALSE & !is.null(out_calculateCov$cov_estimated)) {
      alpha_beta_estimated <- current_params
      cov_estimated <- out_calculateCov$cov_estimated
    } else if (should.normalize.covariate == TRUE) {
      adj <- 1 / as.vector(out_normalizeCovariate$range)
      alpha_beta_estimated <- adj * current_params
      cov_estimated <- NULL
    } else {
        alpha_beta_estimated <- current_params
        cov_estimated <- NULL
    }
      return(list(alpha_beta_estimated = alpha_beta_estimated, cov_estimated = cov_estimated))
  }

  if (outcome.type == 'COMPETINGRISK' & report.boot.conf == FALSE) {
    out_calculateCov <- calculateCov(out_getResults, estimand, prob.bound)
    out_normalizeEstimate <- normalizeEstimate(out_calculateCov, should.normalize.covariate, report.boot.conf, current_params, out_normalizeCovariate)
  } else if ((outcome.type == 'SURVIVAL' | outcome.type == 'BINOMIAL') & report.boot.conf == FALSE) {
    out_calculateCov <- calculateCovSurvival(out_getResults, estimand, prob.bound)
    out_normalizeEstimate <- normalizeEstimate(out_calculateCov, should.normalize.covariate, report.boot.conf, current_params, out_normalizeCovariate)
  } else if (outcome.type == 'PROPORTIONAL' | report.boot.conf == TRUE) {
    out_calculateCov <- NULL
    out_normalizeEstimate <- normalizeEstimate(out_calculateCov, should.normalize.covariate, report.boot.conf, current_params, out_normalizeCovariate)
  }
  alpha_beta_estimated <- out_normalizeEstimate$alpha_beta_estimated
  cov_estimated <- out_normalizeEstimate$cov_estimated

  #######################################################################################################
  # 6. Calculating bootstrap confidence interval (functions: boot, solveEstimatingEquationP)
  #######################################################################################################
  boot.coef     <- rep(NA,2)
  boot.coef_se  <- rep(NA,2)
  boot.p_value  <- rep(NA,2)
  boot.conf_low <- rep(NA,2)
  boot.conf_high<- rep(NA,2)
  if (outcome.type=='POLY-PROPORTIONAL') {
    index_coef    <- c(length(time.point) + 1, 2*length(time.point) + 2)
  } else if (outcome.type=='PROPORTIONAL') {
    index_coef    <- c(length(time.point) + 1)
  }

  if (outcome.type=='PROPORTIONAL' | outcome.type=='POLY-PROPORTIONAL') {
    set.seed(boot.parameter2)
    if (outcome.type=='PROPORTIONAL') {
      boot_function <- function(data, indices) {
        return(solveEstimatingEquationP(nuisance.model=nuisance.model, exposure=exposure, strata=strata,
                                       sorted_data=data[indices, ], estimand=estimand, optim.method=optim.method, data.initial.values=alpha_beta_0
        ))
      }
    } else if (outcome.type=='POLY-PROPORTIONAL') {
      boot_function <- function(data, indices) {
        return(solveEstimatingEquationPP(nuisance.model=nuisance.model, exposure=exposure, strata=strata,
                                        sorted_data=data[indices, ], estimand=estimand, optim.method=optim.method, data.initial.values=alpha_beta_0
        ))
      }
    } else if (outcome.type=='COMPETINGRISK') {
      boot_function <- function(data, indices) {
        coef <- solveEstimatingEquationC(nuisance.model=nuisance.model, exposure=exposure, strata=strata,
                                        sorted_data=data[indices, ], estimand=estimand, optim.method=optim.method, data.initial.values=alpha_beta_0)
        return(coef)
      }
    } else if (outcome.type=='SURVIVAL') {
      boot_function <- function(data, indices) {
        return(solveEstimatingEquationS(nuisance.model=nuisance.model, exposure=exposure, strata=strata,
                                        sorted_data=data[indices, ], estimand=estimand, optim.method=optim.method, data.initial.values=alpha_beta_0
        ))
      }
    } else if (outcome.type=='BINOMIAL') {
      boot_function <- function(data, indices) {
        return(solveEstimatingEquationB(nuisance.model=nuisance.model, exposure=exposure, strata=strata,
                                        sorted_data=data[indices, ], estimand=estimand, optim.method=optim.method, data.initial.values=alpha_beta_0
        ))
      }
    }
    out_boot <- boot(sorted_data, boot_function, R = boot.parameter1)
    for (j in seq_len(length(index_coef))) {
      if (boot.bca == TRUE) {
        out_boot.ci <- boot.ci(out_boot, conf = conf.level, index = index_coef[j], type = c("norm", "bca"))
        boot.coef[j] <- (out_boot.ci$normal[2] + out_boot.ci$normal[3])/2
        ci_range <- out_boot.ci$normal[3] - out_boot.ci$normal[2]
        boot.coef_se[j] <- ci_range/2/qnorm(1 - (1-conf.level)/2)
        boot.p_value[j] <- 2 * (1 - pnorm(abs(boot.coef[j]) / boot.coef_se[j]))
        boot.conf_low[j] <- out_boot.ci$bca[4]
        boot.conf_high[j] <- out_boot.ci$bca[5]
      } else {
        out_boot.ci <- boot.ci(out_boot, conf = conf.level, index = index_coef[j], type = c("norm"))
        boot.coef[j] <- (out_boot.ci$normal[2] + out_boot.ci$normal[3])/2
        ci_range <- out_boot.ci$normal[3] - out_boot.ci$normal[2]
        boot.coef_se[j] <- ci_range/2/qnorm(1 - (1-conf.level)/2)
        boot.p_value[j] <- 2 * (1 - pnorm(abs(boot.coef[j]) / boot.coef_se[j]))
        boot.conf_low[j] <- out_boot.ci$normal[2]
        boot.conf_high[j] <- out_boot.ci$normal[3]
      }
    }
  }
  out_bootstrap <- list(
    boot.coef=boot.coef, boot.coef_se=boot.coef_se, boot.p_value=boot.p_value, boot.conf_low=boot.conf_low, boot.conf_high=boot.conf_high
  )

  #######################################################################################################
  # 7. Output (functions: reportSurvival, reportCompetingRisk, reportPrediction)
  #######################################################################################################
  if (outcome.type == 'PROPORTIONAL') {
    out_summary <- reportProportional(
      nuisance.model, exposure, estimand, out_bootstrap,
      out_getResults, iteration, max_param_diff, sol,
      optim.method$outer.optim.method
    )
    out_data <- NULL
  } else {
    out_summary <- reportEffects (
      outcome.type, report.nuisance.parameter, report.optim.convergence, report.boot.conf, nuisance.model, exposure, estimand, alpha_beta_estimated,
      cov_estimated, out_bootstrap, out_getResults, iteration, max_param_diff, sol,
      conf.level, optim.method$outer.optim.method
    )
    sorted_data$influence.function <- out_calculateCov$influence.function
    sorted_data$ip.weight <- out_getResults$ip.weight
    sorted_data$potential.CIFs <- out_getResults$potential.CIFs
    out_data <- sorted_data
  }
  out <- list(summary = out_summary, coefficient=alpha_beta_estimated, cov=cov_estimated, diagnosis.statistics=out_data)
  return(out)
}
