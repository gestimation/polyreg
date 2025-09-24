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
#' @param outcome.type character Specifies the type of outcome (COMPETING-RISK, SURVIVAL, BINOMIAL, PROPORTIONAL and POLY-PROPORTIONAL).
#' @param conf.level numeric The level for confidence intervals.
#' @param report.nuisance.parameter logical Specifies contents of return. (TRUE = report estimates of nuisance parameters, FALSE = otherwise). Defaults to FALSE.
#' @param report.optim.convergence logical Specifies contents of return. (TRUE = report indicators of convergence of parameter estimation, FALSE = otherwise). Defaults to FALSE.
#' @param report.boot.conf logical Specifies contents of return. (TRUE = report bootstrap confidence intervals, FALSE = otherwise). Defaults to FALSE.
#' @param boot.bca logical Specifies the method of bootstrap confidence intervals (TRUE = BCA method, FALSE = normal approximation).
#' @param boot.parameter1 integer Number of replications for bootstrap confidence intervals. Defaults to 200.
#' @param boot.parameter2 numeric Seed used for bootstrap confidence intervals.
#' @param nleqslv.method character Specifies the method of optimization (nleqslv, Broyden, Newton, optim, BFGS, SANN).
#' @param optim.parameter1 numeric A threshold for determining convergence in outer loop. Defaults to 1e-5.
#' @param optim.parameter2 numeric A threshold for determining convergence in outer loop. Defaults to 1e-5.
#' @param optim.parameter3 numeric Constraint range for parameters. Defaults to 100.
#' @param optim.parameter4 integer Maximum number of iterations in outer loop. Defaults to 50.
#' @param optim.parameter5 integer Maximum number of iterations for nleqslv in outer loop. Defaults to 50.
#' @param optim.parameter6 integer Maximum number of iterations for Levenberg-Marquardt routine. Defaults to 50.
#' @param optim.parameter7 numeric A threshold for determining convergence in Levenberg-Marquardt routine. Defaults to 1e-10.
#' @param optim.parameter8 numeric A threshold for determining convergence in Levenberg-Marquardt routine. Defaults to 1e-5.
#' @param optim.parameter9 numeric Lambda in Levenberg-Marquardt routine. Defaults to 1e-5.
#' @param optim.parameter10 numeric Maximum for lambda in Levenberg-Marquardt routine. Defaults to 40.
#' @param optim.parameter11 numeric Minimum for lambda in Levenberg-Marquardt routine. Defaults to 0.025.
#' @param optim.parameter12 numeric Increment for lambda in Levenberg-Marquardt routine. Defaults to 2.
#' @param optim.parameter13 numeric Decrement for lambda in Levenberg-Marquardt routine. Defaults to 0.5.
#' @param data.initlal.values data.frame A dataset containing initial values. Defaults to NULL.
#' @param should.normalize.covariate logical Indicates whether covariates are normalized (TRUE = normalize, FALSE = otherwise). Defaults to TRUE.
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
    effect.measure1 = "RR",
    effect.measure2 = "RR",
    time.point = NULL,
    outcome.type = "COMPETING-RISK",
    conf.level = 0.95,
    report.nuisance.parameter = FALSE,
    report.optim.convergence = FALSE,
    report.boot.conf = NULL,
    boot.bca = TRUE,
    boot.parameter1 = 200,
    boot.parameter2 = 46,
    nleqslv.method = "nleqslv",
    optim.parameter1 = 1e-6,
    optim.parameter2 = 1e-6,
    optim.parameter3 = 100,
    optim.parameter4 = 50,
    optim.parameter5 = 50,
    optim.parameter6 = 50,
    optim.parameter7 = 1e-10,
    optim.parameter8 = 1e-6,
    optim.parameter9 = 1e-6,
    optim.parameter10 = 40,
    optim.parameter11 = 0.025,
    optim.parameter12 = 2,
    optim.parameter13 = 0.5,
    data.initial.values = NULL,
    should.normalize.covariate = TRUE,
    prob.bound = 1e-5
) {

  #######################################################################################################
  # 1. Pre-processing (function: checkSpell, checkInput, normalizeCovariate, sortByCovariate)
  #######################################################################################################
  computation.time0 <- proc.time()
  cs <- checkSpell(outcome.type, effect.measure1, effect.measure2)
  outcome.type <- cs$outcome.type
  ci <- checkInput(data, nuisance.model, exposure, code.event1, code.event2, code.censoring, code.exposure.ref, outcome.type, conf.level, report.boot.conf, nleqslv.method)
  report.boot.conf <- ci$report.boot.conf

  data <- createAnalysisDataset(formula=nuisance.model, data=data, other.variables.analyzed=c(exposure, strata), subset.condition=subset.condition, na.action=na.action)
  out_normalizeCovariate <- normalizeCovariate(nuisance.model, data, should.normalize.covariate, outcome.type)
  normalized_data <- out_normalizeCovariate$normalized_data
  tp <- read_time.point(nuisance.model, normalized_data, outcome.type, code.censoring, time.point)

  estimand <- list(
    effect.measure1=cs$effect.measure1,
    effect.measure2=cs$effect.measure2,
    time.point=tp,
    code.event1=code.event1,
    code.event2=code.event2,
    code.censoring=code.censoring,
    code.exposure.ref=code.exposure.ref,
    exposure.levels=ci$exposure.levels,
    index.vector=ci$index.vector
  )
  optim.method <- list(
    nleqslv.method = nleqslv.method,
    optim.parameter1 = optim.parameter1,
    optim.parameter2 = optim.parameter2,
    optim.parameter3 = optim.parameter3,
    optim.parameter4 = optim.parameter4,
    optim.parameter5 = optim.parameter5,
    optim.parameter6 = optim.parameter6,
    optim.parameter7 = optim.parameter7,
    optim.parameter8 = optim.parameter8,
    optim.parameter9 = optim.parameter9,
    optim.parameter10 = optim.parameter10,
    optim.parameter11 = optim.parameter11,
    optim.parameter12 = optim.parameter12,
    optim.parameter13 = optim.parameter13
  )

  #######################################################################################################
  # 2. Pre-processing and Calculating initial values alpha_beta_0 (function: calculateInitialValues)
  #######################################################################################################
  if (outcome.type == "COMPETING-RISK" | outcome.type == "SURVIVAL" | outcome.type == "BINOMIAL") {
    alpha_beta_0 <- getInitialValues(
      formula = nuisance.model,
      data = normalized_data,
      exposure = exposure,
      data.initial.values = data.initial.values,
      estimand = estimand,
      specific.time = estimand$time.point,
      outcome.type = outcome.type,
      prob.bound = prob.bound
    )
  } else if (outcome.type == "PROPORTIONAL" || outcome.type == "POLY-PROPORTIONAL") {
    alpha_beta_0 <- getInitialValuesProportional(
      formula = nuisance.model,
      data = normalized_data,
      outcome.type = outcome.type,
      exposure = exposure,
      estimand = estimand,
      data.initial.values = data.initial.values,
      prob.bound = prob.bound,
      out_normalizeCovariate = out_normalizeCovariate
    )
  }
  alpha_beta_0 <- rep(0, estimand$index.vector[8])

  #######################################################################################################
  # 3. Calculating IPCW (function: calculateIPCW, calculateIPCWMatrix)
  #######################################################################################################
  if (outcome.type == "COMPETING-RISK" | outcome.type == "SURVIVAL") {
    ip.weight <- calculateIPCW(nuisance.model, normalized_data, code.censoring, strata, estimand$time.point)
  } else if (outcome.type == "BINOMIAL") {
    ip.weight <- rep(1,nrow(normalized_data))
  } else if (outcome.type == "PROPORTIONAL" | outcome.type == "POLY-PROPORTIONAL") {
    ip.weight.matrix <- calculateIPCWMatrix(nuisance.model, normalized_data, code.censoring, strata, estimand, out_normalizeCovariate)
  }

  #######################################################################################################
  # 4. Parameter estimation (functions: estimating_equation_ipcw, _survival, _proportional)
  #######################################################################################################
  makeObjectiveFunction <- function() {
    out_ipcw <- list()
    initial.CIFs <- NULL
      call_and_capture <- function(fun, ...) {
      out_ipcw <<- do.call(fun, list(...))
      out_ipcw$ret
      }
      estimating_equation_i <- function(p) call_and_capture(
        estimating_equation_ipcw,
        formula = nuisance.model, data = normalized_data, exposure = exposure,
        ip.weight = ip.weight, alpha_beta = p, estimand = estimand,
        optim.method = optim.method, prob.bound = prob.bound,
        initial.CIFs = initial.CIFs
      )
      estimating_equation_s <- function(p) call_and_capture(
      estimating_equation_survival,
      formula = nuisance.model, data = normalized_data, exposure = exposure,
      ip.weight = ip.weight, alpha_beta = p, estimand = estimand,
      prob.bound = prob.bound, initial.CIFs = initial.CIFs
    )
      estimating_equation_p <- function(p) call_and_capture(
      estimating_equation_proportional,
      formula = nuisance.model, data = normalized_data, exposure = exposure,
      ip.weight.matrix = ip.weight.matrix, alpha_beta = p, estimand = estimand,
      optim.method = optim.method, prob.bound = prob.bound,
      initial.CIFs = initial.CIFs
    )
      estimating_equation_pp <- function(p) call_and_capture(
      estimating_equation_pproportional,
      formula = nuisance.model, data = normalized_data, exposure = exposure,
      ip.weight.matrix = ip.weight.matrix, alpha_beta = p, estimand = estimand,
      optim.method = optim.method, prob.bound = prob.bound,
      initial.CIFs = initial.CIFs
    )
    setInitialCIFs <- function(new.CIFs) initial.CIFs <<- new.CIFs
    getResults     <- function() out_ipcw
      list(
      estimating_equation_i = estimating_equation_i,
      estimating_equation_s = estimating_equation_s,
      estimating_equation_p = estimating_equation_p,
      estimating_equation_pp = estimating_equation_pp,
      setInitialCIFs = setInitialCIFs,
      getResults = getResults
    )
  }

  assessConvergence <- function(new_params, current_params, current_obj_value, optim.parameter1, optim.parameter2, optim.parameter3) {
    if (any(abs(new_params) > optim.parameter3)) {
      stop("Estimates are either too large or too small, and convergence might not be achieved.")
    }
    param_diff <- abs(new_params - current_params)
    max.absolute.difference <- max(param_diff)
    relative.difference <- assessRelativeDifference(new_params, current_params)
    obj_value <- get_obj_value(new_params)
    converged <- (relative.difference <= optim.parameter1) || (obj_value <= optim.parameter2) || is_stalled(c(current_obj_value, obj_value))

    criteria1 <- (relative.difference <= optim.parameter1)
    criteria2 <- (obj_value <= optim.parameter2)
    criteria3 <- is_stalled(c(current_obj_value, obj_value))
    converged  <- (criteria1 || criteria2 || criteria3)
    converged.by <- if (!converged) NA_character_
    else if (criteria1)   "Converged in relative difference"
    else if (criteria2) "Converged in objective function"
    else "Stalled"

    list(converged = converged, converged.by=converged.by, relative.difference = relative.difference, max.absolute.difference = max.absolute.difference, obj_value = obj_value)
  }

  assessRelativeDifference <- function(new, old) {
    max(abs(new - old) / pmax(1, abs(old)))
  }

  is_stalled <- function(x, stall_patience=3, stall_eps=1e-3) {
    if (length(x) < stall_patience) return(FALSE)
    recent <- tail(x, stall_patience)
    (diff(range(recent)) / max(1e-12, mean(recent))) <= stall_eps
  }

  choose_estimating_equation <- function(outcome.type, obj) {
    if (outcome.type == "COMPETING-RISK") {
      obj$estimating_equation_i
    } else if (outcome.type == "SURVIVAL" | outcome.type == "BINOMIAL") {
      obj$estimating_equation_s
    } else if (outcome.type == "PROPORTIONAL") {
      obj$estimating_equation_p
    } else if (outcome.type == "POLY-PROPORTIONAL") {
      obj$estimating_equation_pp
    } else {
      stop("Unknown outcome.type: ", outcome.type)
    }
  }

  choose_nleqslv_method <- function(nleqslv.method) {
    if (nleqslv.method %in% c("nleqslv", "Broyden")) {
      "Broyden"
    } else if (nleqslv.method == "Newton") {
      "Newton"
    } else {
      stop("Unsupported nleqslv.method without optim(): ", nleqslv.method)
    }
  }

  get_obj_value <- switch(outcome.type,
                          "COMPETING-RISK"   = function(p) drop(crossprod(obj$estimating_equation_i(p))),
                          "SURVIVAL"         = function(p) drop(crossprod(obj$estimating_equation_s(p))),
                          "BINOMIAL"         = function(p) drop(crossprod(obj$estimating_equation_s(p))),
                          "PROPORTIONAL"     = function(p) drop(crossprod(obj$estimating_equation_p(p))),
                          "POLY-PROPORTIONAL"= function(p) drop(crossprod(obj$estimating_equation_pp(p))),
                          stop("Unknown outcome.type: ", outcome.type)
  )

  obj <- makeObjectiveFunction()
  estimating_fun <- choose_estimating_equation(outcome.type, obj)
  nleqslv_method  <- choose_nleqslv_method(nleqslv.method)
  iteration <- 0L
  max.absolute.difference <- Inf
  out_nleqslv <- NULL
  current_params <- alpha_beta_0
  current_obj_value <- numeric(0)
  trace_df  <- NULL
  store_params <- TRUE

  while ((iteration < optim.parameter4) & (max.absolute.difference > optim.parameter1)) {
    iteration <- iteration + 1
    prev_params <- current_params

    out_nleqslv <- nleqslv(
      prev_params,
      estimating_fun,
      method  = nleqslv_method,
      control = list(maxit = optim.parameter5, allowSingular = FALSE)
    )
    new_params <- out_nleqslv$x

    current_obj_value <- get_obj_value(new_params)
    obj$setInitialCIFs(obj$getResults()$potential.CIFs)
    ac <- assessConvergence(new_params, prev_params, current_obj_value,optim.parameter1, optim.parameter2, optim.parameter3)

    nleqslv.info <- extractOptimizationInfo(out_nleqslv, nleqslv.method)
    computation.time.second <- as.numeric((proc.time() - computation.time0)[3])

    trace_df <- append_trace(
      trace_df,
      iteration = iteration,
      computation.time.second = computation.time.second,
      nleqslv.method = nleqslv.method,
      nleqslv.info = nleqslv.info,
      objective.function = ac$obj_value,
      relative.difference = ac$relative.difference,
      max.absolute.difference = ac$max.absolute.difference,
      converged.by = if (ac$converged) ac$converged.by else FALSE,
      coefficient = if (store_params) new_params else NULL
    )

    current_params <- new_params
    max.absolute.difference <- ac$max.absolute.difference
    if (ac$converged) break
  }
  out_getResults <- obj$getResults()

  #######################################################################################################
  # 5. Calculating variance (functions: calculateCov, calculateCovSurvival)
  #######################################################################################################
  normalizeEstimate <- function(
    outcome.type,
    report.boot.conf,
    should.normalize.covariate,
    current_params,
    out_getResults,
    estimand,
    prob.bound,
    out_normalizeCovariate
  ) {
    exclude_cov <- report.boot.conf == TRUE || outcome.type %in% c("PROPORTIONAL", "POLY-PROPORTIONAL")

    if (exclude_cov) {
      alpha_beta_estimated <- if (should.normalize.covariate) {
        (1 / as.vector(out_normalizeCovariate$range)) * current_params
      } else {
        current_params
      }
      return(list(alpha_beta_estimated = alpha_beta_estimated, cov_estimated = NULL))
    }
    if (should.normalize.covariate) {
      adj <- 1 / as.vector(out_normalizeCovariate$range)
      if (length(adj) != length(current_params))stop("Length of adj (range) must match length of current_params.")
      alpha_beta_estimated <- adj * current_params
      adj_matrix <- diag(adj, length(adj))
      cov_estimated <- adj_matrix %*% out_calculateCov$cov_estimated %*% adj_matrix
      #cov_estimated <- NULL
    } else {
      alpha_beta_estimated <- current_params
      cov_estimated <- out_calculateCov$cov_estimated
    }
    return(list(alpha_beta_estimated = alpha_beta_estimated, cov_estimated = cov_estimated))
  }

  out_calculateCov <- switch(
    outcome.type,
    "COMPETING-RISK" = calculateCov(out_getResults, estimand, prob.bound),
    #    "COMPETING-RISK" = NULL,
    "SURVIVAL" = calculateCovSurvival(out_getResults, estimand, prob.bound),
    "BINOMIAL" = calculateCovSurvival(out_getResults, estimand, prob.bound),
    "PROPORTIONAL" = NULL,
    "POLY-PROPORTIONAL" = NULL,
    stop(sprintf("Unsupported outcome.type for covariance: %s", outcome.type))
  )
  out_normalizeEstimate <- normalizeEstimate(
    outcome.type = outcome.type,
    report.boot.conf = report.boot.conf,
    should.normalize.covariate = should.normalize.covariate,
    current_params = current_params,
    out_getResults = out_getResults,
    estimand = estimand,
    prob.bound = prob.bound,
    out_normalizeCovariate = out_normalizeCovariate
  )

  alpha_beta_estimated <- out_normalizeEstimate$alpha_beta_estimated
  cov_estimated        <- out_normalizeEstimate$cov_estimated

  #######################################################################################################
  # 6. Calculating bootstrap confidence interval (functions: boot, solveEstimatingEquationP)
  #######################################################################################################
  if (isTRUE(report.boot.conf)) {
    set.seed(boot.parameter2)
    boot.coef     <- rep(NA,2)
    boot.coef_se  <- rep(NA,2)
    boot.p_value  <- rep(NA,2)
    boot.conf_low <- rep(NA,2)
    boot.conf_high<- rep(NA,2)

    if (outcome.type=="POLY-PROPORTIONAL") {
      index_coef    <- c(length(estimand$time.point) + 1, 2*length(estimand$time.point) + 2)
    } else if (outcome.type=="PROPORTIONAL") {
      index_coef    <- c(length(estimand$time.point) + 1)
    } else if (outcome.type=="COMPETING-RISK") {
      index_coef <- seq_len(2*out_normalizeCovariate$n_covariate+4)
    } else if (outcome.type=="BINOMIAL" | outcome.type=="SURVIVAL") {
      index_coef <- seq_len(out_normalizeCovariate$n_covariate+2)
    }
    boot_function <- function(data, indices) {
      coef <- solveEstimatingEquation(nuisance.model=nuisance.model, exposure=exposure, strata=strata,
                                      normalized_data = data[indices, , drop = FALSE], outcome.type=outcome.type, estimand=estimand, optim.method=optim.method, out_normalizeCovariate=out_normalizeCovariate, prob.bound=prob.bound, alpha_beta_0=alpha_beta_0)
      return(coef)
    }

    out_boot <- boot(normalized_data, boot_function, R = boot.parameter1)
    for (j in seq_len(length(index_coef))) {
      if (isTRUE(boot.bca)) {
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
    out_bootstrap <- list(
      boot.coef=boot.coef, boot.coef_se=boot.coef_se, boot.p_value=boot.p_value, boot.conf_low=boot.conf_low, boot.conf_high=boot.conf_high
    )
  } else {out_bootstrap <- NULL}

  #######################################################################################################
  # 7. Output (functions: reportSurvival, reportCOMPETING-RISK, reportPrediction)
  #######################################################################################################
  if (outcome.type == "PROPORTIONAL" | outcome.type == "POLY-PROPORTIONAL") {
    out_summary <- reportConstantEffects(
      nuisance.model, exposure, estimand, out_bootstrap, out_getResults, iteration, max.absolute.difference, out_nleqslv, optim.method$nleqslv.method
    )
    out_data <- NULL
  } else {
    out_summary <- reportEffects (
      outcome.type, report.nuisance.parameter, report.optim.convergence, report.boot.conf, nuisance.model, exposure, estimand, alpha_beta_estimated,
      cov_estimated, out_bootstrap, out_getResults, iteration, max.absolute.difference, out_nleqslv,
      conf.level, optim.method$nleqslv.method
    )
    normalized_data$influence.function <- out_calculateCov$influence.function
    normalized_data$ip.weight <- out_getResults$ip.weight
    normalized_data$potential.CIFs <- out_getResults$potential.CIFs
    out_data <- normalized_data
  }
  out <- list(summary=out_summary, coefficient=alpha_beta_estimated, cov=cov_estimated, bootstrap=out_bootstrap, diagnosis.statistics=out_data, optimization.info=trace_df)
  return(out)
}
