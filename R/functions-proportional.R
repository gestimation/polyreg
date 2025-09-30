solveEstimatingEquation <- function(
    nuisance.model,
    exposure,
    strata,
    normalized_data,
    outcome.type,
    estimand,
    optim.method,
    out_normalizeCovariate,
    prob.bound,
    alpha_beta_0
) {
  if (outcome.type == "COMPETING-RISK" | outcome.type == "SURVIVAL") {
    ip.weight <- calculateIPCW(nuisance.model, normalized_data, estimand$code.censoring, strata, estimand$time.point)
  } else if (outcome.type == "BINOMIAL") {
    ip.weight <- rep(1,nrow(normalized_data))
  } else if (outcome.type == "PROPORTIONAL" | outcome.type == "POLY-PROPORTIONAL") {
    ip.weight.matrix <- calculateIPCWMatrix(nuisance.model, normalized_data, estimand$code.censoring, strata, estimand, out_normalizeCovariate)
  }

  # ---- 目的関数/推定式（polyreg() 準拠） ---------------------------------------
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
      optim.method = NULL, prob.bound = prob.bound,
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
      optim.method = NULL, prob.bound = prob.bound,
      initial.CIFs = initial.CIFs
    )

    estimating_equation_pp <- function(p) call_and_capture(
      estimating_equation_pproportional,
      formula = nuisance.model, data = normalized_data, exposure = exposure,
      ip.weight.matrix = ip.weight.matrix, alpha_beta = p, estimand = estimand,
      optim.method = NULL, prob.bound = prob.bound,
      initial.CIFs = initial.CIFs
    )

    setInitialCIFs <- function(new.CIFs) initial.CIFs <<- new.CIFs
    getResults     <- function() out_ipcw

    list(
      estimating_equation_i  = estimating_equation_i,
      estimating_equation_s  = estimating_equation_s,
      estimating_equation_p  = estimating_equation_p,
      estimating_equation_pp = estimating_equation_pp,
      setInitialCIFs = setInitialCIFs,
      getResults     = getResults
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

  # ---- ループ（polyreg() と同じ構造） ------------------------------------------
  obj <- makeObjectiveFunction()
  estimating_fun <- choose_estimating_equation(outcome.type, obj)
  nleqslv_method <- choose_nleqslv_method(optim.method$nleqslv.method)

  iteration <- 0L
  max.absolute.difference <- Inf
  out_nleqslv <- NULL
  current_params <- alpha_beta_0
  current_obj_value <- numeric(0)
  trace_df  <- NULL

  computation.time0 <- proc.time()

  while ((iteration < optim.method$optim.parameter4) & (max.absolute.difference > optim.method$optim.parameter1)) {
    iteration <- iteration + 1L
    prev_params <- current_params

    out_nleqslv <- nleqslv(
      prev_params,
      estimating_fun,
      method  = optim.method$nleqslv_method,
      control = list(maxit = optim.method$optim.parameter5, allowSingular = FALSE)
    )
    new_params <- out_nleqslv$x

    # 目的関数値と内部状態の更新
    current_obj_value <- get_obj_value(new_params)
    obj$setInitialCIFs(obj$getResults()$potential.CIFs)

    # 収束判定（既存の assessConvergence() を利用する前提）
    ac <- assessConvergence(
      new_params, prev_params, current_obj_value,
      optim.method$optim.parameter1, optim.method$optim.parameter2, optim.method$optim.parameter3
    )
    current_params <- new_params
    max.absolute.difference <- ac$max.absolute.difference
    if (ac$converged) break
  }
  return(current_params)
}

solveEstimatingEquationP <- function(
    nuisance.model,
    exposure,
    strata,
    normalized_data,
    estimand,
    out_normalizeCovariate,
    optim.method,
    data.initial.values,
    prob.bound = 1e-5
) {
  outcome.type <- 'PROPORTIONAL'
  ip.weight.matrix <- calculateIPCWMatrix(nuisance.model, normalized_data, estimand$code.censoring, strata, estimand, out_normalizeCovariate)

  makeObjectiveFunction <- function() {
    out_ipcw <- list()
    initial.CIFs <- NULL
    estimating_equation_p <- function(p) {
      out_ipcw <- estimating_equation_proportional(
        formula = nuisance.model,
        data = normalized_data,
        exposure = exposure,
        ip.weight.matrix = ip.weight.matrix,
        alpha_beta = p,
        estimand = estimand,
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
      estimating_equation_p = estimating_equation_p,
      setInitialCIFs = setInitialCIFs,
      getResults = getResults
    )
  }

  iteration <- 0
  max_param_diff  <- Inf
  obj <- makeObjectiveFunction()
  sol_list <- list()
  diff_list <- list()

  current_params <- data.initial.values
  while ((iteration < optim.method$optim.parameter2) & (max_param_diff > optim.method$optim.parameter1)) {
    iteration <- iteration + 1
    if (optim.method$outer.optim.method == "nleqslv" | optim.method$outer.optim.method == "Broyden"){
      sol <- nleqslv(current_params, obj$estimating_equation_p, method="Broyden", control=list(maxit=optim.method$optim.parameter5, allowSingular=FALSE))
      new_params <- sol$x
    } else if (optim.method$optim.method$outer.optim.method == "Newton"){
      sol <- nleqslv(current_params, obj$estimating_equation_p, method="Newton", control=list(maxit=optim.method$optim.parameter5, allowSingular=FALSE))
      new_params <- sol$x
    } else if (optim.method$outer.optim.method == "optim" | optim.method$outer.optim.method == "SANN"){
      sol <- optim(par = current_params,
                   fn = function(params) {
                     sum(obj$estimating_equation_p(params)^2)
                   },
                   method = "SANN",  control = list(maxit=optim.method$optim.parameter5, reltol=optim.method$optim.parameter4)
      )
      new_params <- sol$par
    } else if (optim.method$outer.optim.method == "BFGS"){
      sol <- optim(par = current_params,
                   fn = function(params) {
                     sum(obj$estimating_equation_p(params)^2)
                   },
                   method = "BFGS",  control = list(maxit=optim.method$optim.parameter5, reltol=optim.method$optim.parameter4)
      )
      new_params <- sol$par
    }
    if (any(abs(new_params) > optim.method$optim.parameter3)) {
      stop("Estimates are either too large or too small, and convergence might not be achieved.")
    }
    param_diff <- abs(new_params - current_params)
    max_param_diff   <- max(param_diff)
    current_params <- new_params

    obj$setInitialCIFs(obj$getResults()$potential.CIFs)
    sol_list[[iteration]] <- sol
    diff_list[[iteration]] <- max_param_diff
  }
  return(current_params)
}

solveEstimatingEquationPP <- function(
    nuisance.model,
    exposure,
    strata,
    normalized_data,
    estimand,
    out_normalizeCovariate,
    optim.method,
    data.initial.values,
    prob.bound = 1e-5
) {
  outcome.type <- 'POLY-PROPORTIONAL'
  ip.weight.matrix <- calculateIPCWMatrix(nuisance.model, normalized_data, estimand$code.censoring, strata, estimand, out_normalizeCovariate)

  makeObjectiveFunction <- function() {
    out_ipcw <- list()
    initial.CIFs <- NULL
    estimating_equation_pp <- function(p) {
      out_ipcw <- estimating_equation_pproportional(
        formula = nuisance.model,
        data = normalized_data,
        exposure = exposure,
        ip.weight.matrix = ip.weight.matrix,
        alpha_beta = p,
        estimand = estimand,
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

  current_params <- data.initial.values
  while ((iteration < optim.method$optim.parameter2) & (max_param_diff > optim.method$optim.parameter1)) {
    iteration <- iteration + 1
    if (optim.method$outer.optim.method == "nleqslv" | optim.method$outer.optim.method == "Broyden"){
      sol <- nleqslv(current_params, obj$estimating_equation_pp, method="Broyden", control=list(maxit=optim.method$optim.parameter5, allowSingular=FALSE))
      new_params <- sol$x
    } else if (optim.method$optim.method$outer.optim.method == "Newton"){
      sol <- nleqslv(current_params, obj$estimating_equation_pp, method="Newton", control=list(maxit=optim.method$optim.parameter5, allowSingular=FALSE))
      new_params <- sol$x
    } else if (optim.method$outer.optim.method == "optim" | optim.method$outer.optim.method == "SANN"){
      sol <- optim(par = current_params,
                   fn = function(params) {
                     sum(obj$estimating_equation_pp(params)^2)
                   },
                   method = "SANN",  control = list(maxit=optim.method$optim.parameter5, reltol=optim.method$optim.parameter4)
      )
      new_params <- sol$par
    } else if (optim.method$outer.optim.method == "BFGS"){
      sol <- optim(par = current_params,
                   fn = function(params) {
                     sum(obj$estimating_equation_p(params)^2)
                   },
                   method = "BFGS",  control = list(maxit=optim.method$optim.parameter5, reltol=optim.method$optim.parameter4)
      )
      new_params <- sol$par
    }
    if (any(abs(new_params) > optim.method$optim.parameter3)) {
      stop("Estimates are either too large or too small, and convergence might not be achieved.")
    }
    param_diff <- abs(new_params - current_params)
    max_param_diff   <- max(param_diff)
    current_params <- new_params

    obj$setInitialCIFs(obj$getResults()$potential.CIFs)
    sol_list[[iteration]] <- sol
    diff_list[[iteration]] <- max_param_diff
  }
  return(current_params)
}

solveEstimatingEquationC <- function(
    nuisance.model,
    exposure,
    strata,
    normalized_data,
    estimand = estimand,
    optim.method = optim.method,
    data.initial.values,
    prob.bound = 1e-5
) {
  outcome.type <- 'COMPETINGRISK'
  ip.weight <- calculateIPCW(formula = nuisance.model, data = normalized_data, code.censoring=estimand$code.censoring, strata=strata, specific.time = estimand$time.point)

  makeObjectiveFunction <- function() {
    out_ipcw <- list()
    initial.CIFs <- NULL
    estimating_equation_i <- function(p) {
      out_ipcw <- estimating_equation_ipcw(
        formula = nuisance.model,
        data = normalized_data,
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
    setInitialCIFs <- function(new.CIFs) {
      initial.CIFs <<- new.CIFs
    }
    getResults <- function() {
      out_ipcw
    }
    list(
      estimating_equation_i = estimating_equation_i,
      setInitialCIFs = setInitialCIFs,
      getResults = getResults
    )
  }

  iteration <- 0
  max_param_diff  <- Inf
  obj <- makeObjectiveFunction()
  sol_list <- list()
  diff_list <- list()

  current_params <- data.initial.values
  while ((iteration < optim.method$optim.parameter2) & (max_param_diff > optim.method$optim.parameter1)) {
    iteration <- iteration + 1
    if (optim.method$outer.optim.method == "nleqslv" | optim.method$outer.optim.method == "Broyden"){
      sol <- nleqslv(current_params, obj$estimating_equation_i, method="Broyden", control=list(maxit=optim.method$optim.parameter5, allowSingular=FALSE))
      new_params <- sol$x
    } else if (optim.method$outer.optim.method == "Newton"){
      sol <- nleqslv(current_params, obj$estimating_equation_i, method="Newton", control=list(maxit=optim.method$optim.parameter5, allowSingular=FALSE))
      new_params <- sol$x
    } else if (optim.method$outer.optim.method == "optim" | optim.method$outer.optim.method == "SANN"){
      sol <- optim(par = current_params,
                   fn = function(params) {
                     sum(obj$estimating_equation_i(params)^2)
                   },
                   method = "SANN",  control = list(maxit=optim.method$optim.parameter5, reltol=optim.method$optim.parameter4)
      )
      new_params <- sol$par
    } else if (optim.method$outer.optim.method == "BFGS"){
      sol <- optim(par = current_params,
                   fn = function(params) {
                     sum(obj$estimating_equation_i(params)^2)
                   },
                   method = "BFGS",  control = list(maxit=optim.method$optim.parameter5, reltol=optim.method$optim.parameter4)
      )
      new_params <- sol$par
    }
    if (any(abs(new_params) > optim.method$optim.parameter3)) {
      stop("Estimates are either too large or too small, and convergence might not be achieved.")
    }
    param_diff <- abs(new_params - current_params)
    max_param_diff   <- max(param_diff)
    current_params <- new_params

    obj$setInitialCIFs(obj$getResults()$potential.CIFs)
    sol_list[[iteration]] <- sol
    diff_list[[iteration]] <- max_param_diff
  }
  return(current_params)
}

solveEstimatingEquationS <- function(
    nuisance.model,
    exposure,
    strata,
    normalized_data,
    estimand = estimand,
    optim.method = optim.method,
    data.initial.values,
    prob.bound = 1e-5
) {
  outcome.type <- "SURVIVAL"
  ip.weight <- calculateIPCW(formula = nuisance.model, data = normalized_data, code.censoring=estimand$code.censoring, strata=strata, specific.time = estimand$time.point)

  makeObjectiveFunction <- function() {
    out_ipcw <- list()
    initial.CIFs <- NULL
    estimating_equation_s <- function(p) {
      out_ipcw <- estimating_equation_survival(
        formula = nuisance.model,
        data = normalized_data,
        exposure = exposure,
        ip.weight = ip.weight,
        alpha_beta = p,
        estimand = estimand,
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
      estimating_equation_s = estimating_equation_s,
      setInitialCIFs = setInitialCIFs,
      getResults = getResults
    )
  }

  iteration <- 0
  max_param_diff  <- Inf
  obj <- makeObjectiveFunction()
  sol_list <- list()
  diff_list <- list()

  current_params <- data.initial.values
  while ((iteration < optim.method$optim.parameter2) & (max_param_diff > optim.method$optim.parameter1)) {
    iteration <- iteration + 1
    if (optim.method$outer.optim.method == "nleqslv" | optim.method$outer.optim.method == "Broyden"){
      sol <- nleqslv(current_params, obj$estimating_equation_s, method="Broyden", control=list(maxit=optim.method$optim.parameter5, allowSingular=FALSE))
      new_params <- sol$x
    } else if (optim.method$outer.optim.method == "Newton"){
      sol <- nleqslv(current_params, obj$estimating_equation_s, method="Newton", control=list(maxit=optim.method$optim.parameter5, allowSingular=FALSE))
      new_params <- sol$x
    } else if (optim.method$outer.optim.method == "optim" | optim.method$outer.optim.method == "SANN"){
      sol <- optim(par = current_params,
                   fn = function(params) {
                     sum(obj$estimating_equation_s(params)^2)
                   },
                   method = "SANN",  control = list(maxit=optim.method$optim.parameter5, reltol=optim.method$optim.parameter4)
      )
      new_params <- sol$par
    } else if (optim.method$outer.optim.method == "BFGS"){
      sol <- optim(par = current_params,
                   fn = function(params) {
                     sum(obj$estimating_equation_s(params)^2)
                   },
                   method = "BFGS",  control = list(maxit=optim.method$optim.parameter5, reltol=optim.method$optim.parameter4)
      )
      new_params <- sol$par
    }
    if (any(abs(new_params) > optim.method$optim.parameter3)) {
      stop("Estimates are either too large or too small, and convergence might not be achieved.")
    }
    param_diff     <- abs(new_params - current_params)
    max_param_diff <- max(param_diff)
    current_params <- new_params

    obj$setInitialCIFs(obj$getResults()$potential.CIFs)
    sol_list[[iteration]] <- sol
    diff_list[[iteration]] <- max_param_diff
  }
  return(current_params)
}

solveEstimatingEquationB <- function(
    nuisance.model,
    exposure,
    strata,
    normalized_data,
    estimand = estimand,
    optim.method = optim.method,
    data.initial.values,
    prob.bound = 1e-5
) {
  outcome.type <- 'BINOMIAL'
  ip.weight <- rep(1,nrow(normalized_data))

  makeObjectiveFunction <- function() {
    out_ipcw <- list()
    initial.CIFs <- NULL
    estimating_equation_s <- function(p) {
      out_ipcw <- estimating_equation_ipcw(
        formula = nuisance.model,
        data = normalized_data,
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
    setInitialCIFs <- function(new.CIFs) {
      initial.CIFs <<- new.CIFs
    }
    getResults <- function() {
      out_ipcw
    }
    list(
      estimating_equation_s = estimating_equation_s,
      setInitialCIFs = setInitialCIFs,
      getResults = getResults
    )
  }

  iteration <- 0
  max_param_diff  <- Inf
  obj <- makeObjectiveFunction()
  sol_list <- list()
  diff_list <- list()

  current_params <- data.initial.values
  while ((iteration < optim.method$optim.parameter2) & (max_param_diff > optim.method$optim.parameter1)) {
    iteration <- iteration + 1
    if (optim.method$outer.optim.method == "nleqslv" | optim.method$outer.optim.method == "Broyden"){
      sol <- nleqslv(current_params, obj$estimating_equation_s, method="Broyden", control=list(maxit=optim.method$optim.parameter5, allowSingular=FALSE))
      new_params <- sol$x
    } else if (optim.method$outer.optim.method == "Newton"){
      sol <- nleqslv(current_params, obj$estimating_equation_s, method="Newton", control=list(maxit=optim.method$optim.parameter5, allowSingular=FALSE))
      new_params <- sol$x
    } else if (optim.method$outer.optim.method == "optim" | optim.method$outer.optim.method == "SANN"){
      sol <- optim(par = current_params,
                   fn = function(params) {
                     sum(obj$estimating_equation_s(params)^2)
                   },
                   method = "SANN",  control = list(maxit=optim.method$optim.parameter5, reltol=optim.method$optim.parameter4)
      )
      new_params <- sol$par
    } else if (outer.optim.method == "BFGS"){
      sol <- optim(par = current_params,
                   fn = function(params) {
                     sum(obj$estimating_equation_s(params)^2)
                   },
                   method = "BFGS",  control = list(maxit=optim.method$optim.parameter5, reltol=optim.method$optim.parameter4)
      )
      new_params <- sol$par
    }
    if (any(abs(new_params) > optim.method$optim.parameter3)) {
      stop("Estimates are either too large or too small, and convergence might not be achieved.")
    }
    param_diff <- abs(new_params - current_params)
    max_param_diff   <- max(param_diff)
    current_params <- new_params

    obj$setInitialCIFs(obj$getResults()$potential.CIFs)
    sol_list[[iteration]] <- sol
    diff_list[[iteration]] <- max_param_diff
  }
  return(current_params)
}

