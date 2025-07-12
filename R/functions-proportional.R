solveEstimatingEquationP <- function(
    nuisance.model,
    exposure,
    strata,
    sorted_data,
    estimand = estimand,
    optim.method = optim.method,
    data.initial.values,
    prob.bound = 1e-5
) {
  outcome.type <- 'PROPORTIONAL'
  i_time <- 0
  ip.weight.matrix <- matrix(NA, nrow=nrow(sorted_data), ncol=length(estimand$time.point))
  for (specific.time in estimand$time.point) {
    i_time <- i_time + 1
    ip.weight.matrix[,i_time] <- calculateIPCW(formula = nuisance.model, data = sorted_data, code.censoring=estimand$code.censoring, strata=strata, specific.time = specific.time)
  }

  makeObjectiveFunction <- function() {
    out_ipcw <- list()
    initial.CIFs <- NULL
    estimating_equation_p <- function(p) {
      out_ipcw <- estimating_equation_proportional(
        formula = nuisance.model,
        data = sorted_data,
        exposure = exposure,
        ip.weight.matrix = ip.weight.matrix,
        alpha_beta = p,
        estimand = estimand,
        time.point = estimand$time.point,
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
    sorted_data,
    estimand = estimand,
    optim.method = optim.method,
    data.initial.values,
    prob.bound = 1e-5
) {
  outcome.type <- 'POLY-PROPORTIONAL'
  i_time <- 0
  ip.weight.matrix <- matrix(NA, nrow=nrow(sorted_data), ncol=length(estimand$time.point))
  for (specific.time in estimand$time.point) {
    i_time <- i_time + 1
    ip.weight.matrix[,i_time] <- calculateIPCW(formula = nuisance.model, data = sorted_data, code.censoring=estimand$code.censoring, strata=strata, specific.time = specific.time)
  }

  makeObjectiveFunction <- function() {
    out_ipcw <- list()
    initial.CIFs <- NULL
    estimating_equation_pp <- function(p) {
      out_ipcw <- estimating_equation_pproportional(
        formula = nuisance.model,
        data = sorted_data,
        exposure = exposure,
        ip.weight.matrix = ip.weight.matrix,
        alpha_beta = p,
        estimand = estimand,
        time.point = estimand$time.point,
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
    sorted_data,
    estimand = estimand,
    optim.method = optim.method,
    data.initial.values,
    prob.bound = 1e-5
) {
  outcome.type <- 'COMPETINGRISK'
  ip.weight <- calculateIPCW(formula = nuisance.model, data = sorted_data, code.censoring=estimand$code.censoring, strata=strata, specific.time = estimand$time.point)

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
    sorted_data,
    estimand = estimand,
    optim.method = optim.method,
    data.initial.values,
    prob.bound = 1e-5
) {
  outcome.type <- 'SURVIVAL'
  ip.weight <- calculateIPCW(formula = nuisance.model, data = sorted_data, code.censoring=estimand$code.censoring, strata=strata, specific.time = estimand$time.point)

  makeObjectiveFunction <- function() {
    out_ipcw <- list()
    initial.CIFs <- NULL
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
    sorted_data,
    estimand = estimand,
    optim.method = optim.method,
    data.initial.values,
    prob.bound = 1e-5
) {
  outcome.type <- 'BINOMIAL'
  ip.weight <- rep(1,nrow(sorted_data))

  makeObjectiveFunction <- function() {
    out_ipcw <- list()
    initial.CIFs <- NULL
    estimating_equation_s <- function(p) {
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

