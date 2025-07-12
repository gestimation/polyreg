Surv <- function (time, event)
{
  if (missing(time))
    stop("Must have a time argument")
  if (!is.numeric(time))
    stop("Time variable is not numeric")
  if (!is.na(any(time)))
    warning("Invalid time variable. NA values included")
  #  if (any(time<0))
  #    warning("Invalid time variable. Non-negative values included")
  if (length(event) != length(time))
    stop("Time and event variables are different lengths")
  if (missing(event))
    stop("Must have an event argument")
  if (is.numeric(event)) {
    if (!is.na(any(event)))
      warning("Invalid event variable. NA values included")
    status <- event
  } else if (is.logical(event)) {
    if (!is.na(any(event)))
      warning("Invalid event variable. NA values included")
    status <- as.numeric(event)
    warning("Event variable is logical, converted to numearic")
  } else if (is.factor(event)) {
    status <- as.numeric(as.factor(event)) - 1
    warning("Event variable is a factor, converted to numearic")
  } else stop("Invalid status value, must be logical or numeric")
  if (nlevels(as.factor(event)) > 3)
    warning("Event variable should not have more than three levels")
  ss <- cbind(time = time, status = status)
  type <- "right"

  inputAttributes <- list()
  if (!is.null(attributes(time)))
    inputAttributes$time <- attributes(time)
  cname <- dimnames(ss)[[2]]
  if (length(cname) == 0) {
    cname <- c("time", "status")
  }
  dimnames(ss) <- list(NULL, cname)
  attr(ss, "type") <- type
  if (length(inputAttributes) > 0)
    attr(ss, "inputAttributes") <- inputAttributes
  class(ss) <- "Surv"
  ss
}

Event <- function (time, event)
{
  if (missing(time))
    stop("A time argument is required")
  if (!is.numeric(time))
    stop("Time variable is not numeric")
  if (!is.na(any(time)))
    warning("Invalid time variable. NA values included")
  #  if (any(time<0))
  #    warning("Invalid time variable. Non-negative values included")
  if (length(event) != length(time))
    stop("Time and event variables are different lengths")
  if (missing(event))
    stop("An event argument is required")
  if (is.numeric(event)) {
    if (!is.na(any(event)))
      warning("Invalid event variable. NA values included")
    status <- event
  } else if (is.logical(event)) {
    if (!is.na(any(event)))
      warning("Invalid event variable. NA values included")
    status <- as.numeric(event)
    warning("Event variable is logical, converted to numearic")
  } else if (is.factor(event)) {
    status <- as.numeric(as.factor(event)) - 1
    warning("Event variable is a factor, converted to numearic")
  } else stop("Invalid status value, must be logical or numeric")
  if (nlevels(as.factor(event)) > 3)
    warning("Event variable should not have more than three levels")
  ss <- cbind(time = time, status = status)
  type <- "right"

  inputAttributes <- list()
  if (!is.null(attributes(time)))
    inputAttributes$time <- attributes(time)
  cname <- dimnames(ss)[[2]]
  if (length(cname) == 0) {
    cname <- c("time", "status")
  }
  dimnames(ss) <- list(NULL, cname)
  attr(ss, "type") <- type
  if (length(inputAttributes) > 0)
    attr(ss, "inputAttributes") <- inputAttributes
  class(ss) <- "Event"
  ss
}

createAnalysisDataset <- function(formula, data, other.variables.analyzed=NULL, subset.condition=NULL, na.action=na.pass) {
  if (!is.null(subset.condition)) {
    analysis_dataset <- subset(data, eval(parse(text = subset.condition)))
  } else {
    analysis_dataset <- data
  }
  all_vars <- c(all.vars(formula), other.variables.analyzed)
  analysis_dataset <- analysis_dataset[, all_vars, drop = FALSE]
  return(na.action(analysis_dataset))
}

getInitialValues <- function(
    formula, data, exposure, data.initial.values, estimand, specific.time, outcome.type, prob.bound
) {
  cl <- match.call()
  mf <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "cluster", "offset")
  Terms <- terms(formula, special, data = data)
  mf$formula <- Terms
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Y <- model.extract(mf, "response")

  if (!inherits(Y, c("Event", "Surv"))) {
    if (outcome.type %in% c('COMPETINGRISK', 'SURVIVAL', 'PROPORTIONAL')) {
      stop("Expected a 'Surv' or 'Event'-object when outcome.type is COMPETINGRISK, SURVIVAL or PROPORTIONAL. ")
    } else {
      t <- rep(0, length(Y))
      epsilon <- Y
      if (!all(epsilon %in% c(estimand$code.event1, estimand$code.censoring))) {
        stop("Invalid event codes. Must be 0 or 1, if event codes are not specified. ")
      }
    }
  } else {
    t <- as.numeric(Y[, 1])
    if (any(t<0))
      stop("Invalid time variable. Expected non-negative values. ")
    if (any(is.na(t)))
      stop("Time variable contains NA values")
    if (any(is.na(Y[, 2]))) {
      stop("Event variable contains NA values")
    } else {
      epsilon <- Y[, 2]
    }
    if (!all(epsilon %in% c(estimand$code.event1, estimand$code.censoring)) & (outcome.type == 'SURVIVAL')) {
      stop("Invalid event codes. Must be 0 or 1, with 0 representing censoring, if event codes are not specified. ")
    } else if (!all(epsilon %in% c(estimand$code.event1, estimand$code.event2, estimand$code.censoring))) {
      stop("Invalid event codes. Must be 0, 1 or 2, with 0 representing censoring, if event codes are not specified. ")
    }
  }
  if (!is.null(offsetpos <- attributes(Terms)$specials$offset)) {
    ts <- survival::untangle.specials(Terms, "offset")
    if (length(ts$vars) > 0) {
      Terms <- Terms[-ts$terms]
      offset <- mf[[ts$vars]]
    } else {
      offset <- rep(0, nrow(mf))
    }
  } else {
    offset <- rep(0, nrow(mf))
  }

  a_ <- as.factor(data[[exposure]])
  if (estimand$code.exposure.ref==0) {
    x_a <- as.matrix(model.matrix(~ a_)[, -1])
  } else {
    x_a <- as.matrix(rep(1,length(t)) - model.matrix(~ a_)[, -1])
  }
  a <- as.vector(x_a)
  x_l <- model.matrix(Terms, mf)
  n_para_1 <- ncol(x_l)
  if (any(is.na(data[[exposure]])))
    stop("Exposure variable contains NA values")
  if (any(is.na(x_l)))
    stop("Covariates contain NA values")

  if (!is.null(data.initial.values)) {
    x_l <- model.matrix(Terms, mf)
    if (!(1+ncol(x_l))*2 == length(data.initial.values))
      stop("Invalid initial value dataset. Must contain the same number of initial values as parameters")
    return(data.initial.values)
  }

  binarizeIfContinuous <- function(x) {
    if (outcome.type == 'SURVIVAL' | outcome.type == 'BINOMIAL') {
      if (is.numeric(x) & length(unique(x)) > 2) {
        l <- as.numeric((x >= median(x)) == TRUE)
        out <- calculateInitialValuesSurvival(t, epsilon, a, l, estimand, specific.time, prob.bound)
        return(out)
      } else if (length(unique(x)) == 2) {
        l <- x
        out <- calculateInitialValuesSurvival(t, epsilon, a, l, estimand, specific.time, prob.bound)
        return(out)
      }
    } else {
      if (is.numeric(x) & length(unique(x)) > 2) {
        l <- as.numeric((x >= median(x)) == TRUE)
        out <- calculateInitialValuesCompetingRisk(t, epsilon, a, l, estimand, specific.time, prob.bound)
        return(out)
      } else if (length(unique(x)) == 2) {
        l <- x
        out <- calculateInitialValuesCompetingRisk(t, epsilon, a, l, estimand, specific.time, prob.bound)
        return(out)
      }
    }
  }

  if (all(epsilon %in% c(estimand$code.event1, estimand$code.censoring))) {
    if (n_para_1>1) {
      out_bic_1 <- t(binarizeIfContinuous(x_l[,2])[1,1:2])
      if (n_para_1>2) {
        for (i_para in 3:n_para_1) {
          out_bic_i <- binarizeIfContinuous(x_l[,i_para])
          out_bic_1 <- cbind(out_bic_1, t(out_bic_i[1,2]))
        }
        out_bic_1 <- cbind(out_bic_1, t(out_bic_i[1,3]))
      } else {
        out_bic_1 <- cbind(out_bic_1, t(binarizeIfContinuous(x_l[,2])[1,3]))
      }
      init_vals <- out_bic_1
    } else {
      l <- NULL
      init_vals <- calculateInitialValuesSurvival(t, epsilon, a, l, estimand, specific.time, prob.bound)
    }
  } else {
    if (n_para_1>1) {
      out_bic_1 <- t(binarizeIfContinuous(x_l[,2])[1,1:2])
      out_bic_2 <- t(binarizeIfContinuous(x_l[,2])[1,4:5])
      if (n_para_1>2) {
        for (i_para in 3:n_para_1) {
          out_bic_i <- binarizeIfContinuous(x_l[,i_para])
          out_bic_1 <- cbind(out_bic_1, t(out_bic_i[1,2]))
          out_bic_2 <- cbind(out_bic_2, t(out_bic_i[1,5]))
        }
        out_bic_1 <- cbind(out_bic_1, t(out_bic_i[1,3]))
        out_bic_2 <- cbind(out_bic_2, t(out_bic_i[1,6]))
      } else {
        out_bic_1 <- cbind(out_bic_1, t(binarizeIfContinuous(x_l[,2])[1,3]))
        out_bic_2 <- cbind(out_bic_2, t(binarizeIfContinuous(x_l[,2])[1,6]))
      }
      init_vals <- cbind(out_bic_1, out_bic_2)
    } else {
      l <- NULL
      init_vals <- calculateInitialValuesCompetingRisk(t, epsilon, a, l, estimand, specific.time, prob.bound)
    }
  }
  return(init_vals)
}


######## ここを修正
calculateInitialValuesCompetingRisk <- function(t, epsilon, a, l = NULL, estimand, specific.time, prob.bound) {

  ### この部分は常に計算されないとp_10がなくなる
  #if (is.null(l)) {
  # Equivalent to calc_initial_2_competing_risk
  epsilon0 <- epsilon[a == 0]
  epsilon1 <- epsilon[a == 1]
  t0 <- t[a == 0]
  t1 <- t[a == 1]

  p_10 <- (sum(epsilon0 == estimand$code.event1 & t0 <= specific.time) / length(epsilon0)) + prob.bound
  p_20 <- (sum(epsilon0 == estimand$code.event2 & t0 <= specific.time) / length(epsilon0)) + prob.bound
  p_00 <- 1 - p_10 - p_20
  p_11 <- (sum(epsilon1 == estimand$code.event1 & t1 <= specific.time) / length(epsilon1)) + prob.bound
  p_21 <- (sum(epsilon1 == estimand$code.event2 & t1 <= specific.time) / length(epsilon1)) + prob.bound
  p_01 <- 1 - p_11 - p_21

  alpha_1 <- log((p_10 * p_11) / (p_00 * p_01))
  alpha_2 <- log((p_20 * p_21) / (p_00 * p_01))
  #} else {
  ### 逆にこっちのほうをif分岐させる
  if(!is.null(l)){
    # Equivalent to calc_initial_1_competing_risk
    epsilon00 <- epsilon[a == 0 & l == 0]
    epsilon10 <- epsilon[a == 1 & l == 0]
    epsilon01 <- epsilon[a == 0 & l == 1]
    epsilon11 <- epsilon[a == 1 & l == 1]
    t00 <- t[a == 0 & l == 0]
    t10 <- t[a == 1 & l == 0]
    t01 <- t[a == 0 & l == 1]
    t11 <- t[a == 1 & l == 1]

    p_100 <- (sum(epsilon00 == estimand$code.event1 & t00 <= specific.time) / length(epsilon00)) + prob.bound
    p_200 <- (sum(epsilon00 == estimand$code.event2 & t00 <= specific.time) / length(epsilon00)) + prob.bound
    p_000 <- 1 - p_100 - p_200
    p_110 <- (sum(epsilon10 == estimand$code.event1 & t10 <= specific.time) / length(epsilon10)) + prob.bound
    p_210 <- (sum(epsilon10 == estimand$code.event2 & t10 <= specific.time) / length(epsilon10)) + prob.bound
    p_010 <- 1 - p_110 - p_210
    p_101 <- (sum(epsilon01 == estimand$code.event1 & t01 <= specific.time) / length(epsilon01)) + prob.bound
    p_201 <- (sum(epsilon01 == estimand$code.event2 & t01 <= specific.time) / length(epsilon01)) + prob.bound
    p_001 <- 1 - p_101 - p_201
    p_111 <- (sum(epsilon11 == estimand$code.event1 & t11 <= specific.time) / length(epsilon11)) + prob.bound
    p_211 <- (sum(epsilon11 == estimand$code.event2 & t11 <= specific.time) / length(epsilon11)) + prob.bound
    p_011 <- 1 - p_111 - p_211

    alpha_10 <- log((p_100 * p_110) / (p_000 * p_010))
    alpha_20 <- log((p_200 * p_210) / (p_000 * p_010))
    alpha_11 <- log((p_101 * p_111) / (p_000 * p_011)) - alpha_10
    alpha_21 <- log((p_201 * p_211) / (p_000 * p_011)) - alpha_20
  }

  # Compute beta values
  if (estimand$effect.measure1 == 'RR') {
    beta_1 <- log(p_11 / p_10)
  } else if (estimand$effect.measure1 == 'OR') {
    beta_1 <- log((p_11 / (1 - p_11)) / (p_10 / (1 - p_10)))
  } else if (estimand$effect.measure1 == 'SHR') {
    beta_1 <- log(log(1 - p_11) / log(1 - p_10))
  } else {
    stop("Invalid effect measure code. Must be RR, OR or SHR.")
  }

  if (estimand$effect.measure2 == 'RR') {
    beta_2 <- log(p_21 / p_20)
  } else if (estimand$effect.measure2 == 'OR') {
    beta_2 <- log((p_21 / (1 - p_21)) / (p_20 / (1 - p_20)))
  } else if (estimand$effect.measure2 == 'SHR') {
    beta_2 <- log(log(1 - p_21) / log(1 - p_20))
  } else {
    stop("Invalid effect measure code. Must be RR, OR or SHR.")
  }

  alpha_beta <- if (is.null(l)) {
    cbind(alpha_1, beta_1, alpha_2, beta_2)
  } else {
    cbind(alpha_10, alpha_11, beta_1, alpha_20, alpha_21, beta_2)
  }

  return(alpha_beta)
}

calculateInitialValuesSurvival <- function(t, epsilon, a, l = NULL, estimand, specific.time, prob.bound) {
  if (is.null(l)) {
    # Use calc_initial_2_survival logic
    epsilon0 <- epsilon[a == 0]
    epsilon1 <- epsilon[a == 1]
    t0 <- t[a == 0]
    t1 <- t[a == 1]
    p_10 <- (sum(epsilon0 == estimand$code.event1 & t0 <= specific.time) / length(epsilon0)) + prob.bound
    p_00 <- 1 - p_10
    p_11 <- (sum(epsilon1 == estimand$code.event1 & t1 <= specific.time) / length(epsilon1)) + prob.bound
    p_01 <- 1 - p_11

    if (estimand$effect.measure1 == 'RR') {
      beta_1 <- log(p_11 / p_10)
    } else if (estimand$effect.measure1 == 'OR') {
      beta_1 <- log((p_11 / (1 - p_11)) / (p_10 / (1 - p_10)))
    } else if (estimand$effect.measure1 == 'SHR') {
      beta_1 <- log(log(1 - p_11) / log(1 - p_10))
    } else {
      stop("Invalid effect measure code. Must be RR, OR or SHR.")
    }

    alpha_1 <- log((p_10 * p_11) / (p_00 * p_01))
    return(cbind(alpha_1, beta_1))
  } else {
    # Use calc_initial_1_survival logic
    epsilon00 <- epsilon[a == 0 & l == 0]
    epsilon10 <- epsilon[a == 1 & l == 0]
    epsilon01 <- epsilon[a == 0 & l == 1]
    epsilon11 <- epsilon[a == 1 & l == 1]
    t00 <- t[a == 0 & l == 0]
    t10 <- t[a == 1 & l == 0]
    t01 <- t[a == 0 & l == 1]
    t11 <- t[a == 1 & l == 1]

    p_100 <- (sum(epsilon00 == estimand$code.event1 & t00 <= specific.time) / length(epsilon00)) + prob.bound
    p_000 <- 1 - p_100
    p_110 <- (sum(epsilon10 == estimand$code.event1 & t10 <= specific.time) / length(epsilon10)) + prob.bound
    p_010 <- 1 - p_110
    p_101 <- (sum(epsilon01 == estimand$code.event1 & t01 <= specific.time) / length(epsilon01)) + prob.bound
    p_001 <- 1 - p_101
    p_111 <- (sum(epsilon11 == estimand$code.event1 & t11 <= specific.time) / length(epsilon11)) + prob.bound
    p_011 <- 1 - p_111

    if (any(c(p_100, p_000, p_110, p_010, p_101, p_001, p_111, p_011) == 0, na.rm = TRUE)) {
      stop("Complete separation detected in initial value search")
    }

    alpha_10 <- log((p_100 * p_110) / (p_000 * p_010))
    alpha_11 <- log((p_101 * p_111) / (p_000 * p_011)) - alpha_10

    epsilon0 <- epsilon[a == 0]
    epsilon1 <- epsilon[a == 1]
    t0 <- t[a == 0]
    t1 <- t[a == 1]
    p_10 <- (sum(epsilon0 == estimand$code.event1 & t0 <= specific.time) / length(epsilon0)) + prob.bound
    p_00 <- 1 - p_10
    p_11 <- (sum(epsilon1 == estimand$code.event1 & t1 <= specific.time) / length(epsilon1)) + prob.bound
    p_01 <- 1 - p_11

    if (estimand$effect.measure1 == 'RR') {
      beta_1 <- log(p_11 / p_10)
    } else if (estimand$effect.measure1 == 'OR') {
      beta_1 <- log((p_11 / (1 - p_11)) / (p_10 / (1 - p_10)))
    } else if (estimand$effect.measure1 == 'SHR') {
      beta_1 <- log(log(1 - p_11) / log(1 - p_10))
    } else {
      stop("Invalid effect measure code. Must be RR, OR or SHR.")
    }

    return(cbind(alpha_10, alpha_11, beta_1))
  }
}


normalizeCovariate <- function(formula, data, should.normalize.covariate, outcome.type) {
  mf <- model.frame(formula, data)
  Y <- model.extract(mf, "response")
  response_term <- formula[[2]]
  if (inherits(mf[[1]], "Surv") | inherits(mf[[1]], "Event")) {
    response_vars <- all.vars(response_term)
    covariate_cols <- setdiff(all.vars(formula), response_vars)  # Remove time and event
  } else {
    covariate_cols <- all.vars(formula)[-1]  # Exclude the response variable
  }
  normalized_data <- data
  range_vector <- 1
  if (should.normalize.covariate == TRUE & length(covariate_cols)>0) {
    for (col in covariate_cols) {
      x <- normalized_data[[col]]
      range <- max(x)-min(x)
      normalized_data[[col]] <- x/range
      range_vector <- cbind(range_vector,range)
    }
    if (outcome.type == 'PROPORTIONAL') {
      range_vector <- cbind(range_vector)
    } else if (outcome.type == 'SURVIVAL' | outcome.type == 'BINOMIAL') {
      range_vector <- cbind(range_vector,1)
    } else {
      range_vector <- cbind(range_vector,1,range_vector,1)
    }
  } else {
    if (outcome.type == 'PROPORTIONAL') {
      range_vector <- rep(1, (length(covariate_cols)+1))
    } else if (outcome.type == 'SURVIVAL' | outcome.type == 'BINOMIAL') {
      range_vector <- rep(1, (length(covariate_cols)+2))
    } else {
      range_vector <- rep(1, (2*length(covariate_cols)+4))
    }
  }
  n_covariate <- length(covariate_cols)
  out <- list(normalized_data=normalized_data, range=range_vector, n_covariate=n_covariate, n=nrow(data))
  return(out)
}

sortByCovariate <- function(formula, data, should.sort.data, n_covariate) {
  if (should.sort.data == TRUE & n_covariate>0) {
    terms_obj <- terms(formula)
    covariate_names <- attr(terms_obj, "term.labels")
    missing_vars <- setdiff(covariate_names, names(data))
    if (length(missing_vars) > 0) {
      stop("The following covariates are missing: ", paste(missing_vars, collapse = ", "))
    }
    sorted_data <- data[do.call(order, data[covariate_names]), , drop = FALSE]
    return(sorted_data)
  } else {
    return(data)
  }
}

checkSpell <- function(outcome.type, effect.measure1, effect.measure2) {
  if (requireNamespace("nleqslv", quietly = TRUE)
      & requireNamespace("boot", quietly = TRUE)) {
    suppressWarnings(library(nleqslv))
    suppressWarnings(library(boot))
  } else {
    stop("Required packages 'nleqslv' and/or 'boot' are not installed.")
  }
  if (outcome.type %in% c("COMPETINGRISK", "C", "CR", "COMPETING RISK", "COMPETINGRISKS", "COMPETING RISKS", "Competingrisk", "Competing risk", "Competingrisks", "Competing risks", "competingrisk", "competing risk", "competingrisks", "competing risks")) {
    outcome.type.corrected <<- "COMPETINGRISK"
  } else if (outcome.type %in% c("SURVIVAL", "S", "Survival", "Survival")) {
    outcome.type.corrected <<- "SURVIVAL"
  } else if (outcome.type %in% c("PROPORTIONAL", "P", "Proportional", "proportional")) {
    outcome.type.corrected <<- "PROPORTIONAL"
  } else if (outcome.type %in% c("BINOMIAL", "B", "Binomial", "binomial")) {
    outcome.type.corrected <<- "BINOMIAL"
  } else {
    stop("Invalid input for outcome.type, Choose 'COMPETINGRISK', 'SURVIVAL', 'BINOMIAL', or 'PROPORTIONAL'.")
  }
  if (effect.measure1 %in% c("RR", "rr", "RISK RATIO", "Risk ratio", "risk ratio")) {
    effect.measure1.corrected <<- "RR"
  } else if (effect.measure1 %in% c("OR", "or", "ODDS RATIO", "Odds ratio", "odds ratio")) {
    effect.measure1.corrected <<- "OR"
  } else if (effect.measure1 %in% c("SHR", "shr", "HR", "hr", "SUBDISTRIBUTION HAZARD RATIO",
                                    "Subdistibution hazard ratio", "subdistibution hazard ratio")) {
    effect.measure1.corrected <<- "SHR"
  } else {
    stop("Invalid input for effect.measure1, Choose 'RR', 'OR', or 'SHR'.")
  }
  if (effect.measure2 %in% c("RR", "rr", "RISK RATIO", "Risk ratio", "risk ratio")) {
    effect.measure2.corrected <<- "RR"
  } else if (effect.measure2 %in% c("OR", "or", "ODDS RATIO", "Odds ratio", "odds ratio")) {
    effect.measure2.corrected <<- "OR"
  } else if (effect.measure2 %in% c("SHR", "shr", "HR", "hr", "SUBDISTRIBUTION HAZARD RATIO",
                                 "Subdistibution hazard ratio", "subdistibution hazard ratio")) {
    effect.measure2.corrected <<- "SHR"
  } else {
    stop("Invalid input for effect.measure2, Choose 'RR', 'OR', or 'SHR'.")
  }
}

checkInput <- function(outcome.type, time.point, conf.level, outer.optim.method, inner.optim.method) {
  if ((outcome.type == "COMPETINGRISK" | outcome.type == "SURVIVAL") & length(time.point)>1) {
    time.point.corrected <<- max(time.point)
  } else if (is.null(time.point)) {
    stop("Invalid input for time.point when outcome.type is COMPETINGRISK or SURVIVAL.")
  } else if (min(time.point)<0) {
    stop("time.point should be positive when outcome.type is COMPETINGRISK or SURVIVAL.")
  } else {
    time.point.corrected <<- time.point
  }
  if (outcome.type == "PROPORTIONAL" & length(time.point)==1) {
    outcome.type.corrected <<- "COMPETINGRISK"
  }
  if (outcome.type == "BINOMIAL") {
    time.point.corrected <<- Inf
  }
  if (conf.level <= 0 | conf.level >= 1)
    stop("Confidence level must be between 0 and 1")
  if (outcome.type == "COMPETINGRISK" & !outer.optim.method %in% c("nleqslv","Newton","Broyden","optim","BFGS","SANN","multiroot","partial")) {
    stop("Invalid input for 'optimization'. Choose 'nleqslv', 'Newton', 'Broyden', 'optim', 'BFGS', 'SANN', 'multiroot' or 'partial'.")
  }
  if (outcome.type == "SURVIVAL" & outer.optim.method == "partial") {
    stop("Invalid input for 'optimization'. Choose 'nleqslv', 'Newton', 'Broyden', 'optim', 'BFGS', 'SANN' or 'multiroot'.")
  }
  if (!inner.optim.method %in% c("optim","BFGS","SANN","multiroot")) {
    stop("Invalid input for 'optimization'. Choose 'optim', 'BFGS', 'SANN' or 'multiroot'.")
  }
}
