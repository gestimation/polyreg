Surv <- function(time, event) {
  if (missing(time))
    stop("Must have a time argument")
  if (!is.numeric(time))
    stop("Time variable is not numeric")
  #if (any(is.na(time)))
  #  warning("Invalid time variable. NA values included")
  #if (any(time<0))
  #  warning("Invalid time variable. Non-negative values included")
  if (length(event) != length(time))
    stop("Time and event variables are different lengths")
  if (missing(event))
    stop("Must have an event argument")
  if (is.numeric(event)) {
    if (any(is.na(event)))
      warning("Invalid event variable. NA values included")
    status <- event
  } else if (is.logical(event)) {
    if (any(is.na(event)))
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

Event <- function(time, event) {
  if (missing(time))
    stop("A time argument is required")
  if (!is.numeric(time))
    stop("Time variable is not numeric")
  #if (any(is.na(time)))
  #  warning("Invalid time variable. NA values included")
  #if (any(time<0))
  #  warning("Invalid time variable. Non-negative values included")
  if (length(event) != length(time))
    stop("Time and event variables are different lengths")
  if (missing(event))
    stop("An event argument is required")
  if (is.numeric(event)) {
    if (any(is.na(event)))
      warning("Invalid event variable. NA values included")
    status <- event
  } else if (is.logical(event)) {
    if (any(is.na(event)))
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

defineExposureDesign <- function(data, exposure, code.exposure.ref = NULL, prefix = "a") {
  stopifnot(is.data.frame(data))
  if (!exposure %in% names(data)) {
    stop("exposure = '", exposure, "' is not found in data.")
  }

  a_ <- data[[exposure]]
  a_ <- factor(a_)
  a_ <- base::droplevels(a_)
  lev <- levels(a_)
  K <- length(lev)

  ref_lab <- NULL
  if (!is.null(code.exposure.ref)) {
    ref_lab <- if (is.numeric(code.exposure.ref)) as.character(code.exposure.ref) else code.exposure.ref
    if (length(lev) > 0 && ref_lab %in% lev) {
      a_ <- stats::relevel(a_, ref = ref_lab)
    } else {
      warning("code.exposure.ref = ", ref_lab," is not found among factor levels. The first level is used as reference.")
      ref_lab <- NULL
    }
  }
  if (is.null(ref_lab)) ref_lab <- lev[1L]
  if (K < 1 || K == 1) stop("Exposure has only one level (", lev, ") or no valid levels. Effect estimation is not possible.")
  X <- stats::model.matrix(~ a_)[, -1, drop = FALSE]
  cn <- colnames(X)
  cn <- sub("^a_", paste0(prefix, "_"), cn, perl = TRUE)
  colnames(X) <- cn
  return(list(
    x_a = as.matrix(X),
    exposure.levels = K,
    exposure.labels = lev,
    ref = ref_lab
  ))
}

checkInput <- function(data, formula, exposure, code.event1, code.event2, code.censoring, code.exposure.ref, outcome.type, conf.level, report.sandwich.conf, report.boot.conf, nleqslv.method, should.normalize.covariate) {
  cl <- match.call()
  if (missing(formula)) stop("A formula argument is required")
  mf <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "offset", "cluster")
  out_terms <- terms(formula, special, data = data)
  if (!is.null(attr(out_terms, "specials")$strata))
    stop("strata() cannot appear in formula")
  if (!is.null(attr(out_terms, "specials")$offset))
    stop("offset() cannot appear in formula")
  if (!is.null(attr(out_terms, "specials")$cluster))
    stop("cluster() cannot appear in formula")
  mf$formula <- out_terms
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Y <- model.extract(mf, "response")
  if (outcome.type %in% c("COMPETING-RISK","SURVIVAL","POLY-PROPORTIONAL","POLY-PROPORTIONAL")) {
    if (!inherits(Y, c("Event", "Surv"))) {
      stop("Surv- or Event-object is expected")
    } else {
      t <- as.numeric(Y[, 1])
      if (any(t<0)) stop("Invalid time variable. Expected non-negative values. ")
      if (any(is.na(t))) stop("Time variable contains NA.")

      epsilon <- as.numeric(Y[, 2])
      if (any(is.na(epsilon))) stop("Event variable contains NA.")
      if (outcome.type == "SURVIVAL") {
        if (!all(epsilon %in% c(code.event1, code.censoring)))
          stop("SURVIVAL requires event codes {censoring,event1}.")
      } else {
        if (!all(epsilon %in% c(code.event1, code.event2, code.censoring)))
          stop("COMPETING-RISK requires event codes {censoring,event1,event2}.")
      }
    }
  }

  out_defineExposureDesign <- defineExposureDesign(data, exposure, code.exposure.ref)
  x_a <- out_defineExposureDesign$x_a
  x_l <- model.matrix(out_terms, mf)
  index.vector <- rep(NA, 7)
  index.vector <- calculateIndexForParameter(NA,x_l,x_a)

  if (!is.numeric(conf.level) || length(conf.level) != 1 || conf.level <= 0 || conf.level >= 1)
    stop("conf.level must be a single number between 0 and 1")

  if (outcome.type == "PROPORTIONAL" | outcome.type == "POLY-PROPORTIONAL") {
    should.normalize.covariate.corrected <- FALSE
    report.sandwich.conf.corrected <- FALSE
    if (is.null(report.boot.conf)) {
      report.boot.conf.corrected <- TRUE
    } else {
      report.boot.conf.corrected <- report.boot.conf
    }
  } else {
    should.normalize.covariate.corrected <- should.normalize.covariate
    if (is.null(report.boot.conf)) {
      report.boot.conf.corrected <- FALSE
    } else {
      report.boot.conf.corrected <- report.boot.conf
    }
    if (report.boot.conf == FALSE || is.null(report.boot.conf)) {
      report.sandwich.conf.corrected <- report.sandwich.conf
    } else {
      report.sandwich.conf.corrected <- FALSE
    }
  }

  outer_choices <- c("nleqslv","Newton","Broyden")
  nleqslv.method <- match.arg(nleqslv.method, choices = outer_choices)
  return(list(should.normalize.covariate = should.normalize.covariate.corrected, report.sandwich.conf = report.sandwich.conf.corrected, report.boot.conf = report.boot.conf.corrected, out_defineExposureDesign=out_defineExposureDesign, index.vector=index.vector))
}

read_time.point <- function(formula, data, outcome.type, code.censoring, time.point) {
#  read_time.point <- function(formula, data, outcome.type, exposure, code.censoring, code.exposure.ref, time.point) {
    if (outcome.type %in% c("COMPETING-RISK","SURVIVAL")) {
    if (is.null(time.point) || !length(time.point)) stop("time.point is required when outcome.type is COMPETING-RISK or SURVIVAL.")
    tp <- suppressWarnings(max(time.point, na.rm = TRUE))
    if (!is.finite(tp) || tp < 0) stop("time.point must be non-negative and finite when outcome.type is COMPETING-RISK or SURVIVAL.")
    return(tp)
  } else if (outcome.type == "BINOMIAL") {
    tp <- Inf
    return(tp)
  } else if (outcome.type %in% c("PROPORTIONAL","POLY-PROPORTIONAL") & is.null(time.point)) {
    cl <- match.call()
    mf <- match.call(expand.dots = TRUE)[1:3]
    special <- c("strata", "cluster", "offset")
    Terms <- terms(formula, special, data = data)
    mf$formula <- Terms
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    Y <- model.extract(mf, "response")
    t <- Y[, 1]
    epsilon <- Y[, 2]
    tp <- t[epsilon != code.censoring]
    tp <- sort(unique(tp[is.finite(tp)]))
    return(tp)
  } else {
    return(time.point)
  }
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

clampP <- function(p, eps = 1e-5) pmin(pmax(p, eps), 1 - eps)
clampLogP <- function(x, eps = 1e-5) {
  ex <- exp(x)
  x[ex < eps]         <- log(eps)
  x[(1 - ex) < eps]   <- log(1 - eps)
  x
}

calculateIndexForParameter <- function(i_parameter,x_l,x_a,length.time.point=1) {
  i_parameter[1] <- ncol(x_l)
  i_parameter[2] <- i_parameter[1] + 1
  i_parameter[3] <- i_parameter[1] + ncol(x_a)
  i_parameter[4] <- i_parameter[1] + ncol(x_a) + 1
  i_parameter[5] <- 2 * i_parameter[1] + ncol(x_a)
  i_parameter[6] <- 2 * i_parameter[1] + ncol(x_a) + 1
  i_parameter[7] <- 2 * i_parameter[1] + 2 * ncol(x_a)
  i_parameter[8] <- length.time.point*(2 * i_parameter[1]) + 2 * ncol(x_a)
  return(i_parameter)
}

create_rr_text <- function(coefficient, cov, index, omit.conf.int=TRUE, conf.int=0.95) {
  alpha <- 1 - conf.int
  critical_value <- qnorm(1 - alpha / 2)
  coef <- coefficient[index]
  coef_se <- sqrt(diag(cov)[index])
  conf_low <- coef - critical_value * coef_se
  conf_high <- coef + critical_value * coef_se
  p_value <- floor(2 * (1 - pnorm(abs(coef) / coef_se)))
  if (omit.conf.int==TRUE) {
    if (p_value<0.01) text <- paste0("RR=", round(exp(coef), digit=2), ", p<0.01")
    else text <- paste0("RR=", round(exp(coef), digit=2), ", p=", p_value)
  } else {
    if (p_value<0.01) text <- paste0("RR=", round(exp(coef), digit=2), " (", round(exp(conf_low), digit=2), " to ", round(exp(conf_high), digit=2), ", p<0.01", ")")
    else text <- paste0("RR=", round(exp(coef), digit=2), " (", round(exp(conf_low), digit=2), " to ", round(exp(conf_high), digit=2), ", p=", p_value, ")")
  }
  return(text)
}

createTestData <- function(n, w, first_zero=FALSE, last_zero=FALSE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE) {
  one <- rep(1, n)
  t <- c(1:(n/2), 1:(n/2))
  epsilon <- rep(1, n)
  epsilon[2] <- 2
  epsilon[3] <- 2
  if (first_zero==TRUE) {
    epsilon[1] <- 0
    epsilon[n/2+1] <- 0
  }
  if (last_zero==TRUE) {
    epsilon[n/2] <- 0
    epsilon[n] <- 0
  }
  w <- rep(w, n)
  if (logical_strata==TRUE) {
    strata <- (t %% 2 == 0)
  } else {
    strata <- as.factor((t %% 2 == 0))
  }
  if (na_strata==TRUE) {
    strata[1] <- NA
  }
  subset <- rep(1, n)
  if (subset_present==TRUE) {
    subset[1] <- 0
  }
  d <- as.numeric(epsilon>0)
  return(data.frame(id = 1:n, t = t, epsilon = epsilon, d = d, w = w, strata = strata, subset=subset))
}

normalizeCovariate <- function(formula, data, should.normalize.covariate, outcome.type) {
  mf <- model.frame(formula, data)
  Y <- model.extract(mf, "response")
  response_term <- formula[[2]]
  if (inherits(mf[[1]], "Surv") || inherits(mf[[1]], "Event")) {
    response_vars <- all.vars(response_term)
    covariate_cols <- setdiff(all.vars(formula), response_vars)
  } else {
    covariate_cols <- all.vars(formula)[-1]
  }
  normalized_data <- data
  range_vector <- 1
  if (should.normalize.covariate == TRUE & length(covariate_cols)>0) {
    for (col in covariate_cols) {
      x <- normalized_data[[col]]
      range <- max(x)-min(x)
      normalized_data[[col]] <- x/range
      range_vector <- cbind(range_vector, range)
    }
    if (outcome.type == "PROPORTIONAL" || outcome.type == "POLY-PROPORTIONAL") {
      range_vector <- cbind(range_vector)
    } else if (outcome.type == "SURVIVAL" || outcome.type == "BINOMIAL") {
      range_vector <- cbind(range_vector,1)
    } else {
      range_vector <- cbind(range_vector,1,range_vector,1)
    }
  } else {
    if (outcome.type == "PROPORTIONAL" || outcome.type == "POLY-PROPORTIONAL") {
      range_vector <- rep(1, (length(covariate_cols)+1))
    } else if (outcome.type == "SURVIVAL" || outcome.type == "BINOMIAL") {
      range_vector <- rep(1, (length(covariate_cols)+2))
    } else {
      range_vector <- rep(1, (2*length(covariate_cols)+4))
    }
  }
  n_covariate <- length(covariate_cols)
  out <- list(normalized_data=normalized_data, range=range_vector, n_covariate=n_covariate, n=nrow(data))
  return(out)
}

checkSpell <- function(outcome.type, effect.measure1, effect.measure2) {
  if (outcome.type %in% c("COMPETING-RISK", "COMPETINGRISK", "C", "CR", "COMPETING RISK", "COMPETING-RISKS", "COMPETINGRISKS", "COMPETING RISKS", "Competingrisk", "Competing-risk", "Competing risk", "Competingrisks", "Competing-risks", "Competing risks", "competing-risk", "competingrisk", "competing risk", "competing-risks", "competingrisks", "competing risks")) {
    outcome.type.corrected <- "COMPETING-RISK"
  } else if (outcome.type %in% c("SURVIVAL", "S", "Survival", "Survival")) {
    outcome.type.corrected <- "SURVIVAL"
  } else if (outcome.type %in% c("POLY-PROPORTIONAL", "PP", "Poly-proportional", "poly-proportional")) {
    outcome.type.corrected <- "POLY-PROPORTIONAL"
  } else if (outcome.type %in% c("PROPORTIONAL", "P", "Proportional", "proportional")) {
    outcome.type.corrected <- "PROPORTIONAL"
  } else if (outcome.type %in% c("BINOMIAL", "B", "Binomial", "binomial")) {
    outcome.type.corrected <- "BINOMIAL"
  } else {
    stop("Invalid input for outcome.type, Choose 'COMPETING-RISK', 'SURVIVAL', 'BINOMIAL', 'PROPORTIONAL', or 'POLY-PROPORTIONAL'.")
  }
  if (effect.measure1 %in% c("RR", "rr", "RISK RATIO", "Risk ratio", "risk ratio")) {
    effect.measure1.corrected <- "RR"
  } else if (effect.measure1 %in% c("OR", "or", "ODDS RATIO", "Odds ratio", "odds ratio")) {
    effect.measure1.corrected <- "OR"
  } else if (effect.measure1 %in% c("SHR", "shr", "HR", "hr", "SUBDISTRIBUTION HAZARD RATIO",
                                    "Subdistibution hazard ratio", "subdistibution hazard ratio")) {
    effect.measure1.corrected <- "SHR"
  } else {
    stop("Invalid input for effect.measure1, Choose 'RR', 'OR', or 'SHR'.")
  }
  if (effect.measure2 %in% c("RR", "rr", "RISK RATIO", "Risk ratio", "risk ratio")) {
    effect.measure2.corrected <- "RR"
  } else if (effect.measure2 %in% c("OR", "or", "ODDS RATIO", "Odds ratio", "odds ratio")) {
    effect.measure2.corrected <- "OR"
  } else if (effect.measure2 %in% c("SHR", "shr", "HR", "hr", "SUBDISTRIBUTION HAZARD RATIO",
                                    "Subdistibution hazard ratio", "subdistibution hazard ratio")) {
    effect.measure2.corrected <- "SHR"
  } else {
    stop("Invalid input for effect.measure2, Choose 'RR', 'OR', or 'SHR'.")
  }
  return(list(outcome.type = outcome.type.corrected, effect.measure1 = effect.measure1.corrected, effect.measure2 = effect.measure2.corrected))
}


`%||%` <- function(x, y) if (is.null(x)) y else x
extractOptimizationInfo <- function(sol, method) {
  out <- list(solver = method)
  if (method %in% c("nleqslv","Broyden","Newton")) {
    out$code    <- sol$termcd %||% NA_integer_
    out$message <- sol$message %||% NA_character_
    out$fn.evals <- sol$feval %||% NA_integer_
    out$jac.evals <- sol$jeval %||% NA_integer_
  } else if (method == "multiroot") {
    out$message <- sol$estim.precis %||% NA_character_
    out$iter    <- sol$iter %||% NA_integer_
    out$estim.precis <- sol$estim.precis %||% NA_character_
  } else if (method %in% c("optim","SANN","BFGS")) {
    out$code       <- sol$convergence %||% NA_integer_
    out$message    <- sol$message %||% NA_character_
    if (!is.null(sol$counts)) {
      out$fn.evals <- sol$counts["function"] %||% NA_integer_
      out$gr.evals <- sol$counts["gradient"] %||% NA_integer_
    }
  }
  return(out)
}

append_trace <- function(trace_df, iteration, computation.time.second = NA_real_, nleqslv.method, objective.function, relative.difference,
                         max.absolute.difference, converged.by, nleqslv.info, store_params = FALSE,
                         coefficient = NULL) {
  row <- data.frame(
    iteration = iteration,
    computation.time.second = computation.time.second,
    nleqslv.method = nleqslv.method,
    objective.function = objective.function,
    relative.difference = relative.difference,
    max.absolute.difference = max.absolute.difference,
    converged.by = if (isTRUE(converged.by)) NA_character_ else as.character(converged.by),
    code = nleqslv.info$code %||% NA_integer_,
    msg  = nleqslv.info$message %||% NA_character_,
    fn_evals = nleqslv.info$fn.evals %||% NA_integer_,
    gr_evals = nleqslv.info$gr.evals %||% NA_integer_,
    jac_evals = nleqslv.info$jac.evals %||% NA_integer_,
    stringsAsFactors = FALSE
  )
  if (isTRUE(store_params)) {
    row$coefficient <- list(coefficient)
  }
  if (is.null(trace_df)) return(row)
  rbind(trace_df, row)
}


fcif <- function(cif,rr,type=c("cif","cloglog","logistic")) {
  mcif <- max(cif[,2])
  if (type[1]=="cif") mcif <- mcif*rr
  if (type[1]=="cloglog") mcif <- 1- exp(-mcif*rr)
  if (type[1]=="logistic") mcif <- mcif* rr/(1 + mcif * rr)
  return(mcif)
}

simRR <- function(n,lrr1,lrr2,cens=NULL,type1=c("cif","cloglog","logistic"),type2=c("cif","cloglog","logistic")) {
  A <- rbinom(n,1,0.5)
  L <- rbinom(n,1,0.5)
  rr1 <- exp(cbind(A,L) %*% lrr1)
  rr2 <- exp(cbind(A,L) %*% lrr2)
  f1 <- fcif(cif1,max(rr1),type=type1[1])
  f2 <- fcif(cif2,max(rr2),type=type2[1])
  mmm<- f1+f2
  mcif1 <- fcif(cif1,rr1,type=type1[1])
  mcif2 <- fcif(cif2,rr2,type=type2[1])
  if (mmm>1) warning(" models not satisfying sum <=1\n")
  T1 <- simsubdist(cif1,rr1,type=type1[1])
  T2 <- simsubdist(cif2,rr2,type=type2[1])
  dies <- rbinom(n,1,mcif1+mcif2)
  sel1 <- rbinom(n,1,mcif2/(mcif1+mcif2))+1
  epsilon  <- dies*(sel1)
  T1$epsilon <- epsilon
  T1$A <- A
  T1$L <- L
  T1$time <- T1$timecause
  T1$time2 <- T2$timecause
  T1$status <- epsilon
  T1 <- dtransform(T1,time=100,epsilon==0)
  T1 <- dtransform(T1,status=0,epsilon==0)
  T1 <- dtransform(T1,time=time2,epsilon==2)
  T1 <- dtransform(T1,status=2,epsilon==2)
  data <- T1

  if (!is.null(cens))  {
    cc <- rexp(n)/cens
    data$status <- ifelse(data$time<cc,data$status,0)
    data$time <- pmin(data$time,cc)
  }
  return(data)
}
