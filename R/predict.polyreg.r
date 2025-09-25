#' @title Predicting cumulative incidence probabilities at a specific time point from direct polynomial regression
#'
#' @param formula formula Model formula representing outcome, exposure and covariates
#' @param exposure character Column name representing the exposure (1 = exposed, 0 = not exposed).
#' @param data data.frame Input dataset containing survival data.
#' @param coefficient numeric The coefficient of direct polynomial regression.
#' @param effect.measure1 character Specifies the effect measure for event (RR, OR, SHR).
#' @param effect.measure2 character Specifies the effect measure for competing risk (RR, OR, SHR).
#' @param outcome.type character Specifies the type of outcome (COMPETINGRISK or SURVIVAL).
#' @param inner.optim.method character Specifies the method of optimization (nleqslv, optim, multiroot).
#' @param prob.bound numeric Small threshold for clamping probabilities. Defaults to 1e-5.
#'
#' @return A list of predicted cumulative incidence probabilities according to exposure and cause from direct polynomial regression.
#' @export
#'
#' @examples
#' data(bmt)
#' result <- polyreg(nuisance.model = Event(time, cause)~age+tcell, exposure = 'platelet',
#' cens.model = Event(time,cause==0)~+1, data = bmt, effect.measure1='RR', effect.measure2='RR', time.point=24, outcome.type='COMPETINGRISK')
#' msummary(result$out_summary, statistic = c("conf.int"), exponentiate = TRUE)
#' prediction <- predict.polyreg(formula = Event(time, cause)~age+tcell, exposure = 'platelet', data = bmt, coefficient = result$out_coefficient, effect.measure1='RR', effect.measure2='RR', outcome.type='COMPETINGRISK')


predict.polyreg <- function(
    formula,
    exposure,
    data,
    coefficient,
    effect.measure1,
    effect.measure2,
    outcome.type,
    inner.optim.method = "optim",
    prob.bound = 1e-5
) {
  spell_check(outcome.type, effect.measure1, effect.measure2)
  cl <- match.call()
  mf <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "cluster", "offset")
  Terms <- terms(formula, special, data = data)
  Terms <- delete.response(Terms)
  mf <- model.frame(Terms, data = data)
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
  if (!nrow(data) == nrow(mf))
    stop("Variables contain NA values")
  if (any(is.na(data[[exposure]])))
    stop("Variables contain NA values")

  a_ <- as.factor(data[[exposure]])
  a <- model.matrix(~ a_)[, 2]

  one <- rep(1, nrow(mf))
  effect.modifier <- NULL
  if (!is.null(effect.modifier)) {
    em <- data[[effect.modifier]]
    a_interaction <- as.matrix(cbind(a, a*em))
    x_em <- as.matrix(cbind(one, em))
  } else {
    a_interaction <- a
    x_em <- as.matrix(one)
  }
  x_l <- model.matrix(Terms, mf)
  x_la <- cbind(x_l, a_interaction)

  optim.method <- list(inner.optim.method=inner.optim.method, optim.parameter7=1e-7, optim.parameter8=200)
  estimand <- list(effect.measure1=effect.measure1, effect.measure2=effect.measure2)
  if (outcome.type == 'SURVIVAL') {
    n_para_1 <- ncol(x_l)
    n_para_2 <- n_para_1+1
    if (estimand$effect.measure1 == 'RR') {
      alpha_beta_ <- as.matrix(coefficient)
      if (ncol(alpha_beta_ == 1)) {
        alpha_beta_ <- t(alpha_beta_)
      }
      phi <- x_l %*% alpha_beta_[, 1:n_para_1] + offset
      theta <- one * alpha_beta_[, n_para_2]
      expphi <- exp(phi)
      exptheta <- exp(theta)
      if (all(phi == 0)) {
        p_10 <- one / (one + exptheta)
        p_11 <- exptheta * p_10
      } else {
        denomi_1 <- -(exptheta + one) * expphi
        denomi_2 <- sqrt(exp(2 * phi) * (exptheta + one) * (exptheta + one) + 4 * exp(theta + phi) * (one - expphi))
        denomi <- denomi_1 + denomi_2
        numera <- 2 * exptheta * (one - expphi)
        p_10 <- denomi / numera
        p_11 <- exptheta * p_10
      }
    } else if (estimand$effect.measure1 == 'OR') {
      alpha_beta_ <- as.matrix(coefficient)
      if (ncol(alpha_beta_ == 1)) {
        alpha_beta_ <- t(alpha_beta_)
      }
      phi <- x_l %*% alpha_beta_[, 1:n_para_1] + offset
      theta <- one * alpha_beta_[, n_para_2]
      sqrt1 <- sqrt(exp(-theta - phi))
      sqrt2 <- sqrt(exp(theta - phi))
      if (all(phi == theta)) {
        p_10 <- 0.5 * one
        p_11 <- one/(one + sqrt1)
      } else {
        p_10 <- one/(one + sqrt2)
        p_11 <- one/(one + sqrt1)
      }
    } else {
      stop("Invalid effect_measure. Must be RR or OR.")
    }
    return(cbind(p_10, p_11))
  }

  i_parameter <- rep(NA, 7)
  i_parameter[1] <- ncol(x_l)
  i_parameter[2] <- i_parameter[1] + 1
  i_parameter[3] <- i_parameter[1] + ncol(x_em)
  i_parameter[4] <- i_parameter[1] + ncol(x_em) + 1
  i_parameter[5] <- 2 * i_parameter[1] + ncol(x_em)
  i_parameter[6] <- 2 * i_parameter[1] + ncol(x_em) + 1
  i_parameter[7] <- 2 * i_parameter[1] + 2 * ncol(x_em)
  alpha_1 <- coefficient[1:i_parameter[1]]
  beta_1  <- coefficient[i_parameter[2]:i_parameter[3]]
  alpha_2 <- coefficient[i_parameter[4]:i_parameter[5]]
  beta_2  <- coefficient[i_parameter[6]:i_parameter[7]]
  alpha_tmp_1_vals <- x_l %*% as.matrix(alpha_1) + offset
  alpha_tmp_2_vals <- x_l %*% as.matrix(alpha_2) + offset
  beta_tmp_1_vals  <- x_em %*% as.matrix(beta_1)
  beta_tmp_2_vals  <- x_em %*% as.matrix(beta_2)

  pred_list <- vector("list", nrow(x_l))

  # Initialize previous prediction values
  prev_pred <- NULL
  # Loop through each observation in the sorted x_l
  for (i_x in seq_len(nrow(x_l))) {
    # Skip the calculation if the current observation is the same as the previous one
    if (!is.null(prev_pred) && all(x_l[i_x, ] == x_l[i_x-1, ])) {
      pred_list[[i_x]] <- prev_pred
      next
    }

    log_p0 <- log(c(0.5, 0.5, 0.5, 0.5))
    if (optim.method$inner.optim.method == 'optim' || optim.method$inner.optim.method == 'BFGS') {
      eq_fn <- function(lp) {
        estimating_equation_pred(
          log_p          = lp,
          alpha_tmp_1    = alpha_tmp_1_vals[i_x],
          beta_tmp_1     = beta_tmp_1_vals[i_x],
          alpha_tmp_2    = alpha_tmp_2_vals[i_x],
          beta_tmp_2     = beta_tmp_2_vals[i_x],
          estimand       = estimand,
          optim.method   = optim.method,
          prob.bound     = prob.bound
        )
      }
      sol <- optim(par = log_p0, fn = eq_fn, method = "BFGS", control = list(maxit = optim.method$optim.parameter8, reltol = optim.method$optim.parameter7))
      pred_list[[i_x]] <- exp(sol$par)
      prev_pred <- pred_list[[i_x]]
    } else if (optim.method$inner.optim.method == 'SANN') {
      eq_fn <- function(lp) {
        estimating_equation_pred(
          log_p          = lp,
          alpha_tmp_1    = alpha_tmp_1_vals[i_x],
          beta_tmp_1     = beta_tmp_1_vals[i_x],
          alpha_tmp_2    = alpha_tmp_2_vals[i_x],
          beta_tmp_2     = beta_tmp_2_vals[i_x],
          estimand       = estimand,
          optim.method   = optim.method,
          prob.bound     = prob.bound
        )
      }
      sol <- optim(par = log_p0, fn = eq_fn, method = "SANN", control = list(maxit = optim.method$optim.parameter8, reltol = optim.method$optim.parameter7))
      pred_list[[i_x]] <- exp(sol$par)
      prev_pred <- pred_list[[i_x]]
    } else if (optim.method$inner.optim.method == 'multiroot'){
      eq_fn <- function(lp) {
        estimating_equation_pred(
          log_p          = lp,
          alpha_tmp_1    = alpha_tmp_1_vals[i_x],
          beta_tmp_1     = beta_tmp_1_vals[i_x],
          alpha_tmp_2    = alpha_tmp_2_vals[i_x],
          beta_tmp_2     = beta_tmp_2_vals[i_x],
          estimand = estimand,
          optim.method = optim.method,
          prob.bound = prob.bound
        )
      }
      sol <- multiroot(eq_fn, start = log_p0)
      pred_list[[i_x]] <- exp(sol$root)
      prev_pred <- pred_list[[i_x]]
    }
  }
  out <- do.call(rbind, pred_list)
  return(out)
}
