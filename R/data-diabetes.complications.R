#' Diabetes Complications Dataset
#'
#' A simulated dataset representing event times and covariates for diabetic patients.
#'
#' @format A data frame with 500 rows and 7 variables:
#' \describe{
#'   \item{t}{event or censoring time}
#'   \item{epsilon}{event indicator (1 = event, 0 = censored)}
#'   \item{fruitq1}{binary exposure variable}
#'   \item{x1}{covariate 1}
#'   \item{x2}{covariate 2}
#'   \item{x3}{covariate 3}
#'   \item{offset}{offset term}
#' }
#'
#' @source Simulated data
#' @keywords data
"diabetes.complications"
