#' TimeMetric: Predictive Performance Metrics for Survival Models
#'
#' Unified performance metrics for survival models with right-censoring,
#' competing risks, and two-phase designs such as nested case-control and
#' case-cohort studies.
#'
#' @section Optional backends:
#' Two prediction backends are optional and declared in `Suggests` rather than
#' `Imports`, because each serves a single branch of [tm_predict_cif()]:
#'
#' * `randomForestSRC` supplies the `cr_model` argument (an `rfsrc` fit).
#' * `cmprsk` supplies the `fg_model` argument (a Fine-Gray `crr` fit).
#'
#' Both branches check availability with `requireNamespace()` and raise an
#' informative error naming the package to install.
#'
#' @keywords internal
#' @name TimeMetric-package
#' @aliases TimeMetric
#'
#' @importFrom survival Surv survfit concordancefit
#' @importFrom stats median na.omit pnorm predict uniroot sd runif rnorm
#' @importFrom stats quantile complete.cases reshape
#' @importFrom stats lm as.formula rbinom rweibull
#' @importFrom utils head
"_PACKAGE"

# Non-standard evaluation: these are column names referenced inside dplyr verbs
# and ggplot2 aes() mappings, not undefined globals. Declaring them silences
# "no visible binding for global variable" from R CMD check.
utils::globalVariables(c(
  ".", ".pred", "pred", "surv_obj", "weight_time", "x_var", "times", "status"
))
