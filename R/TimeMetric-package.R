#' TimeMetric: Predictive Performance Metrics for Survival Models
#'
#' Unified performance metrics for survival models with right-censoring,
#' competing risks, and two-phase designs such as nested case-control and
#' case-cohort studies.
#'
#' @section Optional backends:
#' Two prediction backends are optional and declared in `Suggests` rather than
#' `Imports`, because each serves a single branch of [pam.predict_cr()]:
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
#' @importFrom stats quantile complete.cases reshape approx
#' @importFrom stats lm as.formula model.matrix
#' @importFrom utils head
"_PACKAGE"
