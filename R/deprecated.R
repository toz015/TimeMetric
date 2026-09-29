#' Deprecated function names
#'
#' @description
#' The `pam.` prefix carried over from the predecessor PAmeasure package and
#' contained two misspellings (`survial` for `survival`, `surverg` for
#' `survreg`). Every exported function now uses a `tm_` prefix instead.
#'
#' Each old name below still works and forwards to its replacement, emitting a
#' deprecation warning that names the new function. They exist so that existing
#' analysis scripts -- including `paper.code.Rmd` and the reproduction code for
#' Zhuang et al. (2025) -- keep running unchanged.
#'
#' @param ... Passed unchanged to the replacement function.
#'
#' @return Whatever the replacement function returns.
#'
#' @name TimeMetric-deprecated
NULL

#' @rdname TimeMetric-deprecated
#' @details `pam.predicted_survial_eval()` is deprecated; use [tm_survival_eval()].
#' @export
pam.predicted_survial_eval <- function(...) {
  .Deprecated("tm_survival_eval", package = "TimeMetric")
  tm_survival_eval(...)
}

#' @rdname TimeMetric-deprecated
#' @details `pam.predicted_survial_eval_cr()` is deprecated; use [tm_survival_eval_cr()].
#' @export
pam.predicted_survial_eval_cr <- function(...) {
  .Deprecated("tm_survival_eval_cr", package = "TimeMetric")
  tm_survival_eval_cr(...)
}

#' @rdname TimeMetric-deprecated
#' @details `pam.survival_eval()` is deprecated; use [tm_fit_and_eval()].
#' @export
pam.survival_eval <- function(...) {
  .Deprecated("tm_fit_and_eval", package = "TimeMetric")
  tm_fit_and_eval(...)
}

#' @rdname TimeMetric-deprecated
#' @details `pam.coxph_restricted()` is deprecated; use [tm_predict_coxph()].
#' @export
pam.coxph_restricted <- function(...) {
  .Deprecated("tm_predict_coxph", package = "TimeMetric")
  tm_predict_coxph(...)
}

#' @rdname TimeMetric-deprecated
#' @details `pam.surverg_restricted()` is deprecated; use [tm_predict_survreg()].
#' @export
pam.surverg_restricted <- function(...) {
  .Deprecated("tm_predict_survreg", package = "TimeMetric")
  tm_predict_survreg(...)
}

#' @rdname TimeMetric-deprecated
#' @details `pam.predict_cr()` is deprecated; use [tm_predict_cif()].
#' @export
pam.predict_cr <- function(...) {
  .Deprecated("tm_predict_cif", package = "TimeMetric")
  tm_predict_cif(...)
}

#' @rdname TimeMetric-deprecated
#' @details `pam.summary()` is deprecated; use [tm_summarize()].
#' @export
pam.summary <- function(...) {
  .Deprecated("tm_summarize", package = "TimeMetric")
  tm_summarize(...)
}

#' @rdname TimeMetric-deprecated
#' @details `pam.summary_cr()` is deprecated; use [tm_summarize_cr()].
#' @export
pam.summary_cr <- function(...) {
  .Deprecated("tm_summarize_cr", package = "TimeMetric")
  tm_summarize_cr(...)
}

#' @rdname TimeMetric-deprecated
#' @details `pam.sample_design()` is deprecated; use [tm_sample_design()].
#' @export
pam.sample_design <- function(...) {
  .Deprecated("tm_sample_design", package = "TimeMetric")
  tm_sample_design(...)
}

#' @rdname TimeMetric-deprecated
#' @details `cc_weights()` is deprecated; use [tm_case_cohort_weights()].
#' @export
cc_weights <- function(...) {
  .Deprecated("tm_case_cohort_weights", package = "TimeMetric")
  tm_case_cohort_weights(...)
}

#' @rdname TimeMetric-deprecated
#' @details `ncc_weights()` is deprecated; use [tm_nested_case_control_weights()].
#' @export
ncc_weights <- function(...) {
  .Deprecated("tm_nested_case_control_weights", package = "TimeMetric")
  tm_nested_case_control_weights(...)
}

#' @rdname TimeMetric-deprecated
#' @details `plot_pred()` is deprecated; use [tm_plot_pred()].
#' @export
plot_pred <- function(...) {
  .Deprecated("tm_plot_pred", package = "TimeMetric")
  tm_plot_pred(...)
}

#' @rdname TimeMetric-deprecated
#' @details `summary_pred_plot()` is deprecated; use [tm_plot_summary()].
#' @export
summary_pred_plot <- function(...) {
  .Deprecated("tm_plot_summary", package = "TimeMetric")
  tm_plot_summary(...)
}

#' @rdname TimeMetric-deprecated
#' @details `sim_cox_weibull_censored()` is deprecated; use [tm_sim_cox_weibull()].
#' @export
sim_cox_weibull_censored <- function(...) {
  .Deprecated("tm_sim_cox_weibull", package = "TimeMetric")
  tm_sim_cox_weibull(...)
}

#' @rdname TimeMetric-deprecated
#' @details `simulateTwoCauseFineGrayModel()` is deprecated; use [tm_simulate_fine_gray()].
#' @export
simulateTwoCauseFineGrayModel <- function(...) {
  .Deprecated("tm_simulate_fine_gray", package = "TimeMetric")
  tm_simulate_fine_gray(...)
}
