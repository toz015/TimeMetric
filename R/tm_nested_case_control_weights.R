# -------------------------------------------------------------------
# PUBLIC: NCC weights (handles matched if `strata` given)
# -------------------------------------------------------------------
#' Nested case-control (NCC) sampling weights
#'
#' Computes inverse-selection weights for NCC designs. If `strata` is
#' provided, matched NCC weights are computed.
#'
#' @param time Numeric vector of follow-up times.
#' @param status Numeric \{0,1\}, 1 = event, 0 = censored.
#' @param strata Optional factor/character for matching sets (matched NCC).
#' @param m Integer, number of controls per case (required).
#' @return Numeric vector of weights.
#' @export
#' @examples
#' # Nested case-control design with m = 2 matched controls per case.
#' d <- tm_sim_cox_weibull(n = 150, pi_c = 0.3, v = 2,
#'                         beta = c(0.5, -0.5), seed = 2025)
#' set.seed(2025)
#' w <- tm_nested_case_control_weights(time = d$time, status = d$status, m = 2)
#' summary(w)
tm_nested_case_control_weights <- function(time, status, strata = NULL, m = NULL) {
  if (is.null(m)) stop("`m` must be provided for NCC weights.")
  design <- if (is.null(strata)) "ncc" else "matched_ncc"
  weighted_param(time = time, status = status, design = design,
                 subcohort = NULL, strata = strata, m = m)
}

