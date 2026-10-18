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
#' @param strata Optional factor/character defining matching strata in the full
#'   phase-one cohort, not IDs of sampled case-control sets.
#' @param m Integer >= 2, total size of each case-control set: one case plus
#'   m - 1 controls (required).
#' @details Supply full-cohort follow-up before selecting sampled rows. All
#'   cases receive weight 1; non-cases receive the inverse probability of ever
#'   being sampled. Only same-stratum case events contribute in matched designs.
#'   Sampling is without replacement within an occasion, independently across
#'   cases; small risk sets supply all available controls. The eligibility rule
#'   is time > case time. Ties with censoring therefore need a specified ordering
#'   or suitably resolved times. A non-case with no eligible sampling occasion
#'   has weight Inf and cannot be treated as an observed phase-two control.
#' @return Numeric vector of weights, in full-cohort input order.
#' @export
#' @examples
#' # Unmatched NCC: m = 2 means one case and one control per set.
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

