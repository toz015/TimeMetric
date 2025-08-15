# -------------------------------------------------------------------
# PUBLIC: Case–cohort weights (handles stratified if `strata` given)
# -------------------------------------------------------------------
#' Case–cohort sampling weights
#'
#' Computes Prentice-style weights for case–cohort designs. If `strata`
#' is provided, stratified case–cohort weights are computed.
#'
#' @param time Numeric vector of follow-up times.
#' @param status Numeric {0,1}, 1 = event, 0 = censored.
#' @param subcohort Logical vector; TRUE if in subcohort.
#' @param strata Optional factor/character for strata (stratified cc).
#' @return Numeric vector of weights.
#' @export
cc_weights <- function(time, status, subcohort = NULL, strata = NULL) {
  design <- if (is.null(strata)) "casecohort" else "strat_casecohort"
  weighted_param(time = time, status = status, design = design,
                 subcohort = subcohort, strata = strata, m = NULL)
}
