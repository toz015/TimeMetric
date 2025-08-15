# -------------------------------------------------------------------
# PUBLIC: NCC weights (handles matched if `strata` given)
# -------------------------------------------------------------------
#' Nested case–control (NCC) sampling weights
#'
#' Computes inverse-selection weights for NCC designs. If `strata` is
#' provided, matched NCC weights are computed.
#'
#' @param time Numeric vector of follow-up times.
#' @param status Numeric {0,1}, 1 = event, 0 = censored.
#' @param strata Optional factor/character for matching sets (matched NCC).
#' @param m Integer, number of controls per case (required).
#' @return Numeric vector of weights.
#' @export
ncc_weights <- function(time, status, strata = NULL, m = NULL) {
  if (is.null(m)) stop("`m` must be provided for NCC weights.")
  design <- if (is.null(strata)) "ncc" else "matched_ncc"
  weighted_param(time = time, status = status, design = design,
                 subcohort = NULL, strata = strata, m = m)
}

