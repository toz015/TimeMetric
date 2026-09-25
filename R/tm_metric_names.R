# Canonical metric identifiers.
#
# Before standardization the same metric could carry up to five different
# spellings depending on which entry point produced it, several of them
# misspelled ("Pesudo_R", "Psuedo.R") or containing a U+2019 curly apostrophe
# that users had to reproduce exactly for metric selection to match.
#
# The canonical set is ASCII, lower case, and one name per quantity. Legacy
# spellings are still accepted at the argument boundary and resolve with a
# deprecation warning naming the replacement.
#
# Pseudo_R_square and Pseudo_R2_point are genuinely distinct -- an integrated
# measure and a point-in-time estimate, differing numerically on the same data
# -- so both survive under distinct names.
#
# r_e and r_sh were withdrawn before the first CRAN release. Their spellings
# are retained in tm_defunct_metric_spellings() so that a request for one
# reports the withdrawal by name, instead of falling through to a generic
# unknown-metric error.

#' @keywords internal
#' @noRd
tm_metric_aliases <- function() {
  c(
    # pseudo R-squared family
    "Pseudo_R_square"              = "pseudo_r2",
    "Pesudo_R"                     = "pseudo_r2",
    "Psuedo.R"                     = "pseudo_r2",
    "pseudo_r2"                    = "pseudo_r2",
    "Pseudo_R2_point"              = "pseudo_r2_point",
    "pseudo_r2_point"              = "pseudo_r2_point",
    # explained variation components
    "R_square"                     = "r_square",
    "r_square"                     = "r_square",
    "L_square"                     = "l_square",
    "l_square"                     = "l_square",
    "R2_point"                     = "r2_point",
    "r2_point"                     = "r2_point",
    "L2_point"                     = "l2_point",
    "l2_point"                     = "l2_point",
    # concordance
    "Harrells_C"                   = "harrell_c",
    "Harrell\u2019s C"             = "harrell_c",
    "Harrell's C"                  = "harrell_c",
    "harrell_c"                    = "harrell_c",
    "Unos_C"                       = "uno_c",
    "Uno\u2019s C"                 = "uno_c",
    "Uno's C"                      = "uno_c",
    "uno_c"                        = "uno_c",
    "C_index"                      = "c_index",
    "c_index"                      = "c_index",
    # calibration and discrimination over time
    "Brier Score"                  = "brier_score",
    "Brier_Score"                  = "brier_score",
    "brier_score"                  = "brier_score",
    "Time Dependent Auc"           = "td_auc",
    "Time Dependent AUC"           = "td_auc",
    "Time_Dependent_Auc"           = "td_auc",
    "AUC"                          = "td_auc",
    "td_auc"                       = "td_auc"
  )
}

#' Canonical metric names
#'
#' The metric identifiers accepted by the `metrics` argument of the evaluation
#' functions, and emitted in the `Metric` column of their results.
#'
#' @return A character vector of canonical metric names.
#' @examples
#' tm_metric_names()
#' @export
tm_metric_names <- function() {
  sort(unique(unname(tm_metric_aliases())))
}

# Metric spellings withdrawn before the first CRAN release. Checked before
# alias resolution so that both canonical and legacy spellings are reported.
#' @keywords internal
#' @noRd
tm_defunct_metrics <- function() {
  c("r_e", "r_sh")
}

#' @keywords internal
#' @noRd
tm_defunct_metric_spellings <- function() {
  c("r_e", "r_sh", "R_E", "R_sh", "R_sph")
}

# Resolve any accepted spelling to its canonical form. Unknown names are
# returned unchanged so the caller's own validation can report them.
#' @keywords internal
#' @noRd
tm_normalize_metrics <- function(metrics, warn = TRUE) {
  if (is.null(metrics)) return(NULL)

  defunct_key <- tolower(gsub("[ _]", "", metrics))
  if (any(defunct_key %in% tolower(gsub("[ _]", "",
                                        tm_defunct_metric_spellings())))) {
    stop("r_e and r_sh were withdrawn before the first CRAN release; ",
         "see NEWS.md.", call. = FALSE)
  }

  aliases <- tm_metric_aliases()

  # case-insensitive, and treat _ / space / straight or curly apostrophe alike
  flatten <- function(x) {
    x <- gsub("\u2019", "'", x)
    tolower(gsub("[ _]", "", x))
  }
  lookup <- stats::setNames(unname(aliases), flatten(names(aliases)))

  key <- flatten(metrics)
  hit <- key %in% names(lookup)
  out <- metrics
  out[hit] <- unname(lookup[key[hit]])

  if (warn) {
    legacy <- hit & (metrics != out)
    if (any(legacy)) {
      pairs <- paste0("'", metrics[legacy], "' -> '", out[legacy], "'",
                      collapse = ", ")
      warning("Deprecated metric name(s): ", pairs,
              ". The legacy spellings still work but will be removed in a ",
              "future release; see tm_metric_names().",
              call. = FALSE)
    }
  }
  out
}
