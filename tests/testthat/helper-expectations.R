# Reduce values to something small and platform-stable before snapshotting.
# Never snapshot a fitted model or a full prediction matrix.

snap_num <- function(x, digits = 6) round(as.numeric(x), digits)

mat_fingerprint <- function(m, digits = 6) {
  m <- as.matrix(m)
  list(
    dim  = dim(m),
    min  = round(min(m, na.rm = TRUE), digits),
    max  = round(max(m, na.rm = TRUE), digits),
    mean = round(mean(m, na.rm = TRUE), digits),
    n_na = sum(is.na(m))
  )
}

# Structural assertions common to every Metric/Value result table.
expect_metric_table <- function(res) {
  testthat::expect_true(is.data.frame(res))
  testthat::expect_true("Metric" %in% names(res))
  testthat::expect_gt(nrow(res), 0)
  testthat::expect_type(res$Metric, "character")
  testthat::expect_false(any(duplicated(res$Metric)))
  invisible(res)
}
