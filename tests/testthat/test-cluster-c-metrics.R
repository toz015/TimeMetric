test_that("pam.rsh_metric returns D, Dx, R_sh rounded to 4 places", {
  d <- fx_surv()
  pred <- seq(0.9, 0.1, length.out = nrow(d))

  res <- TimeMetric:::pam.rsh_metric(pred, d$time, d$status)

  expect_type(res, "list")
  expect_identical(names(res), c("D", "Dx", "R_sh"))
  expect_length(res$D, 1L)
  expect_length(res$Dx, 1L)
  expect_length(res$R_sh, 1L)
  expect_true(is.finite(res$D))
  expect_true(is.finite(res$Dx))
  expect_gt(res$D, 0)
  # the function rounds its outputs; confirm that is still true
  expect_identical(res$D, round(res$D, 4))
  expect_identical(res$Dx, round(res$Dx, 4))

  expect_snapshot_value(snap_num(c(res$D, res$Dx, res$R_sh)), style = "serialize")
})

test_that("pam.rsh_metric is sensitive to input order (pins FINDING 1)", {
  # The function sorts its internal data frame by survival_time but multiplies
  # by the UNSORTED predicted_data argument. Permuting all three inputs
  # consistently should leave the result unchanged. It does not. Pinned as-is.
  d <- fx_surv()
  pred <- seq(0.9, 0.1, length.out = nrow(d))
  ord <- order(d$time)

  sorted   <- TimeMetric:::pam.rsh_metric(pred[ord], d$time[ord], d$status[ord])
  shuffled <- TimeMetric:::pam.rsh_metric(pred, d$time, d$status)

  # D depends only on time and status, so it is invariant
  expect_equal(sorted$D, shuffled$D, tolerance = 1e-6)
  # Dx depends on the misaligned predictions, so it is not
  expect_false(isTRUE(all.equal(sorted$Dx, shuffled$Dx, tolerance = 1e-6)))

  expect_snapshot_value(
    snap_num(c(sorted$Dx, shuffled$Dx)), style = "serialize"
  )
})

test_that("pam.rsph_metric takes (time, status, risk_score) and returns r2", {
  d <- fx_surv()
  risk <- as.numeric(predict(fx_cox(), newdata = d, type = "lp"))

  res <- TimeMetric:::pam.rsph_metric(d$time, d$status, risk)

  expect_type(res, "list")
  expect_identical(sort(names(res)), c("denominator", "numerator", "r2"))
  expect_length(res$r2, 1L)
  expect_true(is.finite(res$r2))
  expect_gt(res$denominator, 0)
  # r2 is the ratio of the other two components
  expect_equal(res$r2, res$numerator / res$denominator, tolerance = 1e-6)

  expect_snapshot_value(
    snap_num(c(res$r2, res$numerator, res$denominator)), style = "serialize"
  )
})

test_that("pam.rsph_metric rejects mismatched input lengths", {
  d <- fx_surv()
  risk <- as.numeric(predict(fx_cox(), newdata = d, type = "lp"))

  # guarded by stopifnot(length(time) == length(status), ...)
  expect_error(
    TimeMetric:::pam.rsph_metric(d$time[-1], d$status, risk),
    "length"
  )
  expect_error(
    TimeMetric:::pam.rsph_metric(d$time, d$status, risk[-1]),
    "length"
  )
})

test_that("pam.Brier_metric validates its inputs with specific messages", {
  d <- fx_surv()
  pred <- rep(0.5, nrow(d))
  sv <- survival::Surv(d$time, d$status)

  expect_error(
    TimeMetric:::pam.Brier_metric(pred, d$time),
    "must be a survival object created using Surv"
  )
  expect_error(
    TimeMetric:::pam.Brier_metric(pred[-1], sv),
    "Length of suvival_time and predicted_data must match"
  )
  expect_error(
    TimeMetric:::pam.Brier_metric(c(NA_real_, pred[-1]), sv),
    "cannot have NA"
  )
})

test_that("pam.Brier_metric returns a finite scalar in [0, 1]", {
  d <- fx_surv()
  sv <- survival::Surv(d$time, d$status)
  pred <- seq(0.9, 0.1, length.out = nrow(d))

  res <- TimeMetric:::pam.Brier_metric(pred, sv)

  expect_length(res, 1L)
  expect_true(is.finite(res))
  expect_gte(res, 0)
  expect_lte(res, 1)
  expect_snapshot_value(snap_num(res), style = "serialize")
})

test_that("pam.Brier_metric default t_star is the median observed time", {
  d <- fx_surv()
  sv <- survival::Surv(d$time, d$status)
  pred <- seq(0.9, 0.1, length.out = nrow(d))

  # source: if (t_star < 0) t_star <- median(time), where time is the
  # sorted observed time from the Surv object
  expect_equal(
    TimeMetric:::pam.Brier_metric(pred, sv),
    TimeMetric:::pam.Brier_metric(pred, sv, stats::median(d$time)),
    tolerance = 1e-6
  )
})
