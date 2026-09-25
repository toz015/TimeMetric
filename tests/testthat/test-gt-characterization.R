# Characterization of Gt(), the Kaplan-Meier estimate of the censoring
# distribution, pinned BEFORE the survminer::surv_summary() call was replaced.
#
# Every expected value below is a committed literal, printed at 17 significant
# digits from the survminer-based implementation at commit edc606a. The purpose
# is to prove the replacement is behaviour-preserving, so these values must not
# be regenerated to accommodate a difference. If one changes, stop.
#
# Two pre-existing quirks are pinned deliberately, not fixed:
#
#  1. The interpolation branch at R/pam.Ct.R indexes `time`, the RAW input
#     vector, using indices derived from the SORTED summary table. The two are
#     not aligned, so interpolated values can be non-monotonic in the timepoint
#     -- visible in the tied-times block below, where t = 2.5 yields a smaller
#     value than t = 3. This is inherited behaviour and is out of scope for the
#     survminer removal; it is recorded here so that any future fix is a
#     deliberate, visible change rather than an accident.
#
#  2. `na.omit()` runs over the whole summary frame, so the row dropped depends
#     on `upper`/`lower` -- columns Gt() never reads. Any replacement must
#     therefore carry those columns, or it will drop different rows.

gt_ <- function(...) TimeMetric:::Gt(...)

test_that("Gt is pinned on ordinary right-censored data", {
  d  <- fx_surv()
  sv <- survival::Surv(d$time, d$status)
  qs <- unname(stats::quantile(d$time, c(0.1, 0.25, 0.5, 0.75, 0.9)))

  expect_equal(gt_(sv, qs[1]), 0.99489795918367352, tolerance = 1e-15)
  expect_equal(gt_(sv, qs[2]), 0.95172567105939188, tolerance = 1e-15)
  expect_equal(gt_(sv, qs[3]), 0.87289681637206018, tolerance = 1e-15)
  expect_equal(gt_(sv, qs[4]), 0.69578180870964024, tolerance = 1e-15)
  expect_equal(gt_(sv, qs[5]), 0.45373758529021757, tolerance = 1e-15)
})

test_that("Gt is pinned at timepoints falling exactly on an observed time", {
  # exercises the `timepoint %in% res.sum$time` branch, bypassing interpolation
  d   <- fx_surv()
  sv  <- survival::Surv(d$time, d$status)
  st0 <- ifelse(d$status == min(d$status), 1, 0)
  f   <- survival::survfit(survival::Surv(d$time, st0) ~ 1)
  on  <- f$time[c(5, 50, 150)]

  expect_equal(gt_(sv, on[1]), 0.99489795918367352, tolerance = 1e-15)
  expect_equal(gt_(sv, on[2]), 0.94823467122223315, tolerance = 1e-15)
  expect_equal(gt_(sv, on[3]), 0.69578180870964057, tolerance = 1e-15)
})

test_that("Gt is pinned before the first and after the last observed time", {
  d  <- fx_surv()
  sv <- survival::Surv(d$time, d$status)

  expect_equal(gt_(sv, min(d$time) / 2), 1, tolerance = 1e-15)
  expect_equal(gt_(sv, max(d$time) * 2), 0.23852558795853213, tolerance = 1e-15)
})

test_that("Gt is pinned under heavily tied event times", {
  # 40 observations on 5 distinct times, alternating status.
  # Non-monotonicity here is the inherited indexing quirk described above.
  tt <- rep(c(1, 2, 3, 4, 5), each = 8)
  ss <- rep(c(1, 0), 20)
  sv <- survival::Surv(tt, ss)

  expect_equal(gt_(sv, 1),   0.90000000000000002, tolerance = 1e-15)
  expect_equal(gt_(sv, 2.5), 0.24609375,          tolerance = 1e-15)
  expect_equal(gt_(sv, 3),   0.65625,             tolerance = 1e-15)
  expect_equal(gt_(sv, 4.5), 0.24609375,          tolerance = 1e-15)
  expect_equal(gt_(sv, 5),   0.24609375,          tolerance = 1e-15)
})

test_that("Gt is pinned on the uncensored fixture", {
  du  <- fx_surv_uncensored()
  svu <- survival::Surv(du$time, du$status)
  qs  <- unname(stats::quantile(du$time, c(0.25, 0.5, 0.75)))

  expect_equal(gt_(svu, qs[1]), 0.74992019118738285, tolerance = 1e-15)
  expect_equal(gt_(svu, qs[2]), 0.61748817289902214, tolerance = 1e-15)
  expect_equal(gt_(svu, qs[3]), 0.33081244121434905, tolerance = 1e-15)
})

test_that("the summary frame drops exactly one row, on missing confidence limits", {
  # This is the behaviour a replacement must reproduce: na.omit() sees columns
  # Gt() never reads, and the row it removes depends on them.
  d   <- fx_surv()
  st0 <- ifelse(d$status == min(d$status), 1, 0)
  f   <- survival::survfit(survival::Surv(d$time, st0) ~ 1)
  # implementation-agnostic on purpose: this invariant must hold for the
  # survminer version and its replacement alike, so the test is identical
  # before and after the swap rather than rewritten across it.
  sm <- if (exists("gt_surv_summary", envir = asNamespace("TimeMetric"),
                   inherits = FALSE)) {
    TimeMetric:::gt_surv_summary(f)
  } else {
    survminer::surv_summary(f)
  }

  expect_identical(
    names(sm),
    c("time", "n.risk", "n.event", "n.censor", "surv", "std.err", "upper", "lower")
  )
  expect_identical(nrow(sm), 200L)
  expect_identical(nrow(stats::na.omit(sm)), 199L)
  expect_true(all(c("upper", "lower") %in% names(sm)[colSums(is.na(sm)) > 0]))
})

test_that("Gt rejects malformed input", {
  d  <- fx_surv()
  sv <- survival::Surv(d$time, d$status)

  expect_error(gt_(d$time, 1), "not of class Surv")
  expect_error(gt_(sv, -1), "must be positive")
  expect_error(gt_(sv, c(1, 2)), "single time point")
})
