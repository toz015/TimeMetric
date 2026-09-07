# Regression tests pinning R_E against the authors' own reference implementation.
#
# PROVENANCE OF THE EXPECTED VALUES
#   Source : Re.r, the Stare-Perme-Henderson reference implementation
#            https://ibmi3.mf.uni-lj.si/ibmi-english/biostat-center/programje/Re.r
#   Retrieved: 2026-09-06T21:52:16Z
#   SHA-256  : cf6532068a9a86a61d3763b4cad54aff4bb5000a3ce108f64135debab2b69a92
#   Computed : re.coxph(fit)$Re, where
#              fit <- coxph(Surv(time, status) ~ x1 + x2, data = d, x = TRUE, y = TRUE)
#              d   <- tm_sim_cox_weibull(n = 200, pi_c = 0.3, v = 2,
#                                              beta = c(0.5, -0.5), seed = <seed>)
#                       [, c("time", "status", "x1", "x2")]
#
# These literals were produced by running the REFERENCE, not TimeMetric.
# Regenerating them from TimeMetric would make the test circular and worthless.
# Re.r carries no licence statement and is therefore never vendored into this
# repository; to re-derive these numbers, download it yourself and verify the
# hash above.
#
# Full analysis: docs/superpowers/r-e-implementation-audit.md

REF_RE <- c(
  "1001" = 0.3118680658,
  "2002" = 0.3857787756,
  "3003" = 0.4202986593,
  "4004" = 0.3607495364,
  "5005" = 0.2905941385
)

ref_fixture <- function(seed) {
  d <- tm_sim_cox_weibull(n = 200, pi_c = 0.3, v = 2,
                                beta = c(0.5, -0.5), seed = seed)
  d[, c("time", "status", "x1", "x2")]
}

ref_fit <- function(seed) {
  survival::coxph(survival::Surv(time, status) ~ x1 + x2,
                  data = ref_fixture(seed), x = TRUE, y = TRUE)
}

test_that("pam.rsph reproduces the Stare-Perme-Henderson reference exactly", {
  for (sd in names(REF_RE)) {
    fit <- ref_fit(as.integer(sd))

    expect_equal(
      TimeMetric:::pam.rsph(fit)$Re,
      unname(REF_RE[[sd]]),
      tolerance = 1e-9,
      info = paste("seed", sd)
    )
  }
})

test_that("R_E is invariant to a location shift of the linear predictor", {
  # The reference centres the linear predictor before ranking. Both rank() and
  # ebx/sum(ebx) are shift-invariant, so centring cannot change Re -- this test
  # pins that property rather than assuming it.
  d <- ref_fixture(1001)
  fit <- ref_fit(1001)
  shifted <- fit
  shifted$linear.predictors <- fit$linear.predictors + 10

  expect_equal(
    TimeMetric:::pam.rsph(shifted)$Re,
    TimeMetric:::pam.rsph(fit)$Re,
    tolerance = 1e-9
  )
})

test_that("R_E has exactly one implementation", {
  # pam.rsph_metric was a second, incorrect implementation of the same metric
  # (plain ranks instead of inverse-censoring-weighted ranks, 1.3-5.3% error).
  # It was deleted; this guards against a duplicate reappearing.
  ns <- asNamespace("TimeMetric")
  candidates <- ls(ns, all.names = TRUE)
  re_impls <- grep("^pam\\.rsph", candidates, value = TRUE)

  # the generic plus its three registered methods, and nothing else
  expect_setequal(
    re_impls,
    c("pam.rsph", "pam.rsph.coxph", "pam.rsph.survreg", "pam.rsph.aareg")
  )
})

test_that("the third-party reference implementation is not redistributed", {
  # Re.r carries no licence grant, so it must never be committed to this repo.
  root <- skip_without_source_tree()

  expect_false(file.exists(file.path(root, "Re.r")))
  expect_false(file.exists(file.path(root, "inst", "Re.r")))
  expect_false(file.exists(file.path(root, "R", "Re.r")))
})

test_that("summary.rsph derives R_E over time from the same object", {
  d <- ref_fixture(1001)
  obj <- TimeMetric:::pam.rsph(ref_fit(1001), test_data = d)

  res <- TimeMetric:::summary.rsph(obj, times = stats::median(d$time))

  expect_s3_class(res, "data.frame")
  expect_identical(names(res), c("times", "Rti", "dRti"))
  expect_true(all(is.finite(res$Rti)))
})
