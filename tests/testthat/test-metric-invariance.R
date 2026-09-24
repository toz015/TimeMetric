# Numerical invariance gate for the r_sh / r_e removal.
#
# fixtures/metric-baseline.csv was generated from the package BEFORE either
# metric was removed, at commit 7e32e6f, using the deterministic fixtures in
# helper-simdata.R (tm_sim_cox_weibull(n = 200, pi_c = 0.3, v = 2,
# beta = c(0.5, -0.5), seed = 1001)). Every metric that survives the removal
# must keep its exact value.
#
# Three scenarios are pinned, so the gate covers every surviving metric rather
# than only the six in the default set:
#   eval_default     -- tm_survival_eval() with metrics = NULL
#   eval_all         -- tm_survival_eval() with metrics = "all"
#   fit_and_eval_all -- tm_fit_and_eval() with metrics = "all"
#
# If this test fails, the removal perturbed a metric it should not have
# touched. Investigate before regenerating the baseline -- regenerating it to
# make the test pass would destroy the only evidence that the values held.

baseline_rows <- function() {
  path <- testthat::test_path("fixtures", "metric-baseline.csv")
  testthat::skip_if(!file.exists(path), "baseline fixture not present")
  utils::read.csv(path, stringsAsFactors = FALSE)
}

inv_pred <- function() {
  tm_predict_coxph(model = fx_cox(), covs = fx_covs(),
                   new_data = fx_surv(), tau = 10e10)
}

inv_eval <- function(metrics) {
  pred <- inv_pred()
  tm_survival_eval(
    model = fx_cox(), event_time = pred$times,
    predicted_probability = pred$surv_prob, status = pred$status,
    covariates = fx_covs(), new_data = fx_surv(), tau = 10e10,
    metrics = metrics
  )
}

expect_matches_baseline <- function(observed, scenario) {
  baseline <- baseline_rows()
  expected <- baseline[baseline$Scenario == scenario, ]
  testthat::expect_gt(nrow(expected), 0)

  for (i in seq_len(nrow(expected))) {
    m <- expected$Metric[i]
    testthat::expect_true(
      m %in% names(observed),
      info = paste0(scenario, ": metric '", m, "' disappeared from the result")
    )
    testthat::expect_equal(
      unname(observed[[m]]), expected$Value[i], tolerance = 1e-6,
      info = paste0(scenario, ": metric '", m, "' changed value")
    )
  }
  invisible(NULL)
}

as_named <- function(res) {
  stats::setNames(round(as.numeric(res$Value), 6), res$Metric)
}

test_that("tm_survival_eval default metrics are numerically unchanged", {
  expect_matches_baseline(as_named(inv_eval(NULL)), "eval_default")
})

test_that("tm_survival_eval metrics = 'all' is numerically unchanged", {
  expect_matches_baseline(as_named(inv_eval("all")), "eval_all")
})

test_that("tm_fit_and_eval metrics = 'all' is numerically unchanged", {
  res <- tm_fit_and_eval(train_data = fx_surv(), covariates = fx_covs(),
                         models = "coxph", metrics = "all")
  num <- vapply(res, is.numeric, logical(1))
  observed <- stats::setNames(
    round(as.numeric(unlist(res[num])), 6), names(res)[num]
  )
  expect_matches_baseline(observed, "fit_and_eval_all")
})

test_that("the baseline itself records no withdrawn metric", {
  baseline <- baseline_rows()
  expect_false("r_e" %in% baseline$Metric)
  expect_false("r_sh" %in% baseline$Metric)
  # 10 distinct surviving metrics across the three scenarios
  expect_identical(sort(unique(baseline$Metric)), sort(c(
    "brier_score", "harrell_c", "l2_point", "l_square", "pseudo_r2",
    "pseudo_r2_point", "r2_point", "r_square", "td_auc", "uno_c"
  )))
})
