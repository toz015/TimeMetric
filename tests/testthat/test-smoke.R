test_that("package exports the functions this suite characterizes", {
  exports <- getNamespaceExports("TimeMetric")

  expect_type(exports, "character")
  expect_true(all(c(
    "pam.survival_eval", "pam.predicted_survial_eval",
    "pam.predicted_survial_eval_cr", "pam.coxph_restricted",
    "pam.surverg_restricted", "pam.predict_cr", "pam.summary",
    "pam.summary_cr", "pam.sample_design", "cc_weights", "ncc_weights",
    "plot_pred", "summary_pred_plot", "sim_cox_weibull_censored",
    "simulateTwoCauseFineGrayModel"
  ) %in% exports))
})

test_that("functions this suite reaches with ::: are genuinely internal", {
  exports <- getNamespaceExports("TimeMetric")

  expect_false(any(c(
    "Gt", "pam.Brier", "pam.rsh_metric", "pam.rsph_metric",
    "pam.Brier_metric", "pam.rsph", "m_cif",
    "pam.predicted_survial_eval_two_phase"
  ) %in% exports))
})
