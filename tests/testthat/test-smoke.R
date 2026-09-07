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

test_that("the deleted Cluster C duplicates are gone", {
  # pam.Brier_metric, pam.rsh_metric and pam.rsph_metric were removed after the
  # equivalence gate and the R_E audit. See docs/superpowers/ for the evidence.
  ns <- asNamespace("TimeMetric")

  expect_false(exists("pam.Brier_metric", envir = ns, inherits = FALSE))
  expect_false(exists("pam.rsh_metric", envir = ns, inherits = FALSE))
  expect_false(exists("pam.rsph_metric", envir = ns, inherits = FALSE))
})

test_that("functions this suite reaches with ::: are genuinely internal", {
  exports <- getNamespaceExports("TimeMetric")

  expect_false(any(c(
    "Gt", "pam.Brier", "pam.rsph", "m_cif",
    "pam.predicted_survial_eval_two_phase"
  ) %in% exports))
})
