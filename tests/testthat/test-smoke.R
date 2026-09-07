test_that("package exports the tm_ API", {
  exports <- getNamespaceExports("TimeMetric")

  expect_true(all(c(
    "tm_survival_eval", "tm_survival_eval_cr", "tm_evaluate_two_phase",
    "tm_fit_and_eval", "tm_predict_coxph", "tm_predict_survreg",
    "tm_predict_cif", "tm_summarize", "tm_summarize_cr", "tm_sample_design",
    "tm_case_cohort_weights", "tm_nested_case_control_weights",
    "tm_plot_pred", "tm_plot_summary", "tm_sim_cox_weibull",
    "tm_simulate_fine_gray"
  ) %in% exports))
  # every public function carries the prefix, apart from the deprecated aliases
  new_api <- grep("^tm_", exports, value = TRUE)
  expect_length(new_api, 17L)   # 16 renamed functions + tm_metric_names()
  expect_true("tm_metric_names" %in% exports)
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
    "Gt", "pam.Brier", "pam.rsph", "m_cif", "my.survfit"
  ) %in% exports))
  # tm_evaluate_two_phase was promoted from internal to public, so the
  # case-cohort and NCC functionality the paper advertises is now reachable
  expect_true("tm_evaluate_two_phase" %in% exports)
})
