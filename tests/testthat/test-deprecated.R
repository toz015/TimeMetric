# Every pre-rename export must still work, warn, and return the same value as
# its replacement. These wrappers exist so paper.code.Rmd and the Zhuang et al.
# (2025) reproduction code keep running unchanged.

DEPRECATED_MAP <- c(
  "pam.predicted_survial_eval"    = "tm_survival_eval",
  "pam.predicted_survial_eval_cr" = "tm_survival_eval_cr",
  "pam.survival_eval"             = "tm_fit_and_eval",
  "pam.coxph_restricted"          = "tm_predict_coxph",
  "pam.surverg_restricted"        = "tm_predict_survreg",
  "pam.predict_cr"                = "tm_predict_cif",
  "pam.summary"                   = "tm_summarize",
  "pam.summary_cr"                = "tm_summarize_cr",
  "pam.sample_design"             = "tm_sample_design",
  "cc_weights"                    = "tm_case_cohort_weights",
  "ncc_weights"                   = "tm_nested_case_control_weights",
  "plot_pred"                     = "tm_plot_pred",
  "summary_pred_plot"             = "tm_plot_summary",
  "sim_cox_weibull_censored"      = "tm_sim_cox_weibull",
  "simulateTwoCauseFineGrayModel" = "tm_simulate_fine_gray"
)

test_that("every previously-exported name is still exported", {
  exports <- getNamespaceExports("TimeMetric")

  expect_true(all(names(DEPRECATED_MAP) %in% exports))
  expect_true(all(unname(DEPRECATED_MAP) %in% exports))
})

test_that("each deprecated name warns and names its replacement", {
  for (old in names(DEPRECATED_MAP)) {
    new <- DEPRECATED_MAP[[old]]
    fn <- get(old, envir = asNamespace("TimeMetric"))

    # calling with no arguments is enough to reach .Deprecated
    cond <- tryCatch(
      { fn(); NULL },
      deprecated = function(c) c,
      error = function(e) e,
      warning = function(w) w
    )
    expect_false(is.null(cond), info = old)
    expect_match(conditionMessage(cond), new, fixed = TRUE, info = old)
  }
})

test_that("deprecated wrappers forward arguments and return identical results", {
  d <- fx_surv()

  old_res <- suppressWarnings(
    pam.coxph_restricted(model = fx_cox(), covs = fx_covs(),
                         new_data = d, tau = 10e10)
  )
  new_res <- tm_predict_coxph(model = fx_cox(), covs = fx_covs(),
                              new_data = d, tau = 10e10)
  expect_equal(old_res, new_res, tolerance = 1e-12)

  old_w <- suppressWarnings(
    ncc_weights(time = d$time, status = d$status, m = 2)
  )
  new_w <- tm_nested_case_control_weights(time = d$time, status = d$status, m = 2)
  expect_equal(old_w, new_w, tolerance = 1e-12)
})

test_that("the misspelled names are gone from the new API", {
  # survial -> survival, surverg -> survreg
  exports <- getNamespaceExports("TimeMetric")
  new_api <- grep("^tm_", exports, value = TRUE)

  expect_false(any(grepl("survial", new_api)))
  expect_false(any(grepl("surverg", new_api)))
  expect_false(any(grepl("^pam\\.", new_api)))
})
