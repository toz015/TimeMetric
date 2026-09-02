# Shared prediction input for the evaluation tests.
ev_pred <- function() {
  pam.coxph_restricted(model = fx_cox(), covs = fx_covs(),
                       new_data = fx_surv(), tau = 10e10)
}

ev_result <- function() {
  pred <- ev_pred()
  pam.predicted_survial_eval(
    model = fx_cox(),
    event_time = pred$times,
    predicted_probability = pred$surv_prob,
    status = pred$status,
    covariates = fx_covs(),
    new_data = fx_surv(),
    tau = 10e10
  )
}

test_that("pam.predicted_survial_eval returns a Metric/Value data frame", {
  res <- ev_result()

  expect_metric_table(res)
  expect_identical(names(res), c("Metric", "Value"))
  expect_type(res$Value, "double")
  expect_snapshot_value(res$Metric, style = "serialize")
  expect_snapshot_value(snap_num(res$Value), style = "serialize")
})

test_that("the default metric set includes R_sh and R_E in the Metric column", {
  # Metric names live in res$Metric; names(res) is c("Metric", "Value").
  # R_sh reaches rms::cph + pam.schemper; R_E reaches pam.rsph dispatch.
  # Asserting both here proves those two code paths execute by default.
  res <- ev_result()

  expect_true("R_sh" %in% res$Metric)
  expect_true("R_E" %in% res$Metric)
  expect_true("Brier Score" %in% res$Metric)
  expect_true("Pseudo_R_square" %in% res$Metric)
  expect_true("Time Dependent Auc" %in% res$Metric)
  # curly apostrophes, written as unicode escapes so this file stays ASCII
  expect_true("Harrell’s C" %in% res$Metric)
  expect_true("Uno’s C" %in% res$Metric)
  # the ASCII spellings are NOT what the package uses
  expect_false("Harrells_C" %in% res$Metric)
  expect_false("Unos_C" %in% res$Metric)
})

test_that("pam.predicted_survial_eval rejects an unknown metric name", {
  pred <- ev_pred()

  expect_error(
    pam.predicted_survial_eval(
      model = fx_cox(),
      event_time = pred$times,
      predicted_probability = pred$surv_prob,
      status = pred$status,
      covariates = fx_covs(),
      new_data = fx_surv(),
      tau = 10e10,
      metrics = "Harrells_C"        # ASCII spelling is NOT in valid_metrics
    ),
    "Invalid metrics"
  )
})

test_that("concordancefit is called but never imported (FINDING 11)", {
  # R/pam.predicted_survial_eval.R:162,169 call concordancefit() unqualified.
  # survival exports it, but TimeMetric's NAMESPACE does not import it, so the
  # call only resolves when survival happens to be attached to the search path.
  # Loading TimeMetric alone and calling this function fails with
  # "could not find function \"concordancefit\"".
  imported <- unlist(getNamespaceImports("TimeMetric"), use.names = FALSE)

  expect_false("concordancefit" %in% imported)
  # while it genuinely is available to import
  expect_true("concordancefit" %in% getNamespaceExports("survival"))
})

test_that("R_E and pam.rsph_metric's r2 are NOT the same number (FINDING 7)", {
  # Bears directly on the spec's open question of whether R_sph equals R_E.
  # If these were the same quantity the two names could be merged; they are not.
  d <- fx_surv()
  risk <- as.numeric(predict(fx_cox(), newdata = d, type = "lp"))

  r_e <- ev_result()
  r_e_val <- r_e$Value[r_e$Metric == "R_E"]
  rsph_r2 <- TimeMetric:::pam.rsph_metric(d$time, d$status, risk)$r2

  expect_length(r_e_val, 1L)
  expect_false(isTRUE(all.equal(r_e_val, rsph_r2, tolerance = 1e-3)))
  expect_snapshot_value(snap_num(c(r_e_val, rsph_r2)), style = "serialize")
})

test_that("pam.summary pivots one model into a Metric column table", {
  res <- pam.summary(list(value = ev_pred()), tau = 10e10)

  expect_metric_table(res)
  expect_identical(names(res), c("Metric", "value"))
  expect_snapshot_value(res$Metric, style = "serialize")
  expect_snapshot_value(snap_num(res$value), style = "serialize")
})

test_that("pam.summary puts one column per model and rounds to digits", {
  d <- fx_surv()
  p_cox <- pam.coxph_restricted(model = fx_cox(), covs = fx_covs(),
                                new_data = d, tau = 10e10)
  p_reg <- pam.surverg_restricted(model = fx_survreg(), covs = fx_covs(),
                                  new_data = d, tau = 10e10)

  res <- pam.summary(list(cox = p_cox, weibull = p_reg), tau = 10e10)

  expect_metric_table(res)
  expect_identical(names(res), c("Metric", "cox", "weibull"))
  # default digits = 2
  expect_equal(res$cox, round(res$cox, 2), tolerance = 1e-12)
  # R_sh is Cox-only, so the weibull column carries NA for it
  expect_true(is.na(res$weibull[res$Metric == "R_sh"]))
  expect_false(is.na(res$cox[res$Metric == "R_sh"]))

  expect_snapshot_value(snap_num(res$cox), style = "serialize")
  expect_snapshot_value(snap_num(res$weibull), style = "serialize")
})

test_that("pam.summary rejects a non-list or empty models argument", {
  expect_error(pam.summary(list()), "must be a non-empty named list")
  expect_error(pam.summary("not a list"), "must be a non-empty named list")
})

test_that("pam.survival_eval is unconditionally broken (FINDING 12)", {
  # R/pam.survial_eval.R:110,113 call pam.coxph_restricted / pam.surverg_restricted
  # with covariates = and newdata =, but those functions take covs = and
  # new_data =. Every invocation of this exported function fails, whatever the
  # input. Pinned so the eventual fix is visible as a test change.
  d <- fx_surv()

  expect_error(
    pam.survival_eval(train_data = d, covariates = fx_covs(),
                      models = "coxph", metrics = "all"),
    "unused arguments"
  )
})

test_that("pam.survival_eval requires train_data and covariates", {
  expect_error(
    pam.survival_eval(),
    "Please provide 'train_data', 'time_var', 'status_var', and 'covariates'"
  )
})
