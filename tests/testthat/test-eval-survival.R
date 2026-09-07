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
  expect_true("Harrell\u2019s C" %in% res$Metric)
  expect_true("Uno\u2019s C" %in% res$Metric)
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

test_that("concordancefit is imported (FINDING 11 fixed)", {
  # R/pam.predicted_survial_eval.R:162,169 call concordancefit() unqualified.
  # It is now imported explicitly, so the call resolves from the namespace
  # rather than depending on survival happening to be attached.
  imported <- unlist(getNamespaceImports("TimeMetric"), use.names = FALSE)

  expect_true("concordancefit" %in% imported)
})

test_that("a clean session without survival attached now succeeds (FINDING 11)", {
  # The definitive check: a fresh R subprocess that loads TimeMetric and nothing
  # else -- the state a first-time user is in. This previously died with
  # "could not find function concordancefit"; it must now produce metrics.
  skip_on_cran()
  skip_if_covr()
  pkg_root <- skip_without_source_tree()

  script <- sprintf('
    suppressWarnings(pkgload::load_all(%s, quiet = TRUE, attach_testthat = FALSE))
    stopifnot(!"package:survival" %%in%% search())
    d <- sim_cox_weibull_censored(n = 50, pi_c = 0.3, v = 2,
                                  beta = c(0.5, -0.5), seed = 1001)
    d <- d[, c("time", "status", "x1", "x2")]
    m <- survival::coxph(survival::Surv(time, status) ~ x1 + x2,
                         data = d, x = TRUE, y = TRUE)
    p <- pam.coxph_restricted(model = m, covs = c("x1", "x2"),
                              new_data = d, tau = 10e10)
    r <- try(pam.predicted_survial_eval(
      model = m, event_time = p$times, predicted_probability = p$surv_prob,
      status = p$status, covariates = c("x1", "x2"), new_data = d, tau = 10e10
    ), silent = TRUE)
    cat(if (inherits(r, "try-error")) as.character(r) else "NO ERROR")
  ', shQuote(pkg_root))

  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    args = c("--vanilla", "-e", shQuote(script)),
    stdout = TRUE, stderr = TRUE
  ))
  out <- paste(out, collapse = "\n")

  expect_no_match(out, "could not find function")
  expect_no_match(out, "Error")
  expect_match(out, "NO ERROR")   # sentinel printed only on success
})

test_that("R_E comes from the canonical pam.rsph path", {
  # Findings 7 and 13 recorded a gap between R_E and pam.rsph_metric$r2 and
  # framed them as distinct quantities. The R_E audit superseded that: they were
  # one metric implemented twice, and pam.rsph_metric was the wrong one. It has
  # been deleted. R_E now has a single implementation, verified against the
  # authors reference in test-r-e-reference.R.
  d <- fx_surv()
  res <- ev_result()
  r_e <- res$Value[res$Metric == "R_E"]

  expect_length(r_e, 1L)
  expect_true(is.finite(r_e))
  # identical to calling the canonical path directly. The reported R_E is
  # pam.summary.rsph(..., times = tau)$Rti, not pam.rsph(...)$Re -- see
  # R/pam.predicted_survial_eval.R:225-226.
  direct <- TimeMetric:::pam.summary.rsph(
    TimeMetric:::pam.rsph(fx_cox(), test_data = d),
    times = 10e10
  )$Rti
  expect_equal(r_e, direct, tolerance = 1e-8)
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

test_that("pam.survival_eval fits and evaluates from raw data (FINDING 12 fixed)", {
  # R/pam.survial_eval.R:106,112 previously passed covariates=/newdata= to
  # functions taking covs=/new_data=, so every call failed. They now also pass
  # predict = FALSE, which is the branch returning R.squared / L.squared that
  # the caller indexes as r_l_list[1] and [2].
  d <- fx_surv()

  res <- pam.survival_eval(train_data = d, covariates = fx_covs(),
                           models = "coxph", metrics = "all")

  expect_s3_class(res, "data.frame")
  expect_identical(nrow(res), 1L)
  expect_true("Model" %in% names(res))
  expect_identical(res$Model, "coxph")
  # wide layout: one column per metric, unlike the long Metric/Value frame
  # returned by pam.predicted_survial_eval
  expect_true(all(c("Pseudo_R_square", "R_square", "L_square",
                    "Brier Score") %in% names(res)))
  # R_sph and R_sh appear as separate columns, consistent with FINDING 13
  expect_true(all(c("R_sph", "R_sh") %in% names(res)))

  expect_snapshot_value(sort(names(res)), style = "serialize")
  expect_snapshot_value(snap_num(res[["R_square"]]), style = "serialize")
})

test_that("pam.survival_eval honours an explicit model and metric subset", {
  d <- fx_surv()

  res <- pam.survival_eval(
    train_data = d, covariates = fx_covs(),
    models  = c("weibull", "lognormal"),
    metrics = c("R_square", "L_square", "Brier Score")
  )

  expect_s3_class(res, "data.frame")
  expect_identical(nrow(res), 2L)
  expect_setequal(res$Model, c("weibull", "lognormal"))
  expect_false("Harrell\u2019s C" %in% names(res))
})

test_that("pam.survival_eval requires train_data and covariates", {
  expect_error(
    pam.survival_eval(),
    "Please provide 'train_data', 'time_var', 'status_var', and 'covariates'"
  )
})
