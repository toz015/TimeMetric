# Two-phase inputs. tm_case_cohort_weights and tm_nested_case_control_weights both return a bare numeric
# vector, one weight per subject -- no list wrapping.
tp_inputs <- function(design = c("cc", "ncc")) {
  design <- match.arg(design)
  d <- fx_surv()
  pred <- tm_predict_coxph(model = fx_cox(), covs = fx_covs(),
                               new_data = d, tau = 10e10)
  set.seed(2001)
  w <- if (design == "cc") {
    tm_case_cohort_weights(time = d$time, status = d$status,
               subcohort = stats::rbinom(nrow(d), 1, 0.4))
  } else {
    tm_nested_case_control_weights(time = d$time, status = d$status, m = 2)
  }
  list(
    d = d,
    pred = pred,
    weights = w,
    km_cens = survival::survfit(survival::Surv(d$time, 1 - d$status) ~ 1)
  )
}

test_that("tm_case_cohort_weights returns one finite weight per subject", {
  d <- fx_surv()
  set.seed(2001)
  subcohort <- stats::rbinom(nrow(d), 1, 0.4)

  w <- tm_case_cohort_weights(time = d$time, status = d$status, subcohort = subcohort)

  expect_true(is.numeric(w))
  expect_false(is.list(w))
  expect_identical(length(w), 200L)
  expect_true(all(is.finite(w)))
  expect_true(all(w > 0))
  # sampling weights inflate non-sampled subjects above 1
  expect_gt(max(w), 1)
  expect_snapshot_value(snap_num(head(w, 10)), style = "serialize")
  expect_snapshot_value(snap_num(sum(w)), style = "serialize")
})

test_that("tm_nested_case_control_weights requires m and returns one weight per subject", {
  d <- fx_surv()

  expect_error(
    tm_nested_case_control_weights(time = d$time, status = d$status),
    "`m` must be provided for NCC weights"
  )

  w <- tm_nested_case_control_weights(time = d$time, status = d$status, m = 2)

  expect_true(is.numeric(w))
  expect_identical(length(w), 200L)
  expect_true(all(is.finite(w)))
  expect_true(all(w > 0))
  expect_snapshot_value(snap_num(head(w, 10)), style = "serialize")
  expect_snapshot_value(snap_num(sum(w)), style = "serialize")
})

test_that("case-cohort and NCC weighting schemes differ", {
  d <- fx_surv()
  set.seed(2001)
  cc <- tm_case_cohort_weights(time = d$time, status = d$status,
                   subcohort = stats::rbinom(nrow(d), 1, 0.4))
  ncc <- tm_nested_case_control_weights(time = d$time, status = d$status, m = 2)

  expect_false(isTRUE(all.equal(cc, ncc, tolerance = 1e-6)))
})

test_that("two-phase evaluation works for a case-cohort design", {
  inp <- tp_inputs("cc")

  res <- tm_evaluate_two_phase(
    pred_results = inp$pred,
    km_cens_fit = inp$km_cens,
    case_weights = inp$weights
  )

  expect_true(inherits(res, "data.frame"))
  expect_true(all(c("Metric", "Value") %in% names(res)))
  expect_identical(nrow(res), 5L)
  expect_true(all(is.finite(res$Value)))
  expect_snapshot_value(res$Metric, style = "serialize")
  expect_snapshot_value(snap_num(res$Value), style = "serialize")
})

test_that("two-phase evaluation works for a nested case-control design", {
  inp <- tp_inputs("ncc")

  res <- tm_evaluate_two_phase(
    pred_results = inp$pred,
    km_cens_fit = inp$km_cens,
    case_weights = inp$weights
  )

  expect_true(all(c("Metric", "Value") %in% names(res)))
  expect_identical(nrow(res), 5L)
  expect_snapshot_value(res$Metric, style = "serialize")
  expect_snapshot_value(snap_num(res$Value), style = "serialize")
})

test_that("case-cohort and NCC weighting give different metric values", {
  # If these coincided the weights would not be reaching the estimator, which
  # would make the whole two-phase feature a no-op.
  cc_in <- tp_inputs("cc")
  ncc_in <- tp_inputs("ncc")

  cc <- tm_evaluate_two_phase(
    pred_results = cc_in$pred, km_cens_fit = cc_in$km_cens,
    case_weights = cc_in$weights
  )
  ncc <- tm_evaluate_two_phase(
    pred_results = ncc_in$pred, km_cens_fit = ncc_in$km_cens,
    case_weights = ncc_in$weights
  )

  expect_identical(cc$Metric, ncc$Metric)
  expect_false(isTRUE(all.equal(cc$Value, ncc$Value, tolerance = 1e-6)))
})

test_that("two-phase output labels AUC differently from tm_sample_design", {
  # FINDING 15: the evaluator emits "td_auc" while its own summary
  # wrapper emits "td_auc" for the same quantity.
  inp <- tp_inputs("cc")

  direct <- tm_evaluate_two_phase(
    pred_results = inp$pred, km_cens_fit = inp$km_cens,
    case_weights = inp$weights
  )
  summarised <- tm_sample_design(
    models = list(cc = inp$pred), case_weights = inp$weights,
    km_cens = inp$km_cens
  )

  expect_true("td_auc" %in% direct$Metric)
  expect_true("td_auc" %in% summarised$Metric)
  # the two entry points now use identical labels, so results can be joined
  expect_setequal(direct$Metric, summarised$Metric)
})

test_that("tm_sample_design validates its models argument", {
  inp <- tp_inputs("cc")

  expect_error(
    tm_sample_design(models = list(), case_weights = inp$weights,
                      km_cens = inp$km_cens),
    "must be a non-empty named list"
  )
})

test_that("tm_sample_design summarises a two-phase design", {
  inp <- tp_inputs("cc")

  res <- tm_sample_design(
    models = list(cc = inp$pred),
    case_weights = inp$weights,
    km_cens = inp$km_cens
  )

  expect_metric_table(res)
  expect_identical(names(res), c("Metric", "cc"))
  expect_equal(res$cc, round(res$cc, 2), tolerance = 1e-12)
  expect_snapshot_value(res$Metric, style = "serialize")
  expect_snapshot_value(snap_num(res$cc), style = "serialize")
})

test_that("two-phase default metric names are pinned with current spelling", {
  # c("pseudo_r2", "Harrell<u2019>s C", "Uno<u2019>s C", "brier_score",
  #   "td_auc") -- misspelling and curly apostrophes included.
  defaults <- eval(formals(
    tm_evaluate_two_phase
  )$metrics)

  expect_true("pseudo_r2" %in% defaults)
  expect_true("harrell_c" %in% defaults)
  expect_true("uno_c" %in% defaults)
  expect_false("Pesudo_R" %in% defaults)
  expect_true(all(defaults %in% tm_metric_names()))
  expect_snapshot_value(defaults, style = "serialize")
})

test_that("two-phase uses Pesudo_R where the survival path uses Pseudo_R_square", {
  # The same underlying measure carries three different spellings across the
  # package. Recorded so spec step 6's standardization has a pinned starting set.
  inp <- tp_inputs("cc")

  res <- tm_evaluate_two_phase(
    pred_results = inp$pred, km_cens_fit = inp$km_cens,
    case_weights = inp$weights
  )

  expect_true("pseudo_r2" %in% res$Metric)
  expect_false("Pesudo_R" %in% res$Metric)
})
