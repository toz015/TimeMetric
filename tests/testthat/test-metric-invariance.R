# Numerical invariance gate for the r_sh / r_e removal.
#
# fixtures/metric-baseline.csv holds committed literal values generated from
# the package BEFORE either metric was removed, at commit 675c99f, using the
# deterministic fixtures in helper-simdata.R. Every metric that survives the
# removal must keep its exact value on every public evaluator path.
#
# Nine scenarios cover all four public entry points that resolve metric names
# through tm_normalize_metrics(), because the tombstone added by the removal is
# shared by all of them:
#
#   eval_default      tm_survival_eval()      metrics = NULL
#   eval_all          tm_survival_eval()      metrics = "all"
#   fit_and_eval_all  tm_fit_and_eval()       metrics = "all"
#   cr_default        tm_survival_eval_cr()   metrics = NULL
#   cr_all            tm_survival_eval_cr()   metrics = "all"
#   two_phase_cc      tm_evaluate_two_phase() case-cohort weights
#   two_phase_ncc     tm_evaluate_two_phase() nested case-control weights
#   two_phase_unit    tm_evaluate_two_phase() unit weights (reference path)
#   sample_design_cc  tm_sample_design()      case-cohort
#
# Between them the scenarios pin all 11 canonical post-removal metric names.
# The three two-phase weightings are pinned separately because they genuinely
# diverge -- test-two-phase.R already asserts cc and ncc differ -- so a single
# weighting would not detect a change confined to one of them.
#
# If this test fails, the removal perturbed a metric it should not have
# touched. Investigate before regenerating the baseline: regenerating it to
# make the test pass destroys the only evidence that the values held.

baseline_rows <- function() {
  path <- testthat::test_path("fixtures", "metric-baseline.csv")
  testthat::skip_if(!file.exists(path), "baseline fixture not present")
  utils::read.csv(path, stringsAsFactors = FALSE)
}

# --- scenario reconstruction -------------------------------------------------
# Each builder returns a named numeric vector of metric -> value.

inv_named <- function(metric, value) {
  stats::setNames(round(as.numeric(value), 6), as.character(metric))
}

inv_surv_pred <- function() {
  tm_predict_coxph(model = fx_cox(), covs = fx_covs(),
                   new_data = fx_surv(), tau = 10e10)
}

inv_eval <- function(metrics) {
  pred <- inv_surv_pred()
  res <- tm_survival_eval(
    model = fx_cox(), event_time = pred$times,
    predicted_probability = pred$surv_prob, status = pred$status,
    covariates = fx_covs(), new_data = fx_surv(), tau = 10e10,
    metrics = metrics
  )
  inv_named(res$Metric, res$Value)
}

inv_cr_pred <- function() {
  dd <- fx_cr()
  dd$time <- dd$obs.times
  dd$status <- dd$obs.event
  m1 <- survival::coxph(survival::Surv(time, status == 1) ~ X1 + X2,
                        data = dd, x = TRUE, y = TRUE)
  m2 <- survival::coxph(survival::Surv(time, status == 2) ~ X1 + X2,
                        data = dd, x = TRUE, y = TRUE)
  tm_predict_cif(model1 = m1, model2 = m2, newdata = dd,
                 covs = c("X1", "X2"), event.type = 1, tau = max(dd$time))
}

inv_eval_cr <- function(metrics) {
  p <- inv_cr_pred()
  res <- tm_survival_eval_cr(
    pred_cif = p$cif_pred[, -1], event_time = p$times,
    time.cif = p$cif_pred[, 1], status = p$status, event_type = 1,
    metrics = metrics
  )
  inv_named(res$Metric, res$Value)
}

# Mirrors tp_inputs() in test-two-phase.R. set.seed(2001) is required: both
# weight constructors sample, so the weights are only reproducible under it.
inv_tp_inputs <- function(design) {
  d <- fx_surv()
  pr <- tm_predict_coxph(model = fx_cox(), covs = fx_covs(),
                         new_data = d, tau = 10e10)
  set.seed(2001)
  w <- switch(
    design,
    cc   = tm_case_cohort_weights(time = d$time, status = d$status,
                                  subcohort = stats::rbinom(nrow(d), 1, 0.4)),
    ncc  = tm_nested_case_control_weights(time = d$time, status = d$status, m = 2),
    unit = rep(1, nrow(d))
  )
  list(pred = pr, weights = w,
       km_cens = survival::survfit(survival::Surv(d$time, 1 - d$status) ~ 1))
}

inv_two_phase <- function(design) {
  inp <- inv_tp_inputs(design)
  res <- tm_evaluate_two_phase(pred_results = inp$pred,
                               km_cens_fit = inp$km_cens,
                               case_weights = inp$weights)
  inv_named(res$Metric, res$Value)
}

inv_sample_design <- function() {
  inp <- inv_tp_inputs("cc")
  res <- tm_sample_design(models = list(cc = inp$pred),
                          case_weights = inp$weights, km_cens = inp$km_cens)
  inv_named(res$Metric, res[["cc"]])
}

inv_builders <- list(
  eval_default     = function() inv_eval(NULL),
  eval_all         = function() inv_eval("all"),
  fit_and_eval_all = function() {
    res <- tm_fit_and_eval(train_data = fx_surv(), covariates = fx_covs(),
                           models = "coxph", metrics = "all")
    num <- vapply(res, is.numeric, logical(1))
    inv_named(names(res)[num], unlist(res[num]))
  },
  cr_default       = function() inv_eval_cr(NULL),
  cr_all           = function() inv_eval_cr("all"),
  two_phase_cc     = function() inv_two_phase("cc"),
  two_phase_ncc    = function() inv_two_phase("ncc"),
  two_phase_unit   = function() inv_two_phase("unit"),
  sample_design_cc = inv_sample_design
)

# --- the gate ----------------------------------------------------------------

for (scenario_name in names(inv_builders)) {
  local({
    scn <- scenario_name
    test_that(paste0(scn, ": surviving metric values are unchanged"), {
      baseline <- baseline_rows()
      expected <- baseline[baseline$Scenario == scn, ]
      expect_gt(nrow(expected), 0)

      observed <- inv_builders[[scn]]()

      for (i in seq_len(nrow(expected))) {
        m <- expected$Metric[i]
        expect_true(
          m %in% names(observed),
          info = paste0(scn, ": metric '", m, "' disappeared from the result")
        )
        expect_equal(
          unname(observed[[m]]), expected$Value[i], tolerance = 1e-6,
          info = paste0(scn, ": metric '", m, "' changed value")
        )
      }
    })
  })
}

test_that("the baseline covers every canonical metric and no withdrawn one", {
  baseline <- baseline_rows()

  expect_false("r_e" %in% baseline$Metric)
  expect_false("r_sh" %in% baseline$Metric)

  # exactly the 11 canonical names that survive the removal
  expect_identical(sort(unique(baseline$Metric)), sort(c(
    "brier_score", "c_index", "harrell_c", "l2_point", "l_square",
    "pseudo_r2", "pseudo_r2_point", "r2_point", "r_square", "td_auc", "uno_c"
  )))

  # every scenario this file can rebuild is actually present in the fixture
  expect_setequal(unique(baseline$Scenario), names(inv_builders))
})
