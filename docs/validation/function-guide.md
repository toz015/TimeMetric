# TimeMetric function and usage guide

Updated 18 October 2026 for TimeMetric 0.2.1 with the paper branch integration.
The package name is **TimeMetric**, singular. There are **17 current exported
functions**, **15 deprecated compatibility wrappers**, and **11 metric names**.
Internal `pam.*` helpers are implementation details, not supported public APIs.
Read [the validation report](cran-resubmission-2026-10-18.md) before interpreting manuscript results.

## Installation and development

```r
# Install the integrated paper branch.
remotes::install_github("toz015/TimeMetric@paper")
library(TimeMetric)

# To use the revised working tree instead, run from its root:
pkgload::load_all(".")
testthat::test_local()
```

R >= 4.1.0 is declared. `cmprsk` and `randomForestSRC` are optional prediction
backends. The installed `?function_name` help pages contain argument details.

## All current functions

| Function | Inputs and purpose | Return / important defaults |
|---|---|---|
| `tm_survival_eval()` | `event_time`, binary `status`, and a subject-by-time `predicted_probability` matrix; optionally `pred_mean_survival`, `model`, `covariates`, `new_data`. Evaluate right-censored predictions. | `Metric` / numeric `Value` table. Six default metrics; `metrics="all"` gives ten. `t_star=NULL` uses median observed time; `tau=NULL` uses available follow-up for integration. |
| `tm_survival_eval_cr()` | `pred_cif`, `event_time`, `time.cif`, multi-cause `status`; `event_type=1`. Evaluate competing-risk predictions. | Five default metrics; nine with `"all"`. `t_star=NULL` means median observed time; `tau=NULL` means maximum observed time. |
| `tm_evaluate_two_phase()` | Prediction list `pred_results`, aligned `case_weights`, and censoring fit `km_cens_fit`. Evaluate case-cohort or NCC samples. | Five default metrics; seven with `"all"`. `tau=10e10` (that is, 100 billion), `t_star=NULL`. Set an explicit scientifically meaningful horizon. |
| `tm_fit_and_eval()` | `train_data`, `covariates`; optionally independent `predicted_data`, `models`, `metrics`, `tau`, `t_star`. Fit and assess models in one call. | One row per model, `Model` plus metric columns. Default `models="coxph"` fits **only Cox**; `models="all"` fits Cox, exponential, lognormal, Weibull. Seven metrics. Uses different metric estimators from the prediction/evaluation route. |
| `tm_predict_coxph()` | A fitted `coxph` `model`, `covs`, optional `new_data`, `tau`, `predict`. Fit with `x=TRUE, y=TRUE`. | With default `predict=TRUE`: list containing `pred`, `times`, `status`, `surv_prob`, `linear.pred`, `model`, `new_data`, `covs`. Subjects and time-grid columns are sorted by observed time. Probabilities now come from the training fit. |
| `tm_predict_survreg()` | A fitted `survreg` `model`, `covs`, optional `new_data`, `tau`, `predict`. Fit with `x=TRUE, y=TRUE`. | Same prediction structure as Cox. Weibull/exponential and lognormal mean-survival paths are implemented. `pred` is the restricted mean at `tau`; original outcomes are retained. Default `tau=10e10`. |
| `tm_predict_cif()` | `newdata`, `covs`; either cause-specific `model1` / `model2`, `fg_model`, or `cr_model`. `event.type=1`, `tau=NULL`. | List with `times`, `status`, `cif_pred`, `pred`, `linear.pred`. `cif_pred` has time in its first column and one subject per remaining column. Cox/AFT cause-specific pairs, Fine–Gray (`cmprsk`), or competing-risk forests (`randomForestSRC`). |
| `tm_summarize()` | Named list of right-censored prediction objects; `metrics`, `t_star`, `tau`, `digits`. | One metric per row and one column per model. Defaults to six metrics, `digits=2`, `tau=10e10`. |
| `tm_summarize_cr()` | Named list of CIF prediction objects; `event_type`, `metrics`, `t_star`, `tau`, `digits`. | Wide table; five default metrics, `event_type=1`, `digits=2`, `tau=NULL`. |
| `tm_sample_design()` | Named prediction list, aligned `case_weights`, `km_cens`; optional metrics/horizons/digits. | Wide two-phase table; five default metrics, `digits=2`, `tau=NULL` (maximum follow-up). `pred` is required. |
| `tm_case_cohort_weights()` | Full-cohort `time`, binary `status`, logical `subcohort`; optional `strata`. | Numeric weights in input order. All cases receive 1; non-cases receive the ratio of full-cohort to sampled non-case counts, within stratum if supplied. Subset to actual sampled subjects **after** computing weights. |
| `tm_nested_case_control_weights()` | Full-cohort `time`, binary `status`; `m` controls per case, optional matching `strata`. | Numeric inverse-selection weights in input order. Does not select controls or return sample membership. Supply the actual sampled subjects separately. Zero-inclusion non-cases have infinite weights and cannot be included. See its help page for tied-event sampling assumptions. |
| `tm_plot_pred()` | Prediction object with `linear.pred`, `pred`, `times`, `status`; title/labels/shapes, `invert_linear`, `sample_index`, `restrict_time`. | `ggplot` of observed and predicted time versus score. Default `invert_linear=TRUE`. This is a diagnostic plot, not a calibration curve or a survival-curve plot. |
| `tm_plot_summary()` | List of prediction objects, optional `titles`, `plot_fun`, layout and plot settings. | Combined patchwork plot. Defaults: `ncol=2`, `seed=NULL`. `invert_linear` may differ by model; `sample_size` limits plotted subjects. |
| `tm_sim_cox_weibull()` | Required `n`, `pi_c`, Weibull shape `v`, and coefficient vector `beta`; optional `mu`, `sd`, `seed`, `interact`, `nonlinear`. | Data frame with `time`, `status`, covariates, and simulation metadata. Use explicit predictor columns: latent `y_true` and `cens_time` must never enter model formulas. |
| `tm_simulate_fine_gray()` | Required `n`, `v`, `beta1`, `beta2`; optional baseline/censoring parameters, covariate matrix `X`, and `seed=NULL`. | Simulated competing-risk data. Map `obs.times` to `time` and `obs.event` to `status` for prediction helpers. Censoring is controlled by `censor`, not `pi_c`. |
| `tm_metric_names()` | No arguments. | Sorted character vector of all 11 canonical metric identifiers. It is the union across settings, not a promise that every evaluator accepts every identifier. |

For Cox/AFT prediction, `predict=FALSE` invokes a legacy fitted-measure path
returning `R.squared`, `L.squared`, and `Psuedo.R` as formatted strings. Prefer
`predict=TRUE` followed by an explicit evaluator for new work.

## Which metrics belong to each setting?

| Metric | Right-censored | Competing risks | Two-phase | Fit-and-evaluate |
|---|---|---|---|---|
| `pseudo_r2` | default | default | default | default |
| `pseudo_r2_point` | default | default | unavailable | unavailable |
| `r_square`, `l_square` | optional | optional | optional | default |
| `r2_point`, `l2_point` | optional | optional | unavailable | unavailable |
| `harrell_c`, `uno_c` | default | unavailable | default | default |
| `c_index` | unavailable | default | unavailable | unavailable |
| `brier_score`, `td_auc` | default | default | default | default |

Pseudo R-squared combines explained-variation and calibration-related
components. C-indices and AUC measure discrimination; larger values normally
indicate better discrimination. Brier score measures overall probability
prediction error, incorporating calibration and discrimination; smaller is
better. It is not a calibration-only measure.

`r_e`, `r_sh`, and legacy spellings such as `R_E` are **withdrawn**, not
deprecated aliases. They raise an error. Other legacy metric names are
recognized with warnings, but all new examples should use canonical names.

## Prediction contracts and alignment

Right-censored `surv_prob` is **subjects × sorted observed follow-up times**,
including duplicate time columns. Each row follows the order of `times` and
`status`. A probability matrix on an arbitrary unrelated grid cannot be
passed without adaptation: there is no separate grid argument on this path.

CIF orientation is different: `cif_pred` is **time points × (1 + subjects)**.
For the lower-level evaluator, use `pred_cif = p$cif_pred[, -1, drop=FALSE]`
and `time.cif = p$cif_pred[, 1]`. `tm_predict_cif()` uses `event.type`, whereas
the evaluator and summary use `event_type`. `newdata` also differs from the
Cox/AFT helpers' `new_data`.

Cox/AFT helpers sort prediction subjects by follow-up time. A weight vector
returned for the unsorted input data must be reordered by the same permutation.
Use a common set of complete cases for model comparisons. Do not let a model
silently drop subjects while retaining the original weight vector.

## Minimal right-censored example

```r
library(TimeMetric)
d <- tm_sim_cox_weibull(240, pi_c = 0.3, v = 2,
                        beta = c(0.5, -0.5), seed = 6201)
d <- d[, c("time", "status", "x1", "x2")]
train <- d[1:160, ]
test <- d[161:240, ]
fit <- survival::coxph(survival::Surv(time, status) ~ x1 + x2,
                       data = train, x = TRUE, y = TRUE)
tau <- unname(quantile(train$time, 0.6))
p <- tm_predict_coxph(fit, covs = c("x1", "x2"),
                      new_data = test, tau = tau)
tm_summarize(list(Cox = p), tau = tau, t_star = tau,
             metrics = c("pseudo_r2", "harrell_c", "brier_score", "td_auc"))
```

`tau` controls restriction/integration; `t_star` selects point metrics.
For right-censored Brier/AUC, the implementation uses the nearest observed
time, with ties toward the earlier time. CIF and two-phase paths use the
last grid time not exceeding `t_star`. These are not identical conventions.
Choose and report an evaluation horizon within the support of the data.
Use the same horizon for prediction and evaluation; passing a new `tau` to
a summary does not regenerate an already supplied `pred` vector.

## Competing-risk example

`tm_survival_eval_cr()` and `tm_summarize_cr()` accept `censoring="cox"`
(default) or `"km"` for generalized R-squared. Optional numeric
`censor_covariates` are used only with Cox censoring and must be aligned with
the prediction subjects. These metrics weight observed events of any cause at
their original follow-up times; censored subjects receive zero weight, including
those followed beyond the outcome horizon. Censoring support is therefore
needed at original contributing event times, not just at `tau`. Probabilities
are floored at 1e-4. Brier score, AUC and concordance keep their own adjustments.

```r
cr <- tm_simulate_fine_gray(200, v = 2, beta1 = c(0.5, -0.5),
                           beta2 = c(-0.3, 0.3), censor = 0.3, seed = 6202)
cr$time <- cr$obs.times
cr$status <- cr$obs.event
train <- cr[1:140, ]
test <- cr[141:200, ]
m1 <- survival::coxph(survival::Surv(time, status == 1) ~ X1 + X2,
                      data = train, x = TRUE, y = TRUE)
m2 <- survival::coxph(survival::Surv(time, status == 2) ~ X1 + X2,
                      data = train, x = TRUE, y = TRUE)
p <- tm_predict_cif(model1 = m1, model2 = m2, newdata = test,
                    covs = c("X1", "X2"), event.type = 1)
tm_summarize_cr(list(CSH = p), event_type = 1)
# This now also works without requesting a Brier score:
tm_summarize_cr(list(CSH = p), event_type = 1, metrics = "td_auc")
```

## Case-cohort and NCC usage

Compute weights with full-cohort outcomes. For case-cohort, retain subjects
in the subcohort or with an event, then order those rows and their weights
together. Fit a censoring distribution from appropriate cohort data and pass
it as `km_cens` (summary) or `km_cens_fit` (direct evaluator). The executable
[API example](api-smoke.R) shows this sequence, including all 17 exports.

NCC `m` counts **controls per case**. With one case and two eligible controls,
`m=1` gives each control weight 2; `m=2` gives weight 1. The function computes
inclusion probabilities for the documented sampling scheme; it does not
infer which controls were actually sampled or the probabilities for an
arbitrary externally designed matched study.

## Deprecated exported functions

All 15 wrappers forward to a current function and issue a deprecation warning.

| Deprecated name | Replacement |
|---|---|
| `pam.predicted_survial_eval` | `tm_survival_eval` |
| `pam.predicted_survial_eval_cr` | `tm_survival_eval_cr` |
| `pam.survival_eval` | `tm_fit_and_eval` |
| `pam.coxph_restricted` | `tm_predict_coxph` |
| `pam.surverg_restricted` | `tm_predict_survreg` |
| `pam.predict_cr` | `tm_predict_cif` |
| `pam.summary` | `tm_summarize` |
| `pam.summary_cr` | `tm_summarize_cr` |
| `pam.sample_design` | `tm_sample_design` |
| `cc_weights` | `tm_case_cohort_weights` |
| `ncc_weights` | `tm_nested_case_control_weights` |
| `plot_pred` | `tm_plot_pred` |
| `summary_pred_plot` | `tm_plot_summary` |
| `sim_cox_weibull_censored` | `tm_sim_cox_weibull` |
| `simulateTwoCauseFineGrayModel` | `tm_simulate_fine_gray` |

Some paper names were never public exports: `pam.predict_subject_cif`,
`pam.survreg_restricted`, and the separate case-cohort/NCC evaluator names.
Use the replacements above and `tm_evaluate_two_phase()`.

The bundled `moore` dataset is a 48-row Moore's-law transistor-count dataset
(`year`, `time`, `count`), not a survival validation cohort. It is documented
separately with `?moore`.
