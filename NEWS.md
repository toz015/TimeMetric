# TimeMetric 0.2.0

Package-quality release preparing resubmission to JOSS. The evaluation
methodology is unchanged; metric *values* are identical to 0.1.0 except where a
defect is listed below.

## Installation

* **The package can now be installed.** 0.1.0 imported `randomForestSRC` in
  `NAMESPACE` without declaring it anywhere in `DESCRIPTION`, so
  `R CMD INSTALL` failed on any machine that did not already have it and
  `remotes::install_github()` could not succeed.
* All runtime dependencies are declared: `rms`, `survminer`, `expint`, `pec`,
  `magrittr`, `purrr`, `tibble` and `utils` were used but undeclared.
* `randomForestSRC` and `cmprsk` moved to `Suggests` as optional backends of
  `tm_predict_cif()`, each guarded with an informative error naming the package
  to install.
* **`rms` is no longer a dependency.** It was used only to fit a model for the
  Schemper-Henderson estimator and to extract linear predictors; both are now
  `survival` calls, verified numerically identical. Because `rms` requires
  R >= 4.4.0, removing it lowers the declared minimum from R 4.4.0 to
  **R >= 4.1.0**, the floor set by `survival` and `dplyr`.

## Renamed API

Every exported function now carries a `tm_` prefix, correcting the `survial` and
`surverg` misspellings and dropping the `pam.` prefix inherited from PAmeasure.
All previous names still work and forward to their replacement with a
deprecation warning; see `?"TimeMetric-deprecated"`.

* `tm_evaluate_two_phase()` is **newly exported**. The case-cohort and nested
  case-control functionality was previously internal and unreachable by users.
* `tm_metric_names()` returns the canonical metric identifiers.

## Bug fixes

* **`r_sh` (Schemper-Henderson) was computed from a degenerate baseline survival
  curve and its values have changed.** The Cox model feeding the estimator was
  fitted with `rms::cph()` without `surv = TRUE`, so the fit carried no `$surv`
  component; `$surv` resolved to `NULL` and `$time` partial-matched the
  `time.inc` scalar. The interpolation therefore produced a two-valued step
  function instead of the fitted baseline hazard. The metric could return
  negative values, which is impossible for an explained-variation measure.

  The baseline now comes from `survival::survfit()` on a `survival::coxph()`
  fit. Representative changes on the package's own fixtures:

  | scenario | before | after |
  |---|---|---|
  | seed 1001, 30% censoring | **-0.0262** | 0.1271 |
  | seed 2002, 30% censoring | 0.3054 | 0.1677 |
  | seed 3003, 30% censoring | 0.2784 | 0.2144 |
  | uncensored | 0.0808 | 0.1562 |

  The corrected implementation is validated against an independent estimator
  written from the published Schemper-Henderson definition, agreeing to 1e-8 on
  uncensored data, and satisfies the definitional properties: approximately zero
  for a null model, strictly increasing in predictor strength, and bounded in
  [0, 1].

* `tm_fit_and_eval()` (was `pam.survival_eval()`) **never worked**. It passed
  `covariates=`/`newdata=` to functions taking `covs=`/`new_data=`, so every
  call failed with `unused arguments`.
* `tm_plot_pred()` (was `plot_pred()`) returned a **silently empty plot** on its
  documented default call. `sample_index` defaulted to `NULL` and the function
  then subset with `x[NULL]`. `NULL` now means "plot every subject".
* `survival::concordancefit()` was called without being imported, so the primary
  evaluation function failed under `library(TimeMetric)` alone. It worked only
  if the user happened to have attached `survival` separately.
* The `pam.rsph` S3 methods were never registered, so `R_E` could not be
  dispatched from a clean session. They are registered, and `print()` and
  `summary()` now work on the returned object.
* `tm_survival_eval_cr()` returned its `Value` column as **character** while
  `tm_survival_eval()` returned numeric, so `is.finite()` was silently `FALSE`
  for every competing-risks value. Both now return numeric.

## Metric names

Standardised to one ASCII identifier per quantity. Legacy spellings are still
accepted, case-insensitively and tolerant of apostrophe style, with a
deprecation warning:

`brier_score`, `c_index`, `harrell_c`, `l_square`, `l2_point`, `pseudo_r2`,
`pseudo_r2_point`, `r_e`, `r_sh`, `r_square`, `r2_point`, `td_auc`, `uno_c`

* `Harrell's C` and `Uno's C` previously required a U+2019 curly apostrophe to
  match, which was close to undiscoverable.
* The pseudo R-squared family carried four spellings, one of them the misspelled
  `Pesudo_R`. `pseudo_r2` and `pseudo_r2_point` remain distinct measures.
* `R_sph` and `R_E` were two labels for the same Stare, Perme & Henderson
  metric; both now resolve to `r_e`.
* Time-dependent AUC carried different labels in the two-phase evaluator and its
  own summary wrapper, so results could not be joined. Both emit `td_auc`.

## Removed

Thirteen unreachable functions were deleted after a call-graph analysis and, for
the metric implementations, an equivalence audit against the original
Stare-Perme-Henderson reference code:

* `pam.rsph_metric()` was a second implementation of `R_E` that omitted the
  inverse-censoring weighting of ranks and was wrong by 1.3-5.3%. The retained
  `pam.rsph` path reproduces the authors' reference implementation exactly.
* `pam.rsh_metric()` was order-dependent: `R_sh` changed from 0.1049 to 0.0264
  on identical data with rows permuted.
* `pam.Brier_metric()` duplicated the Brier score already computed by a
  weight-aware implementation.
* `pam.prediction_survial_eval()`, `pam.prediction_metrics()`,
  `pam.prediction_metrics_cr()`, `pam.coxph()`, `pam.nlm()`, `pam.survreg()`
  and related helpers were unreachable from any export.

## Testing and infrastructure

* New `testthat` suite: 380 tests, none skipped.
* `R_E` is pinned against the authors' reference implementation with locally
  stored expected values and recorded provenance.
* GitHub Actions: `R-CMD-check` across five platform/version combinations, a
  full-dependency job that fails if any test skips, test coverage, and the
  existing JOSS draft-PDF build, whose path trigger had pointed at a filename
  that never existed.
* `R CMD check --as-cran` reports 0 errors and 0 warnings.
