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
  since-withdrawn `r_sh` estimator and to extract linear predictors, both of
  which are now `survival` calls. Because `rms` requires R >= 4.4.0, removing
  it lowers the declared minimum from R 4.4.0 to **R >= 4.1.0**, the floor set
  by `survival` and `dplyr`.

## Renamed API

Every exported function now carries a `tm_` prefix, correcting the `survial` and
`surverg` misspellings and dropping the `pam.` prefix inherited from PAmeasure.
All previous names still work and forward to their replacement with a
deprecation warning; see `?"TimeMetric-deprecated"`.

* `tm_evaluate_two_phase()` is **newly exported**. The case-cohort and nested
  case-control functionality was previously internal and unreachable by users.
* `tm_metric_names()` returns the canonical metric identifiers.

## Removed metrics

* Removed `r_sh` and `r_e` from the pre-release API following a maintainer
  decision to exclude metrics that are not sufficiently validated for the first
  CRAN release. Requesting either name, or any of the legacy spellings `R_sh`,
  `R_E` and `R_sph`, now raises an error naming the withdrawal.
* `tm_metric_names()` consequently returns 11 canonical names rather than 13.
* Every other metric is numerically unchanged, verified against values recorded
  before the removal across all four public evaluators, including the
  competing-risks and two-phase paths.

## Bug fixes

* **The Brier score's censoring weights were wrong, and `brier_score` from
  `tm_fit_and_eval()` changes on data with tied times.** `Gt()`, the
  Kaplan-Meier estimate of the censoring distribution behind the
  inverse-probability-of-censoring weights, had three independent defects: it
  derived interpolation indices from a sorted table but applied them to the raw
  unsorted input, so G(t) was neither monotone nor invariant to row order; it
  identified censoring events as `status == min(status)`, which inverts on data
  with no censoring and made G decay from 1 instead of staying at 1; and it
  silently replaced a zero censoring survival with the smallest positive value,
  concealing an undefined weight.

  `Gt()` is now a reverse-Kaplan-Meier **step function** with no interpolation
  between jump times, using G(t) for the evaluation-time weight and G(t-) for
  subject-specific event-time weights, and the reverse-Kaplan-Meier tie
  convention. It returns `NA` beyond the last observed time and **errors** when
  the censoring distribution is exhausted, rather than substituting a value.
  Values agree exactly with `pec::ipcw()`.

  Reported Brier scores on tied data at an evaluation time between jump times
  were substantially overstated: on one 200-observation fixture the reported
  value was 0.75 against a correct 0.22. Metrics other than `brier_score` are
  unaffected, as is `tm_survival_eval()`, which computes the Brier score by a
  different route.

* `tm_fit_and_eval()` (was `pam.survival_eval()`) **never worked**. It passed
  `covariates=`/`newdata=` to functions taking `covs=`/`new_data=`, so every
  call failed with `unused arguments`.
* `tm_plot_pred()` (was `plot_pred()`) returned a **silently empty plot** on its
  documented default call. `sample_index` defaulted to `NULL` and the function
  then subset with `x[NULL]`. `NULL` now means "plot every subject".
* `survival::concordancefit()` was called without being imported, so the primary
  evaluation function failed under `library(TimeMetric)` alone. It worked only
  if the user happened to have attached `survival` separately.
* **`tm_survival_eval_cr()` rejected every legacy metric spelling.** It carried
  a duplicate metric-validation block that ran *before* name normalization, so
  `"Harrells_C"`, `"C_index"` and `"Brier Score"` all failed on the
  competing-risks path with `Invalid metrics:` while working on every other
  entry point. Validation now runs after normalization, as it does elsewhere,
  and the documented legacy spellings are accepted there too. Genuinely unknown
  names are still rejected.
* `tm_survival_eval_cr()` returned its `Value` column as **character** while
  `tm_survival_eval()` returned numeric, so `is.finite()` was silently `FALSE`
  for every competing-risks value. Both now return numeric.

## Metric names

Standardized to one ASCII identifier per quantity. Legacy spellings are still
accepted, case-insensitively and tolerant of apostrophe style, with a
deprecation warning:

`brier_score`, `c_index`, `harrell_c`, `l_square`, `l2_point`, `pseudo_r2`,
`pseudo_r2_point`, `r_square`, `r2_point`, `td_auc`, `uno_c`

* `Harrell's C` and `Uno's C` previously required a U+2019 curly apostrophe to
  match, which made them close to impossible to discover.
* The pseudo R-squared family carried four spellings, one of them the misspelled
  `Pesudo_R`. `pseudo_r2` and `pseudo_r2_point` remain distinct measures.
* Time-dependent AUC carried different labels in the two-phase evaluator and its
  own summary wrapper, so results could not be joined. Both emit `td_auc`.

## Removed

Thirteen unreachable functions were deleted after a call-graph analysis and, for
the metric implementations, an equivalence audit against the original
Stare-Perme-Henderson reference code:

* `pam.rsph_metric()` and `pam.rsh_metric()` were duplicate implementations of
  the two metrics that have since been withdrawn entirely (see **Removed
  metrics**).
* `pam.Brier_metric()` duplicated the Brier score already computed by a
  weight-aware implementation.
* `pam.prediction_survial_eval()`, `pam.prediction_metrics()`,
  `pam.prediction_metrics_cr()`, `pam.coxph()`, `pam.nlm()`, `pam.survreg()`
  and related helpers were unreachable from any export.

## Testing and infrastructure

* New `testthat` suite: 498 assertions across 100 tests, none skipped.
* Every surviving metric is pinned to values recorded before the `r_sh` / `r_e`
  withdrawal, across all four public evaluators including the competing-risks
  and two-phase paths.
* GitHub Actions: `R-CMD-check` across five platform/version combinations, a
  full-dependency job that fails if any test skips, test coverage, and the
  existing JOSS draft-PDF build, whose path trigger had pointed at a filename
  that never existed.
* `R CMD check --as-cran` reports 0 errors and 0 warnings.
