# TimeMetric

<!-- badges: start -->
[![R-CMD-check](https://github.com/toz015/TimeMetric/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/toz015/TimeMetric/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/toz015/TimeMetric/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/toz015/TimeMetric/actions/workflows/test-coverage.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://www.r-project.org/Licenses/MIT)
<!-- badges: end -->

`TimeMetric` evaluates the predictive performance of survival models under a
single interface, across right-censored data, competing risks, and two-phase
sampling designs (case-cohort and nested case-control).

It separates **prediction** from **evaluation**: the evaluation functions take
predicted survival probabilities or cumulative incidence functions from any
source, so a model fitted outside the package can be assessed with the same
metrics as one fitted by it.

## Installation

```r
# install.packages("remotes")
remotes::install_github("toz015/TimeMetric")
```

## Quick start

```r
library(TimeMetric)

# simulate right-censored data and fit a model
d   <- tm_sim_cox_weibull(n = 200, pi_c = 0.3, v = 2,
                          beta = c(0.5, -0.5), seed = 1001)
d   <- d[, c("time", "status", "x1", "x2")]
fit <- survival::coxph(survival::Surv(time, status) ~ x1 + x2,
                       data = d, x = TRUE, y = TRUE)

# predicted survival probabilities over the observed time grid
pred <- tm_predict_coxph(model = fit, covs = c("x1", "x2"), new_data = d)

# evaluate
tm_summarize(list(cox = pred))
#>            Metric  cox
#> 1       pseudo_r2 0.24
#> 2 pseudo_r2_point 0.09
#> 3       harrell_c 0.66
#> 4           uno_c 0.66
#> 5     brier_score 0.21
#> 6          td_auc 0.74
```

## Metrics

`tm_metric_names()` returns the canonical identifiers accepted by every
`metrics` argument and emitted in the `Metric` column:

| Name | Measure |
|---|---|
| `pseudo_r2`, `pseudo_r2_point` | Pseudo *R²*, integrated and point-in-time (Li & Wang 2019; Zhuang et al. 2025) |
| `r_square`, `l_square`, `r2_point`, `l2_point` | Explained-variation components |
| `harrell_c`, `uno_c` | Concordance indices (Harrell et al. 1982; Uno et al. 2011) |
| `c_index` | Concordance for competing risks |
| `brier_score` | Brier score (Brier 1950; Graf et al. 1999) |
| `td_auc` | Time-dependent AUC (Heagerty, Lumley & Pepe 2000) |

Not every metric applies to every setting: competing risks use `c_index` rather
than the concordance indices defined for right-censored data.

## Functions

**Evaluation** — take predictions, return a `Metric`/`Value` table:

| Function | Setting |
|---|---|
| `tm_survival_eval()` | Right-censored |
| `tm_survival_eval_cr()` | Competing risks |
| `tm_evaluate_two_phase()` | Case-cohort and nested case-control |
| `tm_fit_and_eval()` | Fits several models from raw data, then evaluates |

**Prediction** — produce the input the evaluators expect:

| Function | Model |
|---|---|
| `tm_predict_coxph()` | Cox proportional hazards |
| `tm_predict_survreg()` | Parametric AFT |
| `tm_predict_cif()` | Cumulative incidence for competing risks |

**Summaries and weights:**

| Function | Purpose |
|---|---|
| `tm_summarize()`, `tm_summarize_cr()` | One column per model, rows as metrics |
| `tm_sample_design()` | Summary for two-phase designs |
| `tm_case_cohort_weights()` | Prentice-style case-cohort weights |
| `tm_nested_case_control_weights()` | Nested case-control weights |

**Plotting and simulation:** `tm_plot_pred()`, `tm_plot_summary()`,
`tm_sim_cox_weibull()`, `tm_simulate_fine_gray()`.

## Two-phase designs

Case-cohort and nested case-control evaluation weights each subject by the
inverse of its sampling probability:

```r
set.seed(2001)
w  <- tm_case_cohort_weights(time = d$time, status = d$status,
                             subcohort = rbinom(nrow(d), 1, 0.4))
km <- survival::survfit(survival::Surv(d$time, 1 - d$status) ~ 1)

tm_evaluate_two_phase(pred_results = pred, km_cens_fit = km, case_weights = w)
```

`tm_nested_case_control_weights(time, status, m = 2)` produces the equivalent
weights for a matched nested case-control sample.

## Optional backends

`tm_predict_cif()` accepts alternative model types through suggested packages:

* `cmprsk` for a Fine-Gray fit, via the `fg_model` argument
* `randomForestSRC` for a survival forest, via the `cr_model` argument

Both are optional; the function reports which package to install if one is
needed and absent.

## Renamed functions

Every exported function previously used a `pam.` prefix inherited from
PAmeasure. Those names still work and forward to their replacements with a
deprecation warning — see `?"TimeMetric-deprecated"` for the full mapping.

## Contributors

Developed by **Li's Lab** (UCLA): Tong Zhu, Zian Zhuang, Wen Su, Xiaowu Dai and
Gang Li.

## References

Brier, G. W. (1950). *Monthly Weather Review* 78(1), 1–3.
Graf, E., Schmoor, C., Sauerbrei, W. & Schumacher, M. (1999). *Statistics in Medicine* 18(17–18), 2529–2545.
Harrell, F. E., Califf, R. M., Pryor, D. B., Lee, K. L. & Rosati, R. A. (1982). *JAMA* 247(18), 2543–2546.
Heagerty, P. J., Lumley, T. & Pepe, M. S. (2000). *Biometrics* 56(2), 337–344.
Li, G. & Wang, X. (2019). *JASA* 114(528), 1815–1825.
Uno, H., Cai, T., Pencina, M. J., D'Agostino, R. B. & Wei, L. J. (2011). *Statistics in Medicine* 30(10), 1105–1117.
Zhuang, Z., Su, W., Kawaguchi, E. & Li, G. (2025). arXiv:2507.15040.
