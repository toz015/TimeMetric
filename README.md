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
#> 2 pseudo_r2_point 0.16
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
in_subcohort <- rbinom(nrow(d), 1, 0.4) == 1
w <- tm_case_cohort_weights(d$time, d$status, subcohort = in_subcohort)
keep <- which(in_subcohort | d$status == 1)
keep <- keep[order(d$time[keep])]
cc_pred <- tm_predict_coxph(fit, covs = c("x1", "x2"), new_data = d[keep, ])
km <- survival::survfit(survival::Surv(d$time, 1 - d$status) ~ 1)
tm_evaluate_two_phase(pred_results = cc_pred, km_cens_fit = km,
                       case_weights = w[keep])
```

`tm_nested_case_control_weights(time, status, m = 2)` computes inverse-selection
weights for two controls per case. Supply matching strata via `strata` when
needed, subset to the subjects actually sampled, and align weights with the
prediction order.

## Optional backends

`tm_predict_cif()` accepts alternative model types through suggested packages:

* `cmprsk` for a Fine-Gray fit, via the `fg_model` argument
* `randomForestSRC` for a survival forest, via the `cr_model` argument

Both are optional; the function reports which package to install if one is
needed and absent.

## Renamed functions

Every exported function previously used a `pam.` prefix inherited from
PAmeasures. Those names still work and forward to their replacements with a
deprecation warning — see `?"TimeMetric-deprecated"` for the full mapping.

## Paper analyses and figures

Run these commands from this checkout. They load its local package code rather
than an older installed TimeMetric version. The analyses also need `pkgload`,
`cmprsk`, `randomForestSRC`, and `riskCommunicator` (the Framingham data source).

```sh
# Recompute all real-data analyses and export their manuscript figures/tables.
Rscript --vanilla figures/run_real_data.R

# Run or resume the full simulation; this starts substantial computation.
Rscript --vanilla figures/run_simulations.R --mode=full

# Redraw completed results only, without refitting or rerunning simulations.
Rscript --vanilla figures/redraw_original.R
# Optional: --section=real or --section=simulations
```

The presentation entry points in `figures/` preserve the selected manuscript
`.updated.pdf` style, including R/L component dots in Simulation 2. The
lower-level drivers in `analysis/` produce diagnostic plots; use the redraw
entry afterward if those drivers were run directly. `paper.code.Rmd` displays
the completed results through the same figure/table outputs; it does not rerun
models or simulations. With `rmarkdown`, `knitr` and Pandoc available, render it
with:

```r
rmarkdown::render("paper.code.Rmd", output_dir = "results/report")
# Optional: params = list(redraw_figures = TRUE)
```

Rendering requires all four full simulation parts and all three real-data
exports; smoke or incomplete results are rejected. See
[the current development record](docs/development-2026-10-18.md) for method
changes and validation. Earlier development/audit records are preserved under
`docs/superpowers/`.

Local outputs, excluded from Git:

| Location | Contents |
|---|---|
| `results/simulations/full/` | Completed simulation, models, checkpoints, raw and summary metrics, figures |
| `results/real_data/runs/<run>/` | Real-data models, predictions and diagnostics |
| `results/real_data/paper/` | Final real-data figures and tables, with their source-run manifest |

The recorded current real-data run is selected through
`results/real_data/paper/source_manifest.csv`; redraw and FHS horizon commands
read that record rather than a historical directory. FHS uses examination 3 as
time zero and a 12-year horizon (4380 days). To re-evaluate its existing training
fit at a selected horizon (changed source/settings create a separate output):

```sh
Rscript --vanilla figures/run_fhs_horizon.R --years=12
```

Keep full simulation models and checkpoints for resumption. Statistical inputs
in `R/`, `analysis/` and `DESCRIPTION` are fingerprinted; changes require a new
simulation output directory rather than silently mixing versions. Smoke runs
are optional, write to a separate directory, and are not manuscript results.
To run package checks without creating paper results:

```sh
Rscript --vanilla -e 'testthat::test_local(".", reporter="summary")'
```

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
