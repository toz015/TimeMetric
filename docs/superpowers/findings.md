# TimeMetric — Findings Log

Behaviour pinned by characterization tests that appears incorrect. Each entry is
resolved deliberately in a later step, never by quietly changing an expected value.

| # | Location | Observed behaviour | Why it looks wrong | Status |
|---|---|---|---|---|
| 1 | `R/pam.rsh_metric.R:48-51` | Sorts `data` by `survival_time`, then computes `Mtx` from the unsorted `predicted_data` argument rather than the reordered column | Predictions are misaligned with the status vector they multiply, so `Dx` and therefore `R_sh` are wrong whenever input is not already time-sorted | Open |
| 2 | `R/pam.rsph.R:79,189,284` | Uses `require(survival)` inside package code | `R CMD check` flags `require()` in package code; imports belong in NAMESPACE | Open |
| 3 | `R/pam.rsh_metric.R` roxygen | Example calls `library(PAmeasure)` | Stale reference to the predecessor package | Open |
| 4 | `R/pam.predicted_survial_eval.R:82-85` | `valid_metrics` and `default_metrics` use U+2019 curly apostrophes in `Harrell's C` / `Uno's C` | Users must type a typographic apostrophe for metric selection to match. Spec section 5 records this for the two-phase function only; it occurs here too | Open |
| 5 | `R/pam.predicted_survival_eval_cr.R:66-76` | Runs an `"all" %in% metrics` / `setdiff` validation block before the `is.null(metrics)` default is applied, then repeats the same block afterwards | Duplicated logic; the first block is a no-op on the `NULL` default | Open |
| 6 | `NAMESPACE` | **FIXED (group 1).** Was: `R CMD INSTALL` fails outright: `lazy loading failed for package 'TimeMetric' ... there is no package called 'randomForestSRC'`. The package imports it in `NAMESPACE` but declares it in no `DESCRIPTION` field. Deleting that one `importFrom` line makes installation succeed | **Most severe finding.** The package cannot be installed on any machine without `randomForestSRC` already present, so the README's `remotes::install_github()` command fails for a reviewer. The import is also stale -- there are no call sites -- so nothing needs it. `R CMD check` aborts at the dependency stage with 1 ERROR before any other check runs | Open |
| 7 | `R/pam.rsph_metric.R` | Returns `list(r2, numerator, denominator)`, where `r2 == numerator/denominator`. No component is named `Re` or `R_E` | Relevant to spec section 5's open question "does `R_sph` equal `R_E`?" -- the component naming suggests an R-squared-type ratio, which must be compared against the `R_E` produced by `pam.summary.rsph(pam.rsph(...))` before the two names are merged | Open |
| 8 | `R/pam.surverg_restricted.R:65` vs `R/pam.coxph_restricted.R:8` | `pam.surverg_restricted` errors with "model must include x = TRUE and y = TRUE"; `pam.coxph_restricted` documents the same requirement in its roxygen but succeeds without it | Asymmetric enforcement of an identical documented contract. Either the coxph path silently relies on something else, or the survreg check is stricter than necessary | Open |
| 9 | `R/pam.predict_cr.R:179,187` | Uses `aft.model$coefficient` (singular) while lines 176 and 190 use `aft.model$coefficients` | Works only because R partial-matches `$` on lists. Would break silently if a survreg object ever gained another `coefficient*` element | Open |
| 10 | `R/pam.predict_cr.R:93` -> `get_CIF_aft` | The `survreg` dispatch path indexes `aft.other` coefficients as `[-1]` to drop an intercept, so `model2` must also be an AFT fit. Passing a `coxph` (no intercept) fails with `non-conformable arguments` from a matrix product | The constraint is undocumented and the error names neither the argument nor the requirement. A user pairing a Cox cause-specific model with a survreg focus model gets an opaque linear-algebra error | Open |
| 11 | `R/pam.predicted_survial_eval.R:162,169` | **FIXED (group 1)** by `importFrom(survival, concordancefit)`. Was: calls `concordancefit()` unqualified. `survival` exports it but `NAMESPACE` never imports it, so `library(TimeMetric)` alone gives `could not find function "concordancefit"` | **Severe.** The package's primary evaluation function only works if the user separately attaches `survival`. It has gone unnoticed because `paper.code.Rmd` calls `library(survival)` first. `R CMD check` reports this as "no visible global function definition" | Open |
| 12 | `R/pam.survial_eval.R:110,113` | Calls `pam.coxph_restricted(..., covariates = covariates, newdata = test_data)` and the same for `pam.surverg_restricted`, but those take `covs =` and `new_data =`. Every call fails with `unused arguments (covariates = ..., newdata = ...)` | **FIXED.** **Severe.** The exported, documented `pam.survival_eval()` cannot run on any input -- it has never worked. Both call sites use two wrong argument names each | Open |
| 13 | `R/pam.predicted_survial_eval.R` | **SUPERSEDED by #30 -- see there.** `R_E` evaluates to 0.311868 on the standard fixture while `pam.rsph_metric()$r2` gives 0.320817 on the same data and risk scores | Direct evidence for spec section 5's open question: `R_sph` and `R_E` are **not** the same quantity and must not be merged into one canonical metric name | Open |
| 14 | `R/pam.predicted_survival_eval_cr.R` vs `R/pam.predicted_survial_eval.R:198` | `pam.predicted_survial_eval_cr()` returns `Value` as a **character** column; `pam.predicted_survial_eval()` returns it as **double** | **FIXED.** The two evaluation entry points disagree on the type of the same column, so a user moving between survival and competing-risks results must call `as.numeric()` in one case and not the other. `is.finite()` silently returns FALSE for every competing-risks value | Open |
| 15 | `R/pam.predicted_survial_eval_two_phase.R` | The evaluator returns the metric labelled `Time Dependent AUC` (capitalised), while its own summary wrapper `pam.sample_design()` returns `Time Dependent Auc` for the identical quantity | A fifth spelling of the time-dependent AUC metric, and the two differ between a function and the wrapper that calls it, so joining results from the two entry points on `Metric` silently fails to match | Open |
| 16 | `R/pam.predicted_survial_eval_two_phase.R:45` | The two-phase path reports `Pesudo_R` where the right-censored path reports `Pseudo_R_square` and `Pseudo_R2_point` | Three spellings of the pseudo R-squared family across the package, one of them misspelled, all user-facing in the `Metric` column | Open |
| 17 | `NAMESPACE`, `R/pam.rsph.R:74` | `NAMESPACE` contains **no `S3method()` directives at all**, so `pam.rsph.coxph` / `.survreg` / `.aareg` are invisible to `UseMethod` outside the package namespace. A clean session calling `TimeMetric:::pam.rsph(fit, ...)` gets `no applicable method for 'pam.rsph' applied to an object of class "coxph"` | **FIXED.** Dispatch works inside the package (which is why `R_E` still computes) and inside testthat, whose test environment inherits the namespace -- so the defect is invisible to ordinary tests. Verified in a clean subprocess. Registration is required before `pam.rsph` or any renamed equivalent can be exported | Open |
| 18 | `R/plot.R:71-78` | `plot_pred()` has `sample_index = NULL` by default and then subsets with `x_var[sample_index]`. `x[NULL]` is a zero-length vector, so the plot's data frame has **0 rows** and the chart is blank | **FIXED.** **Severe and user-facing.** The documented default call `plot_pred(pred)` silently returns an empty ggplot with no warning or error. Passing `sample_index = seq_len(n)` renders correctly, so the default value alone breaks the function | Open |
| 19 | `man/pam.predict_cr.Rd`, `man/pam.survival_eval.Rd` | Both roxygen example blocks are **syntactically invalid R**. `pam.predict_cr` has a malformed call (`... drop = FALSE], failcode = 2)` followed by `data = train_data,`); `pam.survival_eval` has the prose `Use Mayo` sitting inside the code block | **FIXED.** `R CMD check` reports `parse error in file 'lines'` and then **`will not attempt to run examples`** -- so *no* example in the package is ever executed. This is the remaining ERROR after group 1 | Open |
| 20 | `man/plot_pred.Rd` vs `R/plot.R:49` | The `.Rd` documents an argument `sample_size` that the function does not have, while the real `sample_index` argument is undocumented | **FIXED.** Explains how FINDING 18 survived: the argument whose `NULL` default blanks the plot is absent from the documentation, so no reader could see the requirement | Open |
| 21 | `DESCRIPTION` / `LICENSE` | `checking DESCRIPTION meta-information ... NOTE: License stub is invalid DCF` | **FIXED.** The `MIT + file LICENSE` stub is malformed, so the license is not machine-readable | Open |
| 22 | `R/pam.predicted_survial_eval.R`, `_two_phase.R`, `pam.prediction_survial_eval.R`, `pam.survial_eval.R` | `checking code files for non-ASCII characters ... WARNING` -- four source files contain non-ASCII (the curly apostrophes of findings #4 and #15) | **FIXED.** Confirms findings 4/15 at the packaging level: portable packages must use `\uxxxx` escapes in R code | Open |
| 23 | `R/pam.rsph.R:79,189,284` | `checking dependencies in R code ... WARNING: 'library' or 'require' call to 'survival' in package code` | **FIXED.** `R CMD check` confirming FINDING 2 | Open |
| 24 | `data/moore.rda` | `Undocumented code objects: 'moore'` / `Undocumented data sets: 'moore'` | **FIXED.** The shipped dataset has no documentation entry, a CRAN blocker | Open |
| 25 | `R/pam.censor.R:28` and `R/pam.predicted_survival_eval_cr.R:304` | Two different functions are both named `pam.censor`. The competing-risks file's definition is sourced later and silently shadows the other. Their bodies are not identical: the shadowing one wraps R-squared and L-squared in `format(round(...), nsmall = 4)`, returning character | This was the root cause of FINDING 14 -- the character values propagated through `unlist()` and made the whole `Value` column character. Two same-named functions with different behaviour is a latent hazard regardless | Open |
| 26 | `R/pam.coxph_restricted.R`, `R/pam.surverg_restricted.R` | The non-predict branch returns a component named `Psuedo.R` -- "Psuedo", a transposition | A **fourth** spelling of the pseudo R-squared family, alongside `Pesudo_R`, `Pseudo_R_square` and `Pseudo_R2_point` | Open |
| 27 | `R/pam.survial_eval.R` vs `R/pam.predicted_survial_eval.R` | `pam.survival_eval()` returns a **wide** frame (one row per model, one column per metric) while `pam.predicted_survial_eval()` returns a **long** `Metric`/`Value` frame | Two exported evaluation entry points with incompatible output shapes, so results cannot be combined without reshaping | Open |
| 28 | `R/pam.rsph_metric.R` roxygen vs signature | Documents `(predicted_data, survival_time, status)`; the function is `(time, status, risk_score)`. Both the **names and the order** differ, and the documented return component `Re` is actually `r2` | Following the documentation passes predictions where the function expects times. On the standard fixture this errors with `Invalid status value, converted to NA` from an internal `Surv()` -- obscure, and it could misbehave silently on inputs where the transposed arguments happen to look valid | Open |
| 29 | `R/pam.rsph_metric.R` `@references` | Cites Schemper & Henderson (2000) -- the `R_sh` paper -- although the function's own title is "RE Measure of Explained Variation" | Wrong attribution: `R_E` is Stare, Perme & Henderson (2011) in the manuscript. A reader checking the method against the citation is sent to the wrong paper | Open |
| 30 | `pam.rsph_metric` vs `pam.rsph` + `pam.summary.rsph` | **Both are implementations of the same manuscript metric, `R_E`** (both roxygen titles read "RE Measure of Explained Variation"; "sph" = Stare-Perme-Henderson). They disagree numerically: 0.320877 vs 0.311868, 0.392610 vs 0.385779, 0.442465 vs 0.420299 | **Revises finding #13.** The gap is a disagreement between two implementations of one metric, not two distinct quantities. One of them is wrong, or they use different estimators of the same estimand and only one matches the manuscript. Requires the authors' methodological determination before either is deleted or kept | Open |

## Status after the check-cleanup phase

`R CMD check --as-cran`: **0 errors, 0 warnings, 2 NOTEs** (was: 1 ERROR aborting
at the dependency stage).

Remaining NOTEs:

1. **`New submission`** -- inherent to a package not yet on CRAN. Not actionable.
2. ~~Undefined global functions `pam.concordance_metric` and `prediction_metrics`~~
   **RESOLVED** by deleting Cluster B (2026-09-03). The NOTE was caused entirely
   by that dead code calling functions which existed nowhere in the package.

## Status after Cluster B removal

`R CMD check --as-cran`: **0 errors, 0 warnings, 1 NOTE** -- and that NOTE is only
`New submission`, which is inherent to a package not yet on CRAN.

Deleted: `pam.prediction_survial_eval`, `pam.prediction_metrics`,
`pam.prediction_metrics_cr` (408 lines, 3 files). Each was its own file, generated
no `.Rd`, and was referenced nowhere in tests, README, or the paper.

**Cluster A** (`pam.coxph`, `pam.nlm`, `pam.survreg`, `pam.print.rsph`) and
**Cluster C** (`pam.rsh_metric`, `pam.rsph_metric`, `pam.Brier_metric`) remain.
Cluster C is deliberately retained until the two-phase equivalence gate of spec
section 5 is resolved -- it may be the better foundation for
`tm_evaluate_two_phase` than the current model-coupled path.

## Cluster C resolution (2026-09-06)

The equivalence gate and the `R_E` audit resolved all three Cluster C functions.
All are **deleted**.

| Function | Decision | Evidence |
|---|---|---|
| `pam.Brier_metric` | deleted | Equivalence *demonstrated*: matches the two-phase Brier under unit weights across five datasets to within the two-phase's own 4-decimal rounding. No weight argument, so it cannot express case-cohort or NCC estimation. Superseded by `yardstick::brier_survival`. |
| `pam.rsh_metric` | deleted | Order-dependent: `R_sh` moves from 0.1049 to 0.0264 on identical data with rows permuted. Does not match the working `R_sh` path, with a gap that changes sign across datasets. |
| `pam.rsph_metric` | deleted | A second implementation of `R_E`, wrong by 1.3-5.3%. Root cause: plain ranks instead of inverse-censoring-weighted ranks. `pam.rsph` reproduces the authors' reference exactly (0.00e+00 difference, five datasets). |

Findings #1, #7, #13, #28 and #29 are all closed by these deletions -- each
described a defect in code that no longer exists. Finding #13's framing was
wrong and is superseded by #30: `R_sph` and `R_E` were one metric implemented
twice, not two quantities.

`R_E` now has exactly one implementation, pinned against the authors' reference
by `tests/testthat/test-r-e-reference.R` with locally stored expected values.
| 31 | `R/pam.rsph.R:390` (was `pam.print.rsph`) | The reference names this `print.re`, an S3 method on `print` for class `"re"`. TimeMetric renamed it `pam.print.rsph`, which dispatches on a `pam.print` generic that does not exist, so the print method was dead. **FIXED**: restored as `print.rsph` with `S3method(print, rsph)`; `print(obj)` on an `rsph` object now works | The rename broke intended functionality rather than merely being cosmetic. Deleting it would have discarded a working feature; restoring it recovers one | Fixed |
| 32 | `pam.summary` vs `pam.summary.rsph` | `pam.summary` is an exported ordinary function, not a generic, yet `pam.summary.rsph` is named as though it were its S3 method. A user calling `pam.summary()` on an `rsph` object reaches the models-list function and gets an error about a non-empty named list, never the rsph summary | A naming collision that looks like S3 dispatch but is not. Should be resolved during the `tm_` rename -- either make the rsph summary a real `summary.rsph` method or give it a non-colliding name | Open |

## Metric name standardization (2026-09-06)

Canonical set, emitted by every entry point and returned by `tm_metric_names()`:

`brier_score`, `c_index`, `harrell_c`, `l_square`, `l2_point`, `pseudo_r2`,
`pseudo_r2_point`, `r_e`, `r_sh`, `r_square`, `r2_point`, `td_auc`, `uno_c`

Legacy spellings are still accepted at the `metrics` argument and resolve with a
deprecation warning naming the replacement. Matching is case-insensitive and
treats `_`, space, and straight or curly apostrophes as equivalent.

Findings closed by this change:

* **#4, #22** -- the U+2019 curly apostrophes in `Harrell's C` / `Uno's C` are
  gone from the emitted labels; users no longer need to reproduce a typographic
  apostrophe for metric selection to match.
* **#15** -- the evaluator and `tm_sample_design` both emit `td_auc`, so results
  from the two entry points can be joined on `Metric`. Previously one said
  `Time Dependent AUC` and the other `Time Dependent Auc`.
* **#16, #26** -- the four spellings of the pseudo R-squared family
  (`Pesudo_R`, `Pseudo_R_square`, `Psuedo.R`, `Pseudo_R2_point`) collapse to two
  canonical names. `pseudo_r2` and `pseudo_r2_point` remain **distinct**: they
  are an integrated measure and a point-in-time estimate and differ numerically
  on identical data (0.3868 vs 0.1429 on the standard competing-risks fixture).
* **#30** -- `R_sph` and `R_E` both map to `r_e`, following the audit that showed
  them to be one metric.

Metric **values** are unchanged. Every `snap_num(...$Value)` snapshot passed
untouched through this change; only the `Metric` label snapshots moved.
| 33 | `man/tm_predict_cif.Rd` example | Uses `event.type = 1` on `pbc`, where 1 is transplant and 2 is death, so transplant is the event of interest and death the competing risk | A valid competing-risks configuration but an unconventional choice for this dataset; death-as-primary is the usual framing. Deferred as a scientific decision for the authors, not changed unilaterally | Open |
| 34 | `tm_predict_survreg()` vs `tm_predict_coxph()` | `tm_predict_survreg()` returns a **named** `pred` vector; `tm_predict_coxph()` returns unnamed, for the same quantity | Inconsistent return shape between two functions documented as interchangeable inputs to the evaluators. A one-line `unname()` harmonises them, but it changes a return value, so deferred for approval | Open |
| 35 | **FIXED.** `R/pam.schemper.R:104-106`, with the `rms::cph()` calls at `R/tm_survival_eval.R:213` and `R/tm_fit_and_eval.R:144` | **`r_sh` is computed from a degenerate baseline survival curve.** `cph()` is called without `surv = TRUE`, so the fit has no `$surv` component: `train.fit$surv` is `NULL`, and `train.fit$time` partial-matches `time.inc`, a scalar (0.5). `approx(scalar, NULL, ...)` therefore returns `yleft` below 0.5 and `yright` above it -- **2 unique values across 148 event times** instead of the model's baseline hazard | **Severe and methodological.** `surv.tot.cox <- surv0.tot.cox^exp(lin.pred)` is built on a step function with a single step, not the fitted baseline. This explains the anomalous `r_sh` values observed throughout, including negative "explained variation" (-0.026 on the standard fixture). Fixing it changes every `r_sh` value the package reports | Open |
| 36 | **FIXED.** `rms` dependency | `rms` is used for exactly two things, and **neither requires it**: `rms::predictrms(fit, newdata, "lp")` is numerically identical to `predict(coxph_fit, newdata, type = "lp")` (correlation 1, mean difference 8e-18), and the baseline curve from `rms::cph(surv = TRUE)` is identical to `survival::survfit(coxph_fit)` (148/148 values agree to 1e-6) | `rms` is what forces `Depends: R (>= 4.4.0)`. Removing it lowers the floor to about R 4.1 and drops a heavy dependency, while fixing finding 35 at the same time. The substitution is verified numerically, not assumed | Open |
| 37 | `R/pam.rsph.R:485,488` (`my.survfit`) | With heavily tied event times, `r_e` emits repeated warnings: `longer object length is not a multiple of shorter object length`, from `n.cens - n.incom` and `1 - n.event/n.risk`. A result is still produced | Pre-existing, inherited from the `Re.r` port, and **unrelated to the r_sh correction** -- `R/pam.rsph.R` was not modified. Vector-length mismatch under ties means the censoring weights may be silently recycled, so `r_e` under ties is of unverified correctness | Open |
| 38 | `R/pam.rsph.R` | **92% of this file (302 of 329 normalised lines) is verbatim third-party code** from `Re.r`, the Stare-Perme-Henderson reference implementation (SHA-256 `cf653206...`). That file carries **no licence statement**, so it is all-rights-reserved by default. TimeMetric is distributed under MIT | **CRAN blocker.** CRAN policy requires that "the ownership of copyright and intellectual property rights of all components of the package must be clear and unambiguous", and that copied code be acknowledged with its licence respected. MIT-licensing another party's code without permission is not something the authors can do unilaterally. `Authors@R` lists no `cph` role and does not name Stare, Perme or Henderson; the only attribution is a URL in a roxygen comment, and the `@references` cites the wrong paper (Schemper & Henderson, which is `r_sh`, not `r_e`) | Open |
| 39 | **CONFIRMED, source identified.** `R/pam.schemper.R` is **91% verbatim (88 of 97 normalised lines) from `survAUC::schemper()`**, which is on CRAN under **GPL-2**, with an identical `(train.fit, traindata, newdata)` signature. Originally flagged because it contains Italian-language identifiers throughout -- `tempi.eventi`, `num.sogg`, `ind.censura`, `stima.surv`, `f.assegna.surv`, `primo`, `secondo` -- and cites Lusa, Miceli & Mariani (2007), "Estimation of predictive accuracy in survival analysis using R and S-PLUS", whose authors are at the Istituto Nazionale Tumori, Milan and who published code with that paper | Strongly indicates this is a **second port of third-party code**, with the same licensing and attribution questions as finding 38. Provenance should be established before submission; unlike `Re.r` no source URL is recorded, so the origin needs confirming with the authors | Open |
| 40 | `R/pam.schemper.R` vs `survAUC::schemper()` | **`survAUC` is not defective; TimeMetric's calling convention was.** `survAUC::schemper()` requires a fit carrying the baseline, i.e. `rms::cph(..., surv = TRUE)`. Given one it returns 0.513286, matching the independent Schemper-Henderson reference of 0.513287. Given `cph()` without `surv = TRUE` it returns 0.193964 | Refines finding 35. The estimator logic was always correct; TimeMetric never supplied a fit with the baseline stored, so the interpolation fell back to a degenerate curve. The corrected TimeMetric implementation reaches 0.513287 via `survival::survfit()` without needing `rms` at all | Closed |

