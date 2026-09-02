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
| 6 | `NAMESPACE` | `pkgload::load_all()` warns `there is no package called 'randomForestSRC'`, yet the package loads and every function works | Runtime confirmation that `importFrom(randomForestSRC, predict.rfsrc)` is a stale import with no call sites, as the spec's dependency analysis predicted | Open |
| 7 | `R/pam.rsph_metric.R` | Returns `list(r2, numerator, denominator)`, where `r2 == numerator/denominator`. No component is named `Re` or `R_E` | Relevant to spec section 5's open question "does `R_sph` equal `R_E`?" -- the component naming suggests an R-squared-type ratio, which must be compared against the `R_E` produced by `pam.summary.rsph(pam.rsph(...))` before the two names are merged | Open |
| 8 | `R/pam.surverg_restricted.R:65` vs `R/pam.coxph_restricted.R:8` | `pam.surverg_restricted` errors with "model must include x = TRUE and y = TRUE"; `pam.coxph_restricted` documents the same requirement in its roxygen but succeeds without it | Asymmetric enforcement of an identical documented contract. Either the coxph path silently relies on something else, or the survreg check is stricter than necessary | Open |
| 9 | `R/pam.predict_cr.R:179,187` | Uses `aft.model$coefficient` (singular) while lines 176 and 190 use `aft.model$coefficients` | Works only because R partial-matches `$` on lists. Would break silently if a survreg object ever gained another `coefficient*` element | Open |
| 10 | `R/pam.predict_cr.R:93` -> `get_CIF_aft` | The `survreg` dispatch path indexes `aft.other` coefficients as `[-1]` to drop an intercept, so `model2` must also be an AFT fit. Passing a `coxph` (no intercept) fails with `non-conformable arguments` from a matrix product | The constraint is undocumented and the error names neither the argument nor the requirement. A user pairing a Cox cause-specific model with a survreg focus model gets an opaque linear-algebra error | Open |

