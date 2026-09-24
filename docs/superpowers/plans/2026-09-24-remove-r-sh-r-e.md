# Removing `r_sh` and `r_e` — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove the `r_sh` and `r_e` metrics and their two borrowed implementation files from the TimeMetric R package, leaving every other metric numerically unchanged.

**Architecture:** Work inward-out. First pin the surviving metric values to a committed baseline so any numerical drift is caught immediately. Then add a tombstone at the single shared name-resolution chokepoint (`tm_normalize_metrics()`, called by all four public entry points), strip the metrics from the two evaluators that compute them, and only then delete the implementation files — so the package loads and tests pass at every commit. Documentation, test cleanup and full-package verification follow.

**Tech Stack:** R (>= 4.1.0), testthat edition 3 with `expect_snapshot_value`, roxygen2 8.1.0, `R CMD build` / `R CMD check --as-cran`.

**Spec:** `docs/superpowers/specs/2026-09-24-remove-r-sh-r-e-design.md`

## Global Constraints

- **Scope is the R package only.** `paper.md` and `paper.code.Rmd` MUST NOT be modified. `git diff --stat` at the end must not list either file.
- **Do not rewrite git history, do not force-push, do not push to origin, do not submit to CRAN.** Task 9 prepares the history-purge command as a reviewable document only.
- **Tombstone message, exact text:** `r_e and r_sh were withdrawn before the first CRAN release; see NEWS.md.`
- **Never write "removed in 0.2.0"** anywhere. Version 0.2.0 *is* the first CRAN release; nothing was released to remove them from.
- **Make no claim about the metrics' statistical performance** in any file. The rationale is "not sufficiently validated for the first CRAN release" — a maintainer decision, not a finding.
- **NEWS.md wording, exact text:** `Removed `r_sh` and `r_e` from the pre-release API following a maintainer decision to exclude metrics that are not sufficiently validated for the first CRAN release.`
- **R sources must stay ASCII-only** (existing package constraint; `R CMD check` warns otherwise). Use `\uXXXX` escapes if a non-ASCII character is unavoidable.
- Canonical metric count after this work: **11** (`brier_score`, `c_index`, `harrell_c`, `l2_point`, `l_square`, `pseudo_r2`, `pseudo_r2_point`, `r2_point`, `r_square`, `td_auc`, `uno_c`).
- Work on branch `joss-revision`.

## Commit Grouping

**Amended 2026-09-24 on maintainer instruction.** Tasks 1, 6, 7, 8 and 9 each
commit on completion as written. **Tasks 2-5 are one atomic commit group.**

Task 2 deliberately leaves the package in an inconsistent state — name-based
requests rejected while the evaluators still compute both metrics by default —
and Tasks 3 and 4 each leave call sites pointing at files that Task 5 deletes.
Committing mid-group would put a broken package in the history.

Therefore: **run the individual `git commit` steps in Tasks 2, 3, 4 and 5 only
as `git add`**, and make a single commit at the end of Task 5 once the relevant
tests pass again. The red-green cycle still runs inside each task; only the
commit is deferred. The group commit message is given at the end of Task 5.

---

### Task 1: Pin the surviving metric values to a committed baseline

The safety net for the whole plan. This baseline is recorded **before** any code changes, so the invariance test proves the surviving metrics are untouched rather than merely self-consistent. The test must pass now and still pass at Task 8.

**Files:**
- Create: `tests/testthat/fixtures/metric-baseline.csv`
- Create: `tests/testthat/test-metric-invariance.R`

**Interfaces:**
- Consumes: `fx_surv()`, `fx_cox()`, `fx_covs()` from `tests/testthat/helper-simdata.R`; `snap_num()` from `tests/testthat/helper-expectations.R`
- Produces: `tests/testthat/fixtures/metric-baseline.csv` with columns `Metric,Value` — read by `test-metric-invariance.R` in this task and unchanged by every later task.

- [ ] **Step 1: Generate the baseline from the current, unmodified package**

Run this from the package root. It writes the CSV; it is a one-off generator, not a test.

```bash
mkdir -p tests/testthat/fixtures
Rscript -e '
  pkgload::load_all(".", quiet = TRUE)
  source("tests/testthat/helper-simdata.R")
  d   <- fx_surv()
  p   <- tm_predict_coxph(model = fx_cox(), covs = fx_covs(),
                          new_data = d, tau = 10e10)
  res <- tm_survival_eval(
    predicted_probability = p$pred, event_time = p$time,
    status = p$status, tau = 10e10, model = fx_cox(),
    new_data = d, covariates = fx_covs()
  )
  keep <- !res$Metric %in% c("r_e", "r_sh")
  out  <- data.frame(Metric = res$Metric[keep],
                     Value  = round(as.numeric(res$Value[keep]), 6))
  write.csv(out, "tests/testthat/fixtures/metric-baseline.csv", row.names = FALSE)
  print(out)
'
```

Expected: a table printed with **no** `r_e` or `r_sh` row, and a CSV written.

If `tm_survival_eval()` requires different arguments than shown, read `R/tm_survival_eval.R:60-90` and adjust the call — do **not** change what is filtered or rounded.

- [ ] **Step 2: Write the invariance test**

Create `tests/testthat/test-metric-invariance.R`:

```r
# Numerical invariance gate for the r_sh / r_e removal.
#
# The baseline CSV was generated from the package BEFORE the metrics were
# removed. Every surviving metric must keep its exact value through the
# removal. If this test fails, the removal changed a metric it should not
# have touched -- investigate before regenerating the baseline.

test_that("surviving metric values are unchanged by the r_sh / r_e removal", {
  path <- testthat::test_path("fixtures", "metric-baseline.csv")
  skip_if(!file.exists(path), "baseline fixture not present")
  baseline <- utils::read.csv(path, stringsAsFactors = FALSE)

  d   <- fx_surv()
  p   <- tm_predict_coxph(model = fx_cox(), covs = fx_covs(),
                          new_data = d, tau = 10e10)
  res <- tm_survival_eval(
    predicted_probability = p$pred, event_time = p$time,
    status = p$status, tau = 10e10, model = fx_cox(),
    new_data = d, covariates = fx_covs()
  )

  current <- stats::setNames(round(as.numeric(res$Value), 6), res$Metric)

  for (i in seq_len(nrow(baseline))) {
    m <- baseline$Metric[i]
    expect_true(m %in% names(current),
                info = paste0("metric '", m, "' disappeared from the result"))
    expect_equal(unname(current[[m]]), baseline$Value[i],
                 tolerance = 1e-6,
                 info = paste0("metric '", m, "' changed value"))
  }
})
```

- [ ] **Step 3: Run it against the unmodified package to prove the baseline is correct**

Run: `Rscript -e 'devtools::test(filter = "metric-invariance")'`
Expected: **PASS**. A failure here means the baseline generator and the test disagree — fix that now, before any removal.

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/fixtures/metric-baseline.csv tests/testthat/test-metric-invariance.R
git commit -m "test: pin surviving metric values before removing r_sh and r_e"
```

---

### Task 2: Tombstone the metric names

`tm_normalize_metrics()` is the single chokepoint — verified called by `tm_survival_eval.R:87`, `tm_fit_and_eval.R:80`, `tm_survival_eval_cr.R:76`, and `tm_evaluate_two_phase.R:46,182`. One tombstone covers every public entry point.

After this task the evaluators still *compute* both metrics by default (their internal default vectors are untouched until Task 3); only name-based requests are rejected. That is the intended intermediate state — do not assert their absence from results yet.

**Files:**
- Modify: `R/tm_metric_names.R:12-16` (header comment), `:50-54` (alias entries), `:84-112` (`tm_normalize_metrics`)
- Test: `tests/testthat/test-metric-tombstone.R` (create)

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `tm_defunct_metrics()` — internal, no arguments, returns `character(2)`: `c("r_e", "r_sh")`. `tm_normalize_metrics(metrics, warn = TRUE)` keeps its signature and now `stop()`s on any defunct spelling.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-metric-tombstone.R`:

```r
# r_e and r_sh were withdrawn before the first CRAN release. Every spelling
# that used to resolve to them must now raise a message naming the removal,
# rather than a generic "invalid metric" error, so a user with an existing
# script learns why it broke.

tombstone_re <- "withdrawn before the first CRAN release"

test_that("every removed spelling raises the tombstone error", {
  for (nm in c("r_e", "r_sh", "R_E", "R_sh", "R_sph", "R_SH", "r e", "R_sPh")) {
    expect_error(TimeMetric:::tm_normalize_metrics(nm), tombstone_re,
                 info = paste0("spelling: ", nm))
  }
})

test_that("the tombstone fires even when mixed with valid metrics", {
  expect_error(
    TimeMetric:::tm_normalize_metrics(c("harrell_c", "r_e", "brier_score")),
    tombstone_re
  )
})

test_that("surviving metric names still normalise unchanged", {
  expect_identical(TimeMetric:::tm_normalize_metrics("harrell_c"), "harrell_c")
  expect_identical(
    suppressWarnings(TimeMetric:::tm_normalize_metrics("Harrells_C")),
    "harrell_c"
  )
  expect_null(TimeMetric:::tm_normalize_metrics(NULL))
})

test_that("an unrelated unknown name is NOT claimed by the tombstone", {
  # falls through unchanged, for the caller's own validation to report
  expect_identical(TimeMetric:::tm_normalize_metrics("not_a_metric", warn = FALSE),
                   "not_a_metric")
})

test_that("tm_metric_names() no longer offers the removed metrics", {
  nms <- tm_metric_names()
  expect_false("r_e" %in% nms)
  expect_false("r_sh" %in% nms)
  expect_length(nms, 11L)
  expect_identical(nms, sort(c(
    "brier_score", "c_index", "harrell_c", "l2_point", "l_square",
    "pseudo_r2", "pseudo_r2_point", "r2_point", "r_square", "td_auc", "uno_c"
  )))
})
```

- [ ] **Step 2: Run it to verify it fails**

Run: `Rscript -e 'devtools::test(filter = "metric-tombstone")'`
Expected: FAIL — no tombstone exists, so `tm_normalize_metrics("r_e")` returns `"r_e"` silently and `tm_metric_names()` has length 13.

- [ ] **Step 3: Delete the five alias entries**

In `R/tm_metric_names.R`, delete these five lines (currently 50-54) and the `# R-squared type measures` comment heading directly above them:

```r
    # R-squared type measures
    "R_sh"                         = "r_sh",
    "r_sh"                         = "r_sh",
    "R_E"                          = "r_e",
    "R_sph"                        = "r_e",
    "r_e"                          = "r_e",
```

Take care with the comma on the preceding line: `"c_index" = "c_index",` must keep its trailing comma, since the `# calibration and discrimination over time` block still follows.

- [ ] **Step 4: Replace the stale header comment**

Replace this paragraph at `R/tm_metric_names.R:12-16`:

```r
# R_sph and R_E were shown to be two labels for the same Stare-Perme-Henderson
# metric (docs/superpowers/r-e-implementation-audit.md), so both map to r_e.
# Pseudo_R_square and Pseudo_R2_point are genuinely distinct -- an integrated
# measure and a point-in-time estimate, differing numerically on the same data
# -- so both survive under distinct names.
```

with:

```r
# Pseudo_R_square and Pseudo_R2_point are genuinely distinct -- an integrated
# measure and a point-in-time estimate, differing numerically on the same data
# -- so both survive under distinct names.
#
# r_e and r_sh were withdrawn before the first CRAN release. Their spellings
# are retained in tm_defunct_metrics() so that a request for one reports the
# withdrawal by name, instead of falling through to a generic unknown-metric
# error.
```

- [ ] **Step 5: Add the defunct set and the tombstone check**

Insert immediately **before** the `tm_normalize_metrics` definition in `R/tm_metric_names.R`:

```r
# Metric spellings withdrawn before the first CRAN release. Checked before
# alias resolution so that both canonical and legacy spellings are reported.
#' @keywords internal
#' @noRd
tm_defunct_metrics <- function() {
  c("r_e", "r_sh")
}

#' @keywords internal
#' @noRd
tm_defunct_metric_spellings <- function() {
  c("r_e", "r_sh", "R_E", "R_sh", "R_sph")
}
```

Then, inside `tm_normalize_metrics()`, insert the check directly after the `if (is.null(metrics)) return(NULL)` line and **before** `aliases <- tm_metric_aliases()`:

```r
  defunct_key <- tolower(gsub("[ _]", "", metrics))
  if (any(defunct_key %in% tolower(gsub("[ _]", "",
                                        tm_defunct_metric_spellings())))) {
    stop("r_e and r_sh were withdrawn before the first CRAN release; ",
         "see NEWS.md.", call. = FALSE)
  }
```

Two notes. The existing `flatten()` closure is defined further down in the
function, after this insertion point, which is why the check normalises inline
rather than calling it — do not reorder the function to share it. And unlike
`flatten()`, this check needs no curly-apostrophe handling, because no defunct
spelling contains an apostrophe; that keeps the line ASCII-only, as the package
requires.

- [ ] **Step 6: Run the test to verify it passes**

Run: `Rscript -e 'devtools::test(filter = "metric-tombstone")'`
Expected: PASS, all five `test_that` blocks.

- [ ] **Step 7: Confirm nothing else regressed yet**

Run: `Rscript -e 'devtools::test()'`
Expected: `test-metric-invariance` still PASSES. Some existing tests may now fail where they pass `"r_e"`/`"r_sh"` by name — record which, they are fixed in Task 6. Do not fix them here.

- [ ] **Step 8: Stage only — do NOT commit (atomic group, see Commit Grouping)**

```bash
git add R/tm_metric_names.R tests/testthat/test-metric-tombstone.R
```

---

### Task 3: Remove both metrics from `tm_survival_eval()`

**Files:**
- Modify: `R/tm_survival_eval.R` — roxygen at `:8`, `:10`, `:23-24`, `:35-36`; `valid_metrics`/`default_metrics` at `:79-85`; the `r_sh` branch at `:177-221`; the `r_e` branch at `:223-231`; `preferred_order` at `:374-375`

**Interfaces:**
- Consumes: `tm_normalize_metrics()` from Task 2.
- Produces: `tm_survival_eval()` returns a `Metric`/`Value` frame containing neither `r_e` nor `r_sh`. Signature unchanged.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-metric-tombstone.R`:

```r
test_that("tm_survival_eval no longer emits the removed metrics", {
  d   <- fx_surv()
  p   <- tm_predict_coxph(model = fx_cox(), covs = fx_covs(),
                          new_data = d, tau = 10e10)
  res <- tm_survival_eval(
    predicted_probability = p$pred, event_time = p$time,
    status = p$status, tau = 10e10, model = fx_cox(),
    new_data = d, covariates = fx_covs()
  )

  expect_false("r_e" %in% res$Metric)
  expect_false("r_sh" %in% res$Metric)
  # the surviving defaults are all still present
  expect_true(all(c("pseudo_r2", "harrell_c", "uno_c",
                    "brier_score", "td_auc") %in% res$Metric))
})

test_that("tm_survival_eval with metrics = 'all' omits the removed metrics", {
  d   <- fx_surv()
  p   <- tm_predict_coxph(model = fx_cox(), covs = fx_covs(),
                          new_data = d, tau = 10e10)
  res <- tm_survival_eval(
    predicted_probability = p$pred, event_time = p$time,
    status = p$status, tau = 10e10, model = fx_cox(),
    new_data = d, covariates = fx_covs(), metrics = "all"
  )
  expect_false(any(c("r_e", "r_sh") %in% res$Metric))
})
```

- [ ] **Step 2: Run it to verify it fails**

Run: `Rscript -e 'devtools::test(filter = "metric-tombstone")'`
Expected: FAIL — both metrics are still in `default_metrics` and `valid_metrics`, so they appear in the result.

- [ ] **Step 3: Strip the metric vectors**

Replace `R/tm_survival_eval.R:79-85`:

```r
  valid_metrics <- c("pseudo_r2", "r_square", "l_square", 
                     "pseudo_r2_point", "r2_point", "l2_point",
                     "harrell_c", "uno_c",
                      "r_e","r_sh", "brier_score", "td_auc")
  default_metrics <- c("pseudo_r2", "pseudo_r2_point",
                       "harrell_c", "uno_c",
                       "r_e","r_sh", "brier_score", "td_auc")
```

with:

```r
  valid_metrics <- c("pseudo_r2", "r_square", "l_square",
                     "pseudo_r2_point", "r2_point", "l2_point",
                     "harrell_c", "uno_c", "brier_score", "td_auc")
  default_metrics <- c("pseudo_r2", "pseudo_r2_point",
                       "harrell_c", "uno_c", "brier_score", "td_auc")
```

- [ ] **Step 4: Delete both computation branches**

Delete the entire `if ("r_sh" %in% metrics) { ... }` block (currently `:177-221`, ending with the three closing braces before `if ("r_e" ...`) **and** the entire `if ("r_e" %in% metrics) { ... }` block (currently `:223-231`). The `check_factors()` helper is defined inside the `r_sh` block and goes with it — confirm no other reference:

```bash
grep -n 'check_factors' R/tm_survival_eval.R
```

Expected after deletion: no output.

The next surviving statement is `if ("brier_score" %in% metrics) {`.

- [ ] **Step 5: Remove the two entries from `preferred_order`**

In the `preferred_order` vector (was `:368-378`), delete these two lines:

```r
    "r_sh",
    "r_e",
```

- [ ] **Step 6: Strip the roxygen references**

Four sites in the roxygen header of `tm_survival_eval`:

- `:8` — delete the whole `\item` line beginning `For \code{"r_sh"} (Schemper-Henderson), a Cox model fitted with \code{x=TRUE, y=TRUE}`
- `:10` — delete the whole `\item` line beginning `For \code{"r_e"} (rank-based \(R^2\)), a Cox model compatible with \code{pam.rsph()}`
- `:23-24` — in the `@param new_data` text, replace `prediction from \code{model} (e.g., Schemper-Henderson \code{"r_sh"} and rank-based \code{"r_e"}). If supplied, it should contain the variables` with `prediction from \code{model}. If supplied, it should contain the variables`
- `:35-36` — delete both `\item` lines: `"r_sh" - Schemper-Henderson explained variation (R_sh)` and `"r_e" - Rank-based R^2`

After editing, check the enclosing `\itemize{}` blocks still have at least one `\item` each; if one is left empty, delete the whole `\itemize{...}` wrapper.

- [ ] **Step 7: Regenerate documentation and run the tests**

```bash
Rscript -e 'devtools::document()'
Rscript -e 'devtools::test(filter = "metric-tombstone|metric-invariance")'
```

Expected: both filters PASS. **`test-metric-invariance` passing here is the key signal** — it proves removing the two branches did not perturb any surviving metric.

- [ ] **Step 8: Stage only — do NOT commit (atomic group)**

```bash
git add R/tm_survival_eval.R man/tm_survival_eval.Rd tests/testthat/test-metric-tombstone.R
```

---

### Task 4: Remove both metrics from `tm_fit_and_eval()`

**Files:**
- Modify: `R/tm_fit_and_eval.R` — roxygen at `:3`, `:15-16`; metric vector at `:81`; the `r_e` branch at `:137-139`; the `r_sh` branch at `:141-167`

**Interfaces:**
- Consumes: `tm_normalize_metrics()` from Task 2.
- Produces: `tm_fit_and_eval()` returns a wide frame with no `r_e` or `r_sh` column. Signature unchanged.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-metric-tombstone.R`:

```r
test_that("tm_fit_and_eval no longer produces the removed metric columns", {
  d   <- fx_surv()
  res <- tm_fit_and_eval(train_data = d, covariates = fx_covs(),
                         models = "coxph", metrics = "all")

  expect_false("r_e" %in% names(res))
  expect_false("r_sh" %in% names(res))
  expect_true(all(c("pseudo_r2", "r_square", "l_square",
                    "brier_score") %in% names(res)))
})

test_that("tm_fit_and_eval rejects the removed metrics by name", {
  d <- fx_surv()
  expect_error(
    tm_fit_and_eval(train_data = d, covariates = fx_covs(),
                    models = "coxph", metrics = "r_sh"),
    "withdrawn before the first CRAN release"
  )
})
```

- [ ] **Step 2: Run it to verify it fails**

Run: `Rscript -e 'devtools::test(filter = "metric-tombstone")'`
Expected: FAIL on the first block — `metrics = "all"` still expands to a vector containing both.

- [ ] **Step 3: Strip the `"all"` expansion**

Replace `R/tm_fit_and_eval.R:81`:

```r
  metrics <- if (("all" %in% metrics))c("pseudo_r2", "r_square", "l_square", "harrell_c", "uno_c", "r_e", "r_sh", "brier_score", "td_auc") else metrics
```

with:

```r
  metrics <- if (("all" %in% metrics)) c("pseudo_r2", "r_square", "l_square", "harrell_c", "uno_c", "brier_score", "td_auc") else metrics
```

- [ ] **Step 4: Delete both computation branches**

Delete the `r_e` block (currently `:137-139`):

```r
    if ("r_e" %in% metrics) {
      metrics_results[[fit_name]]$r_e <- pam.rsph(fits[[fit_name]])$Re
    }
```

and the entire `if ("r_sh" %in% metrics) { ... }` block (currently `:141-167`), which includes its own local `check_factors()` helper. Verify:

```bash
grep -n 'check_factors\|pam.rsph\|pam.schemper' R/tm_fit_and_eval.R
```

Expected after deletion: no output.

The next surviving statement is `if ("brier_score" %in% metrics) {`.

- [ ] **Step 5: Strip the roxygen references**

- `:3` — in `@description`, replace `It provides metrics such as R_square, L_square, Pseudo_R, Harrell's C, Uno's C, R_sph (distance-based estimator for survival predictive accuracy), R_sh, Brier Score, and Time-dependent AUC.` with `It provides metrics such as R_square, L_square, Pseudo_R, Harrell's C, Uno's C, Brier Score, and Time-dependent AUC.`
- `:15-16` — delete both `\item` lines:

```r
#'     \item "r_e": Explained variation (R_sph).
#'     \item "r_sh": Explained variation (R_sh).
```

- [ ] **Step 6: Regenerate docs and run the tests**

```bash
Rscript -e 'devtools::document()'
Rscript -e 'devtools::test(filter = "metric-tombstone|metric-invariance")'
```

Expected: PASS.

- [ ] **Step 7: Stage only — do NOT commit (atomic group)**

```bash
git add R/tm_fit_and_eval.R man/tm_fit_and_eval.Rd tests/testthat/test-metric-tombstone.R
```

---

### Task 5: Delete the borrowed implementations and their registrations

Both files are now unreachable. This is the task that closes findings 38 and 39.

**Files:**
- Delete: `R/pam.rsph.R`, `R/pam.schemper.R`
- Modify: `R/TimeMetric-package.R` (the `@importFrom stats` block), `NAMESPACE` (regenerated)

**Interfaces:**
- Consumes: Tasks 3 and 4 — both call sites must already be gone.
- Produces: a namespace with no `pam.rsph*`, `print.rsph`, `summary.rsph`, `pam.schemper`, `my.survfit` or `pam.re` object.

- [ ] **Step 1: Prove nothing still references them**

```bash
grep -rn 'pam\.rsph\|pam\.schemper\|summary\.rsph\|print\.rsph\|my\.survfit\|pam\.re\b' R/
```

Expected: matches **only** inside `R/pam.rsph.R` and `R/pam.schemper.R`. If anything else matches, stop — Task 3 or 4 is incomplete.

- [ ] **Step 2: Write the failing test**

Append to `tests/testthat/test-metric-tombstone.R`:

```r
test_that("the borrowed implementations are gone from the namespace", {
  ns <- asNamespace("TimeMetric")
  for (nm in c("pam.rsph", "pam.rsph.coxph", "pam.rsph.aareg",
               "pam.rsph.survreg", "print.rsph", "summary.rsph",
               "pam.schemper", "my.survfit", "pam.re")) {
    expect_false(exists(nm, envir = ns, inherits = FALSE),
                 info = paste0("still present: ", nm))
  }
})
```

- [ ] **Step 3: Run it to verify it fails**

Run: `Rscript -e 'devtools::test(filter = "metric-tombstone")'`
Expected: FAIL — the objects still exist.

- [ ] **Step 4: Delete the files**

```bash
git rm R/pam.rsph.R R/pam.schemper.R
```

- [ ] **Step 5: Drop the two now-unused stats imports**

`stats::approx` and `stats::model.matrix` were used **only** inside the deleted files (verified: 3 and 1 call sites inside, 0 elsewhere). Find and edit the roxygen `@importFrom stats` block:

```bash
grep -rn '@importFrom stats' R/
```

Remove `approx` and `model.matrix` from that list. Leave `uniroot`, `quantile`, `reshape`, `na.omit`, `complete.cases`, `lm`, `median`, and every other entry alone — all are still used.

- [ ] **Step 6: Regenerate NAMESPACE and confirm the S3 registrations are gone**

```bash
Rscript -e 'devtools::document()'
grep -n 'rsph\|schemper\|approx\|model.matrix' NAMESPACE
```

Expected: **no output.** The five `S3method()` lines (`pam.rsph,aareg`, `pam.rsph,coxph`, `pam.rsph,survreg`, `print,rsph`, `summary,rsph`) and the two `importFrom(stats, ...)` entries must all be gone.

- [ ] **Step 7: Verify the package still installs and the invariance gate holds**

```bash
Rscript -e 'devtools::load_all("."); cat("LOADED OK\n")'
Rscript -e 'devtools::test(filter = "metric-tombstone|metric-invariance")'
```

Expected: `LOADED OK`, then PASS. A `could not find function` error here means a call site was missed in Task 3 or 4.

- [ ] **Step 8: Gate the group — the package must be consistent before committing**

```bash
Rscript -e 'devtools::test(filter = "metric-tombstone|metric-invariance")'
```

Expected: PASS. Tests outside those two filters may still fail here; they are
repaired in Task 6. If either filter fails, **do not commit** — the group is not
yet consistent.

- [ ] **Step 9: Commit the whole 2-5 group as one change**

```bash
git add -A R/ NAMESPACE man/ tests/testthat/test-metric-tombstone.R
git commit -m "refactor: remove the r_sh and r_e metrics and their borrowed implementations

Rejects r_e, r_sh and the legacy spellings R_E, R_sh and R_sph at
tm_normalize_metrics(), the shared chokepoint of all four public entry points;
drops both metrics from the defaults and evaluator branches of
tm_survival_eval() and tm_fit_and_eval(); deletes the implementations and their
S3 registrations; drops stats::approx and stats::model.matrix, now unused.

Closes findings 38 and 39: R/pam.rsph.R was 92% verbatim from the unlicensed
Re.r reference implementation, and R/pam.schemper.R 91% verbatim from
survAUC::schemper(), which is GPL-2 and incompatible with TimeMetric's MIT
licence.

Committed as one change because tasks 2-4 each leave the package internally
inconsistent; only after the deletion does it hold together again."
```

---

### Task 6: Clean up the existing test suite

Three test files test only removed behaviour and are deleted. Four others hold assertions that must be updated. Per the spec, these files are deleted from the **working tree only** — they stay in git history, because they contain no borrowed code (verified: zero matches for any `survAUC`-port identifier; `dRti` and `my.survfit` appear only as assertion targets like `expect_identical(names(res), c("times", "Rti", "dRti"))`).

**Files:**
- Delete: `tests/testthat/test-r-e-reference.R`, `tests/testthat/test-r-sh-definition.R`, `tests/testthat/test-rsph-dispatch.R`, `tests/testthat/_snaps/rsph-dispatch.md`
- Modify: `tests/testthat/test-eval-survival.R:30-45`, `:113-133`, `:157-159`, `:188-189`; `tests/testthat/test-smoke.R:32`; `tests/testthat/test-eval-competing-risks.R:120-131`
- Regenerate: `tests/testthat/_snaps/eval-survival.md`

- [ ] **Step 1: Delete the three obsolete test files and the orphaned snapshot**

```bash
git rm tests/testthat/test-r-e-reference.R \
       tests/testthat/test-r-sh-definition.R \
       tests/testthat/test-rsph-dispatch.R \
       tests/testthat/_snaps/rsph-dispatch.md
```

- [ ] **Step 2: Fix the default-metric-set test**

In `tests/testthat/test-eval-survival.R`, replace the `test_that` block starting at line 30 — its title, its comment, and its two `expect_true` lines:

```r
test_that("the default metric set includes R_sh and R_E in the Metric column", {
  # Metric names live in res$Metric; names(res) is c("Metric", "Value").
  # R_sh reaches rms::cph + pam.schemper; R_E reaches pam.rsph dispatch.
  # Asserting both here proves those two code paths execute by default.
  res <- ev_result()

  expect_true("r_sh" %in% res$Metric)
  expect_true("r_e" %in% res$Metric)
  expect_true("brier_score" %in% res$Metric)
```

with:

```r
test_that("the default metric set appears in the Metric column", {
  # Metric names live in res$Metric; names(res) is c("Metric", "Value").
  # r_sh and r_e were withdrawn before the first CRAN release; the assertion
  # that they are absent lives in test-metric-tombstone.R.
  res <- ev_result()

  expect_true("brier_score" %in% res$Metric)
```

Leave the remaining `expect_true` / `expect_false` lines in that block untouched.

- [ ] **Step 3: Delete the `R_E` provenance test**

Delete the entire `test_that("R_E comes from the canonical pam.rsph path", { ... })` block at `tests/testthat/test-eval-survival.R:113-133`. It calls `TimeMetric:::summary.rsph()` and `TimeMetric:::pam.rsph()`, which no longer exist.

- [ ] **Step 4: Fix the two-model summarize test**

At `tests/testthat/test-eval-survival.R:157-159`, delete these three lines:

```r
  # R_sh is Cox-only, so the weibull column carries NA for it
  expect_true(is.na(res$weibull[res$Metric == "r_sh"]))
  expect_false(is.na(res$cox[res$Metric == "r_sh"]))
```

- [ ] **Step 5: Fix the `tm_fit_and_eval` column test**

At `tests/testthat/test-eval-survival.R:188-189`, delete these two lines:

```r
  # R_sph and R_sh appear as separate columns, consistent with FINDING 13
  expect_true(all(c("r_e", "r_sh") %in% names(res)))
```

- [ ] **Step 6: Fix the smoke test's internals list**

At `tests/testthat/test-smoke.R:32`, replace:

```r
    "Gt", "pam.Brier", "pam.rsph", "m_cif", "my.survfit"
```

with:

```r
    "Gt", "pam.Brier", "m_cif"
```

- [ ] **Step 7: Strengthen the competing-risks unknown-metric test**

In `tests/testthat/test-eval-competing-risks.R`, the block `test_that("tm_survival_eval_cr rejects an unknown metric name", ...)` passes `metrics = "r_sh"` and expects `"Invalid metrics"`. That message is now the tombstone. Replace the expected pattern and add a genuinely-unknown case:

```r
test_that("tm_survival_eval_cr rejects a withdrawn metric name", {
  p <- cr_pred()

  expect_error(
    tm_survival_eval_cr(
      pred_cif = p$cif_pred[, -1], event_time = p$times,
      time.cif = p$cif_pred[, 1], status = p$status, event_type = 1,
      metrics = "r_sh"
    ),
    "withdrawn before the first CRAN release"
  )
})

test_that("tm_survival_eval_cr rejects an unknown metric name", {
  p <- cr_pred()

  expect_error(
    tm_survival_eval_cr(
      pred_cif = p$cif_pred[, -1], event_time = p$times,
      time.cif = p$cif_pred[, 1], status = p$status, event_type = 1,
      metrics = "not_a_metric"
    ),
    "Invalid metrics"
  )
})
```

Leave lines 114-115 (`expect_false("r_sh" %in% res$Metric)` and the `r_e` equivalent) exactly as they are — they still assert the right thing.

- [ ] **Step 8: Run the full suite and review every snapshot change before accepting**

```bash
Rscript -e 'devtools::test()'
```

Snapshot mismatches in `_snaps/eval-survival.md` are expected: `res$Metric` and `snap_num(res$Value)` both lose two entries. **Before accepting, confirm the surviving values are unchanged** — `test-metric-invariance` passing is that proof. If it fails, stop and investigate; do not accept the snapshot.

Then:

```bash
Rscript -e 'testthat::snapshot_accept()'
Rscript -e 'devtools::test()'
```

Expected: full suite PASSES, 0 failures, 0 skips outside the documented `skip_if_covr` / `skip_without_source_tree` guards.

- [ ] **Step 9: Commit**

```bash
git add -A tests/
git commit -m "test: drop the r_sh and r_e suites, update the remaining assertions"
```

---

### Task 7: Documentation, NEWS, and the Phase 2 follow-up record

**Files:**
- Modify: `README.md:45-48`, `:64-65`, `:69-70`, `:148-149`
- Modify: `NEWS.md`
- Modify: `docs/superpowers/findings.md`, `docs/superpowers/api-reference-tables.md`, `docs/superpowers/specs/2026-08-30-timemetric-package-quality-design.md`
- Create: `docs/superpowers/phase-2-followups.md`

- [ ] **Step 1: Fix the README example output**

At `README.md:45-48`, the sample output block lists eight rows. Delete the two withdrawn rows and renumber the rest:

```
#> 5            r_sh  0.13
#> 6             r_e  0.31
#> 7     brier_score  0.21
#> 8          td_auc  0.73
```

becomes:

```
#> 5     brier_score  0.21
#> 6          td_auc  0.73
```

- [ ] **Step 2: Fix the README metric table**

Delete these two rows at `README.md:64-65`:

```
| `r_sh` | Schemper & Henderson's *R*<sub>sh</sub> (2000) |
| `r_e` | Stare, Perme & Henderson's *R*<sub>E</sub> (2011) |
```

- [ ] **Step 3: Fix the applicability sentence**

Replace `README.md:69-70`:

```
Not every metric applies to every setting: `r_sh` and `r_e` are defined for
right-censored data, while competing risks use `c_index`.
```

with:

```
Not every metric applies to every setting: competing risks use `c_index` rather
than the concordance indices defined for right-censored data.
```

- [ ] **Step 4: Remove the two orphaned README references**

Delete these two lines at `README.md:148-149`:

```
Schemper, M. & Henderson, R. (2000). *Biometrics* 56(1), 249–255.
Stare, J., Perme, M. P. & Henderson, R. (2011). *Biometrics* 67(3), 750–759.
```

Then verify neither author is cited elsewhere in the file:

```bash
grep -n 'Schemper\|Stare\|Perme' README.md
```

Expected: no output.

- [ ] **Step 5: Rewrite the NEWS entry**

In `NEWS.md` under `# TimeMetric 0.2.0`, find the bug-fix bullet that begins **`**`r_sh` (Schemper-Henderson) was computed from a degenerate baseline survival curve and its values have changed.**`** and delete that entire bullet and its continuation lines. Since `r_sh` never ships in 0.2.0, NEWS must not both fix and remove the same metric.

Add a new section immediately after the `## Renamed API` section:

```markdown
## Removed metrics

* Removed `r_sh` and `r_e` from the pre-release API following a maintainer
  decision to exclude metrics that are not sufficiently validated for the first
  CRAN release. Requesting either name, or any of the legacy spellings `R_sh`,
  `R_E` and `R_sph`, now raises an error naming the withdrawal.
* `tm_metric_names()` consequently returns 11 canonical names rather than 13.
```

Do not add any statement about how the metrics performed.

- [ ] **Step 6: Update the findings log**

In `docs/superpowers/findings.md`, change the `Status` column of the last-column entries:

- Finding 37 (`my.survfit` vector recycling under ties) → `Resolved by removal`
- Finding 38 (unlicensed `Re.r` port) → `Resolved by removal`
- Finding 39 (`survAUC` GPL-2 port) → `Resolved by removal`
- Finding 35 (degenerate `r_sh` baseline) → `Moot -- metric removed`
- Finding 36 (`rms` unnecessary) → `Moot -- metric removed; rms already dropped`

- [ ] **Step 7: Update the API reference table**

In `docs/superpowers/api-reference-tables.md`, delete the `r_sh` and `r_e` rows from the metric table.

- [ ] **Step 8: Amend the original Phase 1 spec**

In `docs/superpowers/specs/2026-08-30-timemetric-package-quality-design.md`, append to the end of the `## Exit criteria` section:

```markdown
**Amended 2026-09-24** by `2026-09-24-remove-r-sh-r-e-design.md`:

* Criterion 9's `R_sph`/`R_E` merge is superseded — both metrics were removed
  rather than merged.
* Criterion 10 ("Deprecated wrappers verified: `paper.code.Rmd` still runs")
  is **deferred to Phase 2**. `paper.code.Rmd:288-293` passes an explicit
  metric vector containing `"R_E"`, so it errors until revised. The file is
  `.Rbuildignore`d and never executed by `R CMD check`, so no packaging gate
  is affected.
```

- [ ] **Step 9: Create the Phase 2 follow-up record**

Create `docs/superpowers/phase-2-followups.md`:

```markdown
# Phase 2 Follow-Ups

Deferred items, each recorded with the evidence already gathered so Phase 2 does
not repeat the investigation.

## From the r_sh / r_e removal (2026-09-24)

Full detail in `specs/2026-09-24-remove-r-sh-r-e-design.md` §4.

### 1. `paper.code.Rmd` is knowingly broken until revised

`paper.code.Rmd:288-293` builds an explicit `metrics` vector containing `"R_E"`
and passes it to `pam.summary(metrics = metrics)` at line 324. That chunk now
raises the tombstone error. The file is `.Rbuildignore`d, so no packaging gate
is affected.

Required edits, with the structural hazard called out:

| Lines | Structure | Required edit |
|---|---|---|
| 92-100 | Positional: 5 metrics <-> 5 labels | Drop `"R_E"` **and** `expression(R[E])` together |
| 288-293 | Explicit `metrics` vector | Drop `"R_E"`. **This is the chunk that errors** |
| 350-358 | Positional: 5 metrics <-> 5 labels | Drop `"R_E"` **and** `expression(R[E])` together |
| 416-425 | Already excludes `R_E` | No change |
| 555-570 | Name-keyed `metric_pretty` | Safe key removal |

`metrics_to_plot` / `metric_levels` are matched element-by-element against
`metric_labels`. Dropping a metric without dropping its label at the same index
silently mislabels every facet, with no error. The authors already did this
paired edit correctly for `R_sh` at 416-425 — follow that precedent.

**Verification constraint:** no `paper.sim*.csv` has ever been committed
(confirmed across all history), so the document cannot run end-to-end without
regenerating three 100-iteration simulation datasets. Use a static structural
check plus a reduced-iteration smoke run; do not claim a full reproduction.

### 2. `paper.md` edits, wording unapproved

`lusa2007estimation`, `schemper2000predictive` and `stare2011measure` are each
cited exactly once, all three in the clause at lines 148-150. Removing it
orphans all three and nothing else.

* Lines 95-96 — delete the `$R_{sh}$` and `$R_E$` rows (right-censored table)
* Lines 105-107 — delete the `$R_E$` row (competing-risks table). **This row was
  already inaccurate**: `tm_survival_eval_cr()` has never emitted `r_e`, as
  `tests/testthat/test-eval-competing-risks.R:115` asserts
* Lines 146-152 — cut the two metrics from the prose sentence
* Lines 222-238 — delete the three orphaned `.csl-entry` blocks

The prose is a scientific claim; the maintainer approves the exact wording.

### 3. No `paper.bib` exists

`paper.md` carries rendered citation text plus inline `.csl-entry` divs — pandoc
output. JOSS expects `paper.md` with `[@key]` citations alongside a `paper.bib`.
Unrelated to the removal, but needed before resubmission.

### 4. Scientific rationale for the withdrawal

The simulation work motivating the removal was done by a collaborator and is not
in this repository. Nothing in the package claims the metrics performed poorly.
Obtain the supporting code or results before writing the JOSS response letter.
```

- [ ] **Step 10: Confirm no manuscript file was touched**

```bash
git status --short paper.md paper.code.Rmd
```

Expected: **no output.** If either appears, revert it — they are out of scope.

- [ ] **Step 11: Commit**

```bash
git add README.md NEWS.md docs/
git commit -m "docs: record the r_sh and r_e removal and the Phase 2 follow-ups"
```

---

### Task 8: Full-package verification

No code changes. This task produces evidence, and its output is what gets reported.

**Files:**
- Create: `docs/superpowers/removal-verification.md`

- [ ] **Step 1: Run the full test suite**

```bash
Rscript -e 'devtools::test()' 2>&1 | tee /tmp/tm-test.log
```

Expected: 0 failures, 0 warnings. Record the pass count.

- [ ] **Step 2: Build the source tarball**

```bash
R CMD build . 2>&1 | tee /tmp/tm-build.log
```

Expected: `TimeMetric_0.2.0.tar.gz` created.

- [ ] **Step 3: Run `R CMD check --as-cran`**

```bash
R CMD check --as-cran TimeMetric_0.2.0.tar.gz 2>&1 | tee /tmp/tm-check.log
tail -30 /tmp/tm-check.log
```

Expected: **0 errors, 0 warnings.** Compare any notes against `docs/superpowers/baseline-check.txt`; no note may be new. If a new note appears, fix it before proceeding — do not tolerate it.

- [ ] **Step 4: Scan the built tarball for borrowed code**

```bash
mkdir -p /tmp/tm-scan && tar xzf TimeMetric_0.2.0.tar.gz -C /tmp/tm-scan
grep -rE 'tempi\.eventi|f\.assegna\.surv|newvare|incom\.sum|Re\.imp|r2nw|dRti' /tmp/tm-scan/TimeMetric/ \
  && echo "FAIL: borrowed fragment found in tarball" \
  || echo "PASS: tarball clean"
```

Expected: `PASS: tarball clean`.

- [ ] **Step 5: Confirm the manuscript is untouched**

```bash
git diff --stat main..HEAD -- paper.md paper.code.Rmd
```

Expected: **no output.**

- [ ] **Step 6: Record the evidence**

Create `docs/superpowers/removal-verification.md` with the actual observed output — the test pass count, the exact check result line, the tarball scan result, and the confirmation that no manuscript file changed. Paste real numbers from the logs, not expected ones.

- [ ] **Step 7: Commit**

```bash
git add docs/superpowers/removal-verification.md
git commit -m "docs: record verification evidence for the r_sh and r_e removal"
```

---

### Task 9: Prepare the history purge — DO NOT RUN IT

Produces a reviewed, executable document. **Nothing in this task rewrites history, force-pushes, or pushes.** The maintainer approves separately.

**Files:**
- Create: `docs/superpowers/history-purge-plan.md`

- [ ] **Step 1: Regenerate the purge path list by content scan**

The path list must be derived fresh, not copied — three paths were undiscoverable from filenames alone.

```bash
git rev-list --all | while read c; do
  git grep -l -E 'tempi\.eventi|f\.assegna\.surv|newvare|incom\.sum|Re\.imp|r2nw' "$c" -- 2>/dev/null
done | awk -F: '{print $2}' | sort | uniq -c | sort -rn
```

Expected paths carrying borrowed **implementation** (purge): `R/pam.rsph.R`, `R/pam.re.R`, `R/pam.rsph.metric.R`, `R/pam.rsph_metric.R`, `R/pam.schemper.R`, `man/pam.rsph.Rd`, `man/pam.re.Rd`, `man/pam.rsph.metric.Rd`, `man/pam.schemper.Rd`.

Expected paths carrying **legitimate discussion** (preserve): `docs/superpowers/r-e-implementation-audit.md`, `docs/superpowers/findings.md`.

Note that `R/pam.rsph.metric.R` and `R/pam.rsph_metric.R` are two distinct historical paths, differing by a dot versus an underscore. Both must be listed.

- [ ] **Step 2: Record the blob hashes**

```bash
for p in R/pam.rsph.R R/pam.re.R R/pam.rsph.metric.R R/pam.rsph_metric.R \
         R/pam.schemper.R man/pam.rsph.Rd man/pam.re.Rd \
         man/pam.rsph.metric.Rd man/pam.schemper.Rd; do
  git log --all --pretty=format:'%H' -- "$p" | while read c; do
    git rev-parse "$c:$p" 2>/dev/null
  done | sort -u | sed "s|^|$p |"
done | tee /tmp/tm-blobs.txt
wc -l /tmp/tm-blobs.txt
```

Expected: 31 lines.

- [ ] **Step 3: Write the purge plan document**

Create `docs/superpowers/history-purge-plan.md` containing: the regenerated path list from Step 1, the 31 blob hashes from Step 2, the `git filter-repo --invert-paths --path ...` command (one `--path` per purge path), and this verification procedure to run **after** any future rewrite:

```bash
# 1. Every recorded blob is unreachable
while read p h; do
  git cat-file -e "$h" 2>/dev/null && echo "STILL REACHABLE: $p $h"
done < /tmp/tm-blobs.txt

# 2. No borrowed fragment survives anywhere under R/ or man/
git rev-list --all | while read c; do
  git grep -l -E 'tempi\.eventi|f\.assegna\.surv|newvare|incom\.sum|Re\.imp|r2nw|dRti' "$c" -- R man 2>/dev/null
done
```

Both must produce no output.

The document must state explicitly, in its own words: the scan is scoped to `R/` and `man/` **deliberately**, because an unscoped scan false-positives on `docs/superpowers/findings.md` and the audit (which quote the borrowed identifiers as evidence) and on the test files (which name `dRti` and `my.survfit` as assertion targets). Those are the legitimate historical discussion this rewrite preserves. `my.survfit` and `n.incom` are excluded as markers because they also appear in `R/pam.rsph_metric.R`, itself on the purge list, so they cannot distinguish success from failure.

It must also record that tests and snapshots are **preserved** in history: they contain no copied implementation, only assertion targets and base64-serialised output.

- [ ] **Step 4: Confirm nothing was executed**

```bash
git log --oneline -1
git status --short
```

Expected: HEAD is the Task 8 commit plus this task's staged document; no rewritten history.

- [ ] **Step 5: Commit**

```bash
git add docs/superpowers/history-purge-plan.md
git commit -m "docs: prepare the borrowed-code history purge (not executed)"
```

- [ ] **Step 6: Stop and report**

Report to the maintainer: the removal is complete and verified, the purge plan is written and unexecuted, and nothing has been pushed. Await separate approval for the history rewrite, any push, and CRAN submission.

---

## Self-Review

**Spec coverage:**

| Spec section | Task |
|---|---|
| §1 Files deleted | 5 (R sources), 6 (tests and snapshot) |
| §2.1 Public API changes | 3, 4 (result contents), 2 (names and aliases) |
| §2.2 Code sites | 2 (`tm_metric_names.R`), 3 (`tm_survival_eval.R`), 4 (`tm_fit_and_eval.R`), 5 (`NAMESPACE`) |
| §2.3 Tombstone error | 2 |
| §3 Dependency audit | 5 Step 5 |
| §4 Manuscript deferred | 7 Step 9 (follow-up record), 7 Step 10 + 8 Step 5 (untouched gates) |
| §5 Other documentation | 7 |
| §6 NEWS.md | 7 Step 5 |
| §7.1 Numerical invariance | 1, re-run as a gate in 3, 4, 5, 6 |
| §7.2 Test edits | 6 |
| §7.3 New tests | 2, 3, 4, 5 |
| §8 History purge | 9 |
| §9 Exit criteria | 8 |

**Type and name consistency:** `tm_defunct_metrics()` and `tm_defunct_metric_spellings()` are defined in Task 2 Step 5 and used nowhere else; `tm_normalize_metrics(metrics, warn = TRUE)` keeps its existing signature throughout; the tombstone regex `"withdrawn before the first CRAN release"` is used identically in Tasks 2, 4 and 6; `metric-baseline.csv` is written in Task 1 Step 1 and read in Task 1 Step 2 at the same path.

**Known ordering dependency:** Task 2 leaves a deliberate intermediate state where the evaluators still compute both metrics by default while name-based requests are rejected. Task 2 Step 7 anticipates that some existing tests fail at that point and defers the fix to Task 6. Do not "fix" them early — that would hide which tests actually depended on the removed metrics.
