# TimeMetric Characterization Test Baseline — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Pin the current numeric behaviour of every reachable function in `TimeMetric` with a deterministic test suite, so that all later renames, dependency changes, and deletions are provably behaviour-preserving.

**Architecture:** `testthat` 3rd edition suite under `tests/testthat/`. Each test combines **hard-coded structural assertions** (component names, types, lengths, ranges, specific error messages) with a **small snapshot of selected numeric values** rounded to 6 decimal places. Whole model objects and large matrices are never snapshotted — they are reduced to compact fingerprints first. Structural assertions catch shape drift; numeric snapshots catch value drift; neither requires transcribing floats into this plan.

**Tech Stack:** R 4.5+, `testthat` (3e), `pkgload`, `roxygen2`, `covr`. (`devtools` is deliberately not used: its `gert`/`ragg` dependencies need system libraries -- libgit2, libpng, freetype, harfbuzz -- that add setup cost for no benefit here. `pkgload::load_all()` and `testthat::test_local()` provide everything this plan needs.) Package runtime deps: `survival`, `rms`, `expint`, `survminer`, `pec`, `tdROC`, `yardstick`, `ggplot2`, `patchwork`, `dplyr`, `magrittr`, `purrr`, `tibble`.

## Global Constraints

- This step **adds tests only**. No renames, no deletions, no refactors, no `DESCRIPTION`/`NAMESPACE` edits beyond adding `testthat` test infrastructure. Those are steps 2-8 of the spec.
- Tests encode behaviour **as it is today, correct or not**. When a test pins something that looks wrong, record it in `docs/superpowers/findings.md` and move on — do not fix it here.
- Every test must be deterministic. Any test touching randomness sets a seed via the shared fixtures in `tests/testthat/helper-simdata.R`.
- **Never snapshot a fitted model object, a full prediction matrix, or any object over ~100 elements.** Reduce to a fingerprint with `mat_fingerprint()` or select specific values with `snap_num()`.
- Numeric snapshots round to 6 decimal places via `snap_num()`; direct comparisons use `tolerance = 1e-6`.
- `expect_error()` must always assert a specific message fragment copied from the source. A bare `expect_error(x)` or `regexp = "."` is not acceptable — it passes on typos and missing arguments alike.
- All new files UTF-8 encoded with ASCII-only content.
- Functions are called by their **current** names (`pam.*`). The `tm_` rename is step 3 of the spec and must not be anticipated here.
- Internal (unexported) functions are reached with `TimeMetric:::`. See the export table below — getting this wrong produces "could not find function" at runtime.
- Target total suite runtime under 60 seconds. Fixtures use `n = 200`, never the `n = 3000`/`n = 10000` of `paper.code.Rmd`.
- Commit after every task. Never run `testthat::snapshot_accept()` after the baseline task without reviewing the diff.
- `pkgload::load_all()` + `testthat::test_local()` are a **local development convenience only**. They load source directly and do not exercise installation, `NAMESPACE` resolution, or `Imports` declarations. CI (spec step 7) must still run `R CMD build` followed by `R CMD check --as-cran` against the built tarball, which is the only thing that validates the installed package. Do not let this substitution propagate into the workflow files.

## Reference: exported vs internal

Calling an internal function without `:::` fails. Verified against `NAMESPACE`.

| Exported (call directly) | Internal (call with `TimeMetric:::`) |
|---|---|
| `pam.survival_eval`, `pam.predicted_survial_eval`, `pam.predicted_survial_eval_cr` | `pam.predicted_survial_eval_two_phase` |
| `pam.coxph_restricted`, `pam.surverg_restricted`, `pam.predict_cr` | `pam.rsh_metric`, `pam.rsph_metric`, `pam.Brier_metric` |
| `pam.summary`, `pam.summary_cr`, `pam.sample_design` | `pam.Brier`, `Gt`, `pam.rsph`, `pam.summary.rsph`, `m_cif` |
| `cc_weights`, `ncc_weights` | `weighted_param`, `integrate_survival` |
| `plot_pred`, `summary_pred_plot` | |
| `sim_cox_weibull_censored`, `simulateTwoCauseFineGrayModel` | |

## Reference: verified signatures and defaults

Copied verbatim from source. Later tasks depend on these being exact.

```r
sim_cox_weibull_censored(n, pi_c, v, beta, mu = NULL, sd = NULL, seed = NULL,
                         interact = FALSE, nonlinear = FALSE)
simulateTwoCauseFineGrayModel(n, v, beta1, beta2, lambda1 = 1, X = NULL, mu = 0,
                              p = 0.7, c_scale = 1, censor = 0, sd.time = NULL,
                              mu.c.time = NULL, independent_c = TRUE,
                              report.mu_and_sd = FALSE, seed = 1234)
pam.coxph_restricted(model, covs, tau = 10e10, new_data = NULL, predict = TRUE)
pam.surverg_restricted(model, covs, tau = 10e10, new_data = NULL, predict = TRUE)
pam.predict_cr(model1 = NULL, model2 = NULL, fg_model = NULL, cr_model = NULL,
               tau = NULL, newdata, event.type = 1, covs)
pam.survival_eval(train_data, covariates, models = "coxph", metrics = "all",
                  predicted_data = NULL, t_star = NULL, tau = NULL)
pam.predicted_survial_eval(model, event_time, predicted_probability,
                           pred_mean_survival = NULL, status, covariates,
                           new_data = NULL, metrics = NULL, t_star = NULL, tau = NULL)
pam.predicted_survial_eval_cr(pred_cif, event_time, time.cif, status,
                              metrics = NULL, t_star = NULL, tau = NULL, event_type = 1)
pam.summary(models, metrics = NULL, t_star = NULL, tau = 10e10, digits = 2)
pam.summary_cr(models, metrics = NULL, t_star = NULL, tau = NULL,
               event_type = 1, digits = 2)
pam.sample_design(models, case_weights, km_cens,
                  metrics = c("Pesudo_R", "Harrell\u2019s C", "Uno\u2019s C",
                              "Brier Score", "Time Dependent Auc"),
                  t_star = NULL, tau = NULL, digits = 2)
cc_weights(time, status, subcohort = NULL, strata = NULL)
ncc_weights(time, status, strata = NULL, m = NULL)
pam.rsh_metric(predicted_data, survival_time, status)
pam.rsph_metric(time, status, risk_score)
pam.Brier_metric(predicted_data, suvival_time, t_star = -1)
pam.Brier(object, pre_sp, t_star = -1)
Gt(object, timepoint)
pam.rsph(fit, ...)
pam.summary.rsph(object, times, band = 5, ...)
plot_pred(data, title = NULL, xlab = "Risk Score", ylab = "Days", ...)
summary_pred_plot(data_list, titles = NULL, plot_fun = plot_pred, ncol = 2, ...)
```

**Return structure of the prediction functions** (`R/pam.coxph_restricted.R:89`) — the component is `times`, **not** `time`:

```r
list(pred = <numeric, restricted mean survival time>,
     times = <numeric, observed times>,
     status = <numeric, event indicator>,
     surv_prob = <matrix, subjects x time grid>)
```

**Return structure of the evaluation functions** — a `data.frame`, not a named list. `names(res)` is `c("Metric", "Value")`; metric names live in `res$Metric` (`R/pam.predicted_survial_eval.R:198`).

**Default metric vectors**, verbatim (note the misspelling and the U+2019 apostrophes):

```r
# pam.predicted_survial_eval  (R/pam.predicted_survial_eval.R:83-85)
c("Pseudo_R_square", "Pseudo_R2_point", "Harrell\u2019s C", "Uno\u2019s C",
  "R_E", "R_sh", "Brier Score", "Time Dependent Auc")

# pam.predicted_survial_eval_cr  (R/pam.predicted_survival_eval_cr.R:64-65)
c("Pseudo_R_square", "Pseudo_R2_point", "C_index", "Brier Score", "Time Dependent Auc")

# pam.predicted_survial_eval_two_phase  (R/pam.predicted_survial_eval_two_phase.R:45)
c("Pesudo_R", "Harrell\u2019s C", "Uno\u2019s C", "Brier Score", "Time Dependent Auc")
```

Use `"’"` escapes in test code so the test files stay ASCII-only while still matching the curly apostrophes in the source.

**Critical gotcha:** `sim_cox_weibull_censored` returns different columns depending on `pi_c`. With `pi_c == 0` it returns `time, status` (all 1), `x1..xp`. With `pi_c > 0` it returns `time, status, x1..xp, y_true, cens_time` plus `mu`/`sd_log` attributes. A naive `coxph(Surv(time, status) ~ ., data = df)` on the censored version silently fits `y_true` and `cens_time` as covariates. Fixtures select covariate columns explicitly.

---

### Task 1: Test infrastructure and local dependency install

**Files:**
- Create: `tests/testthat.R`
- Create: `tests/testthat/test-smoke.R`
- Create: `docs/superpowers/findings.md`
- Modify: `DESCRIPTION` (add `Config/testthat/edition: 3`; give `testthat` a version floor)

**Interfaces:**
- Consumes: nothing
- Produces: a package that loads under `devtools::load_all()` and a working `testthat::test_local()` entry point that every later task calls.

Installing dependencies **locally** so the package loads is distinct from **declaring** them in `DESCRIPTION`, which is spec step 2 and out of scope here.

- [ ] **Step 1: Install the toolchain and every runtime dependency**

```bash
Rscript -e 'install.packages(c("roxygen2","testthat","covr","pkgload","rms","expint","survminer","pec","tdROC","yardstick","patchwork","dplyr","magrittr","purrr","tibble","ggplot2","survival"), repos="https://cloud.r-project.org")'
```

`rms`, `pec`, and `survminer` compile and pull large dependency trees; expect 10-20 minutes. One-time cost.

- [ ] **Step 2: Verify every dependency loads**

```bash
Rscript -e 'for (p in c("survival","rms","expint","survminer","pec","tdROC","yardstick","patchwork","dplyr","magrittr","purrr","tibble","ggplot2")) { ok <- requireNamespace(p, quietly=TRUE); cat(sprintf("%-12s %s\n", p, if (ok) "OK" else "MISSING")) }'
```
Expected: every line `OK`. Do not proceed with any `MISSING`.

- [ ] **Step 3: Verify the package loads**

```bash
Rscript -e 'pkgload::load_all("."); cat("loaded OK\n")'
```
Expected: `loaded OK`. Warnings about undeclared imports are expected here and are fixed in spec step 2.

- [ ] **Step 4: Create the testthat entry point**

Create `tests/testthat.R`:
```r
library(testthat)
library(TimeMetric)

test_check("TimeMetric")
```

- [ ] **Step 5: Enable testthat 3rd edition**

In `DESCRIPTION`, set these two fields, leaving every other field untouched:
```
Suggests: testthat (>= 3.0.0)
Config/testthat/edition: 3
```

- [ ] **Step 6: Write a smoke test proving the harness runs**

Create `tests/testthat/test-smoke.R`:
```r
test_that("package exports the functions this suite characterizes", {
  exports <- getNamespaceExports("TimeMetric")

  expect_type(exports, "character")
  expect_true(all(c(
    "pam.survival_eval", "pam.predicted_survial_eval",
    "pam.predicted_survial_eval_cr", "pam.coxph_restricted",
    "pam.surverg_restricted", "pam.predict_cr", "pam.summary",
    "pam.summary_cr", "pam.sample_design", "cc_weights", "ncc_weights",
    "plot_pred", "summary_pred_plot", "sim_cox_weibull_censored",
    "simulateTwoCauseFineGrayModel"
  ) %in% exports))
})

test_that("functions this suite reaches with ::: are genuinely internal", {
  exports <- getNamespaceExports("TimeMetric")

  expect_false(any(c(
    "Gt", "pam.Brier", "pam.rsh_metric", "pam.rsph_metric",
    "pam.Brier_metric", "pam.rsph", "m_cif",
    "pam.predicted_survial_eval_two_phase"
  ) %in% exports))
})
```

- [ ] **Step 7: Run the suite**

Run: `Rscript -e 'testthat::test_local()'`
Expected: PASS, 0 failures.

If any name in step 6's first block is absent from `exports`, stop: the `NAMESPACE` differs from the spec's analysis. Record the difference in `findings.md` before continuing.

- [ ] **Step 8: Create the findings log**

Create `docs/superpowers/findings.md`:
```markdown
# TimeMetric — Findings Log

Behaviour pinned by characterization tests that appears incorrect. Each entry is
resolved deliberately in a later step, never by quietly changing an expected value.

| # | Location | Observed behaviour | Why it looks wrong | Status |
|---|---|---|---|---|
| 1 | `R/pam.rsh_metric.R:48-51` | Sorts `data` by `survival_time`, then computes `Mtx` from the unsorted `predicted_data` argument rather than the reordered column | Predictions are misaligned with the status vector they multiply, so `Dx` and therefore `R_sh` are wrong whenever input is not already time-sorted | Open |
| 2 | `R/pam.rsph.R:79,189,284` | Uses `require(survival)` inside package code | `R CMD check` flags `require()` in package code; imports belong in NAMESPACE | Open |
| 3 | `R/pam.rsh_metric.R` roxygen | Example calls `library(PAmeasure)` | Stale reference to the predecessor package | Open |
| 4 | `R/pam.predicted_survial_eval.R:82-85` | `valid_metrics` and `default_metrics` use U+2019 curly apostrophes in `Harrell's C` / `Uno's C` | Users must type a typographic apostrophe for metric selection to match. Spec section 5 records this for the two-phase function only; it occurs here too | Open |
| 5 | `R/pam.predicted_survial_eval_cr.R:66-76` | Runs an `"all" %in% metrics` / `setdiff` validation block before the `is.null(metrics)` default is applied, then repeats the same block afterwards | Duplicated logic; the first block is a no-op on the `NULL` default | Open |
```

- [ ] **Step 9: Commit**

```bash
git add tests/testthat.R tests/testthat/test-smoke.R DESCRIPTION docs/superpowers/findings.md
git commit -m "test: add testthat infrastructure and findings log"
```

---

### Task 2: Deterministic fixtures and expectation helpers

**Files:**
- Create: `tests/testthat/helper-simdata.R`
- Create: `tests/testthat/helper-expectations.R`
- Test: `tests/testthat/test-fixtures.R`

**Interfaces:**
- Consumes: `sim_cox_weibull_censored`, `simulateTwoCauseFineGrayModel`
- Produces, used by every later task:
  - `fx_surv()` -> `data.frame(time, status, x1, x2)`, 200 rows, ~30% censored
  - `fx_surv_uncensored()` -> `data.frame(time, status, x1, x2)`, 200 rows, all `status == 1`
  - `fx_cr()` -> `data.frame` with `obs.times`, `obs.event`, plus generated covariate columns
  - `fx_cox()` -> `coxph` fit on `fx_surv()`, `x = TRUE, y = TRUE`
  - `fx_survreg()` -> `survreg` fit on `fx_surv()`
  - `fx_covs()` -> `c("x1", "x2")`
  - `snap_num(x, digits = 6)` -> rounded numeric, safe to snapshot
  - `mat_fingerprint(m, digits = 6)` -> compact `list(dim, min, max, mean, n_na)` for matrices
  - `expect_metric_table(res)` -> structural assertions for `Metric`/`Value` data frames

- [ ] **Step 1: Write the data fixtures**

Create `tests/testthat/helper-simdata.R`:
```r
# Deterministic fixtures shared by all characterization tests.
# Seeds are fixed; never change them without regenerating every snapshot.

fx_covs <- function() c("x1", "x2")

fx_surv <- function() {
  d <- sim_cox_weibull_censored(
    n = 200, pi_c = 0.3, v = 2, beta = c(0.5, -0.5), seed = 1001
  )
  # pi_c > 0 also returns y_true and cens_time; drop them so model
  # formulas built with "." cannot pick them up as covariates.
  d[, c("time", "status", "x1", "x2")]
}

fx_surv_uncensored <- function() {
  sim_cox_weibull_censored(
    n = 200, pi_c = 0, v = 2, beta = c(0.5, -0.5), seed = 1002
  )[, c("time", "status", "x1", "x2")]
}

fx_cr <- function() {
  simulateTwoCauseFineGrayModel(
    n = 200, v = 2, beta1 = c(0.5, -0.5), beta2 = c(-0.3, 0.3),
    censor = 0.3, seed = 1003
  )
}

fx_cox <- function() {
  survival::coxph(
    survival::Surv(time, status) ~ x1 + x2,
    data = fx_surv(), x = TRUE, y = TRUE
  )
}

fx_survreg <- function() {
  survival::survreg(
    survival::Surv(time, status) ~ x1 + x2,
    data = fx_surv(), dist = "weibull"
  )
}
```

- [ ] **Step 2: Write the expectation helpers**

Create `tests/testthat/helper-expectations.R`:
```r
# Reduce values to something small and platform-stable before snapshotting.
# Never snapshot a fitted model or a full prediction matrix.

snap_num <- function(x, digits = 6) round(as.numeric(x), digits)

mat_fingerprint <- function(m, digits = 6) {
  m <- as.matrix(m)
  list(
    dim  = dim(m),
    min  = round(min(m, na.rm = TRUE), digits),
    max  = round(max(m, na.rm = TRUE), digits),
    mean = round(mean(m, na.rm = TRUE), digits),
    n_na = sum(is.na(m))
  )
}

# Structural assertions common to every Metric/Value result table.
expect_metric_table <- function(res) {
  testthat::expect_true(is.data.frame(res))
  testthat::expect_true("Metric" %in% names(res))
  testthat::expect_gt(nrow(res), 0)
  testthat::expect_type(res$Metric, "character")
  testthat::expect_false(any(duplicated(res$Metric)))
  invisible(res)
}
```

- [ ] **Step 3: Write tests pinning the fixtures**

Create `tests/testthat/test-fixtures.R`:
```r
test_that("fx_surv is deterministic and correctly shaped", {
  d1 <- fx_surv()
  d2 <- fx_surv()

  expect_s3_class(d1, "data.frame")
  expect_identical(names(d1), c("time", "status", "x1", "x2"))
  expect_identical(nrow(d1), 200L)
  expect_identical(d1, d2)
  expect_true(all(d1$status %in% c(0, 1)))
  expect_true(all(d1$time > 0))
  expect_gt(sum(d1$status == 0), 0)
})

test_that("fx_surv_uncensored has no censoring", {
  d <- fx_surv_uncensored()

  expect_identical(names(d), c("time", "status", "x1", "x2"))
  expect_identical(nrow(d), 200L)
  expect_true(all(d$status == 1))
})

test_that("fx_cr produces competing-risk codes", {
  d <- fx_cr()

  expect_s3_class(d, "data.frame")
  expect_identical(nrow(d), 200L)
  expect_true(all(c("obs.times", "obs.event") %in% names(d)))
  expect_true(all(d$obs.event %in% c(0, 1, 2)))
  expect_gt(sum(d$obs.event == 2), 0)
  expect_true(all(d$obs.times > 0))
})

test_that("fixture models fit and are reproducible", {
  m <- fx_cox()
  expect_s3_class(m, "coxph")
  expect_identical(names(coef(m)), c("x1", "x2"))
  expect_equal(coef(m), coef(fx_cox()), tolerance = 1e-12)
  expect_false(is.null(m$x))
  expect_false(is.null(m$y))

  s <- fx_survreg()
  expect_s3_class(s, "survreg")
  expect_identical(names(coef(s)), c("(Intercept)", "x1", "x2"))
})

test_that("fixture numbers are snapshot-stable", {
  d <- fx_surv()

  expect_snapshot_value(snap_num(head(d$time, 10)), style = "serialize")
  expect_snapshot_value(sum(d$status), style = "serialize")
  expect_snapshot_value(snap_num(unname(coef(fx_cox()))), style = "serialize")
  expect_snapshot_value(snap_num(unname(coef(fx_survreg()))), style = "serialize")
  expect_snapshot_value(table(fx_cr()$obs.event), style = "serialize")
})
```

- [ ] **Step 4: Run and record snapshots**

Run: `Rscript -e 'testthat::test_local(filter = "fixtures")'`

Expected: PASS. This first run **creates** `tests/testthat/_snaps/fixtures.md`.

If `fx_cr()` errors, read `R/simulateTwoCauseFineGrayModel.R` to find how many covariates it generates when `X = NULL`, adjust `beta1`/`beta2` to that length, and record the signature surprise in `findings.md`. If `obs.event` contains values outside `{0,1,2}`, widen the assertion to the observed set and record that too.

- [ ] **Step 5: Run again to confirm snapshots compare rather than re-record**

Run: `Rscript -e 'testthat::test_local(filter = "fixtures")'`
Expected: PASS with no "adding new snapshot" messages.

- [ ] **Step 6: Commit**

```bash
git add tests/testthat/helper-simdata.R tests/testthat/helper-expectations.R tests/testthat/test-fixtures.R tests/testthat/_snaps/fixtures.md
git commit -m "test: add deterministic fixtures and expectation helpers"
```

---

### Task 3: Characterize Cluster C metric functions

Highest priority: these are the spec's deletion candidates and have zero coverage. All three are internal.

**Files:**
- Test: `tests/testthat/test-cluster-c-metrics.R`

**Interfaces:**
- Consumes: `fx_surv()`, `snap_num()`
- Produces: the regression suite gating deletion of `pam.rsh_metric`, `pam.rsph_metric`, `pam.Brier_metric` in spec step 5

- [ ] **Step 1: Write the characterization tests**

Create `tests/testthat/test-cluster-c-metrics.R`:
```r
test_that("pam.rsh_metric returns D, Dx, R_sh rounded to 4 places", {
  d <- fx_surv()
  pred <- seq(0.9, 0.1, length.out = nrow(d))

  res <- TimeMetric:::pam.rsh_metric(pred, d$time, d$status)

  expect_type(res, "list")
  expect_identical(names(res), c("D", "Dx", "R_sh"))
  expect_length(res$D, 1L)
  expect_length(res$Dx, 1L)
  expect_length(res$R_sh, 1L)
  expect_true(is.finite(res$D))
  expect_true(is.finite(res$Dx))
  expect_gt(res$D, 0)
  # the function rounds its outputs; confirm that is still true
  expect_identical(res$D, round(res$D, 4))
  expect_identical(res$Dx, round(res$Dx, 4))

  expect_snapshot_value(snap_num(c(res$D, res$Dx, res$R_sh)), style = "serialize")
})

test_that("pam.rsh_metric is sensitive to input order (pins FINDING 1)", {
  # The function sorts its internal data frame by survival_time but multiplies
  # by the UNSORTED predicted_data argument. Permuting all three inputs
  # consistently should leave the result unchanged. It does not. Pinned as-is.
  d <- fx_surv()
  pred <- seq(0.9, 0.1, length.out = nrow(d))
  ord <- order(d$time)

  sorted   <- TimeMetric:::pam.rsh_metric(pred[ord], d$time[ord], d$status[ord])
  shuffled <- TimeMetric:::pam.rsh_metric(pred, d$time, d$status)

  # D depends only on time and status, so it is invariant
  expect_equal(sorted$D, shuffled$D, tolerance = 1e-6)
  # Dx depends on the misaligned predictions, so it is not
  expect_false(isTRUE(all.equal(sorted$Dx, shuffled$Dx, tolerance = 1e-6)))

  expect_snapshot_value(
    snap_num(c(sorted$Dx, shuffled$Dx)), style = "serialize"
  )
})

test_that("pam.rsph_metric takes (time, status, risk_score) and returns a list", {
  d <- fx_surv()
  risk <- as.numeric(predict(fx_cox(), newdata = d, type = "lp"))

  res <- TimeMetric:::pam.rsph_metric(d$time, d$status, risk)

  expect_type(res, "list")
  expect_gt(length(res), 0)
  expect_snapshot_value(sort(names(res)), style = "serialize")
  expect_snapshot_value(
    snap_num(unlist(res[vapply(res, function(z) is.numeric(z) && length(z) == 1L,
                               logical(1))])),
    style = "serialize"
  )
})

test_that("pam.rsph_metric rejects mismatched input lengths", {
  d <- fx_surv()
  risk <- as.numeric(predict(fx_cox(), newdata = d, type = "lp"))

  # guarded by stopifnot(length(time) == length(status), ...)
  expect_error(
    TimeMetric:::pam.rsph_metric(d$time[-1], d$status, risk),
    "length"
  )
  expect_error(
    TimeMetric:::pam.rsph_metric(d$time, d$status, risk[-1]),
    "length"
  )
})

test_that("pam.Brier_metric validates its inputs with specific messages", {
  d <- fx_surv()
  pred <- rep(0.5, nrow(d))
  sv <- survival::Surv(d$time, d$status)

  expect_error(
    TimeMetric:::pam.Brier_metric(pred, d$time),
    "must be a survival object created using Surv"
  )
  expect_error(
    TimeMetric:::pam.Brier_metric(pred[-1], sv),
    "Length of suvival_time and predicted_data must match"
  )
  expect_error(
    TimeMetric:::pam.Brier_metric(c(NA_real_, pred[-1]), sv),
    "cannot have NA"
  )
})

test_that("pam.Brier_metric returns a finite scalar in [0, 1]", {
  d <- fx_surv()
  sv <- survival::Surv(d$time, d$status)
  pred <- seq(0.9, 0.1, length.out = nrow(d))

  res <- TimeMetric:::pam.Brier_metric(pred, sv)

  expect_length(res, 1L)
  expect_true(is.finite(res))
  expect_gte(res, 0)
  expect_lte(res, 1)
  expect_snapshot_value(snap_num(res), style = "serialize")
})

test_that("pam.Brier_metric default t_star is the median observed time", {
  d <- fx_surv()
  sv <- survival::Surv(d$time, d$status)
  pred <- seq(0.9, 0.1, length.out = nrow(d))

  # source: if (t_star < 0) t_star <- median(time), where time is the
  # sorted observed time from the Surv object
  expect_equal(
    TimeMetric:::pam.Brier_metric(pred, sv),
    TimeMetric:::pam.Brier_metric(pred, sv, stats::median(d$time)),
    tolerance = 1e-6
  )
})
```

- [ ] **Step 2: Run and record snapshots**

Run: `Rscript -e 'testthat::test_local(filter = "cluster-c-metrics")'`
Expected: PASS, snapshots created.

If `pam.rsph_metric` returns no scalar numeric components, the `vapply` filter yields an empty vector; replace that snapshot with `expect_snapshot_value(snap_num(res$Re), style = "serialize")` using whatever component the names snapshot revealed, and record the real structure in `findings.md`.

- [ ] **Step 3: Run again to confirm stability**

Run: `Rscript -e 'testthat::test_local(filter = "cluster-c-metrics")'`
Expected: PASS, no snapshots added.

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/test-cluster-c-metrics.R tests/testthat/_snaps/cluster-c-metrics.md
git commit -m "test: characterize Cluster C metric functions before deletion decision"
```

---

### Task 4: Characterize the pec dependency chain

Pins `Gt` -> `pam.Brier` -> `predictSurvProb`. This chain decides whether `pec` can leave `Imports`, so lock it before any rewrite. `Gt` and `pam.Brier` are both internal.

**Files:**
- Test: `tests/testthat/test-brier-pec-chain.R`

**Interfaces:**
- Consumes: `fx_surv()`, `fx_cox()`, `snap_num()`, `mat_fingerprint()`
- Produces: the gate for the `pec` demotion decision in spec step 2

- [ ] **Step 1: Write the characterization tests**

Create `tests/testthat/test-brier-pec-chain.R`:
```r
test_that("Gt returns a censoring survival probability in [0, 1]", {
  d <- fx_surv()
  sv <- survival::Surv(d$time, d$status)

  res <- TimeMetric:::Gt(sv, stats::median(d$time))

  expect_length(res, 1L)
  expect_true(is.finite(res))
  expect_gte(res, 0)
  expect_lte(res, 1)
  expect_snapshot_value(snap_num(res), style = "serialize")
})

test_that("Gt is non-increasing in the timepoint", {
  d <- fx_surv()
  sv <- survival::Surv(d$time, d$status)
  qs <- stats::quantile(d$time, c(0.2, 0.4, 0.6, 0.8))

  vals <- vapply(qs, function(tp) TimeMetric:::Gt(sv, tp), numeric(1))

  expect_true(all(diff(vals) <= 1e-8))
  expect_snapshot_value(snap_num(vals), style = "serialize")
})

test_that("Gt rejects non-Surv input", {
  expect_error(TimeMetric:::Gt(1:10, 5), "not of class Surv")
})

test_that("pam.Brier on a coxph fit returns a finite scalar in [0, 1]", {
  # This is the chain that reaches pec::predictSurvProb unqualified.
  d <- fx_surv()
  fit <- fx_cox()

  res <- TimeMetric:::pam.Brier(fit, d, stats::median(d$time))

  expect_length(res, 1L)
  expect_true(is.finite(res))
  expect_gte(res, 0)
  expect_lte(res, 1)
  expect_snapshot_value(snap_num(res), style = "serialize")
})

test_that("pam.Brier default t_star is the median of TRAINING event times", {
  # source R/pam.Brier.R:59-63 --
  #   distime <- sort(unique(as.vector(obj$y[obj$y[, 2] == 1])))
  #   if (t_star0 <= 0) t_star0 <- median(distime)
  # This is the median of the model's own event times, NOT median(d$time).
  d <- fx_surv()
  fit <- fx_cox()

  distime <- sort(unique(as.vector(fit$y[fit$y[, 2] == 1])))
  expected_t <- stats::median(distime)

  expect_equal(
    TimeMetric:::pam.Brier(fit, d),
    TimeMetric:::pam.Brier(fit, d, expected_t),
    tolerance = 1e-6
  )

  # and confirm this genuinely differs from the naive test-data median,
  # so the distinction is protected against a future "simplification"
  expect_false(isTRUE(all.equal(expected_t, stats::median(d$time),
                                tolerance = 1e-6)))
})

test_that("pec::predictSurvProb drives pam.Brier and is reachable", {
  skip_if_not_installed("pec")
  d <- fx_surv()
  fit <- fx_cox()
  t_star <- stats::median(d$time)

  probs <- pec::predictSurvProb(fit, d, t_star)

  expect_identical(NROW(probs), 200L)
  expect_true(all(probs >= 0 & probs <= 1))
  # fingerprint, not the full 200-element matrix
  expect_snapshot_value(mat_fingerprint(probs), style = "serialize")
})
```

- [ ] **Step 2: Run and record snapshots**

Run: `Rscript -e 'testthat::test_local(filter = "brier-pec-chain")'`
Expected: PASS, snapshots created.

If the "differs from the naive median" assertion fails, the two medians coincide for this fixture. Change `fx_surv()`'s seed in that test only by constructing a local dataset with heavier censoring, so the distinction stays meaningful; record that you did so.

- [ ] **Step 3: Run again to confirm stability**

Run: `Rscript -e 'testthat::test_local(filter = "brier-pec-chain")'`
Expected: PASS, no snapshots added.

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/test-brier-pec-chain.R tests/testthat/_snaps/brier-pec-chain.md
git commit -m "test: characterize Gt/pam.Brier/predictSurvProb chain gating pec removal"
```

---

### Task 5: Characterize the prediction module

Covers all three exported prediction functions, including `pam.predict_cr`.

**Files:**
- Test: `tests/testthat/test-prediction-module.R`

**Interfaces:**
- Consumes: `fx_surv()`, `fx_cr()`, `fx_cox()`, `fx_survreg()`, `fx_covs()`, `snap_num()`, `mat_fingerprint()`
- Produces: pinned `list(pred, times, status, surv_prob)` structure that every evaluation function consumes

- [ ] **Step 1: Write the characterization tests**

Create `tests/testthat/test-prediction-module.R`:
```r
test_that("pam.coxph_restricted returns pred, times, status, surv_prob", {
  d <- fx_surv()

  res <- pam.coxph_restricted(model = fx_cox(), covs = fx_covs(),
                              new_data = d, tau = 10e10)

  expect_type(res, "list")
  # component names the evaluation functions depend on; note "times", not "time"
  expect_true(all(c("pred", "times", "status", "surv_prob") %in% names(res)))
  expect_identical(length(res$times), 200L)
  expect_identical(length(res$status), 200L)
  expect_identical(length(res$pred), 200L)
  expect_true(all(res$status %in% c(0, 1)))
  expect_true(all(is.finite(res$times)))
  expect_identical(NROW(res$surv_prob), 200L)
  expect_true(all(res$surv_prob >= 0 & res$surv_prob <= 1, na.rm = TRUE))

  expect_snapshot_value(sort(names(res)), style = "serialize")
  expect_snapshot_value(snap_num(head(res$pred, 10)), style = "serialize")
  expect_snapshot_value(mat_fingerprint(res$surv_prob), style = "serialize")
})

test_that("pam.coxph_restricted requires time and status in new_data", {
  d <- fx_surv()

  expect_error(
    pam.coxph_restricted(model = fx_cox(), covs = fx_covs(),
                         new_data = d[, c("x1", "x2")], tau = 10e10),
    "new_data require time and status columns"
  )
})

test_that("pam.surverg_restricted returns the same component structure", {
  d <- fx_surv()

  res <- pam.surverg_restricted(model = fx_survreg(), covs = fx_covs(),
                                new_data = d, tau = 10e10)

  expect_type(res, "list")
  expect_true(all(c("pred", "times", "status", "surv_prob") %in% names(res)))
  expect_identical(length(res$times), 200L)
  expect_identical(NROW(res$surv_prob), 200L)

  expect_snapshot_value(sort(names(res)), style = "serialize")
  expect_snapshot_value(snap_num(head(res$pred, 10)), style = "serialize")
  expect_snapshot_value(mat_fingerprint(res$surv_prob), style = "serialize")
})

test_that("pam.surverg_restricted validates model class and covariates", {
  d <- fx_surv()

  expect_error(
    pam.surverg_restricted(model = fx_cox(), covs = fx_covs(), new_data = d),
    "model must be an object of class 'survreg'"
  )
  expect_error(
    pam.surverg_restricted(model = fx_survreg(), covs = c("x1", "nope"),
                           new_data = d),
    "All covariates must be present in new_data"
  )
  expect_error(
    pam.surverg_restricted(model = fx_survreg(), covs = fx_covs(),
                           new_data = d[, c("x1", "x2")]),
    "new_data require time and status columns"
  )
})

test_that("pam.predict_cr returns CIF predictions for a competing-risks fit", {
  d <- fx_cr()
  covs <- setdiff(names(d), c("obs.times", "obs.event"))
  dd <- d
  dd$time <- dd$obs.times
  dd$status <- dd$obs.event

  fg <- survival::coxph(
    survival::Surv(time, status == 1) ~ .,
    data = dd[, c("time", "status", covs)], x = TRUE, y = TRUE
  )

  res <- pam.predict_cr(model1 = fg, newdata = dd, covs = covs,
                        event.type = 1, tau = max(dd$time))

  expect_true(is.list(res) || is.matrix(res))
  expect_snapshot_value(
    if (is.list(res)) sort(names(res)) else dim(res),
    style = "serialize"
  )
})

test_that("pam.predict_cr rejects an unrecognised model type", {
  d <- fx_cr()
  covs <- setdiff(names(d), c("obs.times", "obs.event"))
  dd <- d
  dd$time <- dd$obs.times
  dd$status <- dd$obs.event

  expect_error(
    pam.predict_cr(newdata = dd, covs = covs, event.type = 1),
    "Unknown model type"
  )
})
```

- [ ] **Step 2: Run and record snapshots**

Run: `Rscript -e 'testthat::test_local(filter = "prediction-module")'`
Expected: PASS, snapshots created.

`pam.predict_cr` accepts four alternative model arguments (`model1`, `model2`, `fg_model`, `cr_model`). If passing a `coxph` fit as `model1` hits "Unknown model type", read `R/pam.predict_cr.R:60-105` to see which classes each argument dispatches on, use the matching one, and record the accepted classes in `findings.md` — the roxygen does not state them.

- [ ] **Step 3: Run again to confirm stability**

Run: `Rscript -e 'testthat::test_local(filter = "prediction-module")'`
Expected: PASS, no snapshots added.

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/test-prediction-module.R tests/testthat/_snaps/prediction-module.md
git commit -m "test: characterize coxph, survreg and competing-risks prediction functions"
```

---

### Task 6: Characterize the right-censored evaluation path

Covers `pam.predicted_survial_eval`, `pam.summary`, and the exported `pam.survival_eval`.

**Files:**
- Test: `tests/testthat/test-eval-survival.R`

**Interfaces:**
- Consumes: `fx_surv()`, `fx_cox()`, `fx_covs()`, `snap_num()`, `expect_metric_table()`; the `times`/`surv_prob` structure from Task 5
- Produces: pinned metric values and the exact metric-name strings spec step 6 will standardize

- [ ] **Step 1: Write the characterization tests**

Create `tests/testthat/test-eval-survival.R`:
```r
test_that("pam.predicted_survial_eval returns a Metric/Value data frame", {
  d <- fx_surv()
  pred <- pam.coxph_restricted(model = fx_cox(), covs = fx_covs(),
                               new_data = d, tau = 10e10)

  res <- pam.predicted_survial_eval(
    model = fx_cox(),
    event_time = pred$times,
    predicted_probability = pred$surv_prob,
    status = pred$status,
    covariates = fx_covs(),
    new_data = d,
    tau = 10e10
  )

  expect_metric_table(res)
  expect_true("Value" %in% names(res))
  expect_snapshot_value(res$Metric, style = "serialize")
  expect_snapshot_value(snap_num(res$Value), style = "serialize")
})

test_that("the default metric set includes R_sh and R_E in the Metric column", {
  # Metric names live in res$Metric; names(res) is c("Metric", "Value").
  # R_sh reaches rms::cph + pam.schemper; R_E reaches pam.rsph dispatch.
  # Asserting both here proves those two code paths execute by default.
  d <- fx_surv()
  pred <- pam.coxph_restricted(model = fx_cox(), covs = fx_covs(),
                               new_data = d, tau = 10e10)

  res <- pam.predicted_survial_eval(
    model = fx_cox(),
    event_time = pred$times,
    predicted_probability = pred$surv_prob,
    status = pred$status,
    covariates = fx_covs(),
    new_data = d,
    tau = 10e10
  )

  expect_true("R_sh" %in% res$Metric)
  expect_true("R_E" %in% res$Metric)
  expect_true("Brier Score" %in% res$Metric)
  # curly apostrophes, written as an escape so this file stays ASCII
  expect_true("Harrell\u2019s C" %in% res$Metric)
  expect_true("Uno\u2019s C" %in% res$Metric)
})

test_that("pam.predicted_survial_eval rejects an unknown metric name", {
  d <- fx_surv()
  pred <- pam.coxph_restricted(model = fx_cox(), covs = fx_covs(),
                               new_data = d, tau = 10e10)

  expect_error(
    pam.predicted_survial_eval(
      model = fx_cox(),
      event_time = pred$times,
      predicted_probability = pred$surv_prob,
      status = pred$status,
      covariates = fx_covs(),
      new_data = d,
      tau = 10e10,
      metrics = "Harrells_C"        # ASCII spelling is NOT in valid_metrics
    ),
    "Invalid metrics"
  )
})

test_that("pam.summary pivots one model into a Metric column table", {
  d <- fx_surv()
  pred <- pam.coxph_restricted(model = fx_cox(), covs = fx_covs(),
                               new_data = d, tau = 10e10)

  res <- pam.summary(list(value = pred), tau = 10e10)

  expect_metric_table(res)
  expect_true("value" %in% names(res))
  expect_snapshot_value(res$Metric, style = "serialize")
  expect_snapshot_value(snap_num(res$value), style = "serialize")
})

test_that("pam.summary puts one column per model", {
  d <- fx_surv()
  p_cox <- pam.coxph_restricted(model = fx_cox(), covs = fx_covs(),
                                new_data = d, tau = 10e10)
  p_reg <- pam.surverg_restricted(model = fx_survreg(), covs = fx_covs(),
                                  new_data = d, tau = 10e10)

  res <- pam.summary(list(cox = p_cox, weibull = p_reg), tau = 10e10)

  expect_metric_table(res)
  expect_identical(names(res), c("Metric", "cox", "weibull"))
  expect_snapshot_value(snap_num(res$cox), style = "serialize")
  expect_snapshot_value(snap_num(res$weibull), style = "serialize")
})

test_that("pam.summary rejects a non-list or empty models argument", {
  expect_error(pam.summary(list()), "must be a non-empty named list")
  expect_error(pam.summary("not a list"), "must be a non-empty named list")
})

test_that("pam.survival_eval fits and evaluates from raw training data", {
  d <- fx_surv()

  res <- pam.survival_eval(
    train_data = d,
    covariates = fx_covs(),
    models = "coxph",
    metrics = "all"
  )

  expect_true(is.data.frame(res) || is.list(res))
  expect_snapshot_value(
    if (is.data.frame(res)) names(res) else sort(names(res)),
    style = "serialize"
  )
})

test_that("pam.survival_eval requires train_data and covariates", {
  expect_error(
    pam.survival_eval(),
    "Please provide 'train_data', 'time_var', 'status_var', and 'covariates'"
  )
})
```

- [ ] **Step 2: Run and record snapshots**

Run: `Rscript -e 'testthat::test_local(filter = "eval-survival")'`
Expected: PASS, snapshots created.

If `R_sh` comes back `NA` with a message about Cox-only support or factor variables, that is current behaviour — keep the assertion that the metric is *present* and record the `NA` condition in `findings.md`.

If `pam.summary`'s per-model column is not named after the list element, read `R/pam.predicted_survial_eval.R:285-296` for the `reshape` call and correct the expected `names(res)`.

- [ ] **Step 3: Run again to confirm stability**

Run: `Rscript -e 'testthat::test_local(filter = "eval-survival")'`
Expected: PASS, no snapshots added.

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/test-eval-survival.R tests/testthat/_snaps/eval-survival.md
git commit -m "test: characterize right-censored evaluation, summary and metric names"
```

---

### Task 7: Characterize the competing-risks path

Covers `pam.predicted_survial_eval_cr`, the exported `pam.summary_cr`, and the `m_cif` helper.

**Files:**
- Test: `tests/testthat/test-eval-competing-risks.R`

**Interfaces:**
- Consumes: `fx_cr()`, `snap_num()`, `expect_metric_table()`
- Produces: pinned competing-risks metric names, which differ from the survival set (`C_index` rather than Harrell/Uno)

- [ ] **Step 1: Write the characterization tests**

Create `tests/testthat/test-eval-competing-risks.R`:
```r
test_that("m_cif reduces a CIF column to a scalar risk", {
  cif <- seq(0, 0.6, length.out = 50)
  times <- seq(1, 50, length.out = 50)

  res <- TimeMetric:::m_cif(cif, time.cif = times, tau = 50)

  expect_length(res, 1L)
  expect_true(is.finite(res))
  expect_snapshot_value(snap_num(res), style = "serialize")
})

test_that("m_cif is monotone in the CIF magnitude", {
  times <- seq(1, 50, length.out = 50)
  low  <- TimeMetric:::m_cif(seq(0, 0.3, length.out = 50), time.cif = times, tau = 50)
  high <- TimeMetric:::m_cif(seq(0, 0.9, length.out = 50), time.cif = times, tau = 50)

  expect_false(isTRUE(all.equal(low, high, tolerance = 1e-6)))
  expect_snapshot_value(snap_num(c(low, high)), style = "serialize")
})

test_that("pam.predicted_survial_eval_cr returns a Metric/Value table", {
  d <- fx_cr()
  n <- nrow(d)
  time_grid <- sort(unique(d$obs.times))
  k <- length(time_grid)
  # deterministic synthetic CIF, monotone increasing down each subject column
  cif <- outer(
    seq_len(k), seq_len(n),
    function(i, j) pmin(0.99, (i / k) * (0.2 + 0.6 * j / n))
  )

  res <- pam.predicted_survial_eval_cr(
    pred_cif = cif,
    event_time = d$obs.times,
    time.cif = time_grid,
    status = d$obs.event,
    event_type = 1
  )

  expect_metric_table(res)
  expect_snapshot_value(res$Metric, style = "serialize")
  expect_snapshot_value(snap_num(res$Value), style = "serialize")
})

test_that("competing-risks default metrics use C_index, not Harrell/Uno", {
  # default_metrics is c("Pseudo_R_square", "Pseudo_R2_point", "C_index",
  #                      "Brier Score", "Time Dependent Auc")
  # This differs from the right-censored set and must survive standardization.
  d <- fx_cr()
  n <- nrow(d)
  time_grid <- sort(unique(d$obs.times))
  k <- length(time_grid)
  cif <- outer(seq_len(k), seq_len(n),
               function(i, j) pmin(0.99, (i / k) * (0.2 + 0.6 * j / n)))

  res <- pam.predicted_survial_eval_cr(
    pred_cif = cif, event_time = d$obs.times, time.cif = time_grid,
    status = d$obs.event, event_type = 1
  )

  expect_true("C_index" %in% res$Metric)
  expect_false("R_sh" %in% res$Metric)
  expect_false("R_E" %in% res$Metric)
})

test_that("pam.predicted_survial_eval_cr rejects an unknown metric name", {
  d <- fx_cr()
  n <- nrow(d)
  time_grid <- sort(unique(d$obs.times))
  k <- length(time_grid)
  cif <- outer(seq_len(k), seq_len(n),
               function(i, j) pmin(0.99, (i / k) * (0.2 + 0.6 * j / n)))

  expect_error(
    pam.predicted_survial_eval_cr(
      pred_cif = cif, event_time = d$obs.times, time.cif = time_grid,
      status = d$obs.event, event_type = 1, metrics = "R_sh"
    ),
    "Invalid metrics"
  )
})

test_that("pam.summary_cr requires event_time and pred_cif", {
  expect_error(
    pam.summary_cr(list()),
    "must be a non-empty|Please provide"
  )
})

test_that("pam.summary_cr pivots competing-risks results into a Metric table", {
  d <- fx_cr()
  n <- nrow(d)
  time_grid <- sort(unique(d$obs.times))
  k <- length(time_grid)
  cif <- outer(seq_len(k), seq_len(n),
               function(i, j) pmin(0.99, (i / k) * (0.2 + 0.6 * j / n)))

  model_entry <- list(
    pred_cif = cif,
    event_time = d$obs.times,
    time.cif = time_grid,
    status = d$obs.event
  )

  res <- pam.summary_cr(list(fg = model_entry), event_type = 1)

  expect_metric_table(res)
  expect_snapshot_value(res$Metric, style = "serialize")
})

test_that("competing-risks status codes include a third level", {
  # Documents that CR functions accept {0,1,2}, unlike single-event functions.
  # Spec step 6 makes validation function-specific on this basis.
  d <- fx_cr()
  expect_setequal(sort(unique(d$obs.event)), c(0, 1, 2))
})
```

- [ ] **Step 2: Run and record snapshots**

Run: `Rscript -e 'testthat::test_local(filter = "eval-competing-risks")'`
Expected: PASS, snapshots created.

Two shapes to verify against the source rather than assume:

1. `pred_cif` orientation. This plan builds it as `times x subjects`. If `pam.predicted_survial_eval_cr` expects the transpose, `t(cif)` and record the required orientation in `findings.md` — the roxygen does not state it.
2. `pam.summary_cr`'s per-model list element. Read `R/pam.predicted_survival_eval_cr.R:440-480` for the fields it requires from each entry and match them; adjust `model_entry` accordingly.

- [ ] **Step 3: Run again to confirm stability**

Run: `Rscript -e 'testthat::test_local(filter = "eval-competing-risks")'`
Expected: PASS, no snapshots added.

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/test-eval-competing-risks.R tests/testthat/_snaps/eval-competing-risks.md
git commit -m "test: characterize competing-risks evaluation and summary"
```

---

### Task 8: Characterize the two-phase designs

The case-cohort and NCC paths the paper advertises but does not export, plus the exported `pam.sample_design` wrapper. Zero coverage today.

**Files:**
- Test: `tests/testthat/test-two-phase.R`

**Interfaces:**
- Consumes: `fx_surv()`, `fx_cox()`, `fx_covs()`, `snap_num()`, `expect_metric_table()`
- Produces: the regression suite gating export of `tm_evaluate_two_phase` in spec step 3

- [ ] **Step 1: Write the characterization tests**

Create `tests/testthat/test-two-phase.R`:
```r
# Shared local builder: case weights plus a censoring KM fit.
tp_inputs <- function(design = c("cc", "ncc")) {
  design <- match.arg(design)
  d <- fx_surv()
  pred <- pam.coxph_restricted(model = fx_cox(), covs = fx_covs(),
                               new_data = d, tau = 10e10)
  set.seed(2001)
  w <- if (design == "cc") {
    cc_weights(time = d$time, status = d$status,
               subcohort = stats::rbinom(nrow(d), 1, 0.4))
  } else {
    ncc_weights(time = d$time, status = d$status, m = 2)
  }
  list(
    d = d,
    pred = pred,
    weights = if (is.list(w)) w[[1]] else w,
    km_cens = survival::survfit(survival::Surv(d$time, 1 - d$status) ~ 1)
  )
}

test_that("cc_weights returns one finite non-negative weight per subject", {
  d <- fx_surv()
  set.seed(2001)
  subcohort <- stats::rbinom(nrow(d), 1, 0.4)

  w <- cc_weights(time = d$time, status = d$status, subcohort = subcohort)
  wv <- if (is.list(w)) w[[1]] else w

  expect_true(is.numeric(wv))
  expect_identical(length(wv), 200L)
  expect_true(all(is.finite(wv)))
  expect_true(all(wv >= 0))
  expect_snapshot_value(snap_num(head(wv, 10)), style = "serialize")
  expect_snapshot_value(snap_num(sum(wv)), style = "serialize")
})

test_that("ncc_weights requires m and returns one weight per subject", {
  d <- fx_surv()

  expect_error(
    ncc_weights(time = d$time, status = d$status),
    "`m` must be provided for NCC weights"
  )

  w <- ncc_weights(time = d$time, status = d$status, m = 2)
  wv <- if (is.list(w)) w[[1]] else w

  expect_true(is.numeric(wv))
  expect_identical(length(wv), 200L)
  expect_true(all(is.finite(wv)))
  expect_snapshot_value(snap_num(head(wv, 10)), style = "serialize")
  expect_snapshot_value(snap_num(sum(wv)), style = "serialize")
})

test_that("two-phase evaluation works for a case-cohort design", {
  inp <- tp_inputs("cc")

  res <- TimeMetric:::pam.predicted_survial_eval_two_phase(
    pred_results = inp$pred,
    km_cens_fit = inp$km_cens,
    case_weights = inp$weights
  )

  expect_true(inherits(res, "data.frame"))
  expect_true(all(c("Metric", "Value") %in% names(res)))
  expect_gt(nrow(res), 0)
  expect_true(all(is.finite(res$Value) | is.na(res$Value)))
  expect_snapshot_value(res$Metric, style = "serialize")
  expect_snapshot_value(snap_num(res$Value), style = "serialize")
})

test_that("two-phase evaluation works for a nested case-control design", {
  inp <- tp_inputs("ncc")

  res <- TimeMetric:::pam.predicted_survial_eval_two_phase(
    pred_results = inp$pred,
    km_cens_fit = inp$km_cens,
    case_weights = inp$weights
  )

  expect_true(all(c("Metric", "Value") %in% names(res)))
  expect_gt(nrow(res), 0)
  expect_snapshot_value(res$Metric, style = "serialize")
  expect_snapshot_value(snap_num(res$Value), style = "serialize")
})

test_that("case-cohort and NCC weighting give different metric values", {
  # If these coincide the weights are not reaching the estimator, which would
  # make the whole two-phase feature a no-op.
  cc  <- TimeMetric:::pam.predicted_survial_eval_two_phase(
    pred_results = tp_inputs("cc")$pred,
    km_cens_fit  = tp_inputs("cc")$km_cens,
    case_weights = tp_inputs("cc")$weights
  )
  ncc <- TimeMetric:::pam.predicted_survial_eval_two_phase(
    pred_results = tp_inputs("ncc")$pred,
    km_cens_fit  = tp_inputs("ncc")$km_cens,
    case_weights = tp_inputs("ncc")$weights
  )

  expect_identical(cc$Metric, ncc$Metric)
  expect_false(isTRUE(all.equal(cc$Value, ncc$Value, tolerance = 1e-6)))
})

test_that("pam.sample_design validates its models argument", {
  inp <- tp_inputs("cc")

  expect_error(
    pam.sample_design(models = list(), case_weights = inp$weights,
                      km_cens = inp$km_cens),
    "must be a non-empty named list"
  )
})

test_that("pam.sample_design summarises a two-phase design", {
  inp <- tp_inputs("cc")

  res <- pam.sample_design(
    models = list(cc = inp$pred),
    case_weights = inp$weights,
    km_cens = inp$km_cens
  )

  expect_metric_table(res)
  expect_snapshot_value(res$Metric, style = "serialize")
})

test_that("two-phase default metric names are pinned with current spelling", {
  # c("Pesudo_R", "Harrell<U+2019>s C", "Uno<U+2019>s C", "Brier Score",
  #   "Time Dependent Auc") -- misspelling and curly apostrophes included.
  defaults <- eval(formals(
    TimeMetric:::pam.predicted_survial_eval_two_phase
  )$metrics)

  expect_true("Pesudo_R" %in% defaults)
  expect_true("Harrell\u2019s C" %in% defaults)
  expect_false("Pseudo_R" %in% defaults)
  expect_snapshot_value(defaults, style = "serialize")
})
```

- [ ] **Step 2: Run and record snapshots**

Run: `Rscript -e 'testthat::test_local(filter = "two-phase")'`
Expected: PASS, snapshots created.

`cc_weights` and `ncc_weights` delegate to `weighted_param`; the `if (is.list(w)) w[[1]] else w` guard handles either return shape, but record the actual shape in `findings.md` since the roxygen does not document it. If `length(wv)` is not 200, the weights are returned only for a sampled subset — record that and relax the length assertion to match, since it changes how `tm_evaluate_two_phase` must be documented.

If `pam.sample_design` requires additional per-model fields, read `R/pam.predicted_survial_eval_two_phase.R:183-213` for the `sprintf("Model '%s' must contain: %s", ...)` message, which names them exactly.

- [ ] **Step 3: Run again to confirm stability**

Run: `Rscript -e 'testthat::test_local(filter = "two-phase")'`
Expected: PASS, no snapshots added.

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/test-two-phase.R tests/testthat/_snaps/two-phase.md
git commit -m "test: characterize case-cohort and NCC two-phase evaluation"
```

---

### Task 9: Characterize R_E dispatch and plotting

Covers the `UseMethod("pam.rsph")` methods the spec identified as live-but-invisible, plus the two exported plot functions.

**Files:**
- Test: `tests/testthat/test-rsph-dispatch.R`
- Test: `tests/testthat/test-plots.R`

**Interfaces:**
- Consumes: `fx_surv()`, `fx_cox()`, `fx_survreg()`, `fx_covs()`, `snap_num()`
- Produces: proof that `pam.rsph.coxph` and `pam.rsph.survreg` are reachable, protecting them from the deletion sweep in spec step 5

- [ ] **Step 1: Write the dispatch tests**

Create `tests/testthat/test-rsph-dispatch.R`:
```r
test_that("pam.rsph dispatches to the coxph method", {
  d <- fx_surv()

  res <- TimeMetric:::pam.rsph(fx_cox(), test_data = d)

  expect_type(res, "list")
  # components pam.summary.rsph consumes
  expect_true(all(c("meanr", "ranks", "perfr") %in% names(res)))
  expect_true(all(is.finite(res$ranks)))
  expect_snapshot_value(sort(names(res)), style = "serialize")
  expect_snapshot_value(snap_num(head(res$ranks, 10)), style = "serialize")
})

test_that("pam.rsph dispatches to the survreg method", {
  d <- fx_surv()

  res <- TimeMetric:::pam.rsph(fx_survreg(), test_data = d)

  expect_type(res, "list")
  expect_snapshot_value(sort(names(res)), style = "serialize")
})

test_that("pam.rsph dispatch is class-driven, not identical across model types", {
  d <- fx_surv()

  cox_res <- TimeMetric:::pam.rsph(fx_cox(), test_data = d)
  reg_res <- TimeMetric:::pam.rsph(fx_survreg(), test_data = d)

  expect_false(isTRUE(all.equal(cox_res$ranks, reg_res$ranks,
                                tolerance = 1e-6)))
})

test_that("pam.summary.rsph converts an rsph object into R_E over time", {
  d <- fx_surv()
  obj <- TimeMetric:::pam.rsph(fx_cox(), test_data = d)

  res <- TimeMetric:::pam.summary.rsph(obj, times = stats::median(d$time))

  expect_true(is.numeric(res) || is.list(res) || is.data.frame(res))
  expect_snapshot_value(
    if (is.numeric(res)) snap_num(res) else sort(names(res)),
    style = "serialize"
  )
})

test_that("no pam.print generic exists, so pam.print.rsph is unreachable", {
  # Documents FINDING: pam.print.rsph can never be dispatched.
  ns <- asNamespace("TimeMetric")

  expect_false(exists("pam.print", envir = ns, inherits = FALSE))
  expect_true(exists("pam.print.rsph", envir = ns, inherits = FALSE))
})
```

- [ ] **Step 2: Write the plot tests**

Create `tests/testthat/test-plots.R`:
```r
test_that("plot_pred builds a ggplot with a point layer and the given labels", {
  d <- fx_surv()
  pred <- pam.coxph_restricted(model = fx_cox(), covs = fx_covs(),
                               new_data = d, tau = 10e10)

  p <- plot_pred(pred, title = "characterization", xlab = "RS", ylab = "Days")

  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$title, "characterization")
  expect_identical(p$labels$x, "RS")
  expect_identical(p$labels$y, "Days")
  expect_gt(length(p$layers), 0)
  # the plot must actually render, not merely construct
  built <- ggplot2::ggplot_build(p)
  expect_s3_class(built, "ggplot_built")
  expect_gt(nrow(built$data[[1]]), 0)
})

test_that("plot_pred uses the prediction data, not an empty frame", {
  d <- fx_surv()
  pred <- pam.coxph_restricted(model = fx_cox(), covs = fx_covs(),
                               new_data = d, tau = 10e10)

  built <- ggplot2::ggplot_build(plot_pred(pred))

  expect_identical(nrow(built$data[[1]]), 200L)
})

test_that("summary_pred_plot combines one panel per input", {
  d <- fx_surv()
  pred <- pam.coxph_restricted(model = fx_cox(), covs = fx_covs(),
                               new_data = d, tau = 10e10)

  p <- summary_pred_plot(list(a = pred, b = pred), ncol = 2)

  expect_true(inherits(p, "patchwork") || inherits(p, "ggplot"))
  if (inherits(p, "patchwork")) {
    expect_identical(length(p$patches$plots) + 1L, 2L)
  }
  expect_silent(print(p))
})
```

- [ ] **Step 3: Run both files and record snapshots**

Run: `Rscript -e 'testthat::test_local(filter = "rsph-dispatch")'` then `Rscript -e 'testthat::test_local(filter = "plots")'`
Expected: PASS, snapshots created.

If `pam.rsph(fit, test_data = d)` errors because `Gmat` has no default, read `R/pam.rsph.R:76-90` for how it is constructed, replicate that construction in the test, and record the undocumented requirement in `findings.md`.

If `plot_pred` expects a data frame rather than the prediction list, read `R/plot.R:49-90` for the columns it references and build a minimal frame with exactly those, keeping `nrow` at 200 so the second test's assertion stays meaningful.

If `expect_silent(print(p))` fails because ggplot emits a message about removed rows, replace it with `expect_no_error(print(p))`.

- [ ] **Step 4: Run again to confirm stability**

Run: `Rscript -e 'testthat::test_local(filter = "rsph-dispatch")'`
Expected: PASS, no snapshots added.

- [ ] **Step 5: Commit**

```bash
git add tests/testthat/test-rsph-dispatch.R tests/testthat/test-plots.R tests/testthat/_snaps/
git commit -m "test: characterize R_E dispatch and plotting functions"
```

---

### Task 10: Baseline verification and coverage record

**Files:**
- Create: `docs/superpowers/baseline-coverage.md`
- Create: `docs/superpowers/baseline-check.txt`
- Modify: `docs/superpowers/findings.md`

**Interfaces:**
- Consumes: every test file from Tasks 1-9
- Produces: the coverage floor spec step 5 must not drop below

- [ ] **Step 1: Run the entire suite**

Run: `Rscript -e 'testthat::test_local()'`
Expected: all tests PASS, 0 failures, no "adding new snapshot" messages.

- [ ] **Step 2: Confirm every exported function is exercised**

```bash
Rscript -e '
exports <- getNamespaceExports("TimeMetric")
src <- paste(readLines(Sys.glob("tests/testthat/test-*.R")), collapse = "\n")
missing <- exports[!vapply(exports, function(f) grepl(f, src, fixed = TRUE), logical(1))]
if (length(missing)) { cat("UNCOVERED EXPORTS:\n"); cat(paste0("  ", missing, collapse = "\n"), "\n") } else cat("all exports referenced\n")'
```
Expected: `all exports referenced`. Any name listed here needs a test before this task completes.

- [ ] **Step 3: Confirm the suite runs inside the time budget**

Run: `Rscript -e 'system.time(testthat::test_local())'`
Expected: elapsed under 60 seconds. If it exceeds that, reduce fixture `n` from 200 to 100 in `helper-simdata.R`, delete `tests/testthat/_snaps/`, re-run the "record snapshots" step of Tasks 2-9, and update the hard-coded `200L` length assertions to `100L`.

- [ ] **Step 4: Measure and record coverage**

```bash
Rscript -e 'cov <- covr::package_coverage(); print(cov); writeLines(capture.output(print(cov)), "docs/superpowers/baseline-coverage.md")'
```

- [ ] **Step 5: Annotate the coverage record**

Prepend to `docs/superpowers/baseline-coverage.md`:
```
# Baseline Coverage - characterization suite

Measured before any rename, dependency change, or deletion. This is the floor:
spec step 5 (dead-code removal) must not reduce coverage of any function that
survives. Functions reported at 0% here are candidates confirmed unreachable by
the call graph in the spec, section 2.
```

- [ ] **Step 6: Record the starting R CMD check state**

```bash
R CMD build . >/dev/null 2>&1 && R CMD check --no-manual TimeMetric_0.1.0.tar.gz 2>&1 | tail -40 > docs/superpowers/baseline-check.txt
cat docs/superpowers/baseline-check.txt
```
Expected: the pre-existing errors and warnings about undeclared imports remain. This is the recorded starting point, not a passing check — spec step 2 is what clears them.

- [ ] **Step 7: Confirm all findings are logged**

Review `docs/superpowers/findings.md`. It must contain the five seeded entries plus every surprise encountered in Tasks 2-9. Each row needs a location, the observed behaviour, why it looks wrong, and status `Open`.

- [ ] **Step 8: Verify encoding of every new file**

```bash
file -I tests/testthat/*.R tests/testthat.R docs/superpowers/*.md
LC_ALL=C grep -lP '[^\x00-\x7F]' tests/testthat/*.R tests/testthat.R || echo "all test sources ASCII-only"
```
Expected: every file `charset=utf-8` or `charset=us-ascii`, and the grep reporting no non-ASCII test sources. The curly apostrophes are written as `’` escapes, so this must pass.

- [ ] **Step 9: Confirm nothing outside tests and docs changed**

```bash
git diff main --stat
```
Expected: additions under `tests/` and `docs/` only, plus the two `DESCRIPTION` lines from Task 1. Any other modified file means a refactor leaked into this step.

- [ ] **Step 10: Commit**

```bash
git add docs/superpowers/baseline-coverage.md docs/superpowers/baseline-check.txt docs/superpowers/findings.md
git commit -m "test: record baseline coverage and check output before remediation"
```

---

## Done criteria for this step

1. `testthat::test_local()` passes with 0 failures and completes in under 60 seconds
2. Every exported function is referenced by at least one test (Task 10 step 2 reports none uncovered)
3. A second consecutive run adds no new snapshots — proving comparison, not re-recording
4. No snapshot contains a fitted model object or a full prediction matrix
5. Every `expect_error()` asserts a specific message fragment
6. Coverage and `R CMD check` baselines recorded in `docs/superpowers/`
7. `findings.md` lists every suspicious behaviour encountered, all status `Open`
8. `git diff main --stat` shows only additions under `tests/`, `docs/`, and two `DESCRIPTION` lines
