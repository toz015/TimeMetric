# TimeMetric Characterization Test Baseline — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Pin the current numeric behaviour of every reachable function in `TimeMetric` with a deterministic test suite, so that all later renames, dependency changes, and deletions are provably behaviour-preserving.

**Architecture:** `testthat` 3rd edition suite under `tests/testthat/`. Every test pairs **structural assertions** (hard-coded names, types, lengths, ranges — readable and reviewable in the diff) with a **serialized snapshot** of the full return value (generated on first run and committed). Structural assertions catch shape drift even if someone regenerates snapshots; snapshots catch numeric drift to 8 decimal places without anyone hand-transcribing floats into the plan. All randomness is seeded in shared fixtures.

**Tech Stack:** R 4.5+, `testthat` (3e), `devtools`, `roxygen2`, `covr`. Package runtime deps: `survival`, `rms`, `expint`, `survminer`, `pec`, `tdROC`, `yardstick`, `ggplot2`, `patchwork`, `dplyr`, `magrittr`, `purrr`, `tibble`.

## Global Constraints

- This step **adds tests only**. No renames, no deletions, no refactors, no `DESCRIPTION`/`NAMESPACE` edits beyond adding `testthat` test infrastructure. Those are steps 2–8 of the spec.
- Tests encode behaviour **as it is today, correct or not**. When a test pins something that looks wrong, record it in `docs/superpowers/findings.md` and move on — do not fix it here.
- Every test must be deterministic. Any test touching randomness sets a seed via the shared fixtures in `tests/testthat/helper-simdata.R`.
- Snapshot comparisons use `tolerance = 1e-8`.
- All new files UTF-8 encoded with ASCII-only content.
- Functions are called by their **current** names (`pam.*`). The `tm_` rename is step 3 of the spec and must not be anticipated here.
- Target total suite runtime under 60 seconds. Fixtures use `n = 200`, never the `n = 3000`/`n = 10000` of `paper.code.Rmd`.
- Commit after every task. Never use `testthat::snapshot_accept()` after the baseline task without reviewing the diff.

## Reference: verified function signatures

Copied from source. Later tasks depend on these being exact.

```r
sim_cox_weibull_censored(n, pi_c, v, beta, mu = NULL, sd = NULL, seed = NULL,
                         interact = FALSE, nonlinear = FALSE)
simulateTwoCauseFineGrayModel(n, v, beta1, beta2, lambda1 = 1, X = NULL, mu = 0,
                              p = 0.7, c_scale = 1, censor = 0, sd.time = NULL,
                              mu.c.time = NULL, independent_c = TRUE,
                              report.mu_and_sd = FALSE, seed = 1234)
pam.coxph_restricted(model, covs, tau = 10e10, new_data = NULL, predict = T)
pam.surverg_restricted(model, covs, tau = 10e10, new_data = NULL, predict = T)
pam.predict_cr(model1 = NULL, model2 = NULL, fg_model = NULL, cr_model = NULL,
               tau = NULL, newdata, event.type = 1, covs)
pam.predicted_survial_eval(model, event_time, predicted_probability,
                           pred_mean_survival = NULL, status, covariates,
                           new_data = NULL, metrics = NULL, t_star = NULL, tau = NULL)
pam.predicted_survial_eval_cr(pred_cif, event_time, time.cif, status,
                              metrics = NULL, t_star = NULL, tau = NULL, event_type = 1)
pam.predicted_survial_eval_two_phase(pred_results, t_star = NULL, tau = 10e10,
                                     km_cens_fit, case_weights, metrics = c(...))
pam.summary(models, metrics = NULL, t_star = NULL, tau = 10e10, digits = 2)
pam.summary_cr(models, metrics = NULL, t_star = NULL, tau = NULL, event_type = 1, digits = 2)
pam.sample_design(models, case_weights, km_cens, metrics = c(...), t_star = NULL, tau = NULL, digits = 2)
cc_weights(time, status, subcohort = NULL, strata = NULL)
ncc_weights(time, status, strata = NULL, m = NULL)
pam.rsh_metric(predicted_data, survival_time, status)
pam.Brier_metric(predicted_data, suvival_time, t_star = -1)   # note: arg is misspelled in source
pam.Brier(object, pre_sp, t_star = -1)
Gt(object, timepoint)
pam.rsph(fit, ...)                                            # UseMethod generic
pam.summary.rsph(object, times, band = 5, ...)
plot_pred(data, title = NULL, xlab = "Risk Score", ylab = "Days", ...)
summary_pred_plot(data_list, titles = NULL, plot_fun = plot_pred, ncol = 2, ...)
```

**Critical gotcha:** `sim_cox_weibull_censored` returns different columns depending on `pi_c`. With `pi_c == 0` it returns `time, status (all 1), x1..xp`. With `pi_c > 0` it returns `time, status, x1..xp, y_true, cens_time` plus `mu`/`sd_log` attributes. A naive `coxph(Surv(time, status) ~ ., data = df)` on the censored version silently fits `y_true` and `cens_time` as covariates. Fixtures must select covariate columns explicitly.

---

### Task 1: Test infrastructure and local dependency install

**Files:**
- Create: `tests/testthat.R`
- Create: `tests/testthat/test-smoke.R`
- Create: `docs/superpowers/findings.md`
- Modify: `DESCRIPTION` (add `Config/testthat/edition: 3`; move `testthat` to `Suggests` with a version floor)

**Interfaces:**
- Consumes: nothing
- Produces: a loadable package under `devtools::load_all()`, and a working `devtools::test()` entry point that all later tasks call.

Note: installing dependencies **locally** so the package loads is distinct from **declaring** them in `DESCRIPTION`, which is step 2 of the spec and out of scope here. This task only installs.

- [ ] **Step 1: Install the toolchain and every runtime dependency**

Run:
```bash
Rscript -e 'install.packages(c("devtools","roxygen2","testthat","covr","rms","expint","survminer","pec","tdROC","yardstick","patchwork","dplyr","magrittr","purrr","tibble","ggplot2","survival"), repos="https://cloud.r-project.org")'
```

`rms`, `pec`, and `survminer` compile and pull large dependency trees; expect this to take 10-20 minutes. It is a one-time cost.

- [ ] **Step 2: Verify every dependency loads**

Run:
```bash
Rscript -e 'for (p in c("survival","rms","expint","survminer","pec","tdROC","yardstick","patchwork","dplyr","magrittr","purrr","tibble","ggplot2")) { ok <- requireNamespace(p, quietly=TRUE); cat(sprintf("%-12s %s\n", p, if (ok) "OK" else "MISSING")) }'
```
Expected: every line `OK`. Do not proceed with any `MISSING` — later tasks will fail in confusing ways.

- [ ] **Step 3: Verify the package itself loads**

Run:
```bash
Rscript -e 'devtools::load_all("."); cat("loaded OK\n")'
```
Expected: `loaded OK`. Warnings about undeclared imports are expected at this stage and are fixed in spec step 2.

- [ ] **Step 4: Create the testthat entry point**

Create `tests/testthat.R`:
```r
library(testthat)
library(TimeMetric)

test_check("TimeMetric")
```

- [ ] **Step 5: Enable testthat 3rd edition**

In `DESCRIPTION`, ensure these lines are present (edit `Suggests` in place; leave every other field untouched):
```
Suggests: testthat (>= 3.0.0)
Config/testthat/edition: 3
```

- [ ] **Step 6: Write a smoke test that proves the harness runs**

Create `tests/testthat/test-smoke.R`:
```r
test_that("package loads and exports the expected function count", {
  exports <- getNamespaceExports("TimeMetric")
  expect_type(exports, "character")
  expect_true("pam.coxph_restricted" %in% exports)
  expect_true("pam.summary" %in% exports)
})
```

- [ ] **Step 7: Run the suite to verify the harness works**

Run: `Rscript -e 'devtools::test()'`
Expected: PASS, 3 assertions, 0 failures.

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
```

- [ ] **Step 9: Commit**

```bash
git add tests/testthat.R tests/testthat/test-smoke.R DESCRIPTION docs/superpowers/findings.md
git commit -m "test: add testthat infrastructure and findings log"
```

---

### Task 2: Deterministic shared fixtures

**Files:**
- Create: `tests/testthat/helper-simdata.R`
- Test: `tests/testthat/test-fixtures.R`

**Interfaces:**
- Consumes: `sim_cox_weibull_censored`, `simulateTwoCauseFineGrayModel` from Task 1's loaded package
- Produces: fixture functions used by every later task —
  - `fx_surv()` -> `data.frame` with columns `time`, `status`, `x1`, `x2`; 200 rows; ~30% censored
  - `fx_surv_uncensored()` -> `data.frame` with columns `time`, `status`, `x1`, `x2`; 200 rows; all `status == 1`
  - `fx_cr()` -> `data.frame` with columns `obs.times`, `obs.event`, and the generated covariate columns (name them from the recorded snapshot in Task 2 step 3 — `simulateTwoCauseFineGrayModel` builds them via `data.frame(obs.times, obs.event, X)`, so they follow `X`'s `colnames`); 200 rows; competing-risk codes
  - `fx_cox()` -> fitted `coxph` object on `fx_surv()`, `x = TRUE, y = TRUE`
  - `fx_survreg()` -> fitted `survreg` object on `fx_surv()`
  - `fx_covs()` -> `c("x1", "x2")`

- [ ] **Step 1: Write the fixture helper**

Create `tests/testthat/helper-simdata.R`:
```r
# Deterministic fixtures shared by all characterization tests.
# Seeds are fixed; never change them without regenerating every snapshot.

fx_covs <- function() c("x1", "x2")

fx_surv <- function() {
  d <- sim_cox_weibull_censored(
    n = 200, pi_c = 0.3, v = 2, beta = c(0.5, -0.5), seed = 1001
  )
  # pi_c > 0 adds y_true and cens_time; drop them so model formulas stay clean
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

- [ ] **Step 2: Write tests that pin the fixtures themselves**

Create `tests/testthat/test-fixtures.R`:
```r
test_that("fx_surv is deterministic and correctly shaped", {
  d1 <- fx_surv()
  d2 <- fx_surv()

  expect_s3_class(d1, "data.frame")
  expect_identical(names(d1), c("time", "status", "x1", "x2"))
  expect_identical(nrow(d1), 200L)
  expect_identical(d1, d2)                       # same seed, same data
  expect_true(all(d1$status %in% c(0, 1)))
  expect_true(all(d1$time > 0))
  expect_gt(sum(d1$status == 0), 0)              # some censoring occurred
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
  expect_true("obs.times" %in% names(d))
  expect_true("obs.event" %in% names(d))
  expect_true(all(d$obs.event %in% c(0, 1, 2)))
  expect_gt(sum(d$obs.event == 2), 0)            # competing events present
})

test_that("fixture models fit and are reproducible", {
  m <- fx_cox()
  expect_s3_class(m, "coxph")
  expect_identical(names(coef(m)), c("x1", "x2"))
  expect_equal(coef(m), coef(fx_cox()))

  s <- fx_survreg()
  expect_s3_class(s, "survreg")
})

test_that("fixture data and models are snapshot-stable", {
  expect_snapshot_value(fx_surv(), style = "serialize", tolerance = 1e-8)
  expect_snapshot_value(fx_cr(), style = "serialize", tolerance = 1e-8)
  expect_snapshot_value(unname(coef(fx_cox())), style = "serialize", tolerance = 1e-8)
})
```

- [ ] **Step 3: Run the tests and record the snapshots**

Run: `Rscript -e 'devtools::test(filter = "fixtures")'`

Expected: PASS. On this first run testthat **creates** `tests/testthat/_snaps/fixtures.md` and reports the snapshots as added rather than compared. That is correct — this run is what establishes the baseline.

If `fx_cr()` errors, inspect `simulateTwoCauseFineGrayModel`'s argument handling and adjust `beta1`/`beta2` lengths to match the covariate count it generates; record any signature surprise in `findings.md`.

- [ ] **Step 4: Run a second time to verify snapshots now compare**

Run: `Rscript -e 'devtools::test(filter = "fixtures")'`
Expected: PASS with no "adding new snapshot" messages. This proves the snapshots are stable across runs rather than being re-recorded each time.

- [ ] **Step 5: Commit**

```bash
git add tests/testthat/helper-simdata.R tests/testthat/test-fixtures.R tests/testthat/_snaps/fixtures.md
git commit -m "test: add deterministic fixtures and pin their output"
```

---

### Task 3: Characterize Cluster C metric functions

Highest priority: these are the spec's deletion candidates and have zero coverage. `pam.rsh_metric`, `pam.rsph_metric`, and `pam.Brier_metric` are internal, so tests reach them via `:::`.

**Files:**
- Test: `tests/testthat/test-cluster-c-metrics.R`

**Interfaces:**
- Consumes: `fx_surv()` from Task 2
- Produces: the regression suite that gates deletion of `pam.rsh_metric`, `pam.rsph_metric`, `pam.Brier_metric` in spec step 5

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
  expect_true(is.finite(res$D))
  expect_true(is.finite(res$Dx))
  expect_equal(res$D, round(res$D, 4))
  expect_snapshot_value(res, style = "serialize", tolerance = 1e-8)
})

test_that("pam.rsh_metric is sensitive to input order (pins current alignment bug)", {
  # FINDING #1: the function sorts its internal data frame by survival_time but
  # multiplies by the UNSORTED predicted_data argument. Permuting the inputs
  # consistently should not change the result, but it does. Pinned as-is.
  d <- fx_surv()
  pred <- seq(0.9, 0.1, length.out = nrow(d))

  ord <- order(d$time)
  sorted <- TimeMetric:::pam.rsh_metric(pred[ord], d$time[ord], d$status[ord])
  shuffled <- TimeMetric:::pam.rsh_metric(pred, d$time, d$status)

  expect_false(isTRUE(all.equal(sorted$Dx, shuffled$Dx)))
  expect_snapshot_value(list(sorted = sorted, shuffled = shuffled),
                        style = "serialize", tolerance = 1e-8)
})

test_that("pam.Brier_metric requires a Surv object and matching lengths", {
  d <- fx_surv()
  pred <- rep(0.5, nrow(d))
  sv <- survival::Surv(d$time, d$status)

  expect_error(TimeMetric:::pam.Brier_metric(pred, d$time), "Surv")
  expect_error(TimeMetric:::pam.Brier_metric(pred[-1], sv), "match")
  expect_error(TimeMetric:::pam.Brier_metric(c(NA, pred[-1]), sv), "NA")
})

test_that("pam.Brier_metric returns a finite scalar at the median time", {
  d <- fx_surv()
  sv <- survival::Surv(d$time, d$status)
  pred <- seq(0.9, 0.1, length.out = nrow(d))

  res <- TimeMetric:::pam.Brier_metric(pred, sv)

  expect_length(res, 1L)
  expect_true(is.finite(res))
  expect_snapshot_value(res, style = "serialize", tolerance = 1e-8)
})

test_that("pam.rsph_metric returns a list containing Re", {
  d <- fx_surv()
  pred <- seq(0.9, 0.1, length.out = nrow(d))

  res <- TimeMetric:::pam.rsph_metric(pred, d$time, d$status, min(d$time))

  expect_type(res, "list")
  expect_true("Re" %in% names(res))
  expect_snapshot_value(res, style = "serialize", tolerance = 1e-8)
})
```

- [ ] **Step 2: Run and record snapshots**

Run: `Rscript -e 'devtools::test(filter = "cluster-c-metrics")'`
Expected: PASS, snapshots created.

If `pam.rsph_metric`'s fourth argument is not a start time, read its signature in `R/pam.rsph_metric.R` and correct the call; record the actual signature in `findings.md`.

- [ ] **Step 3: Run again to confirm stability**

Run: `Rscript -e 'devtools::test(filter = "cluster-c-metrics")'`
Expected: PASS, no snapshots added.

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/test-cluster-c-metrics.R tests/testthat/_snaps/cluster-c-metrics.md
git commit -m "test: characterize Cluster C metric functions before deletion decision"
```

---

### Task 4: Characterize the pec dependency chain

Pins `Gt` -> `pam.Brier` -> `predictSurvProb`. This chain decides whether `pec` can leave `Imports`, so its behaviour must be locked before any rewrite.

**Files:**
- Test: `tests/testthat/test-brier-pec-chain.R`

**Interfaces:**
- Consumes: `fx_surv()`, `fx_cox()` from Task 2
- Produces: the gate for the `pec` demotion decision in spec step 2

- [ ] **Step 1: Write the characterization tests**

Create `tests/testthat/test-brier-pec-chain.R`:
```r
test_that("Gt returns the censoring survival probability at a timepoint", {
  d <- fx_surv()
  sv <- survival::Surv(d$time, d$status)

  res <- Gt(sv, stats::median(d$time))

  expect_length(res, 1L)
  expect_true(is.finite(res))
  expect_gte(res, 0)
  expect_lte(res, 1)
  expect_snapshot_value(res, style = "serialize", tolerance = 1e-8)
})

test_that("Gt rejects non-Surv input", {
  expect_error(Gt(1:10, 5), "class Surv")
})

test_that("pam.Brier on a coxph fit returns a finite scalar", {
  # This is the chain that reaches pec::predictSurvProb unqualified.
  d <- fx_surv()
  fit <- fx_cox()

  res <- TimeMetric:::pam.Brier(fit, d, stats::median(d$time))

  expect_length(res, 1L)
  expect_true(is.finite(res))
  expect_gte(res, 0)
  expect_snapshot_value(res, style = "serialize", tolerance = 1e-8)
})

test_that("pam.Brier default t_star uses the median observed time", {
  d <- fx_surv()
  fit <- fx_cox()

  default_res <- TimeMetric:::pam.Brier(fit, d)
  median_res  <- TimeMetric:::pam.Brier(fit, d, stats::median(d$time))

  expect_equal(default_res, median_res, tolerance = 1e-8)
})

test_that("pec::predictSurvProb is reachable and drives pam.Brier", {
  skip_if_not_installed("pec")
  d <- fx_surv()
  fit <- fx_cox()

  probs <- pec::predictSurvProb(fit, d, stats::median(d$time))

  expect_true(is.matrix(probs) || is.numeric(probs))
  expect_identical(NROW(probs), 200L)
  expect_true(all(probs >= 0 & probs <= 1))
  expect_snapshot_value(as.numeric(probs), style = "serialize", tolerance = 1e-8)
})
```

- [ ] **Step 2: Run and record snapshots**

Run: `Rscript -e 'devtools::test(filter = "brier-pec-chain")'`
Expected: PASS, snapshots created.

- [ ] **Step 3: Run again to confirm stability**

Run: `Rscript -e 'devtools::test(filter = "brier-pec-chain")'`
Expected: PASS, no snapshots added.

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/test-brier-pec-chain.R tests/testthat/_snaps/brier-pec-chain.md
git commit -m "test: characterize Gt/pam.Brier/predictSurvProb chain gating pec removal"
```

---

### Task 5: Characterize the prediction module

**Files:**
- Test: `tests/testthat/test-prediction-module.R`

**Interfaces:**
- Consumes: `fx_surv()`, `fx_cox()`, `fx_survreg()`, `fx_covs()` from Task 2
- Produces: pinned output shape of `pam.coxph_restricted` and `pam.surverg_restricted`, whose return structure (`surv_prob`, `time`, `status`, `pred`) every evaluation function consumes

- [ ] **Step 1: Write the characterization tests**

Create `tests/testthat/test-prediction-module.R`:
```r
test_that("pam.coxph_restricted returns the documented prediction structure", {
  d <- fx_surv()
  fit <- fx_cox()

  res <- pam.coxph_restricted(model = fit, covs = fx_covs(),
                              new_data = d, tau = 10e10)

  expect_type(res, "list")
  # pin the exact component names downstream evaluators rely on
  expect_snapshot_value(sort(names(res)), style = "serialize")
  expect_snapshot_value(res, style = "serialize", tolerance = 1e-8)
})

test_that("pam.coxph_restricted rejects a non-coxph model", {
  expect_error(
    pam.coxph_restricted(model = fx_survreg(), covs = fx_covs(),
                         new_data = fx_surv()),
    regexp = "."
  )
})

test_that("pam.surverg_restricted returns the documented prediction structure", {
  d <- fx_surv()
  fit <- fx_survreg()

  res <- pam.surverg_restricted(model = fit, covs = fx_covs(),
                                new_data = d, tau = 10e10)

  expect_type(res, "list")
  expect_snapshot_value(sort(names(res)), style = "serialize")
  expect_snapshot_value(res, style = "serialize", tolerance = 1e-8)
})

test_that("pam.surverg_restricted rejects a non-survreg model", {
  expect_error(
    pam.surverg_restricted(model = fx_cox(), covs = fx_covs(),
                           new_data = fx_surv()),
    "survreg"
  )
})

test_that("prediction output length matches the input data", {
  d <- fx_surv()
  res <- pam.coxph_restricted(model = fx_cox(), covs = fx_covs(),
                              new_data = d, tau = 10e10)

  expect_identical(length(res$time), 200L)
  expect_identical(length(res$status), 200L)
})
```

- [ ] **Step 2: Run and record snapshots**

Run: `Rscript -e 'devtools::test(filter = "prediction-module")'`
Expected: PASS, snapshots created.

The `sort(names(res))` snapshot documents the real component names. Read the recorded `_snaps/prediction-module.md` and confirm the final test's `res$time` / `res$status` references match; if the components are named differently, correct that test to the real names and note the discrepancy against the roxygen docs in `findings.md`.

- [ ] **Step 3: Run again to confirm stability**

Run: `Rscript -e 'devtools::test(filter = "prediction-module")'`
Expected: PASS, no snapshots added.

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/test-prediction-module.R tests/testthat/_snaps/prediction-module.md
git commit -m "test: characterize coxph and survreg restricted prediction functions"
```

---

### Task 6: Characterize the right-censored evaluation path

**Files:**
- Test: `tests/testthat/test-eval-survival.R`

**Interfaces:**
- Consumes: `fx_surv()`, `fx_cox()`, `fx_covs()` from Task 2; the prediction structure pinned in Task 5
- Produces: pinned metric values and the exact metric-name strings that spec step 6 will standardize

- [ ] **Step 1: Write the characterization tests**

Create `tests/testthat/test-eval-survival.R`:
```r
test_that("pam.summary returns a data frame of metrics with pinned column names", {
  d <- fx_surv()
  pred <- pam.coxph_restricted(model = fx_cox(), covs = fx_covs(),
                               new_data = d, tau = 10e10)

  res <- pam.summary(list(value = pred), tau = 10e10)

  expect_true(is.data.frame(res))
  expect_gt(nrow(res), 0)
  # Pins the CURRENT metric identifiers, including the "Pesudo" misspelling and
  # any curly apostrophes. Spec step 6 changes these deliberately.
  expect_snapshot_value(res, style = "serialize", tolerance = 1e-8)
})

test_that("pam.summary rejects a non-list or empty models argument", {
  expect_error(pam.summary(list()), "non-empty")
  expect_error(pam.summary("not a list"), "non-empty")
})

test_that("pam.predicted_survial_eval returns metrics for a coxph model", {
  d <- fx_surv()
  pred <- pam.coxph_restricted(model = fx_cox(), covs = fx_covs(),
                               new_data = d, tau = 10e10)

  res <- pam.predicted_survial_eval(
    model = fx_cox(),
    event_time = pred$time,
    predicted_probability = pred$surv_prob,
    status = pred$status,
    covariates = fx_covs(),
    new_data = d,
    tau = 10e10
  )

  expect_type(res, "list")
  expect_snapshot_value(sort(names(res)), style = "serialize")
  expect_snapshot_value(res, style = "serialize", tolerance = 1e-8)
})

test_that("the default metric set is pinned including R_sh and R_E", {
  d <- fx_surv()
  pred <- pam.coxph_restricted(model = fx_cox(), covs = fx_covs(),
                               new_data = d, tau = 10e10)

  res <- pam.predicted_survial_eval(
    model = fx_cox(),
    event_time = pred$time,
    predicted_probability = pred$surv_prob,
    status = pred$status,
    covariates = fx_covs(),
    new_data = d,
    tau = 10e10
  )

  # R_sh reaches rms::cph + pam.schemper; R_E reaches pam.rsph dispatch.
  # Both are default metrics, so this test also proves those paths execute.
  expect_true(any(grepl("R_sh", names(res), fixed = TRUE)))
  expect_true(any(grepl("R_E", names(res), fixed = TRUE)))
})
```

- [ ] **Step 2: Run and record snapshots**

Run: `Rscript -e 'devtools::test(filter = "eval-survival")'`
Expected: PASS, snapshots created.

If `pred$surv_prob` is not the component name, use the name recorded in Task 5's snapshot. If `R_sh` comes back `NA` with a message about Cox-only support or factor variables, that is current behaviour — keep the test and record the condition in `findings.md`.

- [ ] **Step 3: Run again to confirm stability**

Run: `Rscript -e 'devtools::test(filter = "eval-survival")'`
Expected: PASS, no snapshots added.

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/test-eval-survival.R tests/testthat/_snaps/eval-survival.md
git commit -m "test: characterize right-censored evaluation and metric names"
```

---

### Task 7: Characterize the competing-risks path

**Files:**
- Test: `tests/testthat/test-eval-competing-risks.R`

**Interfaces:**
- Consumes: `fx_cr()` from Task 2
- Produces: pinned behaviour of `pam.predict_cr`, `pam.predicted_survial_eval_cr`, `pam.summary_cr`, and the `m_cif` helper

- [ ] **Step 1: Write the characterization tests**

Create `tests/testthat/test-eval-competing-risks.R`:
```r
test_that("m_cif reduces a CIF column to a scalar risk", {
  cif <- seq(0, 0.6, length.out = 50)
  times <- seq(1, 50, length.out = 50)

  res <- TimeMetric:::m_cif(cif, time.cif = times, tau = 50)

  expect_length(res, 1L)
  expect_true(is.finite(res))
  expect_snapshot_value(res, style = "serialize", tolerance = 1e-8)
})

test_that("pam.predicted_survial_eval_cr returns metrics for supplied CIFs", {
  d <- fx_cr()
  n <- nrow(d)
  time_grid <- sort(unique(d$obs.times))
  # deterministic synthetic CIF: monotone increasing per subject
  cif <- outer(seq_len(length(time_grid)), seq_len(n),
               function(i, j) pmin(0.99, i / length(time_grid) * (0.2 + 0.6 * j / n)))

  res <- pam.predicted_survial_eval_cr(
    pred_cif = cif,
    event_time = d$obs.times,
    time.cif = time_grid,
    status = d$obs.event,
    event_type = 1
  )

  expect_type(res, "list")
  expect_snapshot_value(sort(names(res)), style = "serialize")
  expect_snapshot_value(res, style = "serialize", tolerance = 1e-8)
})

test_that("pam.predicted_survial_eval_cr errors when required arguments are missing", {
  d <- fx_cr()
  expect_error(
    pam.predicted_survial_eval_cr(event_time = d$obs.times),
    regexp = "."
  )
})

test_that("competing-risks status codes include a third level", {
  # Documents that CR functions accept {0,1,2}, unlike single-event functions.
  # Spec step 6 makes validation function-specific on this basis.
  d <- fx_cr()
  expect_setequal(sort(unique(d$obs.event)), c(0, 1, 2))
})
```

- [ ] **Step 2: Run and record snapshots**

Run: `Rscript -e 'devtools::test(filter = "eval-competing-risks")'`
Expected: PASS, snapshots created.

If `pam.predicted_survial_eval_cr` requires `pred_cif` with times as rows and subjects as columns in the opposite orientation, transpose `cif` and record the expected orientation in `findings.md` — the roxygen documentation does not currently state it.

- [ ] **Step 3: Run again to confirm stability**

Run: `Rscript -e 'devtools::test(filter = "eval-competing-risks")'`
Expected: PASS, no snapshots added.

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/test-eval-competing-risks.R tests/testthat/_snaps/eval-competing-risks.md
git commit -m "test: characterize competing-risks evaluation path"
```

---

### Task 8: Characterize the two-phase designs

The case-cohort and NCC paths the paper advertises but does not export. Zero coverage today.

**Files:**
- Test: `tests/testthat/test-two-phase.R`

**Interfaces:**
- Consumes: `fx_surv()`, `fx_cox()`, `fx_covs()` from Task 2
- Produces: the regression suite that gates exporting `tm_evaluate_two_phase` in spec step 3

- [ ] **Step 1: Write the characterization tests**

Create `tests/testthat/test-two-phase.R`:
```r
test_that("cc_weights returns one weight per subject", {
  d <- fx_surv()
  set.seed(2001)
  subcohort <- rbinom(nrow(d), 1, 0.4)

  w <- cc_weights(time = d$time, status = d$status, subcohort = subcohort)

  expect_true(is.numeric(w) || is.list(w))
  expect_snapshot_value(w, style = "serialize", tolerance = 1e-8)
})

test_that("ncc_weights requires m and returns weights", {
  d <- fx_surv()

  expect_error(ncc_weights(time = d$time, status = d$status), "`m` must be provided")

  w <- ncc_weights(time = d$time, status = d$status, m = 2)
  expect_true(is.numeric(w) || is.list(w))
  expect_snapshot_value(w, style = "serialize", tolerance = 1e-8)
})

test_that("pam.predicted_survial_eval_two_phase evaluates a case-cohort design", {
  d <- fx_surv()
  pred <- pam.coxph_restricted(model = fx_cox(), covs = fx_covs(),
                               new_data = d, tau = 10e10)
  set.seed(2002)
  subcohort <- rbinom(nrow(d), 1, 0.4)
  w <- cc_weights(time = d$time, status = d$status, subcohort = subcohort)
  km_cens <- survival::survfit(survival::Surv(d$time, 1 - d$status) ~ 1)

  res <- TimeMetric:::pam.predicted_survial_eval_two_phase(
    pred_results = pred,
    km_cens_fit = km_cens,
    case_weights = if (is.list(w)) w[[1]] else w
  )

  expect_true(inherits(res, "data.frame") || inherits(res, "tbl_df"))
  expect_true(all(c("Metric", "Value") %in% names(res)))
  expect_snapshot_value(res, style = "serialize", tolerance = 1e-8)
})

test_that("pam.predicted_survial_eval_two_phase evaluates an NCC design", {
  d <- fx_surv()
  pred <- pam.coxph_restricted(model = fx_cox(), covs = fx_covs(),
                               new_data = d, tau = 10e10)
  w <- ncc_weights(time = d$time, status = d$status, m = 2)
  km_cens <- survival::survfit(survival::Surv(d$time, 1 - d$status) ~ 1)

  res <- TimeMetric:::pam.predicted_survial_eval_two_phase(
    pred_results = pred,
    km_cens_fit = km_cens,
    case_weights = if (is.list(w)) w[[1]] else w
  )

  expect_true(all(c("Metric", "Value") %in% names(res)))
  expect_snapshot_value(res, style = "serialize", tolerance = 1e-8)
})

test_that("two-phase default metric names are pinned with their current spelling", {
  # Default is c("Pesudo_R", "Harrell<U+2019>s C", "Uno<U+2019>s C", "Brier Score",
  # "Time Dependent Auc") -- misspelling and curly apostrophes included.
  # Spec step 6 replaces these; this pins what they are today.
  fml <- formals(TimeMetric:::pam.predicted_survial_eval_two_phase)
  expect_snapshot_value(eval(fml$metrics), style = "serialize")
})
```

- [ ] **Step 2: Run and record snapshots**

Run: `Rscript -e 'devtools::test(filter = "two-phase")'`
Expected: PASS, snapshots created.

`cc_weights` and `ncc_weights` delegate to `weighted_param`; if they return a list rather than a bare vector, the `if (is.list(w)) w[[1]] else w` guard handles it, but record the real return shape in `findings.md` since the roxygen does not document it.

- [ ] **Step 3: Run again to confirm stability**

Run: `Rscript -e 'devtools::test(filter = "two-phase")'`
Expected: PASS, no snapshots added.

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/test-two-phase.R tests/testthat/_snaps/two-phase.md
git commit -m "test: characterize case-cohort and NCC two-phase evaluation"
```

---

### Task 9: Characterize R_E dispatch and plotting

Covers the `UseMethod("pam.rsph")` S3 methods the spec identified as live-but-invisible, plus the two plot functions.

**Files:**
- Test: `tests/testthat/test-rsph-dispatch.R`
- Test: `tests/testthat/test-plots.R`

**Interfaces:**
- Consumes: `fx_surv()`, `fx_cox()`, `fx_survreg()`, `fx_covs()` from Task 2
- Produces: proof that `pam.rsph.coxph` and `pam.rsph.survreg` are reachable, protecting them from the deletion sweep in spec step 5

- [ ] **Step 1: Write the dispatch tests**

Create `tests/testthat/test-rsph-dispatch.R`:
```r
test_that("pam.rsph dispatches to the coxph method", {
  d <- fx_surv()
  fit <- fx_cox()

  res <- TimeMetric:::pam.rsph(fit, test_data = d)

  expect_type(res, "list")
  # components consumed by pam.summary.rsph
  expect_true(all(c("meanr", "ranks", "perfr") %in% names(res)))
  expect_snapshot_value(sort(names(res)), style = "serialize")
  expect_snapshot_value(res, style = "serialize", tolerance = 1e-8)
})

test_that("pam.rsph dispatches to the survreg method", {
  d <- fx_surv()
  fit <- fx_survreg()

  res <- TimeMetric:::pam.rsph(fit, test_data = d)

  expect_type(res, "list")
  expect_snapshot_value(sort(names(res)), style = "serialize")
})

test_that("pam.summary.rsph converts an rsph object into R_E over time", {
  d <- fx_surv()
  obj <- TimeMetric:::pam.rsph(fx_cox(), test_data = d)

  res <- TimeMetric:::pam.summary.rsph(obj, times = stats::median(d$time))

  expect_true(is.numeric(res) || is.list(res) || is.data.frame(res))
  expect_snapshot_value(res, style = "serialize", tolerance = 1e-8)
})

test_that("no pam.print generic exists, so pam.print.rsph is unreachable", {
  # Documents FINDING: pam.print.rsph can never be dispatched.
  expect_false(exists("pam.print", envir = asNamespace("TimeMetric"), inherits = FALSE))
})
```

- [ ] **Step 2: Write the plot tests**

Create `tests/testthat/test-plots.R`:
```r
test_that("plot_pred returns a ggplot object", {
  d <- fx_surv()
  pred <- pam.coxph_restricted(model = fx_cox(), covs = fx_covs(),
                               new_data = d, tau = 10e10)

  p <- plot_pred(pred)

  expect_s3_class(p, "ggplot")
})

test_that("summary_pred_plot combines multiple panels", {
  d <- fx_surv()
  pred <- pam.coxph_restricted(model = fx_cox(), covs = fx_covs(),
                               new_data = d, tau = 10e10)

  p <- summary_pred_plot(list(a = pred, b = pred), ncol = 2)

  expect_true(inherits(p, "patchwork") || inherits(p, "ggplot"))
})
```

- [ ] **Step 3: Run both files and record snapshots**

Run: `Rscript -e 'devtools::test(filter = "rsph-dispatch")'` then `Rscript -e 'devtools::test(filter = "plots")'`
Expected: PASS, snapshots created.

If `pam.rsph(fit, test_data = d)` errors because `Gmat` is a required argument without a default, inspect `R/pam.rsph.R:76` for how callers construct it, replicate that construction in the test, and record the undocumented requirement in `findings.md`.

If `plot_pred` needs a different input shape than the prediction list, read `R/plot.R:49` for the expected columns and build a minimal data frame matching it.

- [ ] **Step 4: Run again to confirm stability**

Run: `Rscript -e 'devtools::test(filter = "rsph-dispatch")'`
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
- Modify: `docs/superpowers/findings.md`

**Interfaces:**
- Consumes: every test file from Tasks 1-9
- Produces: the coverage floor that spec step 5 must not drop below

- [ ] **Step 1: Run the entire suite**

Run: `Rscript -e 'devtools::test()'`
Expected: all tests PASS, 0 failures, 0 warnings about missing snapshots.

- [ ] **Step 2: Confirm the suite runs inside the time budget**

Run: `Rscript -e 'system.time(devtools::test())'`
Expected: elapsed under 60 seconds. If it exceeds that, reduce fixture `n` from 200 to 100 in `helper-simdata.R`, delete `tests/testthat/_snaps/`, and re-run Tasks 2-9 step "record snapshots" to regenerate the baseline against the smaller fixtures.

- [ ] **Step 3: Measure and record coverage**

Run:
```bash
Rscript -e 'cov <- covr::package_coverage(); print(cov); writeLines(capture.output(print(cov)), "docs/superpowers/baseline-coverage.md")'
```

- [ ] **Step 4: Annotate the coverage record**

Prepend to `docs/superpowers/baseline-coverage.md`:
```markdown
# Baseline Coverage — characterization suite

Measured before any rename, dependency change, or deletion. This is the floor:
spec step 5 (dead-code removal) must not reduce coverage of any function that
survives. Functions reported at 0% here are candidates confirmed unreachable by
the call graph in the spec, section 2.

```

- [ ] **Step 5: Verify the package still checks no worse than before**

Run:
```bash
Rscript -e 'devtools::check(document = FALSE, args = c("--no-manual"))' 2>&1 | tail -20
```
Expected: the pre-existing errors and warnings about undeclared imports remain. This is the recorded starting point, not a passing check — the dependency reconciliation in spec step 2 is what clears them. Save the output:
```bash
Rscript -e 'devtools::check(document = FALSE, args = c("--no-manual"))' 2>&1 | tail -40 > docs/superpowers/baseline-check.txt
```

- [ ] **Step 6: Confirm all findings are logged**

Review `docs/superpowers/findings.md`. It must contain the three seeded entries plus every surprise encountered in Tasks 2-9. Each row needs a location, the observed behaviour, why it looks wrong, and status `Open`.

- [ ] **Step 7: Verify encoding of every new file**

Run:
```bash
file -I tests/testthat/*.R tests/testthat.R docs/superpowers/*.md
LC_ALL=C grep -lP '[^\x00-\x7F]' tests/testthat/*.R tests/testthat.R || echo "all test sources ASCII-only"
```
Expected: every file `charset=utf-8` or `charset=us-ascii`, and the grep reporting no non-ASCII test sources.

- [ ] **Step 8: Commit**

```bash
git add docs/superpowers/baseline-coverage.md docs/superpowers/baseline-check.txt docs/superpowers/findings.md
git commit -m "test: record baseline coverage and check output before remediation"
```

---

## Done criteria for this step

1. `devtools::test()` passes with 0 failures and completes in under 60 seconds
2. Snapshots exist and are committed for every characterized function
3. A second consecutive run adds no new snapshots — proving stability, not re-recording
4. Coverage and `R CMD check` baselines recorded in `docs/superpowers/`
5. `findings.md` lists every suspicious behaviour encountered, all status `Open`
6. No renames, deletions, or refactors have occurred — `git diff main --stat` shows only additions under `tests/`, `docs/`, and the two `DESCRIPTION` lines
