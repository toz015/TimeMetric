# CRAN Preparation — Implementation Plan

**Goal:** close the remaining local CRAN-readiness gaps on `joss-revision`, without
touching the manuscript, history, or the remote.

**Scope:** bounded. No design spec; this plan is the whole design.

**Baseline:** commit `d4bb31d`. Suite green at 498 assertions / 0 failures / 0 skips
under `testthat::test_local()`; `R CMD check --as-cran` at 0 errors, 0 warnings,
2 notes.

## Global constraints

- `paper.md` and `paper.code.Rmd` MUST NOT change.
- No history rewrite, push, force-push, tag, release, or CRAN submission.
- `devtools` is not installed: use `testthat::test_local()`,
  `roxygen2::roxygenise()`, `pkgload::load_all()`.
- Every claim reported must be measured, not expected.
- Commit at each checkpoint; stop and report.

---

## Task 1 — Metadata

**Files:** `DESCRIPTION`

- [ ] Add `URL: https://github.com/toz015/TimeMetric`
- [ ] Add `BugReports: https://github.com/toz015/TimeMetric/issues`
- [ ] Add `Language: en-US`
- [ ] Validate: DESCRIPTION parses, fields are well-formed, URLs resolve

```bash
Rscript -e 'd <- read.dcf("DESCRIPTION"); print(colnames(d)); print(d[,c("URL","BugReports","Language")])'
Rscript -e 'urlchecker::url_check(".")'
```

**Checkpoint 1:** report the added fields and the validation output.

---

## Task 2 — Remove `survminer`

Single call site: `R/pam.Ct.R:63`, `survminer::surv_summary(fit)`. Only `$time`
and `$surv` are read, but `na.omit()` runs over the whole frame, so the rows
dropped depend on `upper`/`lower` — columns the code never uses. A replacement
returning only `time`/`surv` would silently change behaviour.

**Order is mandatory: tests first, then the swap, then the dependency removal.**

- [ ] **2a. Characterization tests for `Gt()` BEFORE any change**, in
      `tests/testthat/test-gt-characterization.R`, covering:
      - ordinary right-censored data (`fx_surv()`)
      - the missing-confidence-limit row (`upper`/`lower` NA at the last time)
      - tied event times
      - timepoint exactly on an observed time; between two times (interpolation);
        before the first and after the last time
      - the degenerate no-deaths branch
- [ ] **2b. Run them against unmodified code** — must pass, pinning current values
- [ ] **2c. Replace the call** with an 8-column data frame built from the
      `survfit` object (`time`, `n.risk`, `n.event`, `n.censor`, `surv`,
      `std.err`, `upper`, `lower`), so `na.omit()` drops identically
- [ ] **2d. Gate:** `test-gt-characterization` passes unchanged, and every value
      in `tests/testthat/fixtures/metric-baseline.csv` is unchanged
- [ ] **2e. Only then** drop `survminer` from `DESCRIPTION` Imports; confirm no
      `survminer` reference remains
- [ ] **2f.** Remove the redundant `importFrom(purrr, map2)` — the only call is
      qualified `purrr::map2()`. **Keep `purrr` in Imports.**

**Checkpoint 2:** report the Gt equivalence evidence and the baseline result.

---

## Task 3 — Examples

- [ ] **3a. Re-audit exports now** — do not reuse the earlier count of nine
- [ ] **3b.** Add short, deterministic, scientifically meaningful examples for
      every primary public function lacking one. Use the `moore` dataset or
      `tm_sim_cox_weibull(seed = …)` so examples are reproducible and fast.
- [ ] **3c.** No `\dontrun{}` unless execution is genuinely impossible. Guard
      optional backends with `if (requireNamespace("pkg", quietly = TRUE))`.
- [ ] **3d.** Run and time them:

```bash
Rscript -e 'roxygen2::roxygenise()'
R CMD build . && R CMD check --as-cran --run-donttest TimeMetric_0.2.0.tar.gz
```

**Checkpoint 3:** report which functions gained examples, and the total example
runtime measured from the check log.

---

## Task 4 — Spelling

- [ ] Add `spelling` to `Suggests`
- [ ] Run `spelling::spell_check_package()`
- [ ] Create a **narrowly scoped** `inst/WORDLIST` — only genuine technical terms
- [ ] **Report every accepted word**, so a real misspelling cannot hide in the list

**Checkpoint 4:** the full WORDLIST, word by word, with justification.

---

## Task 5 — Community and citation

- [ ] `CONTRIBUTING.md`: how to report a bug, request a feature, submit a PR, run
      the tests. (A JOSS review-checklist item, and relevant to the desk-rejection
      grounds.)
- [ ] `inst/CITATION` using `bibentry()`. **No JOSS or Zenodo DOI** — neither
      exists yet; cite the package and its GitHub URL only.
- [ ] Verify against the **installed** package:

```bash
R CMD INSTALL --library=/tmp/tmlib .
Rscript -e '.libPaths("/tmp/tmlib"); print(citation("TimeMetric"))'
```

**Checkpoint 5:** the rendered `citation("TimeMetric")` output.

---

## Task 6 — Final CRAN preparation

- [ ] Full suite: `testthat::test_local()`
- [ ] `R CMD build` then `R CMD check --as-cran` on the tarball
- [ ] `spelling::spell_check_package()` clean
- [ ] Confirm `paper.md` / `paper.code.Rmd` byte-identical (SHA-256)
- [ ] Write `cran-comments.md` **from the actual final results**, not anticipated
      ones: test environments, check results, and each NOTE justified

**Checkpoint 6:** full results, then stop.
