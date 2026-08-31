# TimeMetric — Phase 1: Package Quality Remediation

**Date:** 2026-08-30
**Status:** Design approved, pending spec review
**Context:** JOSS desk rejection (editor: arfon, 2026-07-16), issue cited software readiness

## Background

`TimeMetric` was submitted to JOSS and desk-rejected before peer review on three
grounds: no automated tests, a repository not in reviewable state for an R
package, and thin evidence of research impact / open development.

The rejection concerned the *software*, not the manuscript. This spec covers only
the software remediation. Manuscript revision, CRAN submission, and the Zenodo
archive are deliberately deferred to Phase 2 — none of them branch on decisions
made here, and all of them presuppose a package that passes `R CMD check`.

## Goal

A JOSS reviewer can clone the repository, run `devtools::test()`, and see the
suite pass; `R CMD check --as-cran` reports 0 errors and 0 warnings; CI
demonstrates both on every push.

## Non-goals

- Changes to `paper.md` or the JSS manuscript (Phase 2)
- CRAN submission itself (Phase 2; this phase only reaches the quality bar)
- Zenodo DOI and tagged release (Phase 2)
- New statistical methodology or new metrics
- Contributor infrastructure (CONTRIBUTING.md, pkgdown) — deferred, not rejected

## Decisions taken

| Decision | Choice | Rationale |
|---|---|---|
| Target venue | Resubmit to JOSS | Rejection was software-only; manuscript is sound |
| API naming | `tm_` prefix, old names deprecated | Fixes S3-method check warning, typos, and PAmeasure optics without breaking existing analyses |
| Test depth | All four categories incl. reference comparison | Numerical correctness is the package's entire value proposition |
| Sequencing | Package first, paper later | Package work is a prerequisite and does not depend on downstream choices |

---

## 1. Repository hygiene

### Committed key pair — assessed, not an incident

`q` (ed25519 private key, passphrase-protected) and `q.pub` have been committed
since `c873cb9` (2024-11-03). The maintainer has confirmed this pair was
generated locally and **never** registered as a GitHub account key, deploy key,
or in any server's `authorized_keys`. There is therefore no trust relationship to
withdraw and no key requiring revocation. Removal is repository hygiene, not
incident response.

### Files to remove from the working tree

`q`, `q.pub`, `git`, `R/git`, `R/.DS_Store`, `inst/.DS_Store`, `..pdf`, `MD5`,
`demo for cr.R`, `.Rhistory (MacBook-Air.lan's conflicted copy 2025-11-08)`.

`R/` must contain only R source — non-source files there cause `R CMD check` to
fail outright.

### Ignore rules

Extend `.gitignore` to cover these classes (`*.pub`, `.DS_Store` at any depth,
conflicted copies, stray PDFs). Add a committed `.Rbuildignore` — it is currently
listed *in* `.gitignore`, which is backwards; `.Rbuildignore` must be tracked so
it takes effect for other people's builds. It should exclude `paper.md`,
`paper.code.Rmd`, `docs/`, `.github/`, and `TimeMetric.Rproj` from the built
tarball.

### History rewrite — final step of this phase

Purge `q` and `q.pub` from all 57 commits with `git filter-repo`
(`brew install git-filter-repo`), then force-push.

Sequenced **last**, after all other Phase 1 work is committed, so the rewrite and
force-push happen exactly once. The three contributors (Tong Zhu, Zian Zhuang,
WangHD) must re-clone afterward; there are no external forks. GitHub retains
unreferenced blobs by SHA until it garbage-collects — immaterial here given
nothing requires revocation, but GitHub Support can purge on request.

## 2. Metadata correctness

`R CMD check` fails today independently of tests, because `DESCRIPTION` and
`NAMESPACE` disagree about dependencies:

- `NAMESPACE` imports from `magrittr`, `pec`, `purrr`, `randomForestSRC`,
  `tibble` — none declared in `DESCRIPTION`
- `DESCRIPTION` declares `tdROC`, `yardstick` — verify against actual usage

**Method:** scan every `R/*.R` for `pkg::fn` calls and roxygen `@importFrom`
tags, derive the true dependency set, and reconcile both files against it.

**Placement:** heavy or optional backends (`pec`, `randomForestSRC`) move to
`Suggests`, guarded at call sites with `requireNamespace(..., quietly = TRUE)`
and an informative error. Only genuinely required packages stay in `Imports`.

Bump `Version` to `0.2.0`. Convert `Author`/`Maintainer` to `Authors@R` with
ORCIDs, as CRAN prefers and JOSS expects for author identification.

## 3. API rename

R parses `pam.foo` as an S3 method for class `foo` on generic `pam`, so
`R CMD check` emits "apparent S3 methods exported but not registered" — a CRAN
blocker. The prefix also foregrounds the PAmeasure derivation the editor asked
about. Renaming resolves both, and the `survial` / `surverg` typos, in one pass.

| Current export | New export |
|---|---|
| `pam.predicted_survial_eval` | `tm_survival_eval()` |
| `pam.predicted_survial_eval_cr` | `tm_survival_eval_cr()` |
| `pam.predicted_survial_eval_two_phase` | `tm_survival_eval_two_phase()` *(newly exported)* |
| `pam.survival_eval` | `tm_fit_and_eval()` *(pending §3.1)* |
| `pam.coxph_restricted` | `tm_predict_coxph()` |
| `pam.surverg_restricted` | `tm_predict_survreg()` |
| `pam.predict_cr` | `tm_predict_cif()` |
| `pam.summary` | `tm_summarize()` |
| `pam.summary_cr` | `tm_summarize_cr()` |
| `pam.sample_design` | `tm_sample_design()` |
| `cc_weights` | `tm_cc_weights()` |
| `ncc_weights` | `tm_ncc_weights()` |
| `plot_pred` | `tm_plot_pred()` |
| `summary_pred_plot` | `tm_plot_summary()` |
| `sim_cox_weibull_censored` | `tm_sim_cox_weibull()` |
| `simulateTwoCauseFineGrayModel` | `tm_sim_fine_gray()` |

Every old name is retained as a thin wrapper that calls `.Deprecated()` and
forwards to the new one, collected in `R/deprecated.R`. Existing analysis
scripts — including `paper.code.Rmd` and the Zhuang et al. (2025) reproduction —
continue to work with a warning. Source files rename to match their primary
function; roxygen is regenerated.

### 3.1 Open question: three overlapping entry points

`pam.survival_eval(train_data, covariates, models, ...)` fits models internally
then evaluates. `pam.predicted_survial_eval(model, event_time,
predicted_probability, status, ...)` accepts predictions. The unexported
`pam.prediction_survial_eval(object, train_data, predicted_data, ...)` accepts a
fitted object.

A reviewer will ask what distinguishes these. **Recommendation:** keep the
prediction-accepting function as the documented primary interface (it matches the
"Performance Metric Module" the paper describes), keep the fit-and-evaluate
convenience wrapper, and delete `pam.prediction_survial_eval` as superseded — it
is unexported and undocumented. To be confirmed against the maintainer's intent
during implementation.

## 4. Test suite

`tests/testthat/`, one file per metric family, roughly 40–60 tests. Shared
fixtures in `helper-simdata.R` with fixed seeds so every test is deterministic.

**Category 1 — known answers.** Cases with provable values, no reference package
needed: Harrell's C is 1 for a perfectly ordering predictor and 0.5 for a
constant one; Brier score is 0 when predicted survival equals truth; pseudo `R^2`
is 0 for a signal-free predictor.

**Category 2 — reference comparison.** Headline metrics cross-checked against
independent implementations on identical data, within tolerance: Harrell's C vs
`survival::concordance`, Uno's C vs `survC1`, Brier vs `pec` / `SurvMetrics`,
time-dependent AUC vs `timeROC`. This is the most persuasive evidence available
for a metrics package, and it transfers credibility to the metrics that have no
reference implementation (`R_sh`, `R_E`, pseudo `R^2`, and the case-cohort / NCC
variants).

**Category 3 — edge cases.** All observations censored; ties in event times; a
single distinct time point; `t_star` beyond the last event; `tau` defaulting.

**Category 4 — input validation.** `expect_error()` with informative messages on
mismatched vector lengths, `status` outside `{0,1,2}`, missing `time`/`status`
columns, predictions whose dimensions disagree with the data.

Reference packages go in `Suggests`, each test guarded by
`skip_if_not_installed()` so the suite degrades gracefully rather than erroring
where they are absent. Target runtime under 60 seconds.

**Expectation:** reference comparisons are likely to surface genuine numerical
bugs. That is the purpose. Any discrepancy is investigated and either fixed or
documented as a deliberate definitional difference.

## 5. Continuous integration

Replace `main.yml`, which only builds the JOSS draft PDF and runs no checks:

- `R-CMD-check.yaml` — r-lib standard matrix, Ubuntu / macOS / Windows against
  R release and devel
- `test-coverage.yaml` — covr + codecov, badge in README
- `draft-pdf.yaml` — the existing JOSS paper build, preserved unchanged

## 6. Documentation truthfulness

The README currently documents functions that do not exist:
`pam.predicted_survial_eval_casecohort`, `pam.predicted_survial_eval_ncc`, and
`pam.predict_subject_cif`. The case-cohort and NCC functionality is real but
lives in `pam.predicted_survial_eval_two_phase` (unexported) plus
`pam.sample_design` / `cc_weights` / `ncc_weights`; the CIF function is
`pam.predict_cr`.

Rewrite the README against the actual `tm_`-prefixed export list, with a runnable
quick-start example. Every exported function gets a working `@examples` block.

## Exit criteria

1. `R CMD check --as-cran` — 0 errors, 0 warnings
2. `devtools::test()` — all tests pass, under 60 seconds
3. CI green on all platforms in the matrix
4. README and all `.Rd` files reference only functions that exist
5. No non-source files in `R/`; no key material in the working tree or history
6. Deprecated aliases verified: `paper.code.Rmd` still runs

## Phase 2 preview (not in scope)

Manuscript revision, Research Impact section rewrite, `R CMD check` clean for
CRAN submission, tagged release, Zenodo DOI, response letter to the JOSS editor,
and contributor infrastructure.
