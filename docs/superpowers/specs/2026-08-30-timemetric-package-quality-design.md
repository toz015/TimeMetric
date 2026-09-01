# TimeMetric — Phase 1: Package Quality Remediation

**Date:** 2026-08-30
**Status:** Design approved; dependency and dead-code analysis complete
**Context:** JOSS desk rejection (editor: arfon, 2026-07-16) citing software readiness

## Background

`TimeMetric` was submitted to JOSS and desk-rejected before peer review on three
grounds: no automated tests, a repository not in reviewable state for an R
package, and thin evidence of research impact / open development.

The rejection concerned the *software*, not the manuscript. This spec covers only
software remediation. Manuscript revision, CRAN submission, and the Zenodo
archive are deferred to Phase 2 — none branch on decisions made here, and all
presuppose a package that passes `R CMD check`.

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
| API naming | `tm_` prefix, selective deprecation | Fixes S3-method check warning, typos, and PAmeasure optics |
| Test depth | All four categories incl. reference comparison | Numerical correctness is the package's entire value proposition |
| Sequencing | Package first, paper later | Package work is a prerequisite and independent of downstream choices |
| Encoding | UTF-8 files, ASCII-only source | Files already UTF-8; non-ASCII *string literals* must go (see §5) |
| Metric names | ASCII canonical set, legacy accepted with warning | Same metric currently has up to 4 spellings (see §5) |

---

## 1. Repository hygiene

### Committed key pair — assessed, not an incident

`q` (ed25519 private key, passphrase-protected) and `q.pub` have been committed
since `c873cb9` (2024-11-03). The maintainer has confirmed this pair was
generated locally and **never** registered as a GitHub account key, deploy key,
or in any server's `authorized_keys`. There is no trust relationship to withdraw
and no key requiring revocation. Removal is repository hygiene, not incident
response.

### Files to remove from the working tree

`q`, `q.pub`, `git`, `R/git`, `R/.DS_Store`, `inst/.DS_Store`, `..pdf`, `MD5`,
`demo for cr.R`, `.Rhistory (MacBook-Air.lan's conflicted copy 2025-11-08)`.

`R/` must contain only R source — non-source files there cause `R CMD check` to
fail outright.

### Ignore rules

Extend `.gitignore` to cover these classes (`*.pub`, `.DS_Store` at any depth,
conflicted copies, stray PDFs). Add a **committed** `.Rbuildignore` — it is
currently listed *in* `.gitignore`, which is backwards; it must be tracked to
take effect for other people's builds. It should exclude `paper.md`,
`paper.code.Rmd`, `docs/`, `.github/`, and `TimeMetric.Rproj` from the tarball.

### History rewrite — final step of this phase

Purge `q` and `q.pub` from all 57 commits with `git filter-repo`
(`brew install git-filter-repo`), then force-push.

Sequenced **last**, after all other Phase 1 work is committed, so the rewrite and
force-push happen exactly once. The three contributors (Tong Zhu, Zian Zhuang,
WangHD) must re-clone afterward; there are no external forks. GitHub retains
unreferenced blobs by SHA until it garbage-collects — immaterial here, but GitHub
Support can purge on request.

## 2. Dead code removal

Reachability computed as a call graph rooted at the 15 `NAMESPACE` exports plus
`UseMethod("pam.rsph")` dispatch targets. **Ten functions are unreachable**, in
three clusters:

**Cluster A — standalone orphans** (no reference anywhere beyond their own
definition):

| Function | File | Action |
|---|---|---|
| `pam.coxph` | `R/pam.coxph.R` | delete file |
| `pam.nlm` | `R/pam.nlm.R` | delete file |
| `pam.survreg` | `R/pam.survreg.R` | delete file |
| `pam.print.rsph` | `R/pam.rsph.R` | delete function; would dispatch on a `pam.print` generic that does not exist |

**Cluster B — the legacy `prediction_*` chain** (unexported, referenced only from
each other or from roxygen `@examples`):

| Function | File | Action |
|---|---|---|
| `pam.prediction_survial_eval` | `R/pam.prediction_survial_eval.R` | delete file |
| `pam.prediction_metrics` | `R/pam.prediction_metrics.R` | delete file |
| `pam.prediction_metrics_cr` | `R/pam.prediction_metrics_cr.R` | delete file |

**Cluster C — transitively dead**, reachable only from Cluster B:

| Function | File | Note |
|---|---|---|
| `pam.rsh_metric` | `R/pam.rsh_metric.R` | sole caller is `pam.prediction_metrics` |
| `pam.rsph_metric` | `R/pam.rsph_metric.R` | sole caller is `pam.prediction_metrics` |
| `pam.Brier_metric` | `R/pam.Brier_metric.R` | sole caller is `pam.prediction_metrics` |

Cluster C are self-contained metric implementations that take vectors rather than
fitted models. They are dead as the code stands, but before deleting, confirm they
are not a better foundation for `tm_evaluate_two_phase` than the current
model-coupled path. **Delete only after that check.**

### Explicitly NOT dead — verified live

Name-based searching produces false positives for functions passed as *values*
rather than called. Each of these is live:

- `find_mu_c` — `uniroot(find_mu_c, …)`, `simulateTwoCauseFineGrayModel.R:116`
- `handle_error` — `error = handle_error`, `pam.predict_cr.R:120,127`
- `m_cif` — `apply(…, 2, m_cif, …)`, five live call sites
- `integrate_survival` — three live call sites
- `pam.rsph.aareg`, `pam.rsph.coxph`, `pam.rsph.survreg` — reached via
  `UseMethod("pam.rsph")` at `R/pam.rsph.R:74`. These serve the **`R_E`** metric
  (`pam.summary.rsph(pam.rsph(…))`, `pam.predicted_survial_eval.R:226`), *not*
  `R_sh`. Deleting them would remove `R_E` entirely.

Any further pruning must re-run the call graph including function-as-value
references, not grep for `name(`.

### Orphan documentation

`man/find_mu_c.Rd` and `man/integrate_survival.Rd` document unexported internal
helpers, which produces check warnings. Mark both `@keywords internal` / `@noRd`
and drop the generated pages. `man/pam.predicted_survial_eval_two_phase.Rd`
becomes valid once that function is exported (§4).

## 3. Dependency reconciliation

`R CMD check` fails today independently of tests. Determined by scanning both
qualified (`pkg::fn`) and **unqualified** call sites — the latter matter because
`importFrom` makes functions available bare.

### Keep — `pec`

`pec` is a live runtime dependency, not a stale import. `R/pam.Brier.R:66` calls
`predictSurvProb(obj, test_data, t_star0)` **unqualified**, resolved through
`importFrom(pec, predictSurvProb)`. The chain is:

```
tm_fit_and_eval()  ->  pam.Brier()  ->  predictSurvProb()   [pec]
```

`pec` stays in `Imports` until that legacy chain is rewritten or removed. Only
then does it move to `Suggests` as a test reference implementation. **Do not
remove it in this phase unless the chain is rewritten first.**

### Add to Imports — used but undeclared

| Package | Call sites | Why required |
|---|---|---|
| `rms` | `rms::cph` (3×), `rms::predictrms` (1×) | `R_sh` is in the **default** `metrics` vector (`pam.predicted_survial_eval.R:82`), so this is not an optional path |
| `expint` | `expint::gammainc` (`pam.surverg_restricted.R:136`) | required while `tm_predict_survreg()` uses it |

### Remove — genuinely stale

`importFrom(randomForestSRC, predict.rfsrc)` has zero call sites, qualified or
unqualified. Delete the import. Add to `Suggests` only if a test needs it as a
prediction backend.

### Conditional — `survminer`

One call: `survminer::surv_summary()` at `R/pam.Ct.R:63`, inside the legacy
`Gt()` / Brier chain. It is a heavy dependency for a single call.

- If the legacy Brier/`Gt()` chain is removed, `survminer` goes with it.
- Otherwise, replace with base `summary(survfit(...))` **and add a test asserting
  numerical equivalence** before the swap is accepted. Do not swap unverified.

Until one of those happens, `survminer` must be declared in `Imports` — it is
currently used and undeclared.

### Verify and correct

`tdROC` and `yardstick` are genuinely used and correctly in `DESCRIPTION` but
lack `@importFrom` tags. `magrittr`, `purrr`, `tibble` are imported in `NAMESPACE`
but absent from `DESCRIPTION`. Roxygen regeneration plus a corrected `DESCRIPTION`
resolves both directions.

Bump `Version` to `0.2.0`. Convert `Author`/`Maintainer` to `Authors@R` with
ORCIDs, as CRAN prefers and JOSS expects.

## 4. API rename

R parses `pam.foo` as an S3 method for class `foo` on generic `pam`, so
`R CMD check` emits "apparent S3 methods exported but not registered" — a CRAN
blocker. The prefix also foregrounds the PAmeasure derivation the editor asked
about. Renaming resolves both, plus the `survial` / `surverg` typos, in one pass.

| Current export | New export |
|---|---|
| `pam.predicted_survial_eval` | `tm_survival_eval()` |
| `pam.predicted_survial_eval_cr` | `tm_survival_eval_cr()` |
| `pam.predicted_survial_eval_two_phase` | `tm_evaluate_two_phase()` *(newly exported)* |
| `pam.survival_eval` | `tm_fit_and_eval()` |
| `pam.coxph_restricted` | `tm_predict_coxph()` |
| `pam.surverg_restricted` | `tm_predict_survreg()` |
| `pam.predict_cr` | `tm_predict_cif()` |
| `pam.summary` | `tm_summarize()` |
| `pam.summary_cr` | `tm_summarize_cr()` |
| `pam.sample_design` | `tm_sample_design()` |
| `cc_weights` | `tm_case_cohort_weights()` |
| `ncc_weights` | `tm_nested_case_control_weights()` |
| `plot_pred` | `tm_plot_pred()` |
| `summary_pred_plot` | `tm_plot_summary()` |
| `sim_cox_weibull_censored` | `tm_sim_cox_weibull()` |
| `simulateTwoCauseFineGrayModel` | `tm_simulate_fine_gray()` |

### Deprecation policy

Wrappers calling `.Deprecated()` are kept **only** for old public functions that
are demonstrably in use or map unambiguously onto the new API. Usage measured
against `paper.code.Rmd`:

**Wrapper retained** — in active use: `pam.coxph_restricted` (16 calls),
`pam.surverg_restricted` (17), `pam.summary` (11), `sim_cox_weibull_censored` (8),
`summary_pred_plot` (5), `pam.predict_cr` (4), `cc_weights` (2),
`pam.sample_design` (1), `pam.summary_cr` (1).

**No wrapper** — exported but unused in any known analysis, and the rename is a
clean break: `ncc_weights`, `plot_pred`, `simulateTwoCauseFineGrayModel`,
`pam.survival_eval`, `pam.predicted_survial_eval`,
`pam.predicted_survial_eval_cr`.

Wrappers live in `R/deprecated.R`. Source files rename to match their primary
function; roxygen is regenerated.

### `tm_evaluate_two_phase` — promotion to public API

Currently `@keywords internal` and unexported, so the case-cohort and NCC
functionality the paper advertises is **not reachable by users**. This is a
genuine functional gap, not merely a documentation error. Export it, write full
roxygen documentation with runnable examples for *both* designs, and test both
paths (§5).

## 5. Metric name standardization

The same metric currently has up to four different string identifiers depending
on the entry point, so result tables carry different column names for the same
quantity:

| Canonical (new) | Legacy spellings found in `R/` |
|---|---|
| `pseudo_r2` | `Pseudo_R_square`, `Pesudo_R` |
| `r_square` | `R_square` |
| `l_square` | `L_square` |
| `harrell_c` | `Harrells_C`, `Harrell’s C` |
| `uno_c` | `Unos_C`, `Uno’s C`, `Uno's C` |
| `r_sh` | `R_sh` |
| `r_e` | `R_E`, `R_sph`, `rsph` |
| `brier_score` | `Brier Score`, `Brier_Score` |
| `td_auc` | `Time Dependent Auc`, `Time_Dependent_Auc`, `Time Dependent AUC`, `AUC` |

Adopt the ASCII canonical set throughout. Legacy names remain accepted at the
argument boundary via a normalisation lookup that emits a deprecation warning
naming the replacement, so existing scripts keep working for one release cycle.

Two defects this fixes, both of which currently force users to type something
wrong to get a result:

1. **`Pesudo_R`** — misspelling of "Pseudo", 13 occurrences, present in default
   `metrics =` arguments. A user asking for `Pseudo_R` gets nothing back.
2. **Curly apostrophes** — `"Harrell’s C"` / `"Uno’s C"` use U+2019
   (`pam.predicted_survial_eval_two_phase.R:45,177`). Selection only matches if
   the user types a typographic apostrophe. Also a CRAN "found non-ASCII strings"
   NOTE.

Normalisation is case-insensitive and treats `_`, space, `’` and `'` as
equivalent, so all historical spellings resolve.

**To confirm during implementation:** `Pseudo_R2_point` and
`Time Dependent Auc Empirical` may be genuinely distinct variants (point estimate
vs. curve) rather than spelling drift. Verify before folding them into the
canonical names above.

All source files are already UTF-8; this concerns non-ASCII *content* in code,
not file encoding. New and edited files stay UTF-8 with ASCII-only source.

## 6. Test suite

`tests/testthat/`, one file per metric family, roughly 40–60 tests. Shared
fixtures in `helper-simdata.R` with fixed seeds so every test is deterministic.

**Category 1 — known answers.** Provable values, no reference package needed:
Harrell's C is 1 for a perfectly ordering predictor and 0.5 for a constant one;
Brier score is 0 when predicted survival equals truth; pseudo `R^2` is 0 for a
signal-free predictor.

**Category 2 — reference comparison.** Headline metrics cross-checked against
independent implementations on identical data, within tolerance: Harrell's C vs
`survival::concordance`, Uno's C vs `survC1`, Brier vs `pec` / `SurvMetrics`,
time-dependent AUC vs `timeROC`. Most persuasive evidence available for a metrics
package, and it transfers credibility to metrics with no reference
implementation (`R_sh`, `R_E`, pseudo `R^2`, case-cohort / NCC variants).

**Category 3 — edge cases.** All observations censored; ties in event times; a
single distinct time point; `t_star` beyond the last event; `tau` defaulting.

**Category 4 — input validation.** `expect_error()` with informative messages on
mismatched vector lengths, `status` outside `{0,1,2}`, missing `time`/`status`
columns, prediction dimensions disagreeing with the data.

`tm_evaluate_two_phase` gets dedicated coverage for **both** the case-cohort and
nested case-control paths, including weight construction via
`tm_case_cohort_weights` and `tm_nested_case_control_weights`.

Reference packages go in `Suggests`, each test guarded by
`skip_if_not_installed()` so the suite degrades gracefully where they are absent.
Target runtime under 60 seconds.

**Expectation:** reference comparisons are likely to surface genuine numerical
bugs. That is the purpose. Any discrepancy is investigated and either fixed or
documented as a deliberate definitional difference.

## 7. Continuous integration

Replace `main.yml`, which only builds the JOSS draft PDF and runs no checks:

- `R-CMD-check.yaml` — r-lib standard matrix, Ubuntu / macOS / Windows against R
  release and devel
- `test-coverage.yaml` — covr + codecov, badge in README
- `draft-pdf.yaml` — the existing JOSS paper build, preserved unchanged

## 8. Documentation truthfulness

The README documents three functions that do not exist:
`pam.predicted_survial_eval_casecohort`, `pam.predicted_survial_eval_ncc`, and
`pam.predict_subject_cif`. The case-cohort and NCC functionality is real but
lives in the unexported `pam.predicted_survial_eval_two_phase` plus
`pam.sample_design` / `cc_weights` / `ncc_weights`; the CIF function is
`pam.predict_cr`.

Rewrite the README against the actual `tm_` export list, with a runnable
quick-start. Every exported function gets a working `@examples` block.

## Exit criteria

1. `R CMD check --as-cran` — 0 errors, 0 warnings
2. `devtools::test()` — all tests pass, under 60 seconds
3. CI green on all platforms in the matrix
4. README and all `.Rd` files reference only functions that exist
5. No non-source files in `R/`; no key material in the working tree or history
6. No non-ASCII characters in R source; all files UTF-8
7. Metric names canonical everywhere; every legacy spelling still resolves with a warning
8. Retained deprecated wrappers verified: `paper.code.Rmd` still runs

## Phase 2 preview (not in scope)

Manuscript revision, Research Impact section rewrite, CRAN submission, tagged
release, Zenodo DOI, response letter to the JOSS editor, and contributor
infrastructure.
