# Pre-Rewrite Review Report

**Date:** 2026-09-07
**Branch:** `joss-revision` (31 commits ahead of `main`)
**Status:** awaiting approval. **No history rewrite or force-push has been performed.**

---

## 1. Security scan across all refs

The repository has **4 refs and 0 tags**:

| Ref | Type |
|---|---|
| `refs/heads/main` | local branch |
| `refs/heads/joss-revision` | local branch (this work) |
| `refs/remotes/origin/main` | remote-tracking |
| `refs/remotes/origin/HEAD` | symbolic, to `origin/main` |

Two commits touch the key pair:

| Commit | Date | Action |
|---|---|---|
| `c873cb9` | 2024-11-03 | **added** `q` and `q.pub` ("original version") |
| `df546b6` | 2026-09-02 | **removed** them from the working tree (this branch) |

**Refs containing the introducing commit `c873cb9`:**

| Ref | Contains `c873cb9` | Key files at tip | Commits carrying `q` |
|---|---|---|---|
| `refs/heads/main` | **yes** | **2 of 2 present** | **55 of 57** |
| `refs/remotes/origin/main` | **yes** | **2 of 2 present** | **55 of 57** |
| `refs/heads/joss-revision` | yes | 0 of 2 (removed at `df546b6`) | 75 of 88 |

### Consequences for the rewrite

* **`main` is affected and must be rewritten too.** The keys are not merely
  historical there -- they are live in `main`'s current tree, and `origin/main`
  is the published branch. Rewriting only `joss-revision` would leave the key
  fully exposed.
* `git filter-repo` must run over **all refs**, not one branch.
* Both branches will need a force-push, and `origin/HEAD` follows `origin/main`.
* No tags exist, so no tag rewriting or re-signing is required.

### Other secret-like blobs

A scan of every blob in every reachable object found exactly two matches, `q`
and `q.pub`. No `.pem`, `.key`, `id_rsa`, or similar exists anywhere in history.

### Standing assessment

Per the maintainer's determination (findings.md #6 context): this key pair was
generated locally and never registered as a GitHub account key, deploy key, or
in any server's `authorized_keys`. There is no trust relationship to withdraw
and no key requiring revocation. The rewrite is hygiene, not incident response.

---

## 2. Verification results

| Check | Result |
|---|---|
| `R CMD check --as-cran` (on `R CMD build` tarball) | **0 errors, 0 warnings, 1 NOTE** |
| The remaining NOTE | `New submission` -- inherent, not actionable |
| Test suite | **380 passing, 0 failed, 0 skipped** |
| Coverage | **79.50%**, 0 files at 0%, 3 files at 100% |
| Export-coverage gate | every export referenced by at least one test |
| `R CMD INSTALL` | succeeds (failed outright on `main`) |

Baseline for comparison, measured at the start of this work:

| | `main` | `joss-revision` |
|---|---|---|
| `R CMD INSTALL` | **fails** | succeeds |
| `R CMD check` | 1 ERROR, aborts at dependency stage | 0E / 0W / 1N |
| Tests | none | 380 |
| Coverage | n/a (package would not install) | 79.50% |

---

## 3. Diff summary against `main`

```
108 files changed, 6081 insertions(+), 2402 deletions(-)
```

| Kind | Count |
|---|---|
| Added | 45 |
| Deleted | 23 |
| Modified | 14 |
| Renamed | 26 |

---

## 4. Exported API: before and after

**Before: 15 exports. After: 32** (16 renamed functions + `tm_metric_names()` +
15 deprecated aliases). **No function was removed from the public API.**

| Old name (still works, deprecated) | New name |
|---|---|
| `pam.predicted_survial_eval` | `tm_survival_eval` |
| `pam.predicted_survial_eval_cr` | `tm_survival_eval_cr` |
| `pam.survival_eval` | `tm_fit_and_eval` |
| `pam.coxph_restricted` | `tm_predict_coxph` |
| `pam.surverg_restricted` | `tm_predict_survreg` |
| `pam.predict_cr` | `tm_predict_cif` |
| `pam.summary` | `tm_summarize` |
| `pam.summary_cr` | `tm_summarize_cr` |
| `pam.sample_design` | `tm_sample_design` |
| `cc_weights` | `tm_case_cohort_weights` |
| `ncc_weights` | `tm_nested_case_control_weights` |
| `plot_pred` | `tm_plot_pred` |
| `summary_pred_plot` | `tm_plot_summary` |
| `sim_cox_weibull_censored` | `tm_sim_cox_weibull` |
| `simulateTwoCauseFineGrayModel` | `tm_simulate_fine_gray` |

Newly public, with no predecessor:

* **`tm_evaluate_two_phase()`** -- previously internal, so the case-cohort and
  nested case-control functionality the manuscript advertises was unreachable.
* **`tm_metric_names()`** -- the canonical metric identifiers.

All 15 old names are retained in `R/deprecated.R`, each calling `.Deprecated()`
with its replacement and forwarding arguments unchanged. Verified: all 9
TimeMetric functions called by `paper.code.Rmd` still resolve.

### S3 registrations

**Before: 0. After: 5.**

```
S3method(pam.rsph,aareg)     S3method(print,rsph)
S3method(pam.rsph,coxph)     S3method(summary,rsph)
S3method(pam.rsph,survreg)
```

Previously none were registered, so `R_E` could not be dispatched from a clean
session and `print()`/`summary()` on the returned object did not work.

---

## 5. Dependencies: before and after

| | `main` | `joss-revision` |
|---|---|---|
| **Imports** | survival, stats, tdROC, yardstick, ggplot2, patchwork, dplyr | survival, stats, **utils**, tdROC, yardstick, ggplot2, patchwork, dplyr, **magrittr**, **purrr**, **tibble**, **pec**, **rms**, **survminer**, **expint** |
| **Suggests** | testthat | testthat (>= 3.0.0), **withr**, **pkgload**, **covr**, **randomForestSRC**, **cmprsk** |
| **NAMESPACE import of `randomForestSRC`** | present, undeclared -- **blocked installation** | removed |

Eight packages were used at runtime but declared nowhere. `randomForestSRC` and
`cmprsk` are optional `tm_predict_cif()` backends, guarded by
`requireNamespace()` with an error naming the package to install.

---

## 6. Deleted files

**Source -- 9 unreachable functions**, each verified by call graph and, for the
metric implementations, by equivalence audit:

```
R/pam.prediction_survial_eval.R   R/pam.rsh_metric.R
R/pam.prediction_metrics.R        R/pam.rsph_metric.R
R/pam.prediction_metrics_cr.R     R/pam.Brier_metric.R
R/pam.coxph.R                     R/pam.nlm.R
R/pam.survreg.R
```

**Repository hygiene:**

```
q                 q.pub             MD5               ..pdf
git               R/git             R/.DS_Store       inst/.DS_Store
demo for cr.R     .Rhistory (MacBook-Air.lan's conflicted copy 2025-11-08)
```

**Stale documentation** (regenerated or orphaned): `man/cc_weights.Rd`,
`man/find_mu_c.Rd`, `man/integrate_survival.Rd`, `man/pam.survival_eval.Rd`.

---

## 7. Commit list

 1. fd85785  Add Phase 1 design spec for JOSS revision
 2. dc8e4f4  Update Phase 1 spec with dead-code and dependency analysis
 3. 4ff2632  Correct Phase 1 spec after maintainer review
 4. 0de45d9  Apply maintainer review: test-first ordering and S3 correction
 5. a424f4e  Add implementation plan for characterization test baseline
 6. 45db525  Correct characterization plan after maintainer review
 7. 21188a8  test: add testthat infrastructure and findings log
 8. 51d204d  test: add deterministic fixtures and expectation helpers
 9. 4a8e6cc  test: characterize Cluster C metric functions before deletion decision
10. 33b9b23  test: characterize Gt/pam.Brier/predictSurvProb chain gating pec removal
11. ff6fe50  test: characterize coxph, survreg and competing-risks prediction functions
12. 4f4ffcf  test: characterize right-censored evaluation, summary and metric names
13. 86e702a  test: characterize competing-risks evaluation and summary
14. 3406455  test: scope survival attachment and reproduce FINDING 11 in a clean session
15. 2eb0942  test: characterize case-cohort and NCC two-phase evaluation
16. 637f2dd  test: characterize R_E dispatch and plotting
17. 3e2c2d3  test: record baseline coverage and check output before remediation
18. 7814f65  fix: reconcile dependencies and unblock installation (findings 6, 11)
19. b6a2f61  fix: repair roxygen examples so R CMD check can run them (finding 19)
20. cdb0981  fix: register S3 methods, repair pam.survival_eval, make CR Value numeric
21. df546b6  fix: clear all R CMD check errors and warnings
22. 598bb9b  refactor: delete Cluster B dead code
23. daa025b  docs: record the Cluster C equivalence gate
24. 3ca99f4  docs: identify which manuscript metrics Cluster C implements
25. 3a1e7f4  docs: audit the two R_E implementations against the authors' reference
26. 64e0f46  refactor: delete the three Cluster C duplicates, pin R_E to the reference
27. 096b2fd  refactor: delete Cluster A dead code, restore print.rsph as a real S3 method
28. 49366e9  refactor: rename the public API to a tm_ prefix with deprecated aliases
29. a8642ab  refactor: standardize metric names to an ASCII canonical set
30. ef510c0  fix: keep R sources ASCII-only after metric standardization
31. cd7b865  docs: rewrite README against the real API and add CI workflows

---

## 8. What has NOT been done

* **No `git filter-repo` has been run.**
* **No force-push. Nothing has been pushed at all.** `main` and `origin/main`
  are byte-for-byte untouched.
* The working branch `joss-revision` exists only locally.

## 9. Proposed rewrite plan, for approval

1. `brew install git-filter-repo`.
2. Clone a fresh mirror as a safety copy, kept until the rewrite is confirmed
   good.
3. Run, over **all refs**, covering `main` and `joss-revision`:

   ```
   git filter-repo --invert-paths --sensitive-data-removal \
     --path q --path q.pub \
     --path docs/security/key-exposure-report.md
   ```

   The third path is included because the detailed security report was moved out
   of the repository at the maintainer's request and is retained only as a
   private local record. It remains in commit `7eab69b` and must not reach the
   published repository.
4. Verify: no blob named `q` or `q.pub` in any object; the working trees of both
   branches are otherwise unchanged; the suite still passes; `R CMD check` still
   clean.
5. Force-push both branches. `origin/HEAD` follows `origin/main`.
6. All three contributors (Tong Zhu, Zian Zhuang, WangHD) must **re-clone** --
   a pull will not reconcile rewritten history.
7. Optionally ask GitHub Support to garbage-collect unreferenced objects, since
   the old commits remain reachable by SHA until they do.

**This plan will not be executed without explicit approval.**
