# Commit Review — 32 commits on `joss-revision`

**Date:** 2026-09-08
**Purpose:** maintainer review before any destructive history operation.
**Status:** nothing pushed; no rewrite run.

Only **10 of 32** commits touch `R/` and can change package behaviour. The other
22 are specs, tests, and documentation. Review effort is best spent on the ten.

---

## A. Behaviour-changing commits — review these closely

### `7814f65` — reconcile dependencies, unblock installation
`3 R/ files · 12 files · +287 −73`

* Deletes `importFrom(randomForestSRC, predict.rfsrc)`. **This one line was why
  `R CMD INSTALL` failed**, so the package could not be installed at all.
* Adds `importFrom(survival, concordancefit)` — the primary evaluation function
  previously failed under `library(TimeMetric)` alone.
* Adds every undeclared runtime dependency to `Imports`.
* Adds `requireNamespace()` guards for the two optional backends.

**To check:** that `Imports` now lists nothing unused, and that moving
`randomForestSRC`/`cmprsk` to `Suggests` matches your intent for those backends.

### `b6a2f61` — repair roxygen examples
`14 R/ files · +184 −200`

Two example blocks were invalid R, so `R CMD check` refused to run *any*
example. Rewrites both, plus two that called `library(tidyverse)`, and strips
10 stale `library(PAmeasures)` calls.

**To check:** the rewritten examples are scientifically sensible, not merely
syntactically valid. They now execute on every check.

### `cdb0981` — S3 registration, `pam.survival_eval`, CR `Value` type
`3 R/ files · +116 −60`

* Registers the `pam.rsph` methods so `R_E` dispatches from a clean session.
* Fixes `pam.survival_eval`'s two call sites (`covariates=`/`newdata=` →
  `covs=`/`new_data=`, plus `predict = FALSE`). **This exported function had
  never run.**
* Makes the competing-risks `Value` column numeric rather than character.

**To check:** `predict = FALSE` is the branch returning `R.squared`/`L.squared`,
which is what the caller indexes as `r_l_list[1]`/`[2]`. Confirm that reading.

### `df546b6` — clear all check errors and warnings
`20 R/ files · 51 files · +332 −782`

The largest behavioural commit. Fixes `plot_pred`'s blank-by-default plot,
removes three `require(survival)` calls, adds `...` to S3 methods, converts all
non-ASCII in R sources to `\uXXXX` escapes, corrects five `@param`/signature
mismatches, fixes the LICENSE stub and documents the `moore` dataset. Also
deletes `q`, `q.pub` and eight other stray files from the working tree.

**To check:** the non-ASCII conversion preserved metric-name *values* exactly —
`"Harrell’s C"` is byte-identical to the original literal, so metric
matching is unchanged.

### `598bb9b` — delete Cluster B
`3 R/ files · −433`

Removes `pam.prediction_survial_eval`, `pam.prediction_metrics`,
`pam.prediction_metrics_cr`. All unreachable; they called functions that exist
nowhere in the package.

### `64e0f46` — delete Cluster C, pin `R_E`
`3 R/ files · +192 −449`

Removes `pam.Brier_metric`, `pam.rsh_metric`, `pam.rsph_metric` after the
equivalence gate and the `R_E` audit. Adds the reference regression test.

**To check:** this is the commit resting on the audit. `pam.rsph` reproduces the
authors' `Re.r` exactly; `pam.rsph_metric` was wrong by 1.3–5.3% through an
omitted inverse-censoring weighting. Confirm you accept that conclusion.

### `096b2fd` — delete Cluster A, restore `print.rsph`
`4 R/ files · +16 −317`

Deletes `pam.coxph`, `pam.nlm`, `pam.survreg`. **Does not delete**
`pam.print.rsph` — it is restored as `print.rsph`, a working S3 method. The
rename had broken it; deleting would have discarded a feature.

### `49366e9` — the `tm_` rename
`9 R/ files · 54 files · +623 −349`

The largest surface change. 16 functions renamed, 15 deprecated wrappers added,
`tm_evaluate_two_phase` newly exported, `pam.summary.rsph` → `summary.rsph`.

**To check:** the chosen names, and that the deprecation policy (a wrapper for
every old export) matches what you want to support.

### `a8642ab` — metric name standardisation
`5 R/ files · +427 −263`

Canonical ASCII metric identifiers with a legacy-tolerant normalisation layer.

**To check:** metric **values** are unchanged — every `snap_num(...$Value)`
snapshot passed untouched; only `Metric` label snapshots moved. Also confirm
`pseudo_r2` and `pseudo_r2_point` should stay distinct (they differ
numerically: 0.3868 vs 0.1429).

### `ef510c0` — ASCII regression fix
`1 R/ file · +3 −3`

My own regression: the alias map used literal curly apostrophes as keys,
reintroducing the non-ASCII WARNING. Fixed with escapes.

---

## B. Specification and planning — 6 commits

`fd85785`, `dc8e4f4`, `4ff2632`, `0de45d9`, `a424f4e`, `45db525`

The design spec and implementation plan, each revised after your review. Docs
only. Worth skimming for whether the recorded decisions match your intent.

## C. Characterization tests — 11 commits

`21188a8`, `51d204d`, `4a8e6cc`, `33b9b23`, `ff6fe50`, `4f4ffcf`, `86e702a`,
`3406455`, `2eb0942`, `637f2dd`, `3e2c2d3`

The test baseline, written before any source change so later edits were provably
behaviour-preserving. `tests/` and `docs/` only.

## D. Analysis and documentation — 5 commits

`daa025b` (Cluster C equivalence gate), `3ca99f4` (metric identification),
`3a1e7f4` (`R_E` audit), `cd7b865` (README + CI), `a0cc99f` (0.2.0, NEWS,
pre-rewrite report). No `R/` changes.

---

## Suggested review commands

```bash
# the ten behavioural commits, source changes only
git log main..joss-revision --reverse -p -- R/ | less

# any single commit in full
git show 64e0f46

# net effect on the package source, ignoring tests and docs
git diff main joss-revision -- R/ DESCRIPTION NAMESPACE | less

# confirm no test was weakened rather than updated
git diff main joss-revision -- tests/ | less
```
