# Finding 41 — Impact Audit (analysis only)

**Date:** 2026-09-25
**Scope:** analysis only. No package source, test, snapshot, or baseline value
was modified. Corrected implementations were injected into the loaded namespace
in throwaway R sessions and reverted; nothing on disk changed.
**Status of the finding:** still OPEN. This audit informs a decision; it does not
take one.

## 1. The defect

`R/pam.Ct.R`, interpolation branch of `Gt()`. `index1` and `index2` are computed
from `res.sum$time` — the **sorted** censoring-KM summary table — but the
interpolation weights are then built from `time[index1]` and `time[index2]`,
where `time` is the **raw, unsorted input vector** `object[, 1]`. The two index
spaces are unrelated, so the weights come from arbitrary observations.

There are in fact **two** deviations from a correct G(t), and they should not be
conflated:

* **(a) the index misalignment** — an outright bug;
* **(b) linear interpolation of a step function** — Kaplan-Meier is a right-
  continuous step function, so interpolating between adjacent values is a
  modelling choice, not an estimate of G(t). This is present even after (a) is
  repaired.

## 2. Call graph

**Direct callers of `Gt()`** — all in `R/pam.Brier.R`:

| Site | Call | Input vector | Timepoint |
|---|---|---|---|
| `:159` | `Gtstar <- Gt(object, t_star)` | **UNSORTED** (`object` as received) | arbitrary |
| `:163` | `Gti <- Gt(Surv(time, status), time[i])` | sorted (`time` re-sorted at `:152`) | always an **observed** time |
| `:176` | `1 / Gt(Surv(time, status), t_star)` | sorted | arbitrary |

`pam.Brier()` sorts `time`, `status` and `pre_sp` at lines 151-154, so sites
`:163` and `:176` pass an already-sorted vector and the misalignment cannot bite
there. Site `:163` additionally always asks for an observed time, taking the
exact-match branch.

**Transitive callers of `pam.Brier()`**: `R/tm_fit_and_eval.R:137,140` only.

**Instrumented reachability** — `Gt()` was replaced with a counting wrapper and
every public entry point exercised:

| Entry point | `Gt()` calls |
|---|---|
| `tm_fit_and_eval()` | **181** |
| `tm_survival_eval()` | 0 |
| `tm_summarize()` | 0 |
| `tm_survival_eval_cr()` | 0 |
| `tm_summarize_cr()` | 0 |
| `tm_evaluate_two_phase()` | 0 |
| `tm_sample_design()` | 0 |

## 3. Which public outputs can be affected

**Exactly one: `brier_score` as returned by `tm_fit_and_eval()`.**

Nothing else. In particular `tm_survival_eval()`'s `brier_score` is computed by
`tdROC::tdROC()` and never touches `Gt()`, so the two entry points compute the
Brier score by different routes — itself worth noting, though out of scope here.

Exposure requires **both** conditions:

1. the evaluation time `t_star` is **not** an observed time, so the
   interpolation branch is taken; and
2. the censoring-KM has **distinct** survival values around `t_star`.

Condition 2 explains why the effect is invisible on the standard fixture:
`fx_surv` yields 200 KM steps but only **53 distinct survival values**, so most
interpolations sit on a flat stretch where wrong weights interpolate between two
equal numbers and cancel. `fx_surv_uncensored` has 200 distinct values and does
show the error.

**On the default `t_star` path, exposure is a matter of parity.**
`pam.Brier()` sets `t_star <- median(distime)` over the *distinct* event times.
An odd count yields an actual observed time (exact branch, dormant); an even
count yields the average of two adjacent times (interpolating, active). On
`fx_surv` there are **149** distinct event times — odd — which is the sole reason
all 181 calls take the exact branch and the committed baseline is unaffected.
That is luck, not design.

## 4. Row-permutation sensitivity

Identical data, rows permuted (`set.seed(99)`):

| Timepoint set | max abs difference |
|---|---|
| observed times | **0** |
| between observations | 1.11e-16 on `fx_surv`; up to **0.73** on tie-heavy data |
| before first / after last | **0** |

So `Gt()` is **not permutation-invariant** — a property any estimator of a
distribution function must have. The magnitude is entirely data-dependent.

No *public metric* changed under permutation on the committed fixtures, because
of the parity accident in §3.

## 5. Current implementation vs correctly aligned Kaplan-Meier

`cur` = current; `fix` = same algorithm with indices applied to the sorted table;
`G(t)` = right-continuous KM; `G(t-)` = left limit.

| Fixture | Set | max abs cur-fix | max abs cur-G(t) | max abs cur-G(t-) | monotone? |
|---|---|---|---|---|---|
| `fx_surv` | observed | 0 | 0.2385 | 0.0097 | yes |
| | between | 1.1e-16 | 1.1e-16 | 1.1e-16 | yes |
| | before/after | 0 | 0 / 0.2385 | 0 / 0.2385 | yes |
| `fx_surv_uncensored` | observed | 0 | 0.0050 | 0.0050 | yes |
| | between | **0.1594** | 0.1619 | 0.1619 | **NO** |
| tie-heavy 40/5 | observed | 0 | 0 | 0.2461 | yes |
| | between | **0.5977** | 0.6539 | 0.6539 | yes |
| tie-heavy 200/10 | observed | 0 | 0 | 0.1757 | yes |
| | between | **0.7331** | 0.7602 | 0.7602 | yes |

Readings:

* **At observed times the current code is already exact** (`cur-fix = 0`): it
  takes the exact-match branch. The `cur-G(t)` gap of 0.2385 on `fx_surv` comes
  from the separate `na.omit` + `min(surv)` fallback at the final time, not from
  the interpolation.
* **Between observations the current code can be arbitrarily wrong**, up to 0.73.
* **Monotonicity fails** on `fx_surv_uncensored`, confirming the reported
  non-monotonicity is not confined to contrived tie-heavy data.
* **`G(t)` versus `G(t-)` is a second-order question.** The two differ by ~0.01
  on the affected Brier values, against a ~0.18 discrepancy from the bug itself.
  The choice should be made deliberately, but it is not what is driving the
  error. No choice between them is made here.

## 6. Before/after for the affected metric

### Existing baseline fixtures — no change

Every value in `tests/testthat/fixtures/metric-baseline.csv` was recomputed with
the corrected `Gt()`. All 57 values across all 9 scenarios were **unchanged**,
for the parity reason in §3.

### Forced interpolation — `tm_fit_and_eval`, `brier_score`

| Fixture | `t_star` | current | fixed | G(t) | G(t-) | abs diff |
|---|---|---|---|---|---|---|
| `fx_surv` | observed | 0.20 | 0.20 | 0.20 | 0.20 | 0 |
| `fx_surv` | midpoint | 0.21 | 0.21 | 0.21 | 0.21 | 0 |
| tie-heavy 200/10 | 5 (observed) | 0.20 | 0.20 | 0.20 | 0.19 | 0 |
| tie-heavy 200/10 | **5.5** | **0.41** | 0.23 | 0.23 | 0.22 | **0.18** |
| tie-heavy 200/10 | **2.5** | **0.16** | 0.11 | 0.11 | 0.11 | **0.05** |
| tie-heavy 200/10 | **7.3** | **0.44** | 0.26 | 0.25 | 0.24 | **0.18** |

### Unrounded `pam.Brier`, sweeping tie density (n = 200, non-observed `t_star`)

| distinct times | current | fixed | abs | relative |
|---|---|---|---|---|
| 5 | 0.234181 | 0.234181 | 0 | 0% |
| 10 | 0.344385 | 0.220332 | 0.1241 | **56.3%** |
| 20 | 0.254740 | 0.226532 | 0.0282 | 12.5% |
| 50 | **0.752935** | 0.224847 | **0.5281** | **235%** |
| 100 | 0.211257 | 0.211257 | 0 | 0% |
| 200 (no ties) | 0.223067 | 0.223067 | 0 | 0% |

Continuous times with no ties showed **no difference** at the 25th, 50th or 75th
percentile.

The error is **not monotone in tie density** — it depends on where `t_star` falls
relative to the misindexed raw observations. At 50 distinct times the reported
Brier score is **0.75 against a correct 0.22**, a 235% overstatement. A Brier
score above 0.25 for a binary-scale loss is already implausible, so a wrong value
of this size is more likely to be noticed than to mislead silently — but it is
reported without any warning.

## 7. Recommendation

**Correct `Gt()` and intentionally update the affected baselines.**

Reasoning:

* **The current behaviour is not defensible methodologically.** It is not an
  alternative convention; it indexes one vector with another vector's indices.
  It produces a non-monotone G(t) and is not invariant to row order. No
  published definition of the IPCW Brier score has these properties, so option
  (b) — retain with a methodological justification — has no justification
  available.
* **Withdrawal is disproportionate.** Only `tm_fit_and_eval()`'s `brier_score` is
  affected, only through interpolation, and the estimator itself is standard
  (Graf et al. 1999). Withdrawing a correct metric because one helper
  interpolates badly would remove working functionality to avoid a contained
  repair. Withdrawal becomes the right answer only if the authors decide the
  Brier path needs re-derivation rather than repair.
* **The blast radius is small and measured.** The committed baseline does not
  move at all; only forced-interpolation cases do. A correction is therefore
  reviewable, with a diff that shows exactly which values change and why.

Two decisions belong to the authors and are **not** made here:

1. **Interpolation versus step function.** Repairing the alignment leaves linear
   interpolation of a step function in place. Evaluating the KM step function
   directly — `G(t)` or `G(t-)` — is the more standard choice. The measured gap
   between fixed-interpolation and `G(t)` is small (~0.01 on the affected Brier
   values), so this can be settled on methodological grounds rather than under
   numerical pressure.
2. **`G(t)` versus `G(t-)`.** Graf et al.'s IPCW weights are conventionally
   written with the left limit `G(t-)` for subjects with an event at `t`. The
   current code uses neither consistently. Choosing one is a separate, deliberate
   change to the estimator.

**Suggested sequence if the correction is approved:** repair the alignment first
as a pure bug fix, with the baseline change reviewed as its own commit; then take
the interpolation and `G(t)`/`G(t-)` decisions separately, each with its own
equivalence evidence. Do not bundle a bug fix with an estimator change.

Until a decision is taken and acted on, the package should not be described as
ready for CRAN submission.

## Reproduction

The audit scripts are throwaway and were not committed. Each loads the package
with `pkgload::load_all()`, replaces `Gt` via `assignInNamespace()` in a
short-lived session, measures, and restores the original. No file in the
repository was written by any of them.
