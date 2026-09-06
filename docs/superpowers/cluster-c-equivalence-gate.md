# Cluster C Equivalence Gate

**Date:** 2026-09-06
**Question:** are `pam.rsh_metric`, `pam.rsph_metric` and `pam.Brier_metric` a
better foundation for `tm_evaluate_two_phase` than the current model-coupled
path? (spec section 2, Cluster C)

**Answer: no, for all three.** None can accept two-phase weights, and only one
is numerically equivalent to anything the two-phase path computes.

## Gate 1 - can they accept two-phase weights?

| Function | Signature | Weight argument |
|---|---|---|
| `pam.rsh_metric` | `(predicted_data, survival_time, status)` | **none** |
| `pam.rsph_metric` | `(time, status, risk_score)` | **none** |
| `pam.Brier_metric` | `(predicted_data, suvival_time, t_star)` | **none** |

The two-phase path threads `case_weights` into *every* metric it computes:

```
pam.r2_metrics(..., case_weight = case_weights)
survival::concordance(..., weights = dat1$case_weights)                # Harrell
survival::concordance(..., weights = ..., timewt = "n/G2")             # Uno
yardstick::brier_survival(..., case_weights = case_weights)
yardstick::roc_auc_survival(..., case_weights = case_weights)
```

Cluster C has no parameter through which a sampling weight could enter, so it
cannot express case-cohort or nested case-control estimation at all. This alone
disqualifies it as a foundation.

## Gate 2 - metric set overlap

Two-phase offers: `Pesudo_R`, `Harrell's C`, `Uno's C`, `Brier Score`,
`Time Dependent Auc`.

`R_sh` and `R_sph` are **not in that set**. Two of the three Cluster C functions
compute metrics the two-phase path does not offer, so they could not replace
anything there even if they were weighted.

## Gate 3 - numerical comparison, unit weights vs real weights

Brier is the only overlapping metric. Across five independent datasets, with
complete-cohort (unit) weights:

| seed | two-phase Brier | `pam.Brier_metric` | difference |
|---|---|---|---|
| 1001 | 0.205300 | 0.205254 | 0.000046 |
| 2002 | 0.197000 | 0.197031 | 0.000031 |
| 3003 | 0.185300 | 0.185340 | 0.000040 |
| 4004 | 0.186700 | 0.186654 | 0.000046 |
| 5005 | 0.205200 | 0.205167 | 0.000033 |

Every difference is at or below 5e-5, the maximum error of the two-phase's own
`round(brier, 4)`. **These are the same estimator under unit weights.**

Once weights are real, `pam.Brier_metric` cannot follow (single dataset, seed 1001):

| weighting | two-phase Brier | `pam.Brier_metric` |
|---|---|---|
| complete cohort | 0.2053 | 0.205254 |
| case-cohort | 0.2068 | 0.205254 (cannot vary) |
| nested case-control | 0.1811 | 0.205254 (cannot vary) |

## Gate 4 - the two R-squared metrics against the working paths

| seed | main `R_sh` (rms::cph + pam.schemper) | `pam.rsh_metric` | main `R_E` | `pam.rsph_metric$r2` |
|---|---|---|---|---|
| 1001 | -0.026234 | 0.104900 | 0.311868 | 0.320877 |
| 2002 | 0.305381 | 0.136100 | 0.385779 | 0.392610 |
| 3003 | 0.278375 | 0.070600 | 0.420299 | 0.442465 |

`pam.rsh_metric` is not close to the working `R_sh`, and the gap is not a
constant offset -- it changes sign and magnitude across datasets.

`pam.rsph_metric$r2` tracks `R_E` closely but is consistently distinct
(1-5% relative). This reconfirms findings.md #13 across three datasets:
**`R_sph` and `R_E` are different quantities and must stay separate.**

## Gate 5 - ordering, orientation and evaluation time

`pam.rsh_metric` is **order-dependent**, which is finding #1 confirmed exactly
rather than by correlation:

```
R_sh, rows in original order  = 0.104900
R_sh, same data rows permuted = 0.026400
```

A 4x change from row order alone. The cause is at `R/pam.rsh_metric.R:48-51`:
the internal data frame is sorted by `survival_time`, but `Mtx` is computed from
the unsorted `predicted_data` argument.

A methodological note on why this was not obvious: `pam.coxph_restricted` returns
`times` **already sorted**, so any comparison built from its output alone shows no
difference and the bug stays hidden. A genuine permutation is required to see it.

`pam.Brier_metric` is order-independent (it sorts all three inputs together) and
defaults `t_star` to the median observed time.

## Conclusion

Cluster C cannot serve as the foundation for `tm_evaluate_two_phase`. The current
model-coupled path is weight-aware throughout; Cluster C is weight-blind by
construction, and adapting it would mean rewriting each function to thread
sampling weights through its IPCW -- reimplementing what `yardstick` and
`survival::concordance` already provide correctly.
