# R_E Implementation Audit

**Date:** 2026-09-06
**Status:** analysis only -- no code changed
**Question:** which of TimeMetric's two `R_E` implementations is canonical?

## 1. Reference provenance

| | |
|---|---|
| URL | `https://ibmi3.mf.uni-lj.si/ibmi-english/biostat-center/programje/Re.r` |
| Retrieved | 2026-09-06T21:52:16Z |
| SHA-256 | `cf6532068a9a86a61d3763b4cad54aff4bb5000a3ce108f64135debab2b69a92` |
| Size | 18,946 bytes / 499 lines |
| Licence | **None stated anywhere in the file** |

The file is kept in a scratch directory and is **not committed**: with no licence
grant, redistribution is not permitted. Anyone repeating this audit should
re-download it and check the hash above.

`R/pam.rsph.R:6` already cites this URL, so the provenance is the authors' own.

## 2. Reference functions

| Reference | TimeMetric counterpart |
|---|---|
| `re <- function(fit, ...) UseMethod("re")` | `pam.rsph` |
| `re.coxph(fit, Gmat)` | `pam.rsph.coxph` |
| `re.aareg(fit, Gmat)` | `pam.rsph.aareg` |
| `re.survreg(fit, Gmat)` | `pam.rsph.survreg` |
| `summary.re(object, times, band = 5, ...)` | `pam.summary.rsph` |
| `print.re(x, digits = 4, ...)` | `pam.print.rsph` |
| `my.survfit(start, stop, event)` | `my.survfit` (`R/pam.rsph.R:439`) |

`re.coxph` inputs: a fitted `coxph` object supplying `fit$y` and
`fit$linear.predictors`, and an optional `Gmat` of censoring weights. When
`Gmat` is omitted it is built from `my.survfit(...)$surv.i2[n.event != 0]`,
recycled across subjects (times in columns, individuals in rows).

Sorting is `order(Y[,2], -Y[,3])` -- ascending time, events before censorings at
ties. Linear predictors are centred (`bx <- lp[sort.it] - mean(lp)`) and
**`Gmat` is sorted by the same permutation**, so weights stay aligned with
subjects. Risk direction is `rank(-bx)`, so higher risk gives lower rank. Ties
are counted per event time (`k`) and handled through the `[1:k]` slices.
`Re = sum(meanr - rg) / sum(meanr - smin)`.

## 3. Line-level comparison

**`pam.rsph.coxph` is a verbatim port of `re.coxph`.** After normalising names
and whitespace, the computational core is identical line for line: the sort key,
the centring, the `Gmat` construction and re-sorting, the at-risk indicator, the
IPCW rank formula, tie handling, `pi` weights, the variance terms, and the
numerator/denominator. The only substantive addition is optional `test_data`
support, which substitutes out-of-sample linear predictors and rebuilds `Y`.

**`pam.rsph_metric` is an independent reimplementation**, not a port. Compared
with the reference:

| Element | Reference `re.coxph` | `pam.rsph_metric` | Equivalent? |
|---|---|---|---|
| Sort key | `order(Y[,2], -Y[,3])` | `order(time, -status)` | yes |
| Centring of `lp` | `lp - mean(lp)` | none | yes -- `rank()` and `ebx/sum(ebx)` are both shift-invariant |
| At-risk set | `Y[,2] >= ti & Y[,1] < ti` | `time >= t_i` | yes for right-censored data (`Y[,1] == 0`) |
| Censoring weights | `my.survfit(...)$surv.i2[n.event != 0]` | `summary(survfit(Surv(time, 1-status) ~ 1), times = event_times)$surv` | yes -- verified numerically identical on the fixtures |
| **Rank definition** | `rangi <- (rank(-bx[inx]) - .5)/Gmat[inx,it] + .5` | `ranks <- rank(-risk_score[at_risk])` | **NO -- the IPCW adjustment is absent** |
| `r0` | `mean(rangi)` (adjusted ranks) | `mean(ranks)` (plain ranks) | no -- follows from the above |
| `rP`, numerator, denominator | as published | same structure | yes |

## 4. Numerical results, full precision

Fixtures are the characterization-suite generator, `n = 200`, `pi_c = 0.3`,
`v = 2`, `beta = c(0.5, -0.5)`. Identical fitted `coxph` objects, identical
observations and event indicators, full-precision linear predictors, no rounding.

| seed | reference `re.coxph$Re` | `pam.rsph$Re` | abs diff | `pam.rsph_metric$r2` | rel. diff |
|---|---|---|---|---|---|
| 1001 | 0.3118680658 | 0.3118680658 | **0.00e+00** | 0.3208774922 | 2.89% |
| 2002 | 0.3857787756 | 0.3857787756 | **0.00e+00** | 0.3926096232 | 1.77% |
| 3003 | 0.4202986593 | 0.4202986593 | **0.00e+00** | 0.4424652347 | 5.27% |
| 4004 | 0.3607495364 | 0.3607495364 | **0.00e+00** | 0.3654341109 | 1.30% |
| 5005 | 0.2905941385 | 0.2905941385 | **0.00e+00** | 0.3019530292 | 3.91% |

`pam.rsph` reproduces the reference **exactly**, to every digit, on all five
datasets.

## 5. Root cause

A single omitted operation. Substituting the reference's rank line into
`pam.rsph_metric` and changing nothing else:

| seed | reference | `pam.rsph_metric` as written | with IPCW ranks restored |
|---|---|---|---|
| 1001 | 0.3118680658 | 0.3208774922 | **0.3118680658** |
| 2002 | 0.3857787756 | 0.3926096232 | **0.3857787756** |
| 3003 | 0.4202986593 | 0.4424652347 | **0.4202986593** |

The entire 1.3-5.3% disagreement comes from `pam.rsph_metric` using plain ranks

```r
ranks <- rank(-risk_score[at_risk])
```

where Stare, Perme & Henderson weight each rank by the inverse probability of
remaining uncensored:

```r
rangi <- (rank(-bx[inx]) - .5)/Gmat[inx,it] + .5
```

Because the omission is in `ranks`, it propagates into both `r0` (hence
`mean_rank`) and `observed_rank`, but not into `rP`. Sorting, alignment, tie
handling, rank direction, the at-risk set, the `G` estimator, and the
numerator/denominator definitions are all equivalent -- each was checked and none
contributes.

## 6. Inconsistency in the reference itself

`re.coxph` returns `r2nw`, commented `#unweighted measure`, computed as

```r
num   <- sum(meanr - rg)
den   <- sum(meanr - smin)
r2    <- num/den
r2nw  <- sum((meanr - rg))/sum((meanr - smin))     # algebraically identical to r2
```

These are the same expression, and both `meanr`, `rg` and `smin` are built from
`Gmat`-weighted quantities. Verified numerically: `Re == r2nw` exactly. So
`r2nw` is **not** an unweighted measure despite its name and comment. This is a
defect in the reference code, not in TimeMetric, and TimeMetric inherits it
through the port. It affects only the auxiliary `r2nw` field, never `Re`.

## 7. Manuscript impact: none

Both reported quantities come from the canonical path:

* `R/pam.predicted_survial_eval.R:226` -- `R_E` via
  `pam.summary.rsph(pam.rsph(model, test_data = new_data), ...)`
* `R/pam.survial_eval.R:137` -- `R_sph` via `pam.rsph(fits[[fit_name]])$Re`

`pam.rsph_metric` has **no callers anywhere in the package**. Its only caller was
`pam.prediction_metrics`, itself deleted as Cluster B dead code. It is therefore
impossible for any published number, example, or `paper.code.Rmd` output to have
come from it.

**Correcting or deleting `pam.rsph_metric` cannot change any published result.**

## 8. Recommendation

**Canonical: `pam.rsph` / `pam.summary.rsph`.** It reproduces the authors' own
reference implementation exactly and is what the manuscript's numbers already
come from.

**Delete `pam.rsph_metric`.** It is unreachable, numerically wrong by 1.3-5.3%,
documents an argument order it does not have (finding #28), and cites the wrong
paper (finding #29). Repairing it would produce a second implementation that is
merely a slower duplicate of a verified-correct one; the spec's rule against two
undocumented implementations of one metric applies directly.

## 9. Naming

* **`R_sh`** -- Schemper & Henderson (2000). A genuinely distinct metric.
  Unaffected by this audit.
* **`R_sph`** and **`R_E`** -- two labels for the *same* Stare, Perme &
  Henderson (2011) metric. `pam.survival_eval` emits `R_sph`;
  `pam.predicted_survial_eval` emits `R_E` for the identical quantity from the
  identical code path.

Finding #13 is superseded: the gap was never two estimands, it was one estimand
computed correctly and incorrectly. The public label stays unstandardised until
the canonical implementation is adopted.

## 10. Proposed regression tests

Expected values are stored **locally as literals** with provenance comments. No
test downloads or sources the remote file; `R CMD check` stays offline.

```r
# Expected values are the authors' reference implementation, Re.r, retrieved
# 2026-09-06 from
#   https://ibmi3.mf.uni-lj.si/ibmi-english/biostat-center/programje/Re.r
#   SHA-256 cf6532068a9a86a61d3763b4cad54aff4bb5000a3ce108f64135debab2b69a92
# re.coxph(fit)$Re evaluated on sim_cox_weibull_censored(n = 200, pi_c = 0.3,
# v = 2, beta = c(0.5, -0.5), seed = <seed>) with
# coxph(Surv(time, status) ~ x1 + x2, x = TRUE, y = TRUE).
# See docs/superpowers/r-e-implementation-audit.md. Do NOT regenerate these by
# running TimeMetric -- that would make the test circular.
REF_RE <- c(
  "1001" = 0.3118680658,
  "2002" = 0.3857787756,
  "3003" = 0.4202986593,
  "4004" = 0.3607495364,
  "5005" = 0.2905941385
)

test_that("pam.rsph reproduces the Stare-Perme-Henderson reference exactly", {
  for (sd in names(REF_RE)) {
    d <- sim_cox_weibull_censored(n = 200, pi_c = 0.3, v = 2,
                                  beta = c(0.5, -0.5),
                                  seed = as.integer(sd))[, c("time", "status", "x1", "x2")]
    m <- survival::coxph(survival::Surv(time, status) ~ x1 + x2,
                         data = d, x = TRUE, y = TRUE)
    expect_equal(TimeMetric:::pam.rsph(m)$Re, unname(REF_RE[[sd]]),
                 tolerance = 1e-9,
                 info = paste("seed", sd))
  }
})
```

A companion test asserts the reference file is *not* vendored, so the licence
position cannot silently change:

```r
test_that("the third-party reference implementation is not redistributed", {
  root <- skip_without_source_tree()
  expect_false(file.exists(file.path(root, "Re.r")))
  expect_false(file.exists(file.path(root, "inst", "Re.r")))
})
```
