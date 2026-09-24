# API Reference Tables

**Date:** 2026-09-08. Generated from the installed package, not hand-written.

## 1. Function rename table

| # | Old name | New name | Old exported | New exported | Deprecation wrapper | Rd page |
|---|---|---|---|---|---|---|
| 1 | `pam.predicted_survial_eval` | `tm_survival_eval` | yes | yes | `.Deprecated` | yes |
| 2 | `pam.predicted_survial_eval_cr` | `tm_survival_eval_cr` | yes | yes | `.Deprecated` | yes |
| 3 | `pam.survival_eval` | `tm_fit_and_eval` | yes | yes | `.Deprecated` | yes |
| 4 | `pam.coxph_restricted` | `tm_predict_coxph` | yes | yes | `.Deprecated` | yes |
| 5 | `pam.surverg_restricted` | `tm_predict_survreg` | yes | yes | `.Deprecated` | yes |
| 6 | `pam.predict_cr` | `tm_predict_cif` | yes | yes | `.Deprecated` | yes |
| 7 | `pam.summary` | `tm_summarize` | yes | yes | `.Deprecated` | yes |
| 8 | `pam.summary_cr` | `tm_summarize_cr` | yes | yes | `.Deprecated` | yes |
| 9 | `pam.sample_design` | `tm_sample_design` | yes | yes | `.Deprecated` | yes |
| 10 | `cc_weights` | `tm_case_cohort_weights` | yes | yes | `.Deprecated` | yes |
| 11 | `ncc_weights` | `tm_nested_case_control_weights` | yes | yes | `.Deprecated` | yes |
| 12 | `plot_pred` | `tm_plot_pred` | yes | yes | `.Deprecated` | yes |
| 13 | `summary_pred_plot` | `tm_plot_summary` | yes | yes | `.Deprecated` | yes |
| 14 | `sim_cox_weibull_censored` | `tm_sim_cox_weibull` | yes | yes | `.Deprecated` | yes |
| 15 | `simulateTwoCauseFineGrayModel` | `tm_simulate_fine_gray` | yes | yes | `.Deprecated` | yes |
| 16 | *(none, was internal)* | `tm_evaluate_two_phase` | n/a | yes | n/a | yes |

Also newly exported: `tm_metric_names()` (yes).

**Total exports: 32** = 16 renamed functions + `tm_metric_names` + 15 deprecated wrappers.

No function was removed from the public API. Every old name forwards to its
replacement with a warning naming the new function.

## 2. Canonical metric names and legacy aliases

Matching is case-insensitive and treats `_`, space, and straight or curly
apostrophes as equivalent, so the listed spellings are representative rather
than exhaustive.

| Canonical | Accepted legacy spellings |
|---|---|
| `brier_score` | `Brier Score`, `Brier_Score` |
| `c_index` | `C_index` |
| `harrell_c` | `Harrells_C`, `Harrell&rsquo;s C`, `Harrell's C` |
| `l_square` | `L_square` |
| `l2_point` | `L2_point` |
| `pseudo_r2` | `Pseudo_R_square`, `Pesudo_R`, `Psuedo.R` |
| `pseudo_r2_point` | `Pseudo_R2_point` |
| `r_square` | `R_square` |
| `r2_point` | `R2_point` |
| `td_auc` | `Time Dependent Auc`, `Time Dependent AUC`, `Time_Dependent_Auc`, `AUC` |
| `uno_c` | `Unos_C`, `Uno&rsquo;s C`, `Uno's C` |

**13 canonical names, 37 accepted spellings.**

Defaults by entry point:

| Function | Default metrics |
|---|---|
| `tm_survival_eval` | `pseudo_r2`, `pseudo_r2_point`, `harrell_c`, `uno_c`, `brier_score`, `td_auc` |
| `tm_survival_eval_cr` | `pseudo_r2`, `pseudo_r2_point`, `c_index`, `brier_score`, `td_auc` |
| `tm_evaluate_two_phase` | `pseudo_r2`, `harrell_c`, `uno_c`, `brier_score`, `td_auc` |

## 3. Example review

All four rewritten examples execute during `R CMD check`. None is hidden behind
`\dontrun{}`. Runtimes measured on `survival::pbc` (n = 312 complete cases,
40.1% event rate).

| Example | Runs | Time | Output |
|---|---|---|---|
| `tm_predict_coxph` | yes | <0.1s | 104 restricted mean survival times on held-out data |
| `tm_predict_survreg` | yes | <0.1s | 104 restricted mean survival times |
| `tm_predict_cif` | yes | 0.1s | 312 CIF-derived risks, both the cause-specific Cox and Fine-Gray paths |
| `tm_fit_and_eval` | yes | 1.7s | full metric table, then a model and metric subset |

### Scientific assessment

**`tm_predict_coxph`, `tm_predict_survreg`** — sound. `pbc` status is collapsed
with `status == 2` as the event, so transplant is treated as censoring; that is
the conventional handling of this dataset. A 2/3 training split is used and
predictions are made on held-out data, which is the correct way to demonstrate
an evaluation package. Cox and Weibull produce the same subject ordering with
plausible magnitudes (793 vs 679, 1825 vs 1510 days).

**`tm_predict_cif`** — technically correct, but see finding #33 below. It keeps
status as 0/1/2 and evaluates `event.type = 1`. In `pbc`, 1 is *transplant* and
2 is *death*, so the example treats transplant as the event of interest with
death as the competing risk. This is a valid competing-risks setup but an
unusual scientific choice; death-as-primary is the more natural framing. The
Fine-Gray branch is correctly guarded with `requireNamespace("cmprsk")`.

**`tm_fit_and_eval`** — sound, with an internal consistency check that passes.
Harrell's C of 0.84 on `pbc` matches the published literature (~0.83-0.85). The
reported `pseudo_r2` of 0.32 equals `r_square * l_square` = 0.39 x 0.83 = 0.3237,
confirming the pseudo R-squared is the product of its two components as
documented.

### Open items arising from this review

| # | Item | Status |
|---|---|---|
| 33 | `tm_predict_cif`'s example uses transplant as the primary event. Valid, but death-as-primary would be the conventional choice for `pbc`. | **Deferred — maintainer's scientific call.** Not changed unilaterally. |
| 34 | `tm_predict_survreg()` returns a **named** `pred` vector while `tm_predict_coxph()` returns an unnamed one, for the same quantity. | **Deferred.** A safe `unname()` would harmonise them, but it changes a return value, so it is not made without approval. |

Neither affects correctness, `R CMD check`, or the test suite. Both are recorded
so they are not lost.

## Amendment 2026-09-24

`r_sh` and `r_e` were withdrawn before the first CRAN release and removed from
every table above. Their spellings -- `r_e`, `r_sh`, `R_E`, `R_sh`, `R_sph` --
now raise an error naming the withdrawal rather than resolving to a canonical
name. `tm_metric_names()` returns 11 names.
