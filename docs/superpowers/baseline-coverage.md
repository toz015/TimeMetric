# Coverage - after Cluster C removal

Measured 2026-09-06 on the package as committed.

| stage | total |
|---|---|
| characterization baseline (on a patched copy; the package would not install) | 63.09% |
| after Cluster B removal | 75.38% |
| after Cluster C removal | **74.62%** |

The small decrease is arithmetic, not a regression. `pam.Brier_metric`,
`pam.rsh_metric` and `pam.rsph_metric` were each at 100% coverage from their
characterization tests, so deleting them removed covered lines from the
numerator. No surviving file lost coverage.

The three files still at 0% are Cluster A -- `pam.coxph`, `pam.nlm` and
`pam.survreg` -- unreachable code deliberately left in place, outside the scope
of the Cluster B and Cluster C decisions.

```
TOTAL: 74.62%

cc_weights.R                                  100.0%
ncc_weights.R                                 100.0%
pam.concordance.R                             100.0%
pam.schemper.R                                 97.8%
pam.predicted_survial_eval_two_phase.R         94.5%
pam.r2_metrics.R                               92.7%
pam.predicted_survival_eval_cr.R               90.0%
pam.predicted_survial_eval.R                   90.0%
plot.R                                         89.6%
pam.coxph_restricted.R                         89.4%
pam.survial_eval.R                             87.0%
pam.surverg_restricted.R                       86.2%
pam.Ct.R                                       78.4%
simulateTwoCauseFineGrayModel.R                76.5%
pam.weighted_param.R                           75.6%
pam.Brier.R                                    69.2%
pam.sim_data.R                                 66.3%
pam.predict_cr.R                               64.5%
pam.predictSurvProb2survreg.R                  63.2%
pam.rsph.R                                     60.2%
pam.survivalROC.R                              51.6%
pam.coxph.R                                     0.0%
pam.nlm.R                                       0.0%
pam.survreg.R                                   0.0%
```
