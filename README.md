<!-- README.md for TimeMetric (PAmeasure) ---------------------------------------------- -->

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
<!-- Optional: add coverage, pkgdown, DOI badges here -->

The `TimeMetric` (PAmeasure) R package offers a comprehensive framework for evaluating prediction performance in survival models, including scenarios with right-censoring, competing risks, as well as under different designs such as nested case control and case-cohort designs. It provides a comprehensive suites of metrics such as pseudo $R$-squared, concordance indices, Brier Score, and time-dependent AUC, unified under a single platform. This paper presents an overview of the mathematical definitions of these metrics,  implementation details, and application examples of these measures. Demonstrations using simulated and real-world datasets validate the utility and robustness of `TimeMetric`,  making it a valuable tool for researchers and practitioners in survival analysis.


The package focuses on four key categories of performance evaluation:

- **R²-related metrics**: quantify the proportion of variability in survival times explained by the model.
- **Concordance indices (C-indices)**: measure the discriminatory power of the model.
- **Time-dependent AUC**: evaluate model discrimination at different time points.
- **Brier Score**: assess the accuracy of probabilistic survival predictions.


The \CRANpkg{TimeMetrics} package is organized into six main components:

- `pam.coxph_restricted` and `pam.survreg_restricted`: functions for generating predicted survival times from Cox and parametric survival models.
- `pam.predicted_survival_eval`: the primary function for computing survival performance metrics.
- `pam.predict_subject_cif`: calculates cumulative incidence function (CIF) predictions for individual subjects across multiple types of competing risks models.
- `pam.predicted_survival_eval_cr`: extended version for evaluating predictions in the presence of competing risks.
- `pam.predicted_survial_eval_casecohort`: evaluates survival model performance in a case–cohort setting using predicted
survival probabilities and case–cohort sampling weights
- `pam.predicted_survial_eval_ncc`:  evaluates survival model performance under a nested case–control (NCC) de-
sign using predicted survival probabilities and NCC sampling weights.

---

## Installation

```r
# 1. install.packages("remotes")     # if not already installed
remotes::install_github("toz015/PAmeasure", build_vignettes = TRUE)
