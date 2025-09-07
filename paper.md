---
title: 'TimeMetric: A Comprehensive R Package for Evaluating Predictive Performance with Survival Data'
tags:
  - R
  - Survival
authors:
  - name: Tong Zhu
    affiliation: 1
  - name: Zian Zhuang
    affiliation: 1
  - name: Wen Su
    affiliation: 2
  - name: Xiaowu Dai
    affiliation: 1
  - name: Gang Li
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: 1
affiliations:
 - name: Department of Biostatistics, University of California at Los Angeles
   index: 1
 - name: Department of Biostatistics, City University of Hong Kong
   index: 2
date: 6 September 2025
header-includes:
  - \usepackage{graphicx}  
  - \usepackage{booktabs}  
  - \usepackage{array}      
  - \usepackage{amsmath,amssymb}
  - \usepackage{float} 
---

# Summary {#summary .unnumbered}

The evaluation of predictive performance is critical for survival models
in biomedical research, where time-to-event data exhibit unique
characteristics, such as censoring, competing risks, and inherited data
structures from specialized designs like the case-cohort design -- that
requires specialized analytical approaches. Predictive performance
measures enable the evaluation and comparison of the predictive or
prognostic performance of different models or variables. These measures
fall into two fundamental categories: discrimination measures, which
quantify a model's ability to distinguish between high- and low-risk
subjects through outcome ranking (e.g., C-index, time-dependent AUC)
(Harrell et al. 1982; Uno et al. 2011; Heagerty, Lumley, and Pepe 2000),
and calibration measures, which assess the agreement between predicted
probabilities and observed event rates (e.g., Brier score) (Brier 1950;
Graf et al. 1999). Existing packages offer implementations for some of
these metrics, and none provide a comprehensive solution. Table 
\ref{tab:metrics} provides a comparison of predictive
performance measures currently implemented across some existing `R`
packages. `TimeMetric` bridges this gap by delivering a unified
framework for performance evaluation across diverse survival data
scenarios, including right-censored data, competing risks, and special
designs such as nested case-control and case-cohort studies.

Designed for ease of implementation, `TimeMetric` allows researchers to
compute a suite of metrics for evaluating the predictive performance of
a model or algorithm. For each type of survival data, `TimeMetric`
provides two key modules. First, it incorporates a comprehensive suite
of performance evaluation metrics, implemented in the
`Performance Metric Module`, including $R^2$-type measures (pseudo $R^2$
(Li and Wang 2019; Zhuang et al. 2025), Schemper and Henderson's
$R_{\text{sh}}$ (Schemper and Henderson 2000; Lusa, Miceli, and Mariani
2007), and $R_E$ (Stare, Perme, and Henderson 2011)), discrimination
indices (Harrell's C-index (Harrell et al. 1982) and Uno's C-index (Uno
et al. 2011)), calibration tools (Brier score (Brier 1950; Graf et al.
1999)), and time-dependent AUC (Heagerty, Lumley, and Pepe 2000). A
complete summary of available metrics for each data type is presented in
Table \ref{tab:metrics}. `Performance Metric Module` takes
predicted survival probabilities (or cumulative incidence functions for
competing risks data)---regardless of how they are obtained---as input
to compute predictive performance metrics. Second, the
`Prediction Module` provides functions for several common survival
models to calculate predicted survival probabilities (or cumulative
incidence functions for competing risks data) over the observation times
for all individuals in a test dataset.

::: {=latex}
\begin{table}[H]
\centering
\caption{Availability of Performance Evaluation Metrics Across Selected Packages}
\label{tab:metrics}
\resizebox{\textwidth}{!}{%
\begin{tabular}{@{} p{6cm} c c c c c c c c c @{}}
\toprule
 & \multicolumn{9}{c}{\textbf{Packages}} \\ 
\cmidrule(lr){2-10}  
\textbf{Metric} & TimeMetric & PAmeasure & survival & survcomp & survAUC & BP & SurvMetrics & pec & timeROC \\
\midrule
\multicolumn{10}{l}{\emph{Right-Censored Data}} \\
Pseudo $R^2$ & $\times$ & $\times$ &  &  &  &  &  &  &  \\
C-index & $\times$ &  & $\times$ & $\times$ &  &  & $\times$ & $\times$ &  \\
$R_{sh}$ & $\times$ &  &  &  & $\times$ &  &  &  &  \\
$R_E$  & $\times$ &  &  &  &  & $\times$ &  &  &  \\
Brier Score & $\times$ &  & $\times$ & $\times$ &  &  & $\times$ & $\times$ &  \\
Time-Dependent AUC & $\times$ &  &  & $\times$ & $\times$ &  & $\times$ & $\times$ & $\times$ \\
\midrule
\multicolumn{10}{l}{\emph{Right-Censored Competing Risks Data}} \\
Pseudo $R^2$ & $\times$ &  &  &  &  &  &  &  &  \\
C-index & $\times$ &  & $\times$ &  &  &  &  & $\times$ &  \\
$R_E$ & $\times$ &  &  &  &  &  &  &  &  \\
Brier Score & $\times$ &  &  &  &  &  &  & $\times$  &   \\
Time-Dependent AUC & $\times$ &  &  &  &  &  &   & $\times$  & $\times$  \\
\midrule
\multicolumn{10}{l}{\emph{Right-Censored Case-Cohort Data}} \\
C-index & $\times$ &  &  &  &  &  &  &  &  \\
Brier Score & $\times$ &  &  &  &  &  &  &  &   \\
Time-Dependent AUC & $\times$ &  &  &  &  &  &   &  &  \\
\midrule
\multicolumn{10}{l}{\emph{Right-Censored Nested Case-Control Data}} \\
C-index & $\times$ &  &  &  &  &  &  &  &  \\
Brier Score & $\times$ &  &  &  &  &  &  &  &   \\
Time-Dependent AUC & $\times$ &  &  &  &  &  &   &  &  \\
\bottomrule
\end{tabular}
}
\end{table}
:::

# Reference

::::::::::::: {#refs .references .csl-bib-body .hanging-indent entry-spacing="0"}
::: {#ref-brier1950verification .csl-entry}
Brier, Glenn W. 1950. "Verification of Forecasts Expressed in Terms of
Probability." *Monthly Weather Review* 78 (1): 1--3.
:::

::: {#ref-graf1999assessment .csl-entry}
Graf, Erika, Claudia Schmoor, Willi Sauerbrei, and Martin Schumacher.
1999. "Assessment and Comparison of Prognostic Classification Schemes
for Survival Data." *Statistics in Medicine* 18 (17-18): 2529--45.
:::

::: {#ref-harrell1982evaluating .csl-entry}
Harrell, Frank E, Robert M Califf, David B Pryor, Kerry L Lee, and
Robert A Rosati. 1982. "Evaluating the Yield of Medical Tests." *Jama*
247 (18): 2543--46.
:::

::: {#ref-heagerty2000time .csl-entry}
Heagerty, Patrick J, Thomas Lumley, and Margaret S Pepe. 2000.
"Time-Dependent ROC Curves for Censored Survival Data and a Diagnostic
Marker." *Biometrics* 56 (2): 337--44.
:::

::: {#ref-li2019prediction .csl-entry}
Li, Gang, and Xiaoyan Wang. 2019. "Prediction Accuracy Measures for a
Nonlinear Model and for Right-Censored Time-to-Event Data." *Journal of
the American Statistical Association* 114 (528): 1815--25.
:::

::: {#ref-lusa2007estimation .csl-entry}
Lusa, Lara, Rosalba Miceli, and Luigi Mariani. 2007. "Estimation of
Predictive Accuracy in Survival Analysis Using r and s-PLUS." *Computer
Methods and Programs in Biomedicine* 87 (2): 132--37.
:::

::: {#ref-schemper2000predictive .csl-entry}
Schemper, Michael, and Robin Henderson. 2000. "Predictive Accuracy and
Explained Variation in Cox Regression." *Biometrics* 56 (1): 249--55.
:::

::: {#ref-stare2011measure .csl-entry}
Stare, Janez, Maja Pohar Perme, and Robin Henderson. 2011. "A Measure of
Explained Variation for Event History Data." *Biometrics* 67 (3):
750--59.
:::

::: {#ref-uno2011c .csl-entry}
Uno, Hajime, Tianxi Cai, Michael J Pencina, Ralph B D'Agostino, and
Lee-Jen Wei. 2011. "On the c-Statistics for Evaluating Overall Adequacy
of Risk Prediction Procedures with Censored Survival Data." *Statistics
in Medicine* 30 (10): 1105--17.
:::

::: {#ref-zhuang2025time .csl-entry}
Zhuang, Zian, Wen Su, Eric Kawaguchi, and Gang Li. 2025. "Time-Dependent
Pseudo $R^{2}$ for Assessing Predictive Performance in Competing Risks
Data." *arXiv Preprint arXiv:2507.15040*.
:::
:::::::::::::
