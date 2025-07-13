<!-- README.md for PAmeasure ---------------------------------------------- -->

# PAmeasure <img src="man/figures/logo.png" align="right" width="120"/>

[![R-CMD-check](https://github.com/toz015/PAmeasure/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/toz015/PAmeasure/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
<!-- Optional: add coverage, pkgdown, DOI badges here -->

`PAmeasure` provides a tidy, pipe-friendly toolkit for **prediction-accuracy
metrics in competing-risks and survival settings**, including  

* weighted and IPCW-adjusted **C-index**,  
* time-dependent **AUC**,  
* **Li–Wang pseudo-R<sup>2</sup>**, and  
* convenient wrappers for bootstrap confidence intervals.  

The package is aimed at researchers who need a lightweight layer on top of
`survival`, `riskRegression`, or simulation workflows.

---

## Installation

```r
# 1. install.packages("remotes")     # if not already installed
remotes::install_github("toz015/PAmeasure", build_vignettes = TRUE)
