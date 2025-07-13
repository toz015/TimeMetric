<!-- README.md for PAmeasure ---------------------------------------------- -->

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
<!-- Optional: add coverage, pkgdown, DOI badges here -->

`PAmeasure` provides a tidy, pipe-friendly toolkit for **prediction-accuracy
metrics in competing-risks and survival settings**, including  
  
* **Li–Wang pseudo-R<sup>2</sup>**, and  
* **Li–Wang pseudo-R<sup>2</sup>** for competing risk models.
* **C-index**,  
* time-dependent **AUC**,
* Brier score  

The package is aimed at researchers who need to evaluate their model's prediction performace.

---

## Installation

```r
# 1. install.packages("remotes")     # if not already installed
remotes::install_github("toz015/PAmeasure", build_vignettes = TRUE)
