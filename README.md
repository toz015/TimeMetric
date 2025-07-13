<!-- README.md for pkg "statTools" ----------------------------------------- -->

# statTools <img src="man/figures/logo.png" align="right" height="120" />

[![R-CMD-check](https://github.com/yourname/statTools/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yourname/statTools/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
<!-- Add more badges: coverage, pkgdown site, DOI, etc. -->

`statTools` provides tidy, pipe-friendly helpers for
survival analysis and causal inference, including  
* **weighted C-index**,  
* **time-dependent AUC**, and  
* **Li–Wang pseudo-\(R^{2}\)**.  

Built for teaching examples and simulation studies where plain
`survival` objects need extra metrics.

---

## Installation

```r
# 1. Install remotes (or pak) if needed
install.packages("remotes")     # or install.packages("pak")

# 2. Install the development version from GitHub
remotes::install_github("yourname/statTools", build_vignettes = TRUE)
# pak::pkg_install("yourname/statTools")   # alternative
