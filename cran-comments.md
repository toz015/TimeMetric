## Resubmission

This is a resubmission of TimeMetric 0.2.0, now 0.2.1. It addresses every item
in the review by Konstanze Lauseker (CRAN), received 2026-10-04:

* **Software names in Title and Description.** Package names are now in single
  quotes ('PAmeasures'). The predecessor package's name has also been corrected
  from "PAmeasure" to its exact, case-sensitive CRAN name 'PAmeasures'. The Title
  contains no software names.
* **Description must not start with the package name.** It now starts
  "Evaluates the predictive performance of survival models ...".
* **Explain acronyms.** AUC is written out as "area under the receiver
  operating characteristic curve (AUC)".
* **Method references.** The Description now cites the methods implemented, in
  the form `authors (year) <doi:...>`: Li and Wang (2019), Harrell et al. (1982),
  Uno et al. (2011), Heagerty, Lumley and Pepe (2000), Graf et al. (1999),
  Wolbers et al. (2014), Wu and Li (2018), Zhuang et al. (2025) and Lee,
  Zeleniuch-Jacquotte and Liu (2021). The Graf et al. DOI contains `<` and `>`,
  which are URL-encoded as `%3C` and `%3E`.
* **TRUE/FALSE instead of T/F.** `tm_predict_coxph()` and `tm_predict_survreg()`
  now default to `predict = TRUE`, and `T`/`F` no longer appear anywhere in the
  package code or documentation.
* **print()/cat() to the console.** All such output has been removed from
  non-print functions:
  * `R/pam.survivalROC.R` has been removed entirely (see below).
  * `R/tm_predict_survreg.R`: a failed weighted least-squares fit now raises
    `warning()` instead of using `cat()`.
  * `R/tm_predict_cif.R`: a failed model fit now raises `warning()` instead of
    using `print()`.
  * An unused debugging `print()` in an internal competing-risks helper was
    removed.

We also made two related changes while reviewing the package against CRAN
policy:

* `R/pam.survivalROC.R` was a copy of `survivalROC::survivalROC()`, which is
  licensed GPL (>= 2). Rather than ship that code under this package's MIT
  license, TimeMetric now imports 'survivalROC' and calls it directly. The test
  suite confirms that every metric value is unchanged.
* `tm_simulate_fine_gray()` and `tm_plot_summary()` no longer call `set.seed()`
  with a fixed default. Their `seed` argument now defaults to `NULL`, so the
  user's random number stream is changed only when they ask for it.

## Test environments

* local: Ubuntu 24.04, x86_64-pc-linux-gnu, R 4.3.3
* GitHub Actions: Ubuntu (R release, devel and oldrel-1), Windows (R release)
  and macOS (R release)
* win-builder: R-devel and R-release (TODO: run on the final tarball and
  record the results here before submitting)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

The CRAN incoming spell check may list author names from the Description's
references (for example Heagerty, Lumley, Wolbers and Zeleniuch-Jacquotte) and
"PAmeasures" as possibly misspelled. These are surnames and an exact package
name, and they are spelled correctly.

## Downstream dependencies

None. This is a new package with no reverse dependencies.
