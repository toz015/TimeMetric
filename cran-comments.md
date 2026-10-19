## Resubmission

This is a resubmission of TimeMetric 0.2.0 as version 0.2.1, addressing the
review received on 2026-10-04.

* The Description begins "Evaluates the predictive performance ...".
* Software names are quoted and case-sensitive ('PAmeasures'); the Title
  contains no software names.
* AUC is expanded as "area under the receiver operating characteristic curve".
* Nine method references use the form authors (year) <doi:...>. Embedded angle
  brackets in the Graf et al. DOI are encoded as %3C and %3E.
* Prediction defaults and logical expressions use TRUE/FALSE in source and
  generated documentation.
* Model-fit diagnostics use warning(). The copied survivalROC implementation
  has been removed in favor of the declared 'survivalROC' dependency.
  Incomplete training/evaluation rows are handled before calling its backend,
  with a suppressible warning; regression tests verify no console output.

The default two-phase AUC requires 'tidyr', now declared in Imports. An
additional CI job exercises this path with required dependencies only.
Simulation and plot seed arguments default to NULL: they do not reset the
random stream, but random draws advance it normally.

Independent validation also identified numerical errors, corrected in this
resubmission: Cox validation predictions now use the fitted training baseline;
finite-horizon AFT means include the survival tail and preserve observed
outcomes; predictions and weights remain aligned during sorting; point
pseudo R-squared uses event probabilities; AUC retains below-chance values;
and nested case-control weights use the documented number of controls per
case and matching strata. Metric-subset labels and standalone competing-risk
AUC are also corrected. NEWS.md describes the changed results and tests.
The obsolete 'expint' dependency has been removed.

## Local validation

R 4.5.3, aarch64-apple-darwin25.3.0, macOS Tahoe 26.5.1:

* Source test suite with NOT_CRAN=true: 661 assertions across 141 tests;
  0 failures, 0 warnings, 0 skips, including optional prediction backends.
* Independent numerical checks: 11/11 pass against fitted-model survival
  predictions, numerical integration, and hand-computable sampling/ranking
  references.
* All 17 current exported functions exercised by a public/simulated-data smoke
  test. Default two-phase metrics also pass with only required dependencies.

## R CMD check results

On the built 0.2.1 source tarball, with network access:
`R CMD check --as-cran --run-donttest` reports **0 errors, 0 warnings, 2 NOTEs**.
Examples, tests, incoming URL checks and the PDF manual pass.

1. CRAN incoming feasibility: "New submission".
2. HTML manual: the local HTML Tidy is too old, and 'V8' is unavailable, so
   HTML validation and math-rendering checks are skipped on this machine.
   These are local toolchain limitations, not claimed successful checks.

## Platform checks

Fresh GitHub Actions checks on this PR cover Linux release/devel/oldrel-1,
Windows release/devel, macOS release, and a separate full-dependency test job.
The local results and checked tarball identity are recorded in
`docs/validation/cran-resubmission-2026-10-18.md`.

## Downstream dependencies

None; this is a new package submission. Earlier manuscript numerical results
must be regenerated because of the corrections. Manuscript sources are
excluded from the package tarball.
