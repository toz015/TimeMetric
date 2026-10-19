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

The integrated paper revision also documents and tests observed-event
inverse-probability-of-censoring weights for competing-risk R-squared, with
Cox or Kaplan-Meier censoring estimation. Cox prediction now handles training
offset centering with case weights and strata. Nested case-control `m` remains
the number of controls per case. Independent reference tests cover these
changes and retain the earlier numerical regression checks.

The final review also corrected censoring weights at tied follow-up times:
all subjects tied at a time now share the same censoring survival left limit.
Point R-squared is tested when censoring occurs exactly at the evaluation
horizon. The competing-risk Cox censoring model explicitly uses Breslow ties,
matching the documented no-covariate Nelson-Aalen reference. Cox prediction
supports strata-by-covariate interactions when evaluation data are supplied.
Independent regression tests reproduce these defects on the previous revision
and pass after the corrections.

## Local validation

R 4.5.3, aarch64-apple-darwin25.3.0, macOS Tahoe 26.5.1:

* Source test suite with NOT_CRAN=true: 766 assertions across 165 tests;
  0 failures, 0 warnings, 0 skips, including optional prediction backends.
* Independent numerical checks: 11/11 pass against fitted-model survival
  predictions, numerical integration, and hand-computable sampling/ranking
  references.
* Additional checks cover 12 combinations of Cox interactions, tie methods
  and offsets, including unequal case weights where supported. Probabilities
  agree with 'survival' reference predictions within 2.8e-16. Cox and AFT
  metric results are invariant to permutations of tied observations;
  conditional Breslow censoring probabilities match reference curves.
* All 17 current exported functions exercised by a public/simulated-data smoke
  test. Default two-phase metrics also pass with only required dependencies.
* All 30 small paper simulation jobs returned finite results without job
  failures. These are workflow checks, not the full manuscript simulation.

## R CMD check results

On the built 0.2.1 source tarball, with network access:
`R CMD check --as-cran --run-donttest` reports **0 errors, 0 warnings, 2 NOTEs**.
Examples, tests, incoming URL checks and the PDF manual pass.

1. CRAN incoming feasibility: "New submission".
2. HTML manual: the local HTML Tidy is too old, and 'V8' is unavailable, so
   HTML validation and math-rendering checks are skipped on this machine.
   These are local toolchain limitations, not claimed successful checks.

## Platform checks

All nine GitHub Actions checks passed on the reviewed implementation
(`fa9183d`): Linux release/devel/oldrel-1, Windows release/devel, macOS release,
full suggested dependencies, required dependencies only, and test coverage.
The full-dependency source test job skips no tests; its installed-package
check reports `Status: OK`. Source-only checks are intentionally skipped
when their analysis/development inputs are excluded from the installed
package. The release-preparation changes only update excluded CRAN notes
and the validation report; the packaged implementation is unchanged.
See `docs/validation/paper-integration-2026-10-18.md` for reproduction details.

## Downstream dependencies

None; this is a new package submission. Earlier manuscript numerical results
must be regenerated because of the corrections. Manuscript sources are
excluded from the package tarball.
