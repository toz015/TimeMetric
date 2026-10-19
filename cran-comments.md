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

The final API smoke check also found an unused R-squared calculation when
only Brier score or concordance was requested. Metric subsets now avoid that
calculation. If R-squared is requested with undefined censoring weights, the
Cox/AFT helpers return missing values and one suppressible diagnostic rather
than recycling shorter calibration vectors. Tests cover all four model types;
individual survival predictions remain available.

A manuscript cross-check additionally found that trapezoidal AFT CIF
integration incorrectly assumed zero density at time zero. Weibull models
with scale greater than one lost substantial initial event probability.
Adaptive quadrature on the target-cause probability scale corrects this;
closed-form Weibull identities (including scale five), independent time-scale
integration, grid-invariance tests, and lognormal/loglogistic identities cover
it. The documented time summary integrates one minus the CIF.

The integrated branch also includes the competing-risk C-index correction:
comparable pairs use reciprocal censoring-survival products rather than the
products themselves. Independent pair enumeration covers both pair types,
observations at the horizon, unused zero-survival rows and public evaluation
wrappers. Previous competing-risk C-index results require reevaluation.

## Calibration integration

This candidate adds ordinary and sampling-weighted survival calibration and
unweighted competing-risk calibration. ICI is included by default; E50 and
E90 are optional. The runtime dependencies 'Hmisc', 'pmcalibration' and
'prodlim' are declared, with references and limitations in the help pages.

Three- and four-valued ordinary predictions use linear Cox calibration,
matching unit-weighted results and independent 'survival::survfit' predictions.
Regression tests require silent successful computation, preventing the spline
backend error and unsuppressible frequency-table output found during review.
The manuscript upgrade scripts defer publication until all three datasets
contain ICI, and honor explicit input paths without requiring an existing
published manifest. An end-to-end synthetic pre-ICI migration test preserves
old metrics and input files and runs in the full-dependency CI job.

## Local validation of the integrated corrections

R 4.5.3, aarch64-apple-darwin25.3.0, macOS Tahoe 26.5.1:

* Source suite: 917 assertions across 192 tests; 0 failures, 0 warnings,
  0 skips, including optional prediction backends.
* Separate manuscript migration regression: 17 assertions pass using simulated
  saved predictions, with no clinical inputs or downloads required.
* Prior main-branch manuscript validation (before adding ICI): all nine R listings execute. PBC and GBSG examples reproduce all
  35 reported table cells on public data. The competing-risk example uses an
  explicitly labelled simulated fixture because the public Framingham teaching
  data are unavailable locally; no claim is made to reproduce its clinical
  results. These must be regenerated after the AFT correction.
* Four Cox fits agree with independent 'survival::survfit' predictions to
  1.3e-15. The manuscript's fixed-scale AFT fixture is checked against the
  analytic competing-risk Weibull CIF at all 87,300 subject/time combinations
  (maximum absolute error 9.99e-16).
* The analysis scripts now consistently use the manuscript's 12-year FHS
  horizon. The real-data runner loads only dependencies of selected datasets.
  Full manuscript simulation results have not been regenerated.

## R CMD check results

On the integrated 0.2.1 source tarball, the network-enabled
`R CMD check --as-cran --run-donttest` run on 2026-10-19 reports
**0 errors, 0 warnings, 3 NOTEs**. Examples, installed-package tests, incoming
DESCRIPTION/URL checks and the PDF manual pass.

1. CRAN incoming feasibility: "New submission".
2. Future file timestamps: unable to verify current time using the external
   clock service, despite network access for the other checks.
3. HTML manual: local HTML Tidy is too old and 'V8' is unavailable, so
   HTML validation and math-rendering checks are skipped.

## Platform checks

The calibration integration will be merged only after all nine checks pass
on its final PR head: six platform checks, the full-dependency suite (including
the new manuscript migration regression), required-dependency smoke test and
coverage job. The PR description records the final head and check results.
The earlier numerical corrections were already merged through PR #3.

## Downstream dependencies

None; this is a new package submission. Earlier manuscript numerical results
must be regenerated because of the corrections. Manuscript sources are
excluded from the package tarball.
