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

## Local validation of the manuscript correction

R 4.5.3, aarch64-apple-darwin25.3.0, macOS Tahoe 26.5.1:

* Source suite: 795 assertions across 171 tests; 0 failures, 0 warnings,
  0 skips, including optional prediction backends.
* All nine manuscript R listings execute. PBC and GBSG examples reproduce all
  35 reported table cells on public data. The competing-risk example uses an
  explicitly labelled simulated fixture because the public Framingham teaching
  data are unavailable locally; no claim is made to reproduce its clinical
  results. These must be regenerated after the AFT correction.
* Four Cox fits agree with independent 'survival::survfit' predictions to
  1.3e-15. The manuscript's fixed-scale AFT fixture is checked against the
  analytic competing-risk Weibull CIF at every subject/time combination.
* The analysis scripts now consistently use the manuscript's 12-year FHS
  horizon. The real-data runner loads only dependencies of selected datasets.
  Full manuscript simulation results have not been regenerated.

## R CMD check results

On the updated 0.2.1 source tarball,
`R CMD check --as-cran --run-donttest` reports **0 errors, 0 warnings, 2 NOTEs**.
Examples, installed-package tests and the PDF manual pass.
`_R_CHECK_CRAN_INCOMING_REMOTE_=false` was set because network access was
unavailable; remote incoming checks have not been rerun for this correction.

1. Future file timestamps: unable to verify current time (network unavailable).
2. HTML manual: local HTML Tidy is too old and 'V8' is unavailable, so
   HTML validation and math-rendering checks are skipped.

## Platform checks

The earlier c09e9fa release candidate passed nine GitHub Actions checks.
Those checks predate the new AFT CIF integration correction. Fresh remote
platform checks and review are still required for the updated patch; no CI
success is claimed for it. PR #3 remains unmerged because the automatic
approval service could not execute the authorized publication step.

## Downstream dependencies

None; this is a new package submission. Earlier manuscript numerical results
must be regenerated because of the corrections. Manuscript sources are
excluded from the package tarball.
