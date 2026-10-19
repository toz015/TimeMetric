# Paper branch integration — 18 October 2026

PR #2 was merged into `main` as `d299cff`. Its latest tested head was
`13e26ed`. The paper revision starts at Zian Zhuang's `66b24bb`; the integration
preserves both histories and all earlier numerical and CRAN fixes.

## Resolved behavior

- `tm_predict_coxph()` retains the training-baseline, blockwise prediction and
  strata support from the paper revision. For offset models, it constructs
  relative risks from the fitted design matrix and offsets centered at their
  weighted training mean, matching `survival::survfit()`. Independent tests
  cover Efron, Breslow and exact ties; fitted/held-out predictions; unequal
  training weights; and stratified fits. The earlier offset discrepancy was
  approximately 0.812 in probability units.
- NCC `m` means **number of controls per case**, with `m=1` valid. Matching
  restricts both cases and controls to the same stratum. Each simultaneous
  case samples independently from eligible controls, excluding all cases at
  that time. Events precede censoring at ties. Small risk sets are saturated.
  The paper's examples and tests now follow this contract. Invalid sizes and
  malformed matching vectors produce errors.
- Competing-risk generalized R-squared retains the paper revision's original-
  time observed-event IPCW method and its Cox/KM options. It assigns zero weight
  to censored rows, including rows beyond the outcome horizon, and floors
  censoring probabilities at 1e-4. Independent moment and censoring-reference
  tests specify the changed estimator; other metrics retain their adjustments.
- All PR #2 corrections remain: AFT restricted means/outcome preservation,
  single-covariate prediction, prediction/weight ordering, event-probability
  point R-squared, AUC direction, metric selection/labels, standalone CR AUC,
  required dependencies, missing-value diagnostics, and seed behavior.

Seven textual conflicts were reconciled in source, documentation, ignore rules
and snapshots. Generated Rd help was refreshed. Only the competing-risk
snapshots differ from merged main; numerical tolerances were not relaxed.

The independent follow-up review in `fa9183d` also corrects censoring left
limits at tied times, including point R-squared at a censored horizon. The
competing-risk censoring Cox fit uses Breslow ties to agree with its documented
Nelson-Aalen reference. Cox models with strata-by-covariate interactions use
subject-specific `survfit()` curves and require explicit `new_data`.
The new tests reproduce failures on `1eaa36f` and pass on `fa9183d`.

## Validation

Environment: R 4.5.3 on arm64 macOS Tahoe 26.5.1.

| Check | Result |
|---|---|
| Full source tests, optional dependencies available, `NOT_CRAN=true` | 165 tests; 766 assertions; zero failures, warnings or skips |
| Independent scientific audit | 11/11 pass |
| Additional interaction, tie-method and offset checks | 12 combinations agree with `survfit()` within 2.8e-16; Cox/AFT tied-row permutation checks pass |
| Conditional Breslow censoring with tied times | Agrees with individual reference curves within 1.2e-16 |
| Public API smoke test | All 17 current exports exercised |
| Installed-package runtime dependency isolation | Default two-phase metrics pass with required dependencies only |
| Matched NCC paper script | Hand-calculated one-control example passes |
| Small paper simulations | 30 jobs; zero failures; zero non-finite jobs |
| Spelling | No findings |
| Network-enabled `R CMD check --as-cran --run-donttest` | 0 errors, 0 warnings, 2 NOTEs |

The local NOTEs are “New submission” and unavailable HTML validation/math
rendering checks (old HTML Tidy; missing V8). Examples, installed-package tests,
URL checks and the PDF manual pass. Source-only tests are intentionally skipped
when the installed tarball omits the analysis/development directories; the
full source run above executes them.

All nine GitHub checks passed on `fa9183d`, covering Linux
release/devel/oldrel-1, macOS release, Windows release/devel, all suggested
dependencies, runtime dependencies only, and coverage. The full-dependency
source test job skips no tests; its installed-package check reports
`Status: OK`. The final release-preparation commit updates excluded notes and
this report only; packaged source remains identical to the reviewed revision.
The [PR checks](https://github.com/toz015/TimeMetric/pull/3/checks) record the
subsequent checks on that documentation update.

## Reproduction and limits

From the integrated source checkout:

```sh
NOT_CRAN=true Rscript -e 'testthat::test_local(reporter="summary")'
Rscript docs/validation/scientific-audit.R
Rscript docs/validation/api-smoke.R
Rscript analysis/verify_matched_ncc.R
Rscript analysis/run_simulations.R --mode=smoke --output=/tmp/timemetric-paper-smoke
R CMD build .
R CMD check --as-cran --run-donttest TimeMetric_0.2.1.tar.gz
```

After installing the built tarball, the runtime dependency check is
`Rscript .github/scripts/minimal-dependencies.R` in a fresh R session.

This integration does not rerun the full 2,800-job manuscript simulation,
real-data analyses, or the rendered paper. The earlier paper development record
describes pre-integration calculations. Recompute affected manuscript outputs
with this tested source and retain fresh manifests before publication. The
analysis scripts intentionally reject checkpoint directories with old source
fingerprints. The package tarball excludes manuscript drivers and results.
