# TimeMetric 0.2.1 CRAN resubmission validation

This revision combines the CRAN-response branch in PR #2 with the numerical
corrections developed in the earlier local audit. The manuscript workspace
is preserved separately. Numerical results in existing papers must be
regenerated against the corrected package.

## Corrections and independent checks

| Contract | Before correction | Corrected result |
|---|---:|---:|
| Cox curves invariant to changed validation event labels | max difference 0.545968 | 0 |
| Cox curves agree with predictions from the fitted training model | max difference 0.329778 | 0 |
| Weibull restricted means agree with survival quadrature | max error 0.766558 | below 2e-11 |
| Lognormal restricted means agree with survival quadrature | max error 0.982042 | below 2e-9 |
| AFT prediction preserves original outcomes | 17 changed labels per model | 0 |
| NCC control weight: one control drawn from two eligible controls | Inf | 2 |
| NCC control weight: both eligible controls sampled | 2 | 1 |
| Explained variation invariant to joint subject permutation | 0.0363 vs 0.2794 | both 0.2794 |
| Point pseudo R-squared matches event-probability regression | 0.0934 | 0.5862, matching reference |
| Reversed held-out AUC matches empirical case/control ranking | 0.94 | 0.06, matching reference |

Eleven checks pass (the AFT outcome check runs for both distributions).
`scientific-audit.R` reproduces this table and exits nonzero on disagreement.
The testthat suite additionally covers single-predictor AFT data frames,
metric-subset names/order, standalone competing-risk AUC, matched/tied NCC
risk sets, and two-phase metric selection. Only seven affected numerical
fixture rows and corresponding snapshots changed; tolerances were not relaxed.

## CRAN-specific regressions

- Missing observations are excluded consistently before model fitting and
  evaluation; the warning is suppressible and the AUC backend writes no
  console output for this case.
- WLS and CIF fit failures emit suppressible warnings.
- Fine-Gray simulation and plot subsampling respect the caller's RNG stream;
  explicit seeds remain reproducible. Help now distinguishes resetting the
  stream from advancing it through random draws.
- `tidyr` is a required, namespace-imported dependency of the default
  two-phase AUC. `.github/scripts/minimal-dependencies.R` installs no software
  and restricts a fresh R process to the declared hard-dependency closure
  before exercising all five default two-phase metrics.
- Executable R source contains no T/F symbols or print()/cat() calls.
- Spelling checks report no findings; generated Rd usage matches the source.

## Local validation

Environment: R 4.5.3, arm64 macOS Tahoe 26.5.1.

- 661 assertions across 141 test blocks, zero failures/warnings/skips, with
  `NOT_CRAN=true` and optional backends installed.
- All 11 independent numerical checks pass.
- All 17 current exports exercised with public/simulated data.
- Default two-phase metrics pass in a library restricted to required packages.
- Line coverage: 87.10% in the measured run.

The built source tarball passed a network-enabled
`R CMD check --as-cran --run-donttest`: **0 errors, 0 warnings, 2 NOTEs**.
Incoming feasibility reports "New submission"; local HTML validation and math
rendering are skipped because the installed Tidy is old and V8 is absent.
URL checks, examples, tests and the PDF manual passed.

Checked artifact: `TimeMetric_0.2.1.tar.gz`.
SHA-256: `0cb26e62fc34b4052e99a9b0a1ec99ad08a0c85499650a061b6c76bf2e8434a9`.
The report and submission comments are excluded from the package tarball.

## Reproduction

```sh
NOT_CRAN=true Rscript -e 'testthat::test_local()'
Rscript docs/validation/scientific-audit.R
Rscript docs/validation/api-smoke.R
R CMD build .
R CMD check --as-cran --run-donttest TimeMetric_0.2.1.tar.gz
R CMD INSTALL TimeMetric_0.2.1.tar.gz
Rscript .github/scripts/minimal-dependencies.R
```

The GitHub matrix covers Linux R release/devel/oldrel-1, Windows R
release/devel, and macOS R release, plus full-dependency and minimal-dependency
jobs. Results for the current revision are available in
[PR #2 checks](https://github.com/toz015/TimeMetric/pull/2/checks).
Local checks are not a substitute for the current revision's platform results.

This validation concerns package behavior and the tested numerical contracts.
It does not reproduce the manuscript's full simulation study or real-data
results, or certify every statistical estimator/boundary case. The function
inventory is in [function-guide.md](function-guide.md).
