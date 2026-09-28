# cran-comments.md — DRAFT

**Status: DRAFT. Not submitted.** The package now passes `R CMD check --as-cran`
on five platform/version combinations in CI with no errors and no warnings, but
the release-coordination items in the checklist below are not finished. Nothing
in this file may be sent to CRAN until they are.

Every result below was measured on 2026-09-28 at commit `951812b`. Nothing is
predicted.

## Submission

New submission. This package has never been on CRAN and has no CRAN check
history.

## Test environments

**Completed — all `Status: OK`**

| Environment | Platform | R |
|---|---|---|
| GitHub Actions | `x86_64-pc-linux-gnu` | release (4.6.1) |
| GitHub Actions | `x86_64-pc-linux-gnu` | devel |
| GitHub Actions | `x86_64-pc-linux-gnu` | oldrel-1 (4.5.3) |
| GitHub Actions | `x86_64-w64-mingw32` | release (4.6.1) |
| GitHub Actions | `aarch64-apple-darwin23` | release (4.6.1) |
| local | `aarch64-apple-darwin25.3.0` | 4.5.3 |

The CI runs cover both `R CMD check --as-cran` and, in a separate
full-dependency job, the test suite with every suggested package installed and
`NOT_CRAN=true`.

**Not yet run — no result claimed**

* win-builder (`devel`, `release`, `oldrel`)
* macOS builder
* R-hub

CI covers Windows, Linux and macOS including R-devel, so these are confirmation
rather than new coverage, but they have not been run and are not claimed.

## Source tarball checked

```
filename : TimeMetric_0.2.0.tar.gz
size     : 83751 bytes
sha256   : 19268afb0ad22a942728dfe8427d3cafe17525e49217d7cfc9be45b776319549
commit   : 951812bec949c5ffe183704853a8848f581e8a66
```

## R CMD check results

**0 errors, 0 warnings.** All five CI platforms report `Status: OK`.

Locally the check reports **2 notes**, neither of which appeared on any CI
platform:

### NOTE 1 — expected for a first submission

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Tong Zhu <toz015@ucla.edu>'

New submission
```

### NOTE 2 — local toolchain only

```
* checking HTML version of manual ... NOTE
Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.
Please obtain a recent version of HTML Tidy by downloading a binary
release or compiling the source code from <https://www.html-tidy.org/>.
Skipping checking math rendering: package 'V8' unavailable
```

This reports that two checks **were skipped** on the local machine — an old
system HTML Tidy and no `V8` — not that anything failed. It is absent from every
CI platform, which confirms it is a property of that machine and not of the
package.

## Test results

* **607 assertions across 124 tests, 0 failures, 0 warnings, 0 skips** under
  `testthat::test_local()`.
* The full-dependency CI job installs every suggested package, sets
  `NOT_CRAN=true`, and **fails if any test skips**. It reports
  `No tests skipped.` and passes.
* Inside `R CMD check` the suite reports 548 passed and 28 skipped. Those skips
  are structural, not coverage gaps: 26 are `expect_snapshot*` assertions, which
  testthat skips unless `NOT_CRAN=true`, and 3 guard tests that read
  `DESCRIPTION`/`NAMESPACE` from the package source, which does not exist when
  tests run against an installed package. The full-dependency job above is what
  exercises all of them.
* Examples: `checking examples ... OK`; 17 examples in 2.35s; no `\dontrun{}`.
* `spelling::spell_check_package()`: zero findings.

## Notable changes since the previous internal build

Two corrections changed reported metric values, both documented in `NEWS.md`:

* **Inverse-probability-of-censoring weights.** `Gt()` is now a reverse
  Kaplan-Meier step function whose values agree exactly with `pec::ipcw()`. It
  previously interpolated between jump times using misaligned indices, inverted
  the censoring indicator on uncensored data, and substituted a positive value
  for a zero censoring survival. `brier_score` from `tm_fit_and_eval()` changes
  on tied data evaluated between jump times.
* **Evaluation-time selection.** The nearest observed time to `t_star` was
  chosen with `which.min()`, which is ambiguous when `t_star` is the median of
  an even-length vector: the two middle times are mathematically equidistant, so
  the winner was decided by a rounding difference in the last bit — and that
  falls differently on platforms with 80-bit `long double` than on those where
  it is plain double. The same data therefore gave different `pseudo_r2_point`,
  `r2_point`, `l2_point`, `brier_score` and `td_auc` on different machines. Ties
  are now broken toward the earlier event time, which is also independent of
  input row order, and all metrics agree bit-for-bit across arm64 macOS, x86_64
  Linux and x86_64 Windows.

No test tolerance was loosened for either correction.

## Downstream dependencies

None. New package, no reverse dependencies.

## Remaining work before submission

No known correctness defect is outstanding, and no `R CMD check` finding
requires action beyond the two notes above.

1. **Confirmation builds.** Run win-builder, the macOS builder and R-hub. CI
   already covers the same platforms, so this is corroboration.
2. **Release coordination.** The remediation work is on `joss-revision` and is
   open as draft PR #1 against `main`. It is not merged. Decide whether to merge
   before submitting, since the CRAN tarball should be built from the intended
   release commit.
3. **Version and date.** Confirm `Version: 0.2.0` and `Date:` in DESCRIPTION are
   what should be released; `Date` currently reads 2026-09-07 and predates this
   work.
4. **Workflow registration.** The CI workflows exist only on `joss-revision`.
   Until that branch reaches the default branch, pushes to `main` trigger no
   checks and `workflow_dispatch` is unavailable there.

## Not part of this submission

`paper.md` and `paper.code.Rmd` still describe two metrics that were withdrawn
from the package, and are tracked in `docs/superpowers/phase-2-followups.md`.
Both are excluded from the tarball by `.Rbuildignore`, so neither affects the
CRAN submission. They are a manuscript task for the JOSS resubmission, **not a
CRAN blocker**.
