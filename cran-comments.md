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

| Service / environment | Platform | R |
|---|---|---|
| **macOS builder** (mac.r-project.org) | `aarch64-apple-darwin23`, macOS 26.6, Apple M1 | 4.6.1 Patched (2026-07-27 r90311) |
| GitHub Actions | `x86_64-pc-linux-gnu` | release (4.6.1) |
| GitHub Actions | `x86_64-pc-linux-gnu` | devel |
| GitHub Actions | `x86_64-pc-linux-gnu` | oldrel-1 (4.5.3) |
| GitHub Actions | `x86_64-w64-mingw32` | release (4.6.1) |
| GitHub Actions | `aarch64-apple-darwin23` | release (4.6.1) |
| local | `aarch64-apple-darwin25.3.0` | 4.5.3 |

The macOS builder reported `Status: OK` with **no notes at all**, and
`checking tests ... OK`. Result:
<https://mac.R-project.org/macbuilder/results/1790624251-de440352373c50a3/>

The GitHub Actions runs cover both `R CMD check --as-cran` and, in a separate
full-dependency job, the test suite with every suggested package installed and
`NOT_CRAN=true`.

**Submitted, results not yet seen**

* win-builder R-release — uploaded successfully (FTP 226)
* win-builder R-devel — uploaded successfully (FTP 226)

win-builder mails its results to the maintainer address and does not list them
publicly, so the outcome of these two runs has **not been read and is not
claimed here**. The maintainer should check the mail sent to
`toz015@ucla.edu` and record the result before submission.

**Not run**

* R-hub — unavailable in the preparation environment: the `rhub` package cannot
  be installed because its dependency `gert` requires the system `libgit2`
  library, which is absent. No R-hub result is claimed.

## Source tarball checked

The same tarball was used for every check below.

```
filename : TimeMetric_0.2.0.tar.gz
size     : 83756 bytes
sha256   : 0a2adf3abbddf0d9ad532c9958f17300321fe9c57ed6890c0a638863b65974f0
commit   : 8a6c879a04d6638db7bc7ede07df2cc4f7de5db0
```

## R CMD check results

**0 errors, 0 warnings** on every environment checked. All five CI platforms and
the macOS builder report `Status: OK`; the macOS builder reports no notes at all.

Locally the check reports **2 notes**, neither of which appeared on the macOS
builder or on any CI platform:

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
system HTML Tidy and no `V8` — not that anything failed. It is absent from the
macOS builder and from every CI platform, which confirms it is a property of
that machine and not of the package.

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

1. **Read the win-builder results.** Both flavours were submitted successfully;
   their results were mailed to the maintainer and have not been read. R-hub was
   not run (see Test environments).
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
