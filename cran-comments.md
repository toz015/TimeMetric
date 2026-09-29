# cran-comments.md — DRAFT

**Status: DRAFT. Not submitted.** The package now passes `R CMD check --as-cran`
on five platform/version combinations in CI with no errors and no warnings, but
the release-coordination items in the checklist below are not finished. Nothing
in this file may be sent to CRAN until they are.

Every result below was measured on the release commit
`f54b6f9d130f25b526b7a1680814fb25b70e53ca`, prepared from the merged `main`
branch on 2026-09-29, or is attributed to the earlier commit it was measured on.
Nothing is predicted.

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

**win-builder** — `x86_64-w64-mingw32`, Windows Server 2022 x64 (build 20348),
gcc 14.3.0.

*Before* the `README.md` link fix, on the tarball from commit `8a6c879`:

* **R-release 4.6.1 (2026-06-24 ucrt)** — `Status: 1 NOTE` (run 19:41:09 UTC)
* **R-devel (2026-09-25 r90590 ucrt)** — `Status: 1 NOTE` (run 19:51:01 UTC)

*After* the fix, on the tarball from commit `145d488`:

* **R-release 4.6.1 (2026-06-24 ucrt)** — `Status: 1 NOTE` (run 20:55:13 UTC).
  The invalid-file-URI item is **gone**; the note now contains only
  "New submission" and the three possible misspellings.
* **R-devel** — submitted twice and processed both times (the upload directory
  was empty on each follow-up check), but **no result e-mail was delivered** on
  either occasion, the second after more than 19 hours. No R-devel result is
  therefore claimed for the post-fix tarball.

  This gap carries very little information, for a reason that can be checked
  directly rather than assumed. The item in question is
  `Found the following (possibly) invalid file URI`, which reports a relative
  link in `README.md` whose target is missing from the built package. The
  packaged `README.md` now contains **no relative links at all** -- its six links
  are all absolute `https://` URLs -- so there is no file URI left for that check
  to report, on any R version. The check is a file-existence test on the package
  contents, not R-version-dependent behaviour, and pre-fix the two flavours
  produced byte-identical note text. `ubuntu-latest (devel)` also passes in CI.

All four logs read so far report **0 errors and 0 warnings**, with every check
other than CRAN incoming feasibility reporting `OK` — including `checking tests`,
`checking examples`, `checking PDF version of manual` and `checking HTML version
of manual`.

That last one matters: the HTML-manual note seen locally (NOTE 2 below) does
**not** appear on win-builder, on the macOS builder, or on any CI platform, which
establishes it as a property of the local machine rather than of the package.

**Not run**

* R-hub — unavailable in the preparation environment: the `rhub` package cannot
  be installed because its dependency `gert` requires the system `libgit2`
  library, which is absent. No R-hub result is claimed.

## Source tarball checked

This is the tarball intended for submission, built from the release commit.

```
filename : TimeMetric_0.2.0.tar.gz
size     : 83666 bytes
sha256   : 0e34791c476efeac6898630ca0d61e43f60086c077ef33a9af8852efd8fe29b4
commit   : f54b6f9d130f25b526b7a1680814fb25b70e53ca
Version  : 0.2.0
Date     : 2026-09-29
```

The external checks recorded above were run on the immediately preceding
content, before `DESCRIPTION`'s `Date` field was updated for the release. The
only difference between those tarballs and this one is that field, plus the
build-time `Packaged:` stamp.

## R CMD check results

`R CMD check --as-cran` on the tarball above, macOS `aarch64-apple-darwin25.3.0`,
R 4.5.3: **0 errors, 0 warnings, 2 notes.**

### NOTE 1 — CRAN incoming feasibility

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Tong Zhu <toz015@ucla.edu>'

New submission
```

Expected for a first submission.

On win-builder this note additionally listed three possible misspellings. It
does not list them locally, and the reason is environmental rather than a
difference in the package: the spell-check portion of this check needs a system
spell-checker, and neither `aspell` nor `hunspell` is installed on the machine
used here (`utils::aspell()` reports "No suitable spell-checker program found").
The local run therefore omits that portion silently. The words are:

* **Brier** — a surname. The Brier score is named after Glenn W. Brier, who
  introduced it in Brier (1950), *Verification of forecasts expressed in terms
  of probability*, Monthly Weather Review 78(1), 1-3. It is the standard name of
  the metric in the survival literature.
* **PAmeasure** — the name of the predecessor R package that this package
  extends and generalises, cited in the Description field for that reason.
* **TimeMetric** — the name of this package.

All three are intentional and correct. They were reported at positions
`Brier (16:38)`, `PAmeasure (20:9)` and `TimeMetric (13:18, 17:5)`, all of which
fall inside the `Description` field.

The `Found the following (possibly) invalid file URI: LICENSE.md` item that
appeared in the pre-fix win-builder logs is gone; see "Fixed during preparation"
below.

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
macOS builder, from win-builder (where `checking HTML version of manual` reports
`OK`), and from every CI platform. It is a property of this machine, not of the
package, and is not expected on CRAN's systems.

### A third note appeared once, transiently

One run of the check on this tarball reported a third note:

```
* checking for future file timestamps ... NOTE
unable to verify current time
```

This check contacts an external time service to confirm that no file timestamp
lies in the future. The message means that call did not succeed, not that any
timestamp is wrong. An immediate re-run of the same check on the same tarball
reported only the two notes above. It is recorded here because it is network
dependent and may recur on any machine, including CRAN's.

### Fixed during preparation

win-builder reported `Found the following (possibly) invalid file URI:
URI: LICENSE.md / From: README.md` on both flavours. This was real: the MIT badge
in `README.md` linked to `LICENSE.md` by relative path, but `LICENSE.md` is
listed in `.Rbuildignore` and so is absent from the built package. The badge now
points to <https://www.r-project.org/Licenses/MIT>. The packaged `README.md`
contains no relative links at all, and the item is absent from the post-fix
win-builder R-release log.

## Test results

Measured on the release commit.

* **607 assertions across 124 tests, 0 failures, 0 warnings, 0 skips** under
  `testthat::test_local()`.
* The full-dependency CI job installs every suggested package, sets
  `NOT_CRAN=true`, and **fails if any test skips**. On the merged `main` it
  reports `No tests skipped.` and passes.
* Inside `R CMD check` the suite reports 548 passed and 28 skipped. Those skips
  are structural, not coverage gaps: 26 are `expect_snapshot*` assertions, which
  testthat skips unless `NOT_CRAN=true`, and 3 guard tests that read
  `DESCRIPTION`/`NAMESPACE` from the package source, which does not exist when
  tests run against an installed package. The full-dependency job exercises all
  of them.
* Examples: `checking examples ... OK`; 17 examples in 2.08s; no `\dontrun{}`.
* `spelling::spell_check_package()`: zero findings against `inst/WORDLIST`.

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

1. **win-builder R-devel post-fix was never delivered** (see Test environments).
   R-release confirms the file-URI item is gone, and the packaged `README.md`
   contains no relative links for that check to flag, so this is recorded as an
   undelivered result rather than an outstanding risk. R-hub was not run.
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
