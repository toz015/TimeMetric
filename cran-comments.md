# cran-comments.md — DRAFT

**Status: DRAFT. Not submitted, and not yet ready to submit.** The package passes
its local checks, but the platform checks listed under "Pending" have not been
run. Nothing in this file may be sent to CRAN until those are complete.

Every result below was measured locally on 2026-09-26 at commit `8ae6d95`.
Results not yet obtained are listed as pending, never predicted.

## Submission

New submission. This package has never been on CRAN and has no CRAN check
history.

## Test environments

**Completed**

* local: macOS (aarch64-apple-darwin25.3.0), R 4.5.3 (2026-03-11)
  * `R CMD check --as-cran` — 0 errors, 0 warnings, 2 notes
  * `R CMD check --as-cran --run-donttest` — 0 errors, 0 warnings, 2 notes

**Pending — not run, no result claimed**

* win-builder: `devel`, `release`, `oldrelease`
* macOS builder
* R-hub / R-devel: Windows, Linux, R-devel
* GitHub Actions matrix, including the full-dependency job with `NOT_CRAN=true`

## Source tarball checked

```
filename : TimeMetric_0.2.0.tar.gz
size     : 81725 bytes
sha256   : 3d9405e3b34a378db0ebd7fef58f0a82b23e44fc87d1bb78084c041059a82403
commit   : 8ae6d952caa59ba8e9d619c87784dc2ef509b658
```

## R CMD check results

Both invocations: **0 errors, 0 warnings, 2 notes.**

### NOTE 1 — expected for a first submission

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Tong Zhu <toz015@ucla.edu>'

New submission
```

Unavoidable and needs no action.

### NOTE 2 — local toolchain, not a package finding

```
* checking HTML version of manual ... NOTE
Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.
Please obtain a recent version of HTML Tidy by downloading a binary
release or compiling the source code from <https://www.html-tidy.org/>.
Skipping checking math rendering: package 'V8' unavailable
```

This reports that two checks **were skipped** on the check machine, not that
anything failed: the macOS system HTML Tidy is older than R expects and `V8` is
not installed locally. It says nothing about the package and is expected to be
absent on CRAN's machines and on win-builder. That must be confirmed there
rather than assumed.

## Local results

* **Tests: 562 assertions across 116 tests, 0 failures, 0 warnings, 0 skips**
  under `testthat::test_local()`.
* Under `R CMD check` the same suite reports 503 passed and 28 skipped. The skips
  are structural, not coverage gaps: 26 are `expect_snapshot*` assertions, which
  testthat skips unless `NOT_CRAN=true`, and 3 guard tests that read
  `DESCRIPTION`/`NAMESPACE` from the package source, which does not exist when
  tests run against an installed package. Every test runs in at least one
  environment. The `NOT_CRAN=true` CI job is what exercises the 26.
* Examples: `checking examples ... OK`. 17 examples, 2.03s total, no
  `\dontrun{}` anywhere.
* `spelling::spell_check_package()`: zero findings, against a 72-entry
  `inst/WORDLIST` of technical terms, cited author names and package
  identifiers.

## Notable changes since the last internal build

The inverse-probability-of-censoring weights behind the Brier score were
corrected. `Gt()`, the Kaplan-Meier estimate of the censoring distribution, is
now a reverse-Kaplan-Meier step function whose values agree exactly with
`pec::ipcw()`; it previously interpolated between jump times using misaligned
indices, inverted the censoring indicator on uncensored data, and substituted a
positive value for a zero censoring survival. `brier_score` from
`tm_fit_and_eval()` changes on tied data evaluated between jump times. See
NEWS.md.

## Downstream dependencies

None. New package, no reverse dependencies.

## Unresolved blockers

**None at the package level.** No known correctness defect is outstanding, and
no `R CMD check` finding requires action beyond the two notes above.

Remaining work before submission is verification and release coordination, not
package defects — see the checklist below.

## Pre-submission checklist

1. **Platform checks.** Run win-builder (`devel`, `release`, `oldrelease`),
   macOS builder, and R-hub/R-devel. Confirm in particular that NOTE 2
   disappears, since it is a property of this machine.
2. **CI.** Confirm the GitHub Actions matrix is green, including the
   full-dependency job with `NOT_CRAN=true`, which runs the 26 snapshot
   assertions a plain `R CMD check` skips.
3. **Release coordination.** A combined purge of borrowed third-party
   implementation blobs and committed key material is prepared and dry-run
   verified in `docs/superpowers/history-purge-plan.md`, and has **not** been
   run. It rewrites published history and requires a force-push, so it must be
   sequenced with collaborators before any release. This is a repository
   coordination task; it is not an `R CMD check` finding and does not affect the
   tarball, which excludes the development documentation.
4. **Version and date.** Confirm `Version` and `Date` in DESCRIPTION are what
   should be released.

## Not part of this submission

`paper.md` and `paper.code.Rmd` still describe two metrics that were withdrawn
from the package, and are tracked in `docs/superpowers/phase-2-followups.md`.
Both are excluded from the tarball by `.Rbuildignore`, so neither affects the
CRAN submission. They are a manuscript task for the JOSS resubmission, **not a
CRAN blocker**.
