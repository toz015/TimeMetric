# cran-comments.md — DRAFT

**Status: DRAFT. This package is NOT ready for submission.** See "Outstanding
blockers" below. Nothing in this file may be sent to CRAN until those are
resolved and the pending checks have actually been run.

Every result below was produced locally on 2026-09-25 and is reported as
measured. Results not yet obtained are listed as pending, not predicted.

## Test environments

**Completed**

* local: macOS (aarch64-apple-darwin25.3.0), R 4.5.3 (2026-03-11)
  * `R CMD check --as-cran`
  * `R CMD check --as-cran --run-donttest`

**Pending — NOT yet run, no result claimed**

* win-builder (`devel`, `release`, `oldrelease`)
* macOS builder
* R-hub: Windows, Linux, R-devel
* GitHub Actions matrix, including the full-dependency job with `NOT_CRAN=true`
* R-devel on any platform

The package has never been submitted to CRAN and has no CRAN check history.

## R CMD check results

Both runs: **0 errors, 0 warnings, 2 notes.**

### NOTE 1 — package-related, expected

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Tong Zhu <toz015@ucla.edu>'

New submission
```

This is a first submission, so the note is unavoidable and needs no action.

### NOTE 2 — local toolchain, not a package finding

```
* checking HTML version of manual ... NOTE
Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.
Please obtain a recent version of HTML Tidy by downloading a binary
release or compiling the source code from <https://www.html-tidy.org/>.
Skipping checking math rendering: package 'V8' unavailable
```

This reports that two checks **were skipped** on the check machine, not that
anything failed. The macOS system HTML Tidy is older than the version R expects,
and `V8` is not installed locally. It says nothing about the package and is
expected to disappear on CRAN's machines and on win-builder. It must be
confirmed absent there before submission, rather than assumed.

## Local results

* Test suite: **525 assertions across 107 tests, 0 failures, 0 warnings,
  0 skips** under `testthat::test_local()`.
* Under `R CMD check` the same suite reports 466 passed and 28 skipped. The
  skips are structural, not gaps: 26 are `expect_snapshot*` assertions, which
  testthat skips unless `NOT_CRAN=true`, and 3 are guards on tests that read
  `DESCRIPTION`/`NAMESPACE` from the package source, which does not exist when
  tests run against an installed package. Every test runs in at least one
  environment.
* Examples: `checking examples ... OK` in both runs. 17 examples, 2.44s total,
  no `\dontrun{}` anywhere.
* `spelling::spell_check_package()`: zero findings, against a 72-entry
  `inst/WORDLIST` of technical terms, cited author names, and package
  identifiers.

## Downstream dependencies

None. This is a new package with no reverse dependencies.

## Outstanding blockers — resolve before submitting

1. **Finding 41 — `Gt()` interpolation is non-monotone. UNRESOLVED.**
   `R/pam.Ct.R` derives interpolation indices from the sorted summary table but
   builds the weights from the raw, unsorted input vector. G(t) is therefore not
   monotone non-increasing: on a tied fixture, G(3) = 0.65625 exceeds
   G(2.5) = 0.246094, which is impossible for a survival function. Values at
   exactly observed times are unaffected.

   `Gt()` supplies the inverse-probability-of-censoring weighting used by the
   Brier score. **Whether any reported public metric value is affected has not
   been determined.** The decision — correct the interpolation and re-baseline,
   or withdraw any affected public metric — is deliberately deferred and must be
   made before submission. See `docs/superpowers/findings.md`.

2. **Pending platform checks.** None of win-builder, macOS builder, R-hub,
   R-devel, or the GitHub Actions matrix has been run. The full-dependency CI
   job with `NOT_CRAN=true` is the run that exercises the 26 snapshot assertions
   skipped by a plain check, and must be confirmed green.

3. **History rewrite prepared but not executed.** A combined purge of borrowed
   third-party implementation blobs and committed key material is written and
   dry-run verified in `docs/superpowers/history-purge-plan.md`, and has not been
   run. Nothing has been pushed.

4. **Manuscript not updated.** `paper.md` and `paper.code.Rmd` still reference
   two metrics withdrawn from the package, and `paper.code.Rmd` errors as a
   result. Both are deliberately out of scope for the package work and are
   tracked in `docs/superpowers/phase-2-followups.md`.

Items 3 and 4 do not affect the tarball, which excludes both the manuscript and
the development documentation. Item 1 is a package-correctness question and is
the blocking one.
