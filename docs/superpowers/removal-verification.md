# Verification — `r_sh` / `r_e` removal

**Date:** 2026-09-24
**Commit verified:** `c7aaac36bccd26af41966e4b2f2d0e144e17d5cd`
**Method:** fresh `git clone` of the repository into a scratch directory,
checked out at that commit, verified there rather than in the working directory.
**Platform:** R 4.5.3 (2026-03-11), aarch64-apple-darwin25.3.0.

## Test results

Three environments, because they legitimately differ. Every number below is
measured, not expected.

| Environment | Pass | Fail | Warn | Skip |
|---|---|---|---|---|
| `testthat::test_local()` (source tree present) | **498** | 0 | 0 | **0** |
| `R CMD check --as-cran` (default) | 441 | 0 | 0 | 28 |
| `R CMD check --as-cran` with `NOT_CRAN=true` | 486 | 0 | 0 | 3 |

100 test blocks across 12 files.

**The skips under `R CMD check` are not a gap in coverage.** They have two
structural causes, both pre-existing and unrelated to this change:

* **26 "On CRAN" skips (default run only).** testthat skips `expect_snapshot*`
  assertions unless `NOT_CRAN=true`. Setting it removes all 26, which is what
  the third row shows and what the full-dependency CI job does.
* **3 "no package source tree" skips (both check runs).** These are
  `skip_without_source_tree()` guards on tests that read `DESCRIPTION`,
  `NAMESPACE` or `man/` from the package *source*, which does not exist when
  the tests run against an installed package inside `R CMD check`. Those three
  run normally under `test_local()`, which is why row 1 shows zero skips.

Every test therefore executes in at least one environment; none is never-run.

**Optional backends genuinely exercised, not skipped away.** `randomForestSRC`
3.7.0 and `cmprsk` 2.2.12 are both installed, so the `skip_if_not_installed()`
guards in `test-optional-backends.R` passed through: 28 assertions, 0 skips
under `test_local()`.

## Source tarball

```
filename : TimeMetric_0.2.0.tar.gz
size     : 75008 bytes
sha256   : 2c08af9babaed4fc748911e43d0ae690ba694e8ba18cca210a52e40e40ff3fc1
```

Digest computed with `shasum -a 256` over the exact artifact that `R CMD check
--as-cran` was run against, and independently recomputed over the same file
afterwards with the same result. It is **64 lowercase hex characters**, matching
`^[0-9a-f]{64}$`. This is not a rebuilt artifact: the checked tarball was
retained and re-hashed in place, so the digest belongs to the file that was
actually checked.

## `R CMD check --as-cran`

`Status: 2 NOTEs` — **0 errors, 0 warnings, 2 notes**. 54 check stages ran.

### NOTE 1 — CRAN incoming feasibility

```
Maintainer: 'Tong Zhu <toz015@ucla.edu>'

New submission
```

Expected and unavoidable: the package has never been on CRAN. This is the note
the spec's exit criterion 1 anticipates as individually justified.

### NOTE 2 — HTML version of manual

```
Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.
Please obtain a recent version of HTML Tidy by downloading a binary
release or compiling the source code from <https://www.html-tidy.org/>.
Skipping checking math rendering: package 'V8' unavailable
```

A local toolchain limitation, not a package defect: macOS ships an old HTML Tidy
and `V8` is not installed. The note reports that two *checks were skipped*, not
that anything failed. It would not appear on CRAN's own machines.

Key stages, all `OK`: package dependencies, installed package size, DESCRIPTION
meta-information, code files for non-ASCII characters, dependencies in R code,
S3 generic/method consistency, R code for possible problems, Rd files, Rd
cross-references, unstated dependencies in examples and in tests, **examples**,
**tests**.

## Comparison with `baseline-check.txt`

The baseline is **not a like-for-like comparison, and saying "no new notes"
would be misleading.**

| | Baseline (2026-09-02) | Now |
|---|---|---|
| Version | 0.1.0 | 0.2.0 |
| Command | `R CMD check --no-manual` | `R CMD check --as-cran` |
| Status | **1 ERROR** | **2 NOTEs** |
| Check stages reached | **4** | **54** |

The baseline **aborted at the dependency stage**:

```
* checking package dependencies ... ERROR
Namespace dependencies missing from DESCRIPTION Imports/Depends entries:
  'magrittr', 'pec', 'purrr', 'randomForestSRC', 'tibble'
```

Because it aborted after 4 stages, it recorded **no notes at all** — not because
the package was clean, but because 50 later stages never ran. The two notes now
present therefore cannot be classified as "new relative to baseline"; there was
no baseline for those stages to be measured against. The baseline file says so
itself: *"This is the recorded STARTING POINT, not a passing check."*

What the comparison does establish: the ERROR that made the package
uninstallable is gone, and the check now runs to completion with nothing worse
than two notes, one of which is unavoidable for a first submission and the other
an artifact of the local toolchain.

## Tarball scans

| Scan | Result |
|---|---|
| `Rplots.pdf` | **ABSENT** |
| Deleted + historical source filenames (13 checked: `pam.rsph.R`, `pam.re.R`, `pam.rsph.metric.R`, `pam.rsph_metric.R`, `pam.schemper.R`, their four `.Rd` pages, three deleted test files, one snapshot) | **all ABSENT** |
| `paper.md`, `paper.code.Rmd` | **ABSENT** (`.Rbuildignore`d) |
| Git metadata and dev files (`.git`, `.github`, `docs/`, `.Rproj`, `.Rproj.user`, `.gitignore`, `Screenshot*`, `.Rcheck`, `superpowers`) | **all ABSENT** |
| Borrowed-code fragments across `R/` and `man/` (9 markers) | **all clean** |

### One scan hit, investigated and dismissed

A whole-tarball content scan for 14 fragments returned a single match:
`my.survfit` in `tests/testthat/test-metric-tombstone.R`. In context:

```r
test_that("the borrowed implementations are gone from the namespace", {
  ns <- asNamespace("TimeMetric")
  for (nm in c("pam.rsph", "pam.rsph.coxph", "pam.rsph.aareg",
               "pam.rsph.survreg", "print.rsph", "summary.rsph",
               "pam.schemper", "my.survfit", "pam.re")) {
    expect_false(exists(nm, envir = ns, inherits = FALSE),
```

The name appears inside an assertion that the function is **absent**. Naming a
function in order to prove it is gone is not reproducing it. This is exactly the
false-positive class the spec anticipated when it scoped fragment verification
to `R/` and `man/`; restricted to that code surface, all 9 markers are clean.

## Working tree

Clean at `c7aaac3` after verification; the fresh clone was used so the working
directory was never touched by the build or check.

## Required before CRAN submission

**The full-dependency CI job with `NOT_CRAN=true` must be confirmed green after
the eventual push and before any CRAN submission.** Local verification covers
this machine only. A plain `R CMD check` exercises 441 of 498 assertions, because
testthat skips snapshot comparisons unless `NOT_CRAN=true`; the CI job is what
closes that gap across the platform matrix. A green local run is not a
substitute, and the push has not happened yet.

## Not done

No history rewrite, no push, no force-push, no tag, no release, no CRAN
submission. Task 9 (history purge) remains prepared-only and awaits separate
maintainer approval.
