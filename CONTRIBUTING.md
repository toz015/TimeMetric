# Contributing to TimeMetric

Thanks for your interest in TimeMetric. This document covers how to report
problems, how to propose changes, and how to run the checks a change is expected
to pass.

TimeMetric is maintained by a small academic group. Contributions are reviewed as
time allows, and no response time or level of support is promised.

## Reporting a bug

Open an issue at
<https://github.com/toz015/TimeMetric/issues>.

A report is most useful when it includes:

* a minimal reproducible example — the smallest code that shows the problem,
  using simulated data (`tm_sim_cox_weibull()` or `tm_simulate_fine_gray()`) or a
  public dataset such as `survival::pbc`, so no private data is needed;
* what you expected and what happened instead, with the exact error or the
  metric values you obtained;
* the output of `sessionInfo()`.

For anything numerical, please say which metric and which evaluator you used.
Most metrics have several entry points, and the values differ legitimately
between them.

## Reporting a security issue

**Do not open a public issue for a security problem.** That includes credentials
or private data committed to a repository, and any defect that could expose data
belonging to someone else.

Email the maintainer directly at <toz015@ucla.edu> with "TimeMetric security" in
the subject line, and include enough detail to reproduce the problem. Please
allow time for the issue to be assessed before discussing it publicly.

## Proposing a change

1. Open an issue first for anything beyond a typo, so the approach can be agreed
   before you spend time on it.
2. Fork the repository and create a branch off `main`.
3. Make the change, with tests (see below).
4. Run the full check suite locally.
5. Open a pull request describing what changed and why, and link the issue.

Keep pull requests focused. A PR that fixes one defect is far easier to review
than one that also reformats unrelated code.

## Local development setup

TimeMetric requires R >= 4.1.0.

```r
# runtime and check dependencies
install.packages(c(
  "survival", "stats", "utils", "tdROC", "yardstick", "ggplot2",
  "patchwork", "dplyr", "magrittr", "purrr", "tibble", "pec", "expint"
))

# development and test dependencies
install.packages(c(
  "testthat", "withr", "pkgload", "covr", "roxygen2", "spelling"
))

# optional backends, exercised by tests when present
install.packages(c("randomForestSRC", "cmprsk"))
```

Then load the package from source:

```r
pkgload::load_all(".")
```

`randomForestSRC` and `cmprsk` are optional. Tests that need them skip cleanly
when they are absent, so a missing optional backend is visible as a skip rather
than a silent gap.

## Running the checks

```r
# tests
testthat::test_local()

# regenerate NAMESPACE and man/ after changing any roxygen block
roxygen2::roxygenise()

# spelling; add genuine technical terms to inst/WORDLIST, never a misspelling
spelling::spell_check_package(".")
```

```bash
# the full check, against a built tarball rather than the source directory
R CMD build .
R CMD check --as-cran TimeMetric_*.tar.gz

# run examples that are skipped by default
R CMD check --as-cran --run-donttest TimeMetric_*.tar.gz
```

A change is expected to leave `R CMD check` with no errors and no warnings, and
to introduce no new notes.

## Tests are expected with behaviour changes

Any change to what the package computes, returns, or accepts should come with
tests.

* **Fixing a defect**: add a test that fails before the fix and passes after, so
  the defect cannot return unnoticed.
* **Changing existing behaviour**: characterize the current behaviour first, in a
  separate commit, then change it. The diff should show the expectation moving,
  not a value quietly appearing.
* **Refactoring with no intended behaviour change**: the existing tests must pass
  untouched. If you find yourself editing expected values, the refactor changed
  behaviour and needs to be treated as such.

Fixtures are deterministic and seeded (`tests/testthat/helper-simdata.R`). Please
reuse them rather than introducing unseeded random data, and do not change a seed
without regenerating every affected snapshot deliberately.

`tests/testthat/fixtures/metric-baseline.csv` pins the value of every metric on
every public evaluator path. It is a record of values measured at a known commit.
**Do not regenerate it to make a test pass.** If a value moves, that is the
finding.

## Statistical and methodological changes

Changes to an estimator, a weighting scheme, or a metric definition need more
than passing tests. Please include, in the pull request or a linked document:

* the definition you are implementing, with a citation to the source;
* why the current implementation is wrong or incomplete, as evidence rather than
  assertion;
* a reproducible comparison — a script with a fixed seed showing old and new
  values on the same data, and agreement with an independent implementation or a
  published reference value where one exists;
* a note on which reported metric values change, and by how much.

Numerical differences, however small, should be reported rather than absorbed
into a tolerance.

## Code style

Follow the surrounding code. Keep R sources ASCII-only, using `\uXXXX` escapes
where a non-ASCII character is unavoidable — `R CMD check` warns otherwise.
Document exported functions with roxygen, including `@examples` that run.

## Licence

TimeMetric is released under the MIT licence. By contributing, you agree that
your contribution is licensed under the same terms. Do not contribute code
copied from another project unless its licence permits it and you say so
explicitly in the pull request, so that provenance and licensing can be
recorded.
