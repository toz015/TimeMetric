# Baseline Coverage - characterization suite

Measured 2026-09-02, before any rename, dependency change, or deletion. This is
the floor: spec step 5 (dead-code removal) must not reduce coverage of any
function that survives. Files at 0% are the deletion candidates the spec's call
graph identified as unreachable.

**Caveat.** Coverage could not be measured on the package as committed, because
`R CMD INSTALL` fails outright (findings.md #6: `NAMESPACE` imports
`randomForestSRC`, which is declared in no `DESCRIPTION` field and is not
installed). The measurement below was taken on an identical copy with that one
`importFrom` line removed, which is sufficient to make installation succeed.
Re-measure on the real package once that finding is fixed.

The two clean-subprocess reproductions (findings #11 and #17) skip under covr,
which instruments the source into a temporary library and breaks a plain
`pkgload::load_all()` in a child process. They run normally in `testthat::test_local()`
and in CI.

```
```
