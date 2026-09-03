# Baseline Coverage - after Cluster B removal

Measured 2026-09-03 on the package **as committed** -- unlike the first baseline,
which had to be taken on a patched copy because findings.md #6 made the real
package impossible to install.

Coverage rose from 63.09% to 75.38%. The package gained no tests in between; the
increase comes from deleting the 408 lines of unreachable Cluster B code that
were dragging the denominator down.

This is the floor for any further removal: deleting Cluster C or Cluster A must
not reduce coverage of a function that survives.

The two clean-subprocess reproductions and the Rd-reading test skip under covr,
which installs an instrumented copy to a temporary library where a child-process
pkgload::load_all() fails and man/ does not exist. They run normally under
testthat::test_local() and in CI.

```
```
