## Resubmission

This is TimeMetric 0.2.2, resubmitted after the review of 0.2.1 reported two
invalid local file URIs in README.md:

* docs/development-2026-10-18.md
* docs/simulation-ici.md

Both links now use absolute HTTPS URLs to their existing documents in the
public GitHub repository. The docs directory is excluded from the CRAN
archive; the former relative links therefore failed outside the checkout.

To prevent recurrence, CI now scans local link targets in the built source
archive, including Markdown links/images, HTML links/images and Rd URLs.
It separately checks remote URLs with 'urlchecker'. The local-file audit
reproduces both reported failures on the submitted 0.2.1 archive and passes
on 0.2.2 (99 document links inspected). The dedicated remote URL check also
passes on the exact 0.2.2 tarball. It is important to run both checks: the
remote URL checker alone did not flag these missing local files.

No statistical functions, numerical behavior or runtime dependencies changed
from 0.2.1. The previous DESCRIPTION, logical-default and suppressible-console
output corrections remain included. Maintainer details are unchanged.

## Validation

* 917 passing assertions across 192 package tests; no failures, warnings or
  skips, including optional prediction backends.
* All local document targets in the source tarball pass the archive audit.
* All remote URLs extracted by 'urlchecker' from that same archive pass.
* Network-enabled R CMD check --as-cran --run-donttest on 2026-10-20:
  0 errors, 0 warnings, 3 NOTEs. Examples, installed-package tests,
  DESCRIPTION checks and PDF manual pass.
* R 4.5.3, aarch64-apple-darwin25.3.0, macOS Tahoe 26.5.1.

The NOTEs are:
1. New submission (the package is not yet published on CRAN).
2. Unable to verify the current time using the external clock service.
3. Local HTML Tidy too old and V8 unavailable, so HTML validation and
   math-rendering checks are skipped. The separate local/remote link audits
   pass independently of that HTML validation step.

The final integration PR must pass the nine existing GitHub jobs, including
six platform checks and the full-dependency job with the new archive-link
checks, before this tarball is resubmitted. Final results are recorded in the
PR description and submission receipt.

## Downstream dependencies

None; this is a new package submission.
