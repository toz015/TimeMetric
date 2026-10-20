## Resubmission

This is TimeMetric 0.2.2, prepared in response to the review of 0.2.1
reporting two invalid local file URIs in README.md:

* docs/development-2026-10-18.md
* docs/simulation-ici.md

Both now use absolute HTTPS URLs to their existing documents in the public
GitHub repository. The docs directory is excluded from the CRAN archive;
the former relative links therefore failed outside the checkout.

The built-archive audit now checks Markdown links/images, HTML links/images,
Rd URLs, DESCRIPTION URLs and CITATION URLs. It rejects missing local targets,
Windows drive paths, other host-specific paths and links outside the package.
Regression fixtures reproduce both original failures and exercise these edge
cases. External URLs are checked separately with 'urlchecker'; incoming
checks are also enabled in the local R CMD check --as-cran run.

Independent review identified a dependency compatibility issue: 'tdROC'
uses Surv() in a formula without importing it, which fails with older
'survival' versions unless users attach 'survival'. We now require
survival >= 3.8-3, and test competing-risk AUC and Brier score in a fresh
R session without attaching survival. The minimum-dependency CI job uses
survival 3.8-3 explicitly. Statistical function implementations are unchanged.
The CITATION version now derives from DESCRIPTION.

The previous software-name, acronym, reference, TRUE/FALSE and suppressible
console-output corrections remain included. Maintainer details are unchanged.

## Validation

* 193 package tests, 922 passing assertions, no failures, warnings or skips,
  with both survival 3.8-11 and the minimum supported survival 3.8-3.
* 50 passing archive-link regression assertions on macOS; an additional
  case-mismatch assertion runs on case-sensitive filesystems.
* 111 unique document links inspected in the actual source archive; all
  local targets pass. The submitted 0.2.1 archive still reproduces exactly
  the two reported invalid file URIs.
* All remote URLs checked by 'urlchecker' in the same archive pass.
* Network-enabled R CMD check --as-cran --run-donttest: 0 errors,
  0 warnings, 2 NOTEs (new submission; local HTML Tidy too old / V8
  unavailable). Examples, installed-package tests and PDF manual pass.
* R 4.5.3, aarch64-apple-darwin25.3.0, macOS Tahoe 26.5.1.

## Downstream dependencies

None; this is a new package submission.
