# TimeMetric — Removing `r_sh` and `r_e`

**Date:** 2026-09-24
**Branch:** `joss-revision`
**Status:** design approved in chat; implementation not started
**Supersedes:** findings 35, 36 (moot), 37, 38, 39 (resolved by removal)

## Background

`r_sh` (Schemper-Henderson explained variation) and `r_e` (Stare-Perme-Henderson
rank-based explained variation) are withdrawn from the package on the
maintainer's decision, ahead of the first CRAN release.

Two independent problems converge on the same removal:

1. **Validation.** The metrics are not sufficiently validated for a first CRAN
   release. Simulation work bearing on their behaviour was carried out by a
   collaborator and is not available in this repository, so **this spec makes no
   claim about how the metrics perform.** The detailed scientific rationale is
   deferred to Phase 2, when the supporting code or results are supplied.
2. **Provenance.** The two implementations are overwhelmingly third-party code
   (findings 38, 39). `R/pam.rsph.R` is 92% verbatim from the unlicensed `Re.r`
   reference implementation; `R/pam.schemper.R` is 91% verbatim from
   `survAUC::schemper()`, which is GPL-2 and therefore incompatible with
   TimeMetric's MIT licence. Deleting the files is the only remedy that removes
   the CRAN blocker, because CRAN's concern is what a package *distributes*, not
   what it exports.

The implementations are deleted outright. They are not archived, retained on
another branch, or published in a companion repository.

## Goal

`r_sh` and `r_e` cease to exist **within the R package** — as code, as accepted
argument values, and as documented features.

## Scope of this phase: the R package only

This phase changes the package and nothing else. The manuscript and its analysis
code are **not** touched: `paper.md` and `paper.code.Rmd` keep their current
contents, and their required edits are recorded as Phase 2 follow-up items in
§4. The package API is stabilised first; the manuscript is revised against a
settled API afterwards, with the maintainer.

## Non-goals

* Any claim about the statistical performance of either metric.
* **Any edit to `paper.md` or `paper.code.Rmd`.** Deferred to Phase 2 (§4).
* Removing the *concepts* from historical discussion. NEWS entries, findings
  logs, and design documents may continue to name `r_e`, `r_sh` and `rsph`.
  Only borrowed implementation code and borrowed documentation are purged.
* Rewriting history, pushing, or submitting to CRAN. Each requires separate
  maintainer approval.

## Decisions taken

| Decision | Rationale |
|---|---|
| Delete `R/pam.rsph.R` and `R/pam.schemper.R` outright | Only route that closes findings 38 and 39 |
| No archive, branch, or companion repo | Maintainer instruction, 2026-09-24 |
| Tombstone error naming the metrics, not a generic "invalid metric" | A user with an existing script learns why it broke |
| Tombstone wording avoids "removed in 0.2.0" | 0.2.0 *is* the first CRAN release; nothing was released to remove it from |
| NEWS states a maintainer decision, not a performance finding | No verifiable evidence in hand |
| Manuscript deferred entirely to Phase 2 | Revise the paper against a stable API, not a moving one |
| History purge scoped to implementation blobs only | Preserve unrelated development history |

---

## 1. Files deleted

| Path | Lines | Reason |
|---|---|---|
| `R/pam.rsph.R` | 504 | Unlicensed `Re.r` port (finding 38) |
| `R/pam.schemper.R` | 141 | GPL-2 `survAUC` port (finding 39) |
| `tests/testthat/test-r-e-reference.R` | 103 | Tests a removed metric |
| `tests/testthat/test-r-sh-definition.R` | 116 | Tests a removed metric |
| `tests/testthat/test-rsph-dispatch.R` | 128 | Tests removed internals |
| `tests/testthat/_snaps/rsph-dispatch.md` | — | Orphaned by the above |

**645 lines of borrowed code leave the distribution.**

Verified before writing this spec:

* Both R files are tagged `@keywords internal` / `@noRd`, so **no `man/*.Rd`
  page needs deleting**. Only `man/tm_survival_eval.Rd` and
  `man/tm_fit_and_eval.Rd` are regenerated.
* Nothing outside the two files calls into them except `tm_survival_eval()` and
  `tm_fit_and_eval()`. `pam.rsph.R` calls only its own internal helpers
  (`my.survfit`, `pam.re`); `pam.schemper.R` calls nothing package-local.
* No shared helper goes dead. `pam.censor()` is still used by
  `tm_survival_eval_cr.R:226`; `Gt()` by `pam.Brier.R` (3 sites);
  `weighted_param()` by both two-phase weight functions.
* The `rms::` string at `pam.schemper.R:95` is inside a comment, not a call.
  No dependency implication.

## 2. Public API changes

**No exported function is removed.** Neither `pam.rsph` nor `pam.schemper` was
ever in `export()`; they reached users only through the two evaluators.

### 2.1 What a caller sees

1. `tm_metric_names()` returns **11 names instead of 13**.
2. `metrics = "r_e"` or `"r_sh"` raises a tombstone error (§2.3).
3. The legacy spellings `R_sh`, `R_E`, `R_sph` no longer resolve — five entries
   leave `tm_metric_aliases()`.
4. `tm_survival_eval()` and `tm_fit_and_eval()` return **two fewer rows** in the
   `Metric` column; `tm_summarize()` returns two fewer columns in wide form.
5. `NAMESPACE` loses five S3 registrations — `pam.rsph.{coxph,aareg,survreg}`,
   `print.rsph`, `summary.rsph`. All were internal-only.

`tm_survival_eval_cr()` is unaffected; it never computed either metric.

### 2.2 Code sites to change

| File | Sites |
|---|---|
| `R/tm_metric_names.R` | Delete 5 alias entries (lines 50–54); update the header comment, which currently explains the `R_sph`/`R_E` merge |
| `R/tm_survival_eval.R` | `valid_metrics` (79–85, two lists), the `r_sh` branch (177–221), the `r_e` branch (223–231), the ordering vector (374–375), roxygen at 8, 10, 23–24, 35–36 |
| `R/tm_fit_and_eval.R` | Default metric vector (81), `r_e` branch (137–140), `r_sh` branch (142–167), roxygen at 3, 15, 16 |
| `NAMESPACE` / roxygen | Drop the 5 S3 registrations; drop `@importFrom stats approx` and `@importFrom stats model.matrix` (§3) |

### 2.3 Tombstone error

`tm_normalize_metrics()` gains a defunct set checked **before** alias
resolution, so both canonical and legacy spellings hit it:

```
r_e and r_sh were withdrawn before the first CRAN release; see NEWS.md.
```

Triggered by `r_e`, `r_sh`, `R_E`, `R_sh`, `R_sph`, under the same
case- and separator-insensitive matching the alias table already uses. Any other
unknown name keeps falling through to the existing `Invalid metrics: ...` error
from `tm_survival_eval.R:95–97`.

## 3. Dependency audit

Counted call sites inside the deleted files against the rest of `R/`:

| Import | In deleted files | Elsewhere | Action |
|---|---|---|---|
| `stats::approx` | 3 | **0** | **Remove** |
| `stats::model.matrix` | 1 | **0** | **Remove** |
| `survival::survfit` | 9 | 10 | Keep |
| `stats::uniroot`, `quantile`, `reshape`, `na.omit`, `complete.cases`, `lm`, `median` | 0 | 1–6 each | Keep |

**No `DESCRIPTION` change.** Neither deleted file uses a package that is not
still required elsewhere. (`rms` was already dropped under finding 36.)

## 4. Manuscript — deferred to Phase 2, recorded here

**Nothing in this section is applied in this phase.** `paper.md` and
`paper.code.Rmd` are left exactly as they are. The analysis below was completed
so that Phase 2 starts from verified facts rather than a fresh investigation.

### 4.0 Consequence of deferring, stated plainly

`paper.code.Rmd:288–293` builds an explicit `metrics` vector containing `"R_E"`
and passes it to `pam.summary(metrics = metrics)` at line 324. Once `r_e` is
removed, that chunk raises the tombstone error. **Leaving the document untouched
therefore leaves it knowingly broken until Phase 2.**

This is an accepted, deliberate trade: the API settles first. Two things follow
from it, and neither is hidden:

* **No gate is affected.** `.Rbuildignore` excludes both `paper.code.Rmd` and
  `paper.md` from the tarball, so `R CMD check` never executes them. The exit
  criteria in §9 are unaffected.
* **An earlier exit criterion is amended.** Criterion 10 of
  `2026-08-30-timemetric-package-quality-design.md` reads "Deprecated wrappers
  verified: `paper.code.Rmd` still runs." That criterion **cannot be met in this
  phase** and is deferred to Phase 2 along with the document itself. It is not
  quietly dropped; §5 records the amendment in the original spec.

### 4.1 `paper.code.Rmd` — required Phase 2 edits

The file uses **positionally paired parallel lists**: `metrics_to_plot` /
`metric_levels` are matched element-by-element against `metric_labels`.
Dropping a metric from one list without dropping its label at the same index
silently mislabels every downstream facet — no error, wrong figure. The authors
already performed this paired edit correctly for `R_sh` at lines 416–425, which
is the precedent to follow.

| Lines | Structure | Required edit |
|---|---|---|
| 92–100 | Positional: 5 metrics ↔ 5 labels | Drop `"R_E"` **and** `expression(R[E])` together → 4 ↔ 4 |
| 288–293 | `metrics` vector passed to `pam.summary(metrics = metrics)` at 324 | Drop `"R_E"`. **This is the chunk that errors if left unedited** |
| 350–358 | Positional: 5 metrics ↔ 5 labels | Drop `"R_E"` **and** `expression(R[E])` together → 4 ↔ 4 |
| 416–425 | Already excludes `R_E` from both lists | No change; optionally clear the commented fragments |
| 555–570 | **Name-keyed** `metric_pretty` list | Safe key removal: drop `"R_E"` from `metrics_to_plot` and the `"R_E" =` entry |

Six already-commented `#"R_sh"` / `#expression(R[sh])` fragments can be cleared
at the same time.

The `pam.summary()` calls at 67, 250 and 504 pass no `metrics=` argument and
default to the full set. They need no edit; their CSVs simply lose two rows.

**Verification constraint for Phase 2.** No `paper.sim*.csv` has ever been
committed — confirmed against all history. The document cannot be executed
end-to-end without first regenerating three simulation datasets, each a
100-iteration loop. Phase 2 verification should therefore be a static
structural check (every positional metric/label pair equal in length, every
name-keyed lookup resolving, every `factor(levels=)` covering the values
present) plus a reduced-iteration smoke run. A full 100-iteration reproduction
should not be claimed casually.

### 4.2 `paper.md` — required Phase 2 edits, wording unapproved

Verified: `lusa2007estimation`, `schemper2000predictive` and `stare2011measure`
are each cited **exactly once**, all three within the same clause at lines
148–150. Removing that clause orphans all three and nothing else.

**Edit 1 — metrics table, right-censored section (lines 95–96).** Delete both:

```diff
-$R_{sh}$ & \times &  &  &  & \times &  &  &  &  \\
-$R_E$  & \times &  &  &  &  & \times &  &  &  \\
```

**Edit 2 — metrics table, competing-risks section (lines 105–107).** Delete:

```diff
-% $R_{sh}$ &  &  &  &  &  &  &  &  &  \\
-% no R_sh,
-$R_E$ & \times &  &  &  &  &  &  &  &  \\
```

> **Pre-existing error, for the maintainer's attention.** This row asserts that
> TimeMetric supplies $R_E$ for competing-risks data. `tm_survival_eval_cr()`
> has never emitted `r_e` — `tests/testthat/test-eval-competing-risks.R:115`
> asserts its absence. The row was already inaccurate before this change.

**Edit 3 — prose (lines 146–152).**

Before:

```
`Performance Metric Module`, including $R^2$-type measures (pseudo $R^2$
(Li and Wang 2019; Zhuang et al. 2025), Schemper and Henderson's
$R_{\text{sh}}$ (Schemper and Henderson 2000; Lusa, Miceli, and Mariani
2007), and $R_E$ (Stare, Perme, and Henderson 2011)), discrimination
indices (Harrell's C-index (Harrell et al. 1982) and Uno's C-index (Uno
et al. 2011)), calibration tools (Brier score (Brier 1950; Graf et al.
1999)), and time-dependent AUC (Heagerty, Lumley, and Pepe 2000).
```

After (draft — **not approved**):

```
`Performance Metric Module`, including $R^2$-type measures (pseudo $R^2$
(Li and Wang 2019; Zhuang et al. 2025)), discrimination
indices (Harrell's C-index (Harrell et al. 1982) and Uno's C-index (Uno
et al. 2011)), calibration tools (Brier score (Brier 1950; Graf et al.
1999)), and time-dependent AUC (Heagerty, Lumley, and Pepe 2000).
```

**Edit 4 — references (lines 222–238).** Delete the three orphaned
`.csl-entry` blocks: `ref-lusa2007estimation`, `ref-schemper2000predictive`,
`ref-stare2011measure`.

### 4.3 Separate Phase 2 observation

`paper.md` carries rendered citation text plus inline `.csl-entry` divs, and no
`paper.bib` exists anywhere in the repository. That is pandoc *output*, whereas
JOSS expects `paper.md` with `[@key]` citations alongside a `paper.bib`.
Unrelated to this removal, but it will need addressing before resubmission.

## 5. Other documentation

| File | Change |
|---|---|
| `README.md` | Example output rows 47–48; metric table rows 64–65; the applicability sentence at 69; the two orphaned references at 148–149 |
| `NEWS.md` | §6 |
| `docs/superpowers/findings.md` | 37, 38, 39 → Resolved by removal; 35, 36 → Moot |
| `docs/superpowers/api-reference-tables.md` | Drop the two metric rows |
| `docs/superpowers/specs/2026-08-30-…-design.md` | Two amendments: §5's `R_sph`/`R_E` merge is superseded, and **exit criterion 10 ("`paper.code.Rmd` still runs") is deferred to Phase 2** per §4.0 |
| `docs/superpowers/phase-2-followups.md` | **New.** A standing list of deferred items, seeded with §4.1, §4.2, §4.3 and the amended criterion 10, so nothing depends on this spec being re-read |

Historical documents (`r-e-implementation-audit.md`, `cluster-c-equivalence-gate.md`,
`commit-review.md`, `pre-rewrite-review.md`, the characterization plan) are
**left intact**. They are the record of how the decision was reached.

## 6. NEWS.md

The unreleased 0.2.0 entry currently contains a long paragraph explaining that
`r_sh` *was fixed* from a degenerate baseline curve. Since 0.2.0 is the first
CRAN release and `r_sh` never ships in it, that paragraph is moot and is
replaced rather than supplemented — NEWS must not both fix and remove the same
metric.

Replacement text, as specified by the maintainer:

> Removed `r_sh` and `r_e` from the pre-release API following a maintainer
> decision to exclude metrics that are not sufficiently validated for the first
> CRAN release.

No performance claim is made.

## 7. Tests

### 7.1 Numerical invariance is the gate

Every surviving metric must produce an identical value. The `eval-survival.md`
snapshot loses two rows; the remaining values are compared **before** the new
snapshot is accepted, never after.

### 7.2 Edits

| File | Sites |
|---|---|
| `test-eval-survival.R` | 30–37, 113–132, 157–159, 188–189 |
| `test-smoke.R` | Internals list at 32; the comment at 19–25 |
| `test-eval-competing-risks.R` | 114–115 already assert absence — they strengthen. Extend the unknown-metric test at 126 to cover `r_e` |

### 7.3 New tests

* `tm_survival_eval()` and `tm_fit_and_eval()` reject `r_e` and `r_sh` with the
  tombstone message.
* Every legacy spelling (`R_E`, `R_sh`, `R_sph`) hits the same tombstone.
* `tm_metric_names()` contains neither, and has length 11.
* The package namespace no longer holds `pam.rsph`, `pam.rsph.coxph`,
  `pam.rsph.aareg`, `pam.rsph.survreg`, `print.rsph`, `summary.rsph`,
  `pam.schemper`, `my.survfit`, or `pam.re`.

## 8. History rewrite — prepared, not executed

Merges with the key-material rewrite already planned as step 8 of the Phase 1
spec. **Not run, not pushed.** Awaits separate maintainer approval.

### 8.1 Scope — implementation code only

The rewrite purges *borrowed implementation code and borrowed documentation*.
It does **not** purge discussion. The strings `r_e`, `r_sh` and `rsph` appearing
in NEWS entries, findings logs, design documents, commit messages and test files
are legitimate development history and are preserved.

**Paths purged** (31 distinct blobs, counted across all history):

| Path | Distinct blobs |
|---|---|
| `R/pam.rsph.R` | 10 |
| `R/pam.re.R` | 4 |
| `R/pam.rsph.metric.R` | 1 |
| `R/pam.rsph_metric.R` | 4 |
| `R/pam.schemper.R` | 5 |
| `man/pam.rsph.Rd` | 1 |
| `man/pam.re.Rd` | 3 |
| `man/pam.rsph.metric.Rd` | 1 |
| `man/pam.schemper.Rd` | 2 |

Three of these were **not** discoverable by reasoning about filenames and were
found only by scanning historical blob *content*:

* `R/pam.rsph.R` was renamed from `R/pam.re.R` at 97% similarity, so the `Re.r`
  port is reachable under both paths.
* `R/pam.rsph.metric.R` and `R/pam.rsph_metric.R` are **two distinct historical
  paths**, differing by a dot versus an underscore. Purging one leaves the other.
* `man/pam.rsph.Rd` carries the borrowed function's generated documentation
  (`Re.fix`, `Re.imp`, `r2nw`, `sen0` all appear in it).

The path list must therefore be regenerated by content scan immediately before
the rewrite, not copied from this table, in case later commits add another.

The `man/*.Rd` pages are included because they are generated *from* the borrowed
sources and reproduce their signatures and argument documentation, which falls
under copied documentation.

**Paths explicitly preserved.** `tests/testthat/test-rsph-dispatch.R`,
`test-r-e-reference.R`, `test-r-sh-definition.R`, `_snaps/rsph-dispatch.md`,
`docs/superpowers/findings.md` and `docs/superpowers/r-e-implementation-audit.md`
stay in history. They were checked for borrowed content directly:

* **No copied implementation.** Zero matches for the Italian identifiers of the
  `survAUC` port.
* **Identifiers appear only as assertion targets, not as code.** The tests do
  name `dRti` and `my.survfit`, but only inside
  `expect_identical(names(res), c("times", "Rti", "dRti"))` — asserting the
  column names of a return value — and inside a namespace inventory list at
  `test-smoke.R:32`. Naming a function is not reproducing it.
* **The snapshot holds output, not source** — base64-serialised R value and name
  vectors.
* **The findings log and audit quote the borrowed identifiers as evidence**,
  which is the record of how the licensing problem was diagnosed.

All of this is original work or legitimate analysis and is preserved.

### 8.2 Verification — blobs and fragments, not words

A word search for `r_e` or `rsph` is the wrong test. Verification confirms:

1. **Blob unreachability.** Each of the 31 blob hashes is recorded before the
   rewrite, then confirmed unreachable from any ref afterwards — absent from
   `git rev-list --objects --all`, with `git cat-file -e` failing.
2. **Fragment absence, scoped to the code surface.** A content scan over all
   reachable blobs **under `R/` and `man/` only** must return nothing for:
   * `survAUC` port — `tempi.eventi`, `f.assegna.surv`
   * `Re.r` port — `newvare`, `incom.sum`, `Re.imp`, `r2nw`, `dRti`
3. **Tarball scan.** The same fragment set is absent from the built `.tar.gz`.

**Why the scan is scoped, and why these seven fragments.** The scoping and the
fragment set were both chosen by testing them against real history, not by
inspection. An unscoped scan produces false positives on every candidate
fragment but one: `docs/superpowers/findings.md` and
`r-e-implementation-audit.md` quote the Italian identifiers as evidence, and
`test-smoke.R` and `test-r-e-reference.R` name internals like `my.survfit` and
`dRti`. Those are the legitimate historical discussion this rewrite preserves,
so flagging them would be a defect in the check, not a finding.

Scoped to `R/` and `man/`, all seven fragments return clean outside the purged
paths — verified against every commit in the repository. Two candidates are
excluded deliberately: `my.survfit` and `n.incom` also appear in
`R/pam.rsph_metric.R`, which is itself on the purge list, so they cannot
distinguish a successful rewrite from a failed one.

Pre-rewrite file hashes, for the record:

```
3092e89f71a11ea3567e3b77010b0e34d06b4530a25131097e4f648029a2e858  R/pam.rsph.R
1e0074171d39acfaf11bef7314dcba8e83fc58dca4f065412ecf410584c869aa  R/pam.schemper.R
```

## 9. Exit criteria

1. `devtools::test()` — full suite passes, none skipped.
2. `R CMD build` then `R CMD check --as-cran` on the tarball: 0 errors,
   0 warnings, no new notes against the pre-change baseline.
3. Every surviving metric's value byte-identical to the pre-change snapshot.
4. `tm_metric_names()` has length 11; the tombstone fires for all five spellings.
5. Fragment scan of the built tarball is clean (§8.2.3).
6. `git diff --stat` touches **no** manuscript file. `paper.md` and
   `paper.code.Rmd` are byte-identical to their pre-change state.
7. The Phase 2 follow-up list exists and carries every deferred item from §4.
8. History rewrite prepared and reviewed, **not executed**.

## Approval gates

| Action | Status |
|---|---|
| Implement §1–§3, §5–§7 (package only) | Approved in chat, 2026-09-24 |
| Edit `paper.md` or `paper.code.Rmd` | **Out of scope this phase** — Phase 2, §4 |
| Run the history rewrite | **Blocked** — separate approval |
| Push to origin | **Blocked** — separate approval |
| Submit to CRAN | **Blocked** — separate approval |
