# Phase 2 Follow-Ups

Deferred items, each recorded with the evidence already gathered so Phase 2 does
not repeat the investigation.

## From the `r_sh` / `r_e` removal (2026-09-24)

The package-side removal is complete. Everything below is **out of scope for
that change** and was deliberately left untouched. Full detail in
`specs/2026-09-24-remove-r-sh-r-e-design.md` §4.

### 1. `paper.code.Rmd` is knowingly broken until revised

`paper.code.Rmd:288-293` builds an explicit `metrics` vector containing `"R_E"`
and passes it to `pam.summary(metrics = metrics)` at line 324. That chunk now
raises:

```
r_e and r_sh were withdrawn before the first CRAN release; see NEWS.md.
```

The file is `.Rbuildignore`d, so no packaging gate is affected. It was not
edited because the manuscript is revised against a settled API, with the
maintainer, not mid-removal.

Required edits, with the structural hazard called out:

| Lines | Structure | Required edit |
|---|---|---|
| 92-100 | Positional: 5 metrics <-> 5 labels | Drop `"R_E"` **and** `expression(R[E])` together |
| 288-293 | Explicit `metrics` vector | Drop `"R_E"`. **This is the chunk that errors** |
| 350-358 | Positional: 5 metrics <-> 5 labels | Drop `"R_E"` **and** `expression(R[E])` together |
| 416-425 | Already excludes `R_E` from both lists | No change needed |
| 555-570 | Name-keyed `metric_pretty` list | Safe key removal |

`metrics_to_plot` / `metric_levels` are matched element-by-element against
`metric_labels`. Dropping a metric without dropping its label at the same index
silently mislabels every downstream facet, with no error raised. The authors
already performed this paired edit correctly for `R_sh` at 416-425 -- follow
that precedent.

The `pam.summary()` calls at 67, 250 and 504 pass no `metrics=` argument and so
default to the full set. They need no edit; their CSVs simply lose two rows.

**Verification constraint.** No `paper.sim*.csv` has ever been committed --
confirmed across all history. The document cannot be executed end-to-end
without first regenerating three simulation datasets, each a 100-iteration
loop. Use a static structural check (every positional metric/label pair equal
in length, every name-keyed lookup resolving, every `factor(levels=)` covering
the values present) plus a reduced-iteration smoke run. Do not claim a full
100-iteration reproduction casually.

### 2. `paper.md` edits, wording unapproved

Verified: `lusa2007estimation`, `schemper2000predictive` and `stare2011measure`
are each cited **exactly once**, all three within the same clause at lines
148-150. Removing that clause orphans all three and nothing else.

* **Lines 95-96** -- delete the `$R_{sh}$` and `$R_E$` rows from the
  right-censored section of the metrics table.
* **Lines 105-107** -- delete the `$R_E$` row from the competing-risks section.
  **This row was already inaccurate**: `tm_survival_eval_cr()` has never emitted
  `r_e`, as `tests/testthat/test-eval-competing-risks.R` asserts. Removing it
  corrects a pre-existing error rather than only reflecting the withdrawal.
* **Lines 146-152** -- cut the two metrics from the prose sentence. Draft:

  ```
  `Performance Metric Module`, including $R^2$-type measures (pseudo $R^2$
  (Li and Wang 2019; Zhuang et al. 2025)), discrimination
  indices (Harrell's C-index (Harrell et al. 1982) and Uno's C-index (Uno
  et al. 2011)), calibration tools (Brier score (Brier 1950; Graf et al.
  1999)), and time-dependent AUC (Heagerty, Lumley, and Pepe 2000).
  ```

* **Lines 222-238** -- delete the three orphaned `.csl-entry` blocks:
  `ref-lusa2007estimation`, `ref-schemper2000predictive`, `ref-stare2011measure`.

The prose is a scientific claim. **The maintainer approves the exact wording**;
the draft above is not approved.

### 3. No `paper.bib` exists

`paper.md` carries rendered citation text plus inline `.csl-entry` divs -- that
is pandoc *output*. JOSS expects `paper.md` with `[@key]` citations alongside a
`paper.bib`. Unrelated to the removal, but it needs addressing before
resubmission.

### 4. Scientific rationale for the withdrawal

The simulation work motivating the removal was carried out by a collaborator
and is not in this repository. Nothing in the package, NEWS, or the findings log
claims the metrics performed poorly -- only that they were not sufficiently
validated for the first CRAN release. Obtain the supporting code or results
before writing the JOSS response letter, and state the rationale there.

### 5. History purge, prepared but not executed

`history-purge-plan.md` records the borrowed-code blobs and the verification
procedure. It requires separate maintainer approval and has **not** been run.
Nothing has been pushed.

## Amended exit criterion

`specs/2026-08-30-timemetric-package-quality-design.md` criterion 10
("Deprecated wrappers verified: `paper.code.Rmd` still runs") cannot be met
until item 1 above is done. It is deferred to Phase 2, not dropped.
