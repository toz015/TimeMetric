# Paper development record — 2026-10-18

This records the current paper checkout. The earlier plans and audits in
`docs/superpowers/` are preserved unchanged; their status labels describe those
historical checkpoints, not necessarily the current implementation.

## Method changes

- Cox predictions use the training fit's baseline cumulative hazard with
  matching predictor centering, weights, ties, strata and factor coding.
  Existing evaluation-grid integration is retained. This is not a new claim
  that unrestricted mean survival is identified beyond observed follow-up.
- Matched nested case-control inclusion probabilities count sampling occasions
  within matching strata, including tied cases and small risk sets. `m` denotes
  total set size (one case and `m - 1` controls); inputs are phase-one cohorts.
- Competing-risks R-squared uses observed-event IPCW at original follow-up
  times, including events beyond the evaluation horizon. Cox (default,
  optionally conditional) and KM censoring estimates use a probability floor
  of 1e-4. Censored rows have zero R-squared weight. This changes CR R-squared
  values relative to the earlier method.
- Point R-squared includes events at the evaluation time (`T <= t`). Ordinary
  survival compares the event indicator with `1 - S`. Independent moment checks
  cover the corrected definitions; old point-metric outputs are not reused.

## Completed paper runs

- Full simulation: 2,800 successful finite jobs, with 100 repetitions per cell:
  Simulation 1, 1,600; Simulation 2 linear/nonlinear, 100 each; Simulation 3,
  1,200. Models, checkpoints, raw metrics and Monte Carlo summaries are retained
  locally under `results/simulations/full/` and excluded from Git.
- PBC, Framingham and GBSG2 were rerun with the updated methods. The currently
  published source is selected by `results/real_data/paper/source_manifest.csv`;
  run-specific models and diagnostics remain under `results/real_data/runs/`.
- Framingham uses actual examination 3 as time zero and a 12-year horizon
  (`tau = t_star = 4380` days), retaining the four model specifications and
  original seed/split. Original-follow-up IPCW weights are saved with the run.
- Simulation 3 errors retain the new saved uncensored finite-sample reference
  means and their Monte Carlo uncertainty; they are not oracle errors.

## Presentation and reproduction

- `figures/run_simulations.R` and `figures/run_real_data.R` are the preferred
  complete entry points. Computation lives in `analysis/`; presentation lives
  in `figures/`. `figures/redraw_original.R` redraws without model refitting.
- Simulation graphics retain the selected prior `.updated.pdf` style, including
  R/L component dots, panel ordering, dimensions and gray/white quantile boxes.
  Real-data graphics use the package's native plotting functions.
- `paper.code.Rmd` now reports the completed outputs rather than maintaining a
  second implementation. It checks simulation completion and real-data export
  checksums, uses canonical metric names, and displays the seven figures and
  three real-data tables. Rendering does not start computations; optional
  redraw delegates to the presentation entry point. README documents rendering.
- Obsolete smoke/verification runs and superseded real-data outputs were
  removed. Current source-run lookup replaces obsolete default paths. The
  previous tutorial and accumulated workflow notes remain in a local archive
  outside the repository; no model or formal result was discarded.

## Validation and limits

- The full package test suite passed after the method changes and directory
  cleanup, including independent prediction, inclusion-probability, IPCW and
  point-moment checks. This is not a new CRAN submission or a cross-platform
  `R CMD check --as-cran` claim.
- All seven paper figures were redrawn successfully after cleanup. Checksums
  confirmed that 3,108 retained calculation-source/statistical files were
  unchanged; presentation regeneration did not rerun the simulations.
- The updated report was rendered against the completed outputs and checked
  for all seven embedded figures and three real-data tables. Rendering and
  figure redraw use the existing statistical outputs, not a fresh 2,800-job run.
- GitHub publication includes source, documentation and regression tests on
  the `paper` branch. Models, checkpoints, generated reports and run logs
  remain local; the package main branch is not a publication target here.
- `paper.md` (the package/JOSS paper) and the earlier history-purge plans have
  not been revised or executed in this work. Their remaining publication tasks
  require separate review; this record does not mark them complete.
