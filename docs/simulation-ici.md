# Run the simulations with ICI

The code includes ICI in Simulation 1, both parts of Simulation 2 and
Simulation 3. These are the existing right-censored designs. The data-generating
mechanisms, sample sizes and evaluation horizons are unchanged. No
competing-risks or sampling-design simulation is added here.

Simulation 1 uses common random numbers. Every scenario fits its Cox model to
training data generated with seed 1234 (the Simulation 2 training seed), and
repetition r of every scenario uses test seed r. The scenarios therefore differ
only through v and beta. Because v enters only through a monotone transform of
event time, the Cox fit, Brier score, AUC, C-index and ICI are the same across v
for a given beta, while R-squared changes. The Simulation 1 Cox training fits
use `timefix = FALSE`: with 1/v = 2 the default tie correction merged
thousands of distinct times. Simulation 3 uses the Simulation 1 beta = 1, v = 1
fit, which is now the seed-1234 model. Small differences across v can remain
for beta = 1.5 at 1/v = 2, where the C-index and calibration backends apply
their own default tie correction to test times.

ICI uses `tm_calibration()` with `nk = 5`. Its outcome horizon is the existing
`t_star` for each job. Predicted risks are requested directly from the fitted
model at exactly that horizon, rather than selected from the test-follow-up
grid. Cox probabilities use the fitted training baseline; AFT probabilities
use the fitted distribution, and RF uses its prediction adapter at `t_star`.
ICI values are saved without rounding. Other metrics
retain their previous rounding. E50 and E90 are not added to simulation outputs.

This time-alignment correction changes the computation version. Preserve
`ici_full` and use a new output directory, such as `ici_crn`, for a new
run. To compare the old and corrected ICI on the saved uncensored Simulation 3
replicates without refitting prediction models or recomputing other metrics:

```sh
Rscript --vanilla analysis/check_sim3_ici_time_alignment.R --input=results/simulations/ici_full --output=results/simulations/ici_time_alignment_check
```

The comparison retains the reference definition: each version uses its own
mean across the 100 uncensored `n=3000` repetitions at the same horizon.

## Terminal

Run from the `TimeMetric-paper` directory:

```sh
# Print the job plan only. No model fitting or evaluation occurs.
Rscript --vanilla figures/run_simulations.R --mode=plan

# Full run with 100 repetitions per cell, including ICI.
Rscript --vanilla figures/run_simulations.R --mode=full --output=results/simulations/ici_crn
```

On the current development machine, `pmcalibration` is available in
`/tmp/timemetric-calibration-lib`. If that library is still being used, the
full-run command is:

```sh
R_LIBS=/tmp/timemetric-calibration-lib Rscript --vanilla figures/run_simulations.R --mode=full --output=results/simulations/ici_crn
```

If the temporary library is no longer available, install `pmcalibration` in your
usual R library before launching. The remaining analysis dependencies are
listed in the repository README.

Repeat the same full-run command to resume. Source files, settings and package
versions must match that run. Historical results in `results/simulations/full/`
are preserved. A new calculation version requires a different output directory.

Optional selectors can be appended to the full-run command:

```sh
# Only the two Simulation 2 parts.
--part=sim2_linear,sim2_nonlinear

# Only repetitions 1 through 10.
--repeats=1:10
```

## RStudio console

Set the working directory to `TimeMetric-paper`, then run:

```r
# Include the development dependency library if it exists.
calibration_library <- "/tmp/timemetric-calibration-lib"
if (dir.exists(calibration_library)) {
  .libPaths(c(calibration_library, .libPaths()))
}
stopifnot(requireNamespace("pmcalibration", quietly = TRUE))

# Let the child R process use the same package libraries.
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
status <- system2(
  file.path(R.home("bin"), "Rscript"),
  c("--vanilla", "figures/run_simulations.R", "--mode=full",
    "--output=results/simulations/ici_crn")
)
stopifnot(status == 0L)
```

## Outputs and later report rendering

The new output directory contains:

- `tables/metrics_raw.csv` and `tables/metrics_summary.csv`, including ICI.
- `tables/calibration_diagnostics.csv`, with per-job and per-model horizons,
  prediction-grid times, requested knots, calibration methods and errors.
- `tables/sim3_references.csv`, `sim3_errors.csv` and `sim3_accuracy.csv`,
  including ICI. The reference remains the mean from the uncensored tests at
  the larger sample size. It is not assumed to equal zero.
- `run_status.csv`, `completion.csv`, models and resumable checkpoints.
- Simulation figures with an ICI panel. Simulation 2 retains its R/L dots.

Calibration failures retain the other metrics and leave a missing ICI with a
recorded warning. Such jobs are not considered successfully complete and are
retried when the run is resumed. The report requires complete ICI results.

The runner generates figures automatically. To redraw completed results later:

```sh
Rscript --vanilla figures/redraw_original.R --section=simulations --sim-output=results/simulations/ici_crn
```

Both commands use the manuscript presentation: gray bars without numeric bar
labels, the original metric/model order, R/L component dots in Simulation 2,
and gray/white horizon boxplots in three censoring rows for Simulation 3.
ICI is appended in a right-hand column. The existing four bar-chart panels
keep their original 2-by-2 positions and physical dimensions; Simulation 3
keeps its original five columns and adds ICI as the sixth.
The lower-level `analysis/run_simulations.R` produces generic diagnostics;
the figure entry points preserve these under `figures/diagnostics/analysis_format/`.

An incomplete simulation part does not prevent completed parts from receiving
manuscript figures. Incomplete parts are drawn under
`figures/original_format_preview/`, with an explicit preview caption and
`replication_counts.csv`; they are not published as final manuscript figures.
One exception applies to Simulation 3: ICI is not estimable in a test sample
with no subject followed to the horizon (recorded in
`calibration_diagnostics.csv` with `at_risk = 0`). If those are the only
missing values, the final Figure 6 is drawn from the finite repetitions
without on-figure labels. The affected cells and finite counts are written to
`figures/paper.fig6.partial_repetitions.csv` for the manuscript caption
(currently n = 300, 70% censoring, 90% horizon: 52/100).
The computation command still returns a failure status for incomplete results.
Redrawing uses saved tables only and does not rerun simulations or require
loading the package's calibration dependencies.

After all four parts finish, render the report using:

```r
rmarkdown::render(
  "paper.code.Rmd", output_dir = "results/report",
  params = list(simulation_output = "results/simulations/ici_crn")
)
```

This code update does not run the simulations, render the report, or replace
the manuscript figures. Update the manuscript's numerical interpretation and
captions after inspecting the completed results.
