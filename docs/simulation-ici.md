# Run the simulations with ICI

The code includes ICI in Simulation 1, both parts of Simulation 2 and
Simulation 3. These are the existing right-censored designs. The data-generating
mechanisms, seeds, sample sizes, fitted prediction models and evaluation horizons
are unchanged. No competing-risks or sampling-design simulation is added here.

ICI uses `tm_calibration()` with `nk = 5`. Its outcome horizon is the existing
`t_star` for each job. Predicted risks come from the last supplied grid point at
or before that horizon. ICI values are saved without rounding. Other metrics
retain their previous rounding. E50 and E90 are not added to simulation outputs.

## Terminal

Run from the `TimeMetric-paper` directory:

```sh
# Print the job plan only. No model fitting or evaluation occurs.
Rscript --vanilla figures/run_simulations.R --mode=plan

# Full run with 100 repetitions per cell, including ICI.
Rscript --vanilla figures/run_simulations.R --mode=full --output=results/simulations/ici_full
```

On the current development machine, `pmcalibration` is available in
`/tmp/timemetric-calibration-lib`. If that library is still being used, the
full-run command is:

```sh
R_LIBS=/tmp/timemetric-calibration-lib Rscript --vanilla figures/run_simulations.R --mode=full --output=results/simulations/ici_full
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
    "--output=results/simulations/ici_full")
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
Rscript --vanilla figures/redraw_original.R --section=simulations --sim-output=results/simulations/ici_full
```

After all four parts finish, render the report using:

```r
rmarkdown::render(
  "paper.code.Rmd", output_dir = "results/report",
  params = list(simulation_output = "results/simulations/ici_full")
)
```

This code update does not run the simulations, render the report, or replace
the manuscript figures. Update the manuscript's numerical interpretation and
captions after inspecting the completed results.
