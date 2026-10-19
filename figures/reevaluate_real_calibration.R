#!/usr/bin/env Rscript
# Add cohort ICI to an existing real-data run without refitting its models.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
if (!file.exists(script)) script <- gsub("~+~", " ", script, fixed = TRUE)
root <- normalizePath(file.path(dirname(script), ".."))
source(file.path(root, "analysis", "common.R"))
source(file.path(root, "figures", "original_format.R"))
source(file.path(root, "figures", "paper_outputs.R"))
input <- paper_current_real_run(root)
output <- NULL
for (arg in commandArgs(trailingOnly = TRUE)) {
  if (startsWith(arg, "--input=")) input <- sub("^--input=", "", arg) else
  if (startsWith(arg, "--output=")) output <- sub("^--output=", "", arg) else
    stop("Use --input=/path/to/source_run and --output=/path/to/new_run.")
}
input <- normalizePath(input, mustWork = TRUE)
if (is.null(output)) output <- file.path(root, "results", "real_data", "runs",
  paste0(format(Sys.time(), "%Y-%m-%d_%H%M%S", tz = "America/Los_Angeles"), "_ici"))
if (file.exists(output)) stop("Choose a new output directory; existing runs are preserved.")
paper_load(root)
dir.create(output, recursive = TRUE)
output <- normalizePath(output)
stopifnot(all(file.copy(list.files(input, full.names = TRUE, all.files = TRUE,
  no.. = TRUE), output, recursive = TRUE)))
provenance <- list(source_run = input, output_run = output, models_refitted = FALSE,
  evaluated_at = format(Sys.time(), tz = "America/Los_Angeles", usetz = TRUE),
  source_signature = paper_signature(root, list(calibration = c("ici"))),
  script_md5 = unname(tools::md5sum(script)))

write_metrics <- function(table, stem) {
  utils::write.csv(table, paste0(stem, ".csv"), row.names = FALSE)
  paper_write_md(table, paste0(stem, ".md"))
  tex <- c(paste0("\\begin{tabular}{l", strrep("r", ncol(table)-1), "}"),
    "\\toprule", paste0(paste(gsub("_", "\\_", names(table), fixed = TRUE),
      collapse = " & "), " \\\\"), "\\midrule")
  for (i in seq_len(nrow(table))) tex <- c(tex, paste0(paste(c(
    gsub("_", "\\_", table$Metric[i], fixed = TRUE),
    formatC(as.numeric(table[i, -1]), format = "f", digits = 4)),
    collapse = " & "), " \\\\"))
  writeLines(c(tex, "\\bottomrule", "\\end{tabular}"), paste0(stem, ".tex"))
}
unchanged_old_metrics <- function(updated, previous) {
  actual <- updated[match(previous$Metric, updated$Metric), names(previous), drop = FALSE]
  rownames(actual) <- rownames(previous) <- NULL
  stopifnot(identical(actual, previous))
}
replace_ici <- function(previous, ici) {
  # A published run may already contain ICI when this script is run again.
  retained <- previous[previous$Metric != "ici", , drop = FALSE]
  table <- rbind(retained, ici[, names(previous), drop = FALSE])
  unchanged_old_metrics(table, retained)
  stopifnot(sum(table$Metric == "ici") == 1L)
  table
}
diagnostics <- list()
record <- function(cal, dataset, model, prediction_time) {
  error <- abs(cal$curve$predicted - cal$curve$calibrated)
  stopifnot(isTRUE(all.equal(unname(cal$metrics["ici"]), mean(error))),
            nrow(cal$curve) == length(error), all(is.finite(error)))
  diagnostics[[length(diagnostics)+1L]] <<- data.frame(
    Dataset = dataset, Model = model, Evaluation_time_days = cal$t_star,
    Prediction_grid_time_days = prediction_time, N = length(error),
    Method = cal$method, ICI = mean(error),
    Boundary_fitting_probabilities = sum(cal$curve$predicted <= sqrt(.Machine$double.eps) |
      cal$curve$predicted >= 1-sqrt(.Machine$double.eps)))
}

# PBC keeps its original default horizon: median test follow-up (1829 days),
# with predictions carried from the last knot at or before it (1827 days).
pbc_file <- file.path(input, "pbc", "pbc_comparison.rds")
pbc <- readRDS(pbc_file)
pbc_original <- pbc
calibration <- list()
for (version in c("before", "after")) {
  previous <- utils::read.csv(file.path(input, "pbc", paste0("metrics_", version, ".csv")),
                             check.names = FALSE)
  ici <- tm_summarize(pbc[[version]], metrics = "ici", tau = 10e10, digits = 4)
  table <- replace_ici(previous, ici)
  write_metrics(table, file.path(output, "pbc", paste0("metrics_", version)))
  if (version == "after") {
    write_metrics(table, file.path(output, "tables", "pbc_metrics"))
    for (name in names(pbc$after)) {
      p <- pbc$after[[name]]
      t <- median(p$times)
      k <- findInterval(t, sort(p$times))
      cal <- tm_calibration(1-p$surv_prob[, k], p$times, p$status, t)
      stopifnot(table[[name]][table$Metric == "ici"] == round(cal$metrics["ici"], 4))
      calibration[[name]] <- cal
      record(cal, "PBC", name, sort(p$times)[k])
    }
  }
}
pbc$calibration <- calibration
pbc$calibration_reevaluation <- c(provenance, list(source_predictions = pbc_file))
stopifnot(identical(pbc$models, pbc_original$models),
          identical(pbc$after, pbc_original$after), identical(pbc$before, pbc_original$before))
saveRDS(pbc, file.path(output, "pbc", "pbc_comparison.rds"))
pbc_settings <- utils::read.csv(file.path(output, "pbc", "analysis_settings.csv"))
pbc_settings$calibration_time <- diagnostics[[1]]$Evaluation_time_days
pbc_settings$calibration_prediction_time <- diagnostics[[1]]$Prediction_grid_time_days
pbc_settings$calibration_method <- "pmcalibration restricted cubic spline; five knots"
utils::write.csv(pbc_settings, file.path(output, "pbc", "analysis_settings.csv"), row.names = FALSE)

# Update both saved FHS horizons, retaining the selected 12-year manuscript view.
files <- list.files(file.path(input, "framingham"),
  pattern = "^framingham_visit3_.*y\\.rds$", full.names = TRUE, recursive = TRUE)
if (!length(files)) stop("No saved Framingham predictions found.")
for (file in files) {
  original <- readRDS(file)
  result <- original
  tau <- result$settings$tau
  ici <- tm_summarize_cr(result$predictions,
    metrics = "ici", event_type = 1,
    tau = tau, t_star = tau, digits = 4, censoring = "cox")
  table <- replace_ici(original$metrics, ici)
  calibration <- list()
  for (name in names(result$predictions)) {
    p <- result$predictions[[name]]
    k <- findInterval(tau, p$cif_pred[, 1])
    cal <- tm_calibration_cr(as.numeric(p$cif_pred[k, -1]), p$times, p$status,
                              t_star = tau, event_type = 1)
    stopifnot(table[[name]][table$Metric == "ici"] == round(cal$metrics["ici"], 4))
    calibration[[name]] <- cal
    record(cal, paste0("FHS_", tau/365, "y"), name, p$cif_pred[k, 1])
  }
  result$metrics <- table
  result$calibration <- calibration
  result$calibration_reevaluation <- c(provenance, list(source_predictions = file))
  stopifnot(identical(result$models, original$models),
            identical(result$predictions, original$predictions),
            identical(result$cohort, original$cohort))
  target <- file.path(output, substring(file, nchar(input)+2L))
  saveRDS(result, target)
  write_metrics(table, file.path(dirname(target), "metrics"))
}
active_file <- file.path(output, "framingham", "active_horizon.csv")
active <- utils::read.csv(active_file, stringsAsFactors = FALSE)
stopifnot(nrow(active) == 1L, active$Horizon_years == 12,
          startsWith(active$Source_predictions, paste0(input, "/")))
active$Source_predictions <- file.path(output,
  substring(active$Source_predictions, nchar(input)+2L))
active$Fit_source_md5 <- unname(tools::md5sum(file.path(output, "framingham",
  "framingham_visit3_10y.rds")))
utils::write.csv(active, active_file, row.names = FALSE)
write_metrics(readRDS(active$Source_predictions)$metrics,
              file.path(output, "tables", "framingham_metrics"))

# The entire two-phase analysis and every source figure must remain identical.
retained <- c(list.files(file.path(input, "gbsg"), recursive = TRUE, full.names = TRUE),
  list.files(input, pattern = "paper\\.fig[123].*", recursive = TRUE, full.names = TRUE),
  file.path(input, "tables", paste0("gbsg_metrics.", c("csv", "md", "tex"))))
copied <- file.path(output, substring(retained, nchar(input)+2L))
stopifnot(identical(unname(tools::md5sum(retained)), unname(tools::md5sum(copied))))
diagnostics <- do.call(rbind, diagnostics)
utils::write.csv(diagnostics, file.path(output, "calibration_diagnostics.csv"), row.names = FALSE)
saveRDS(provenance, file.path(output, "calibration_reevaluation.rds"))
capture.output(sessionInfo(), file = file.path(output, "calibration_sessionInfo.txt"))
writeLines(c("# Cohort calibration update", "", paste("Source run:", input),
  "Saved models, splits, predictions and figures are unchanged; no model was refitted.",
  "ICI was added for PBC (1829 days) and FHS (10 and 12 years; manuscript uses 12).",
  "Previously reported metrics are retained exactly; this update recomputes ICI only.",
  "PBC: pmcalibration restricted cubic spline, five knots, no confidence interval computation.",
  "FHS: prodlim nearest-neighborhood Aalen-Johansen calibration, automatic bandwidth.",
  "Every ICI was checked against mean absolute error on the subject-level calibrated probabilities.",
  "Boundary clipping for fitting retains every subject and the original probabilities in metrics.",
  "GBSG sample-design inputs, metrics, saved results, and all figures are byte-identical.",
  "Historical configuration files describe original fits; calibration_reevaluation.rds records this update."),
  file.path(output, "CALIBRATION_UPDATE.md"))
paper_publish_real_outputs(root, output)
print(diagnostics, row.names = FALSE)
cat("Cohort ICI update PASS. Output:", output, "\n")
