#!/usr/bin/env Rscript
# Reevaluate saved CR predictions after the reciprocal-IPCW correction.
# Keep the original run intact and publish a separate recorded evaluation.
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
  paste0(format(Sys.time(), "%Y-%m-%d_%H%M%S", tz = "America/Los_Angeles"), "_c_index_ipcw"))
if (file.exists(output)) stop("Choose a new output directory; existing runs are preserved.")
dir.create(output, recursive = TRUE)
output <- normalizePath(output)
stopifnot(all(file.copy(list.files(input, full.names = TRUE, all.files = TRUE,
  no.. = TRUE), output, recursive = TRUE)))
paper_load(root)

# Independent comparable-pair summation, separate from the matrix estimator.
reference <- function(time, status, predicted, tau) {
  fit <- survival::survfit(survival::Surv(time, status == 0) ~ 1)
  left <- vapply(time, function(t) min(c(1, fit$surv[fit$time < t])), 0.)
  right <- vapply(time, function(t) min(c(1, fit$surv[fit$time <= t])), 0.)
  numerator <- denominator <- 0
  for (i in which(status == 1 & time <= tau)) {
    for (j in seq_along(time)) {
      probability <- if (time[i] < time[j]) left[i] * right[i] else
        if (status[j] == 2) left[i] * left[j] else NA_real_
      if (is.na(probability)) next
      stopifnot(is.finite(probability), probability > 0)
      denominator <- denominator + 1 / probability
      numerator <- numerator + (predicted[i] < predicted[j]) / probability
    }
  }
  numerator / denominator
}
write_metrics <- function(table, stem) {
  utils::write.csv(table, paste0(stem, ".csv"), row.names = FALSE)
  paper_write_md(table, paste0(stem, ".md"))
  header <- gsub("_", "\\_", names(table), fixed = TRUE)
  tex <- c(paste0("\\begin{tabular}{l", strrep("r", ncol(table)-1), "}"),
    "\\toprule", paste0(paste(header, collapse = " & "), " \\\\"), "\\midrule")
  for (i in seq_len(nrow(table))) {
    row <- c(gsub("_", "\\_", table$Metric[i], fixed = TRUE),
      formatC(as.numeric(table[i, -1]), format = "f", digits = 4))
    tex <- c(tex, paste0(paste(row, collapse = " & "), " \\\\"))
  }
  writeLines(c(tex, "\\bottomrule", "\\end{tabular}"), paste0(stem, ".tex"))
}
provenance <- list(method = "Wolbers et al. (2014), equation (3.4): reciprocal IPCW",
  source_run = input, output_run = output,
  evaluated_at = format(Sys.time(), tz = "America/Los_Angeles", usetz = TRUE),
  models_refitted = FALSE, source_signature = paper_signature(root, list()),
  reevaluation_script_md5 = unname(tools::md5sum(script)))
files <- list.files(file.path(input, "framingham"),
  pattern = "^framingham_visit3_.*y\\.rds$", full.names = TRUE, recursive = TRUE)
if (!length(files)) stop("No saved Framingham horizon predictions found.")
changes <- list()
for (file in files) {
  original <- readRDS(file)
  result <- original
  tau <- result$settings$tau
  updated <- tm_summarize_cr(result$predictions, metrics = "c_index",
    event_type = 1, tau = tau, t_star = tau, digits = 4)
  row <- match("c_index", result$metrics$Metric)
  stopifnot(!is.na(row))
  for (model in names(result$predictions)) {
    p <- result$predictions[[model]]
    prediction <- apply(p$cif_pred[, -1, drop = FALSE], 2,
      getFromNamespace("m_cif", "TimeMetric"), time.cif = p$cif_pred[, 1], tau = tau)
    independent <- round(reference(p$times, p$status, prediction, tau), 4)
    stopifnot(updated[[model]] == independent)
    changes[[length(changes) + 1L]] <- data.frame(Horizon_days = tau, Model = model,
      Previous = result$metrics[[model]][row], Corrected = independent)
    result$metrics[[model]][row] <- independent
  }
  stopifnot(identical(result$predictions, original$predictions),
    identical(result$models, original$models),
    identical(result$metrics[-row, ], original$metrics[-row, ]))
  result$concordance_reevaluation <- c(provenance, list(
    source_result = file, source_result_md5 = unname(tools::md5sum(file))))
  relative <- substring(file, nchar(input) + 2L)
  target <- file.path(output, relative)
  saveRDS(result, target)
  write_metrics(result$metrics, file.path(dirname(target), "metrics"))
}
active_file <- file.path(output, "framingham", "active_horizon.csv")
active <- utils::read.csv(active_file, stringsAsFactors = FALSE)
stopifnot(nrow(active) == 1L, active$Horizon_years == 12,
          startsWith(active$Source_predictions, paste0(input, "/")))
active$Source_predictions <- file.path(output,
  substring(active$Source_predictions, nchar(input) + 2L))
utils::write.csv(active, active_file, row.names = FALSE)
selected <- readRDS(active$Source_predictions)
write_metrics(selected$metrics, file.path(output, "tables", "framingham_metrics"))
changes <- do.call(rbind, changes)
utils::write.csv(changes, file.path(output, "c_index_changes.csv"), row.names = FALSE)
saveRDS(provenance, file.path(output, "c_index_reevaluation.rds"))
capture.output(sessionInfo(), file = file.path(output, "c_index_sessionInfo.txt"))
writeLines(c("# C-index reevaluation", "",
  paste("Source run:", input),
  "Saved models, predictions, figures and non-C-index metrics are retained.",
  "Framingham C-index values at all saved horizons were reevaluated with reciprocal IPCW.",
  "Each value was checked against an independent comparable-pair summation.",
  "Existing configuration.rds files describe the original fits/calculations.",
  "c_index_reevaluation.rds records the current evaluation sources and provenance.",
  "No model was refitted; PBC, GBSG and simulations were not recomputed."),
  file.path(output, "C_INDEX_REEVALUATION.md"))
paper_publish_real_outputs(root, output)
print(changes, row.names = FALSE)
cat("Updated manuscript outputs:", file.path(root, "results", "real_data", "paper"), "\n")
