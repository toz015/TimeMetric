#!/usr/bin/env Rscript
# Reevaluate ICI only, with the saved prediction model and original test seeds.
f <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
if (!file.exists(f)) f <- gsub("~+~", " ", f, fixed = TRUE)
root <- normalizePath(file.path(dirname(f), ".."))
source(file.path(root, "analysis", "common.R"))
source(file.path(root, "analysis", "simulation_calibration.R"))
input <- file.path(root, "results", "simulations", "ici_full")
output <- file.path(root, "results", "simulations", "ici_time_alignment_check")
for (arg in commandArgs(trailingOnly = TRUE)) {
  if (startsWith(arg, "--input=")) input <- sub("^--input=", "", arg) else
    if (startsWith(arg, "--output=")) output <- sub("^--output=", "", arg) else
      stop("Use --input=old-run and --output=comparison-directory.")
}
input <- normalizePath(input, mustWork = TRUE)
output <- normalizePath(output, mustWork = FALSE)
if (identical(input, output)) stop("Comparison output must differ from the original run.")
paper_load(root)
raw_path <- file.path(input, "tables", "metrics_raw.csv")
model_path <- file.path(input, "models", "sim1_10.rds")
original_hashes <- tools::md5sum(c(raw_path, model_path))
old <- utils::read.csv(raw_path)
old <- old[old$part == "sim3" & old$censor == 0 & old$Metric == "ici", ]
if (!setequal(old$n, c(300, 3000)) || !setequal(old$quantile, c(50, 90)))
  stop("The comparison requires both full-run sample sizes and horizons.")
cells <- split(old, interaction(old$n, old$quantile, drop = TRUE))
if (length(cells) != 4L || !all(vapply(cells, function(d)
    nrow(d) == 100L && setequal(d$rep, 1:100) && all(is.finite(d$Value)), logical(1))))
  stop("The original run must contain 100 finite ICI values per uncensored cell.")
model <- readRDS(model_path)$Cox
if (!inherits(model, "coxph")) stop("The original training Cox model is missing.")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
rows <- vector("list", nrow(old))
for (i in seq_len(nrow(old))) {
  r <- old[i, ]
  d <- tm_sim_cox_weibull(r$n, pi_c = 0, v = r$v,
                         beta = r$beta*c(1, -1), seed = r$seed)
  d <- d[, c("time", "status", "x1", "x2")]
  ordered <- d[order(d$time), ]
  # Only outcomes are needed here: exact-time risk is requested from the model.
  pred <- list(times = ordered$time, status = ordered$status)
  warnings <- character()
  result <- withCallingHandlers(
    paper_simulation_calibration(pred, r$t_star, list(nk = 5L),
                                   model = model, new_data = d),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  if (!is.finite(result$ici)) stop("Corrected ICI failed for ", r$id,
                                  ": ", result$diagnostics$error)
  rows[[i]] <- data.frame(id = r$id, n = r$n, quantile = r$quantile,
    rep = r$rep, seed = r$seed, t_star = r$t_star,
    old_probability_time = max(d$time[d$time <= r$t_star]),
    new_probability_time = result$diagnostics$prediction_grid_time,
    old_ici = r$Value, new_ici = result$ici,
    ici_change = result$ici-r$Value,
    warnings = paste(unique(warnings), collapse = " | "))
  if (i %% 100L == 0L) message("Compared ", i, "/", nrow(old), " saved evaluations.")
}
paired <- do.call(rbind, rows)
means <- stats::aggregate(cbind(old_ici, new_ici, ici_change) ~ n + quantile,
                           paired, mean)
comparison <- do.call(rbind, lapply(sort(unique(paired$quantile)), function(q) {
  small <- paired[paired$n == 300 & paired$quantile == q, ]
  large <- paired[paired$n == 3000 & paired$quantile == q, ]
  old_reference <- mean(large$old_ici)
  new_reference <- mean(large$new_ici)
  old_gap <- mean(small$old_ici) - old_reference
  new_gap <- mean(small$new_ici) - new_reference
  data.frame(quantile = q, repetitions_per_cell = 100L,
    old_mean_n300 = mean(small$old_ici), new_mean_n300 = mean(small$new_ici),
    old_reference_n3000 = old_reference, new_reference_n3000 = new_reference,
    old_gap = old_gap, new_gap = new_gap, gap_change = new_gap-old_gap,
    gap_reduction_percent = 100*(1-abs(new_gap)/abs(old_gap)),
    paired_change_mcse = sqrt(stats::var(small$ici_change)/nrow(small) +
                               stats::var(large$ici_change)/nrow(large)),
    new_gap_using_old_reference = mean(small$new_ici)-old_reference)
}))
utils::write.csv(paired, file.path(output, "paired_ici.csv"), row.names = FALSE)
utils::write.csv(means, file.path(output, "cell_means.csv"), row.names = FALSE)
utils::write.csv(comparison, file.path(output, "gap_comparison.csv"), row.names = FALSE)
saveRDS(list(input = input, original_file_md5 = original_hashes,
  changed_source_md5 = tools::md5sum(file.path(root, "analysis",
    c("simulation_calibration.R", "check_sim3_ici_time_alignment.R"))),
  nk = 5L, censoring = 0, refit_prediction_model = FALSE,
  reference = "Each version's mean across its 100 uncensored n=3000 evaluations"),
  file.path(output, "comparison_settings.rds"))
stopifnot(identical(original_hashes, tools::md5sum(names(original_hashes))),
          all(paired$new_probability_time == paired$t_star))
print(comparison, row.names = FALSE)
cat("Saved ICI time-alignment comparison to:", output, "\n")
