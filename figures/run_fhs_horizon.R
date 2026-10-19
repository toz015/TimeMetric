#!/usr/bin/env Rscript
# Re-evaluate an existing third-examination fit at a new horizon. No refitting,
# resampling or simulation-source changes are needed to change tau/t_star.
f <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
if (!file.exists(f)) f <- gsub("~+~", " ", f, fixed = TRUE)
root <- normalizePath(file.path(dirname(f), ".."))
source(file.path(root, "analysis", "common.R"))
source(file.path(root, "analysis", "framingham.R"))
source(file.path(root, "figures", "original_format.R"))
years <- 12
input <- NULL
out <- NULL
for (a in commandArgs(trailingOnly = TRUE)) {
  if (startsWith(a, "--years=")) years <- as.numeric(sub("^--years=", "", a)) else
  if (startsWith(a, "--input=")) input <- sub("^--input=", "", a) else
  if (startsWith(a, "--output=")) out <- sub("^--output=", "", a) else
    stop("Use --years=12, --input=/path/to/fit.rds, --output=/path.")
}
if (is.null(input)) input <- file.path(paper_current_real_run(root), "framingham",
                                      "framingham_visit3_10y.rds")
stopifnot(length(years) == 1L, is.finite(years), years > 0, file.exists(input))
tag <- paste0(format(years, trim = TRUE, scientific = FALSE), "y")
explicit_output <- !is.null(out)
if (!explicit_output) out <- file.path(dirname(input), paste0("horizon_", tag))
paper_load(root)
old <- readRDS(input)
stopifnot(old$settings$examination == 3L, old$settings$days_per_year == 365)
tau <- years * 365
cohort <- old$cohort
e <- new.env()
utils::data("framingham", package = "riskCommunicator", envir = e)
prepared <- paper_framingham_cohort(e$framingham)
stopifnot(identical(cohort$RANDID, prepared$cohort$RANDID),
          identical(cohort$time, cohort$endpoint_time - cohort$landmark_time))
covs <- prepared$covs
reduced_covs <- prepared$reduced_covs
test <- cohort[cohort$Split == "Test", c("time", "status", covs)]
stopifnot(identical(test$time, old$predictions$FG$times),
          identical(test$status, old$predictions$FG$status))

# Inspect empirical support before calculating metrics.
support <- do.call(rbind, lapply(unique(c(10, 11, 12, 12.5, years)), function(y) {
  do.call(rbind, lapply(c("Train", "Test"), function(s) {
    d <- cohort[cohort$Split == s, ]
    km <- survival::survfit(survival::Surv(time, status == 0) ~ 1, data = d)
    data.frame(Split = s, N = nrow(d), Horizon_years = y, tau_days = y * 365,
      At_risk = sum(d$time >= y * 365), Maximum_followup_days = max(d$time),
      G_horizon = if (y * 365 <= max(d$time))
        summary(km, times = y * 365, extend = TRUE)$surv else NA_real_)
  }))
}))
target_support <- support[support$Horizon_years == years, ]
stopifnot(all(target_support$At_risk > 0), all(is.finite(target_support$G_horizon)),
          all(target_support$G_horizon > 0))
settings <- old$settings
settings$tau <- settings$t_star <- tau
settings$horizon_years <- years
settings$cr_r2_weighting <- "Observed events; original follow-up; Cox without covariates; G floor 1e-4"
settings$fit_source <- normalizePath(input)
settings$fit_source_md5 <- unname(tools::md5sum(input))
settings$horizon_runner_md5 <- unname(tools::md5sum(f))
# A refit can serialize formula environments differently despite an identical
# training fit. Fingerprint the numerical fit and cohort for repeatable reruns;
# retain the full input-file checksum separately as provenance.
fit_file <- tempfile(fileext = ".rds")
saveRDS(list(cohort = cohort, models = lapply(old$models, function(model) {
  fields <- intersect(c("coef", "coefficients", "var", "x", "y", "scale",
                        "uftime", "bfitj"), names(model))
  model[fields]
})), fit_file, version = 3)
settings$fit_numeric_md5 <- unname(tools::md5sum(fit_file))
unlink(fit_file)
settings$refitted <- FALSE
settings$tail_support <- if (min(target_support$At_risk) < 20 ||
                            min(target_support$G_horizon) < 0.1)
  "Sparse: exploratory horizon; inspect followup_support.csv" else "See followup_support.csv"
manifest_settings <- settings
manifest_settings$fit_source_md5 <- NULL
# Preserve the recorded result if this runner or its computation sources changed.
# An explicitly selected output keeps the strict source/settings guard.
config_file <- file.path(out, "configuration.rds")
if (!explicit_output && file.exists(config_file) &&
    !identical(readRDS(config_file), paper_signature(root, manifest_settings))) {
  base <- paste0(out, "_", format(Sys.time(), "%Y-%m-%d_%H%M%S", tz = "America/Los_Angeles"))
  out <- base
  index <- 1L
  while (dir.exists(out)) {
    out <- paste0(base, "_", index)
    index <- index + 1L
  }
  message("Source/settings changed; preserving the existing horizon result and writing to ", out)
}
paper_manifest(root, out, manifest_settings)
utils::write.csv(support, file.path(out, "followup_support.csv"), row.names = FALSE)
m <- old$models
predictions <- list(
  FG = tm_predict_cif(fg_model = m$FG, newdata = test, covs = covs, tau = tau),
  FG.Reduced = tm_predict_cif(fg_model = m$FG.Reduced, newdata = test,
                             covs = reduced_covs, tau = tau),
  CSH.AFT = tm_predict_cif(model1 = m$AFT.Hypertension, model2 = m$AFT.Death,
                          newdata = test, covs = covs, tau = tau),
  CSH.AFT.FIX = tm_predict_cif(model1 = m$AFT.Fixed.Hypertension,
    model2 = m$AFT.Fixed.Death, newdata = test, covs = covs, tau = tau))
checks <- do.call(rbind, lapply(names(predictions), function(n) {
  p <- predictions[[n]]
  stopifnot(identical(p$times, test$time), identical(p$status, test$status),
    identical(p$cif_pred, old$predictions[[n]]$cif_pred),
    identical(p$linear.pred, old$predictions[[n]]$linear.pred),
    all(is.finite(p$pred)), all(p$pred > 0 & p$pred <= tau + 1e-8))
  data.frame(Model = n, Same_CIF_and_risk_scores = TRUE,
    tau_days = tau, CIF_knot_at_horizon = max(p$cif_pred[p$cif_pred[, 1] <= tau, 1]),
    Grid_maximum_days = max(p$cif_pred[, 1]))
}))
metrics <- tm_summarize_cr(predictions,
  metrics = c("r_square", "l_square", "pseudo_r2", "pseudo_r2_point",
              "c_index", "brier_score", "td_auc"),
  event_type = 1, tau = tau, t_star = tau, digits = 4, censoring = "cox")
stopifnot(all(is.finite(as.matrix(metrics[, -1]))))
utils::write.csv(metrics, file.path(out, "metrics.csv"), row.names = FALSE)
paper_write_md(metrics, file.path(out, "metrics.md"))
utils::write.csv(checks, file.path(out, "prediction_checks.csv"), row.names = FALSE)
utils::write.csv(cohort, file.path(out, "cohort_and_split.csv"), row.names = FALSE)
ipcw <- getFromNamespace("tm_cr_event_ipcw", "TimeMetric")(
  test$time, test$status, censoring = "cox")
r2_weights <- data.frame(RANDID = cohort$RANDID[cohort$Split == "Test"],
  time = test$time, status = test$status,
  censoring_survival_at_original_time = ipcw$probability,
  normalized_event_ipcw = ipcw$weight)
stopifnot(all(r2_weights$normalized_event_ipcw[test$status == 0] == 0),
          abs(sum(r2_weights$normalized_event_ipcw) - 1) < 1e-12)
utils::write.csv(r2_weights, file.path(out, "r2_original_followup_weights.csv"), row.names = FALSE)
point_predictions <- data.frame(RANDID = cohort$RANDID[cohort$Split == "Test"],
                               time = test$time, status = test$status)
for (n in names(predictions)) {
  p <- predictions[[n]]
  k <- max(which(p$cif_pred[, 1] <= tau))
  point_predictions[[paste0(n, "_CIF_at_horizon")]] <- p$cif_pred[k, -1]
  point_predictions[[paste0(n, "_restricted_mean_days")]] <- p$pred
}
utils::write.csv(point_predictions, file.path(out, "test_predictions.csv"), row.names = FALSE)
plot <- paper_original_real_plot("framingham", predictions, tau)
paper_original_save(plot, out, "paper.fig2", 9, 7)
plot_data <- do.call(rbind, lapply(seq_along(predictions), function(i) {
  cbind(Model = names(predictions)[i], plot[[i]]$data)
}))
stopifnot(all(plot_data$times <= tau))
utils::write.csv(plot_data, file.path(out, "paper.fig2.plot_data.csv"), row.names = FALSE)
result <- list(settings = settings, cohort = cohort, models = m, predictions = predictions,
               metrics = metrics, support = support, r2_weights = r2_weights,
               prediction_checks = checks, plot = plot)
saveRDS(result, file.path(out, paste0("framingham_visit3_", tag, ".rds")))
# Twelve years is the selected manuscript horizon. Keep its main figure and
# consolidated table in sync while retaining the earlier horizon directories.
if (years == 12) {
  parent <- dirname(input)
  for (ext in c("pdf", "png")) {
    stopifnot(file.copy(file.path(out, paste0("paper.fig2.", ext)),
      file.path(parent, paste0("paper.fig2.", ext)), overwrite = TRUE))
  }
  saveRDS(plot, file.path(parent, "paper.fig2.plot.rds"))
  utils::write.csv(plot_data, file.path(parent, "paper.fig2.plot_data.csv"), row.names = FALSE)
  table_dir <- file.path(dirname(parent), "tables")
  if (basename(parent) == "framingham" && dir.exists(table_dir)) {
    utils::write.csv(metrics, file.path(table_dir, "framingham_metrics.csv"), row.names = FALSE)
    paper_write_md(metrics, file.path(table_dir, "framingham_metrics.md"))
    header <- gsub("_", "\\_", names(metrics), fixed = TRUE)
    tex <- c(paste0("\\begin{tabular}{l", paste(rep("r", ncol(metrics)-1), collapse = ""), "}"),
      "\\toprule", paste0(paste(header, collapse = " & "), " \\\\"), "\\midrule")
    for (i in seq_len(nrow(metrics))) {
      row <- c(gsub("_", "\\_", metrics$Metric[i], fixed = TRUE),
        formatC(as.numeric(metrics[i, -1]), format = "f", digits = 4))
      tex <- c(tex, paste0(paste(row, collapse = " & "), " \\\\"))
    }
    writeLines(c(tex, "\\bottomrule", "\\end{tabular}"),
      file.path(table_dir, "framingham_metrics.tex"))
  }
  utils::write.csv(data.frame(Horizon_years = years, tau_days = tau, t_star_days = tau,
    Examination = 3L, Source_predictions = normalizePath(file.path(out,
      paste0("framingham_visit3_", tag, ".rds"))),
    Fit_source_md5 = settings$fit_source_md5), file.path(parent, "active_horizon.csv"),
    row.names = FALSE)
}
writeLines(c(
  paste("Horizon:", years, "years =", tau, "days. Seed and split unchanged; no model refit."),
  "Time zero: each participant's actual examination-3 date (PERIOD == 3).",
  "Observed follow-up = min(TIMEHYP, TIMEDTH) from baseline minus examination-3 TIME.",
  "Native package plot: points show follow-up capped at tau; line shows integral of 1 - hypertension CIF through tau.",
  "Competing deaths remain a distinct outcome; this line is not an unrestricted event-time mean.",
  "The original model-specific grids, step interpolation and quadrature are retained.",
  "R2 uses original-follow-up event-observation IPCW, a Cox censoring model without covariates and G floor 1e-4.",
  "Censored rows have zero R2 weight even beyond the horizon; weights can depend on follow-up after tau.",
  "R2 original-time probabilities and normalized weights are in r2_original_followup_weights.csv.",
  "Fine-Gray CIF is constant after its last training-event knot; see prediction_checks.csv.",
  settings$tail_support,
  paste("Training/test numbers at risk:", paste(target_support$At_risk, collapse = "/")),
  "The original 10-year results are preserved. This directory contains only the requested new horizon."),
  file.path(out, "README.txt"))
print(target_support, row.names = FALSE)
print(metrics, row.names = FALSE)
cat("FHS", tag, "horizon calculation PASS. Support:", settings$tail_support, "\n")
cat("Outputs:", normalizePath(out), "\n")
