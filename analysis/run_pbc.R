#!/usr/bin/env Rscript
arg <- grep("^--file=", commandArgs(), value = TRUE)
script <- sub("^--file=", "", arg[1])
if (!file.exists(script)) script <- gsub("~+~", " ", script, fixed = TRUE)
source(file.path(dirname(normalizePath(script)), "common.R"))
root <- paper_root(script)
out <- file.path(root, "results", "pbc")
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  if (length(args) != 1 || !startsWith(args, "--output=")) stop("Use --output=/path.")
  out <- sub("^--output=", "", args)
}
paper_load(root)
paper_manifest(root, out, list(seed = 1234567L, analysis = "PBC baseline-only comparison"))

set.seed(1234567)
e <- new.env(); utils::data("pbc", package = "survival", envir = e)
df <- e$pbc[!is.na(e$pbc$trt), ]
df$log_albumin <- log(df$albumin); df$log_bili <- log(df$bili)
df$log_protime <- log(df$protime); df$status <- as.integer(df$status == 2)
covs <- c("age", "log_albumin", "log_bili", "log_protime", "edema")
df <- df[, c("time", "status", covs)]
idx <- sample(seq_len(nrow(df)), round(2/3 * nrow(df)))
train <- df[idx, ]; test <- df[-idx, ]
fits <- list(
  Cox = survival::coxph(survival::Surv(time, status) ~ ., train, x = TRUE, y = TRUE),
  Cox.Reduce = survival::coxph(survival::Surv(time, status) ~ age + edema, train, x = TRUE, y = TRUE),
  Weibull.AFT = survival::survreg(survival::Surv(time, status) ~ ., train,
                                 dist = "weibull", x = TRUE, y = TRUE))
old <- paper_legacy_cox(root)
make_predictions <- function(cox_fn) list(
  Cox = cox_fn(fits$Cox, covs, new_data = test),
  Cox.Reduce = cox_fn(fits$Cox.Reduce, c("age", "edema"), new_data = test),
  Weibull.AFT = tm_predict_survreg(fits$Weibull.AFT, covs, new_data = test))
before <- make_predictions(old)
after <- make_predictions(tm_predict_coxph)
tables <- lapply(list(before = before, after = after), function(p) {
  tm_summarize(p, metrics = paper_metrics, tau = 10e10, digits = 4)
})
for (name in names(tables)) {
  utils::write.csv(tables[[name]], file.path(out, paste0("metrics_", name, ".csv")), row.names = FALSE)
  paper_write_md(tables[[name]], file.path(out, paste0("metrics_", name, ".md")))
}
comparison <- do.call(rbind, lapply(names(fits), function(m) {
  data.frame(Model = m, Metric = tables$after$Metric, Before = tables$before[[m]],
             After = tables$after[[m]], Difference = tables$after[[m]] - tables$before[[m]])
}))
utils::write.csv(comparison, file.path(out, "comparison.csv"), row.names = FALSE)
paper_write_md(comparison, file.path(out, "comparison.md"))

# Verify the actual PBC held-out call, including the reduced one-column pathway.
changed <- test
changed$status[which(test$status == 1 & test$time <= median(test$time))[1]] <- 0
checks <- do.call(rbind, lapply(c("Cox", "Cox.Reduce"), function(m) {
  cv <- if (m == "Cox") covs else c("age", "edema")
  perturbed <- tm_predict_coxph(fits[[m]], cv, new_data = changed)
  ref <- survival::survfit(fits[[m]], newdata = test[order(test$time), ], se.fit = FALSE)
  expected <- t(summary(ref, times = sort(test$time), extend = TRUE)$surv)
  data.frame(Model = m, Probability_change_when_labels_changed = max(abs(after[[m]]$surv_prob - perturbed$surv_prob)),
             Mean_change_when_labels_changed = max(abs(after[[m]]$pred - perturbed$pred)),
             Max_difference_from_survfit = max(abs(after[[m]]$surv_prob - expected)))
}))
stopifnot(all(checks$Probability_change_when_labels_changed == 0),
           all(checks$Mean_change_when_labels_changed == 0),
           all(checks$Max_difference_from_survfit < 1e-12),
           all(is.finite(comparison$After)))
utils::write.csv(checks, file.path(out, "prediction_checks.csv"), row.names = FALSE)
saveRDS(list(train = train, test = test, models = fits, before = before, after = after),
         file.path(out, "pbc_comparison.rds"))
settings <- paper_support(fits$Cox, test, 10e10)
settings$probability_time_requested <- median(test$time)
settings$probability_time_used <- sort(test$time)[TimeMetric:::tm_nearest_time_index(sort(test$time), median(test$time))]
utils::write.csv(settings, file.path(out, "analysis_settings.csv"), row.names = FALSE)

plot_data <- do.call(rbind, lapply(c("before", "after"), function(version) {
  ps <- if (version == "before") before else after
  do.call(rbind, lapply(names(ps), function(m) {
    p <- ps[[m]]
    data.frame(Version = version, Model = m, Observed = p$times, Status = p$status, Predicted = p$pred)
  }))
}))
utils::write.csv(plot_data, file.path(out, "plot_data.csv"), row.names = FALSE)
p <- ggplot2::ggplot(plot_data, ggplot2::aes(Observed, Predicted, colour = factor(Status))) +
  ggplot2::geom_point(alpha = 0.6, size = 1.2) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  ggplot2::facet_grid(Version ~ Model) + ggplot2::theme_bw(base_size = 11) +
  ggplot2::labs(x = "Observed follow-up (days)", y = "Model-derived mean-time summary (days)",
                colour = "Event", caption = "Descriptive display; not a calibration curve. Legacy integration grid retained.")
ggplot2::ggsave(file.path(out, "pbc_before_after.png"), p, width = 11, height = 7, dpi = 150)
ggplot2::ggsave(file.path(out, "pbc_before_after.pdf"), p, width = 11, height = 7)
writeLines(c("PBC check passed: 208 training and 104 test participants; seed 1234567.",
  "Before and after both use upstream 0.2.1 evaluation; only the Cox prediction adapter differs.",
  "Legacy time-grid integration and tail carry-forward are retained to isolate the baseline repair.",
  "This full-domain comparison is a numerical diagnostic, not a claim that an unrestricted mean is identified.",
  "Use analysis_settings.csv to inspect the actual integration endpoint and probability horizon."),
  file.path(out, "README.txt"))
print(checks); print(comparison)
cat("PBC PASS. Results:", out, "\n")
