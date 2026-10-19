#!/usr/bin/env Rscript
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
script <- sub("^--file=", "", file_arg[1])
if (!file.exists(script)) script <- gsub("~+~", " ", script, fixed = TRUE)
source(file.path(dirname(normalizePath(script)), "common.R"))
root <- paper_root(script)
out <- file.path(root, "results", "gbsg2")
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  if (length(args) != 1L || !startsWith(args, "--output=")) stop("Use --output=/path.")
  out <- sub("^--output=", "", args)
  if (!nzchar(out)) stop("The output path cannot be empty.")
}
paper_load(root)
paper_manifest(root, out, list(analysis = "GBSG2 case-cohort", seed = 2025L,
  subcohort_fraction = 0.10, evaluation_weights = "Aligned to sorted prediction subjects"))

set.seed(2025)
e <- new.env(); utils::data("GBSG2", package = "pec", envir = e)
df <- e$GBSG2
df$status <- df$cens; df$id <- seq_len(nrow(df))
covs <- names(df)[match("horTh", names(df)):match("estrec", names(df))]
df <- stats::na.omit(df[df$time != 0, c("time", "status", "id", covs)])
idx <- sample(seq_len(nrow(df)), round(2/3 * nrow(df)))
train <- df[idx, ]; test <- df[-idx, ]
sub_train <- sample(train$id, round(nrow(train) * 0.10))
sub_test <- sample(test$id, round(nrow(test) * 0.10))
km <- survival::survfit(survival::Surv(time, 1 - status) ~ 1, data = train)

# Obtain inclusion weights on the complete phase-one data, before selecting rows.
train$subco <- train$id %in% sub_train
train$weight <- tm_case_cohort_weights(train$time, train$status, train$subco)
train$selected <- train$subco | train$status == 1
cc_train <- train[train$selected, ]
random_train <- train[train$id %in% sample(train$id, nrow(cc_train)), ]
test$subco <- test$id %in% sub_test
test$weight <- tm_case_cohort_weights(test$time, test$status, test$subco)
test$selected <- test$subco | test$status == 1
cc_test <- test[test$selected, ]
# Cox adapters sort subjects by time. Sort IDs, outcomes, covariates and weights
# together before prediction so every evaluation weight refers to its subject.
cc_test <- cc_test[order(cc_test$time), ]
weights <- cc_test$weight
model_cols <- c("time", "status", covs)
fit1 <- survival::coxph(survival::Surv(time, status) ~ .,
                        data = cc_train[, model_cols], weights = cc_train$weight,
                        x = TRUE, y = TRUE)
fit2 <- survival::coxph(survival::Surv(time, status) ~ .,
                        data = random_train[, model_cols], x = TRUE, y = TRUE)
eval <- cc_test[, model_cols]
models <- list(Cox.cc = fit1, Cox.random = fit2)
predictions <- lapply(models, tm_predict_coxph, covs = covs, new_data = eval)
for (p in predictions) {
  stopifnot(identical(p$times, eval$time), identical(as.numeric(p$status), as.numeric(eval$status)),
            all(is.finite(p$surv_prob)), all(is.finite(p$pred)))
}
stopifnot(all(is.finite(weights)), all(weights > 0))
metrics <- tm_sample_design(predictions, case_weights = weights, km_cens = km,
                            metrics = paper_metrics, tau = 10e10, digits = 4)
stopifnot(all(is.finite(as.matrix(metrics[, -1]))))
utils::write.csv(metrics, file.path(out, "metrics.csv"), row.names = FALSE)
paper_write_md(metrics, file.path(out, "metrics.md"))
phase_one <- rbind(cbind(Split = "Train", train), cbind(Split = "Test", test))
utils::write.csv(phase_one, file.path(out, "phase_one_sampling.csv"), row.names = FALSE)
utils::write.csv(cc_test[, c("id", "time", "status", "subco", "weight")],
                 file.path(out, "evaluation_subjects_and_weights.csv"), row.names = FALSE)
checks <- do.call(rbind, lapply(names(models), function(m) {
  changed <- eval; changed$status[1] <- 1 - changed$status[1]
  a <- predictions[[m]]
  b <- tm_predict_coxph(models[[m]], covs, new_data = changed)
  data.frame(Model = m,
    Probability_change_when_labels_changed = max(abs(a$surv_prob - b$surv_prob)),
    Mean_change_when_labels_changed = max(abs(a$pred - b$pred)),
    Weights_aligned = identical(a$times, cc_test$time))
}))
stopifnot(all(checks$Probability_change_when_labels_changed == 0),
           all(checks$Mean_change_when_labels_changed == 0), all(checks$Weights_aligned))
utils::write.csv(checks, file.path(out, "prediction_checks.csv"), row.names = FALSE)
settings <- data.frame(Full_cohort_n = nrow(df), Phase_one_train_n = nrow(train),
  Phase_one_test_n = nrow(test), Case_cohort_train_n = nrow(cc_train),
  Random_train_n = nrow(random_train), Case_cohort_test_n = nrow(cc_test),
  Train_subcohort_n = length(sub_train), Test_subcohort_n = length(sub_test),
  Probability_time_requested = median(eval$time),
  Probability_time_used = max(eval$time[eval$time <= median(eval$time)]),
  Time_horizon = 10e10, Actual_integration_end = max(eval$time))
utils::write.csv(settings, file.path(out, "analysis_settings.csv"), row.names = FALSE)
support <- do.call(rbind, lapply(names(models), function(m) {
  cbind(Model = m, paper_support(models[[m]], eval, median(eval$time)))
}))
utils::write.csv(support, file.path(out, "support_diagnostics.csv"), row.names = FALSE)
plot_data <- do.call(rbind, lapply(names(predictions), function(m) {
  p <- predictions[[m]]
  data.frame(Model = m, id = cc_test$id, Risk_score = p$linear.pred,
              Observed = p$times, Predicted = p$pred, Status = p$status)
}))
utils::write.csv(plot_data, file.path(out, "plot_data.csv"), row.names = FALSE)
p <- ggplot2::ggplot(plot_data, ggplot2::aes(Risk_score)) +
  ggplot2::geom_line(ggplot2::aes(y = Predicted)) +
  ggplot2::geom_point(ggplot2::aes(y = Observed, shape = factor(Status)), alpha = 0.65) +
  ggplot2::facet_wrap(~Model) + ggplot2::theme_bw() +
  ggplot2::labs(x = "Model linear predictor", y = "Days", shape = "Event",
    title = "GBSG2: recurrence-free survival",
    caption = "Line: model-derived mean-time summary; points: observed follow-up. Descriptive, not a calibration curve.")
ggplot2::ggsave(file.path(out, "gbsg2.png"), p, width = 10, height = 5, dpi = 150)
ggplot2::ggsave(file.path(out, "gbsg2.pdf"), p, width = 10, height = 5)
saveRDS(list(train = train, test = test, cc_train = cc_train, random_train = random_train,
              cc_test = cc_test, models = models, predictions = predictions,
              weights = weights, km = km, metrics = metrics), file.path(out, "gbsg2.rds"))
writeLines(c("Seed 2025; 2:1 split; 10% unstratified case-cohort sampling within each split.",
  "Weights calculated from full phase-one data before selecting phase-two rows.",
  "Both Cox models use their training baseline hazards.",
  "Evaluation weights are explicitly aligned with the time-sorted prediction subjects.",
  "This alignment also corrects the older analysis's unsorted-weight handoff; changes are not attributable to the baseline repair alone.",
  "Original marginal training censoring KM, horizons, and legacy time-grid mean integration retained.",
  "GBSG2 is a case-cohort illustration; the matched NCC repair is tested separately."),
  file.path(out, "README.txt"))
print(settings); print(metrics); print(checks)
cat("GBSG2 analysis PASS. Results:", out, "\n")
