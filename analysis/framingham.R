# Third-examination illustration. Source common.R and load this checkout first.
# No examination is selected by looking for the last visit before a future event.
paper_framingham_cohort <- function(framingham) {
  covs <- c("SEX", "TOTCHOL", "AGE", "SYSBP", "DIABP", "CURSMOKE",
            "CIGPDAY", "BMI", "DIABETES", "HEARTRTE", "GLUCOSE", "educ")
  required <- c("RANDID", "PERIOD", "TIME", "PREVHYP", "DEATH", "TIMEDTH",
                "HYPERTEN", "TIMEHYP", covs)
  if (!all(required %in% names(framingham))) stop("Missing Framingham columns.")
  baseline <- framingham[framingham$PERIOD == 1, ]
  third <- framingham[framingham$PERIOD == 3, ]
  stopifnot(!anyDuplicated(baseline$RANDID), !anyDuplicated(third$RANDID))

  # Preserve the original baseline exclusions and competing-event tie rule.
  # TIMEHYP/TIMEDTH are measured from the first examination, not the third.
  baseline$endpoint_time <- pmin(baseline$TIMEDTH, baseline$TIMEHYP)
  baseline$status <- ifelse(baseline$endpoint_time == baseline$TIMEDTH,
                            baseline$DEATH * 2L, baseline$HYPERTEN)
  eligible <- baseline[stats::complete.cases(baseline[, c("RANDID", "endpoint_time", "status")]) &
                         !is.na(baseline$PREVHYP) & baseline$PREVHYP == 0 &
                         baseline$endpoint_time > 0, c("RANDID", "endpoint_time", "status")]
  at_third <- merge(eligible, third[, c("RANDID", "TIME", "PREVHYP", covs)],
                    by = "RANDID", sort = TRUE)
  event_free <- at_third[!is.na(at_third$TIME) & !is.na(at_third$PREVHYP) &
                          at_third$PREVHYP == 0 & at_third$endpoint_time > at_third$TIME, ]
  cohort <- event_free[stats::complete.cases(event_free[, covs]), ]
  names(cohort)[names(cohort) == "TIME"] <- "landmark_time"
  cohort$time <- cohort$endpoint_time - cohort$landmark_time
  cohort$PREVHYP <- NULL
  cohort <- cohort[, c("RANDID", "landmark_time", "endpoint_time", "time", "status", covs)]
  rownames(cohort) <- NULL
  stopifnot(!anyDuplicated(cohort$RANDID), all(cohort$time > 0),
            all(cohort$status %in% 0:2))
  flow <- data.frame(
    Step = c("Baseline participants", "Baseline eligible with positive endpoint time",
             "Third-examination attendees", "Baseline eligible attending third examination",
             "Event-free at third examination with positive remaining follow-up",
             "Complete third-examination covariates"),
    N = c(nrow(baseline), nrow(eligible), nrow(third), nrow(at_third),
          nrow(event_free), nrow(cohort)))
  list(cohort = cohort, flow = flow, covs = covs,
       reduced_covs = c("SEX", "TOTCHOL", "AGE", "DIABP"))
}

paper_run_framingham <- function(root, out = file.path(root, "results", "framingham_visit3_12y")) {
  for (p in c("riskCommunicator", "cmprsk")) {
    if (!requireNamespace(p, quietly = TRUE)) stop("Install ", p, " before running this illustration.")
  }
  seed <- 1234567L
  tau <- 12 * 365 # Original paper's convention: 365 days per year.
  e <- new.env()
  utils::data("framingham", package = "riskCommunicator", envir = e)
  prepared <- paper_framingham_cohort(e$framingham)
  cohort <- prepared$cohort
  covs <- prepared$covs
  reduced_covs <- prepared$reduced_covs
  settings <- list(analysis = "Framingham third-examination illustration", seed = seed,
                   examination = 3L, days_per_year = 365, tau = tau, t_star = tau,
                   baseline_prevalent_hypertension_excluded = TRUE,
                   complete_cases = "All full-model covariates at examination 3",
                   training_followup = "Full observed follow-up after examination 3",
                   cr_r2_weighting = "Observed events; original follow-up; Cox without covariates; G floor 1e-4",
                   riskCommunicator = as.character(utils::packageVersion("riskCommunicator")),
                   cmprsk = as.character(utils::packageVersion("cmprsk")))
  paper_manifest(root, out, settings)
  utils::write.csv(prepared$flow, file.path(out, "cohort_flow.csv"), row.names = FALSE)

  set.seed(seed)
  idx <- sample(seq_len(nrow(cohort)), round(2/3 * nrow(cohort)))
  df <- cohort[, c("time", "status", covs)]
  train <- df[idx, ]; test <- df[-idx, ]
  membership <- cohort
  membership$Split <- ifelse(seq_len(nrow(cohort)) %in% idx, "Train", "Test")
  utils::write.csv(membership, file.path(out, "cohort_and_split.csv"), row.names = FALSE)

  outcome_summary <- do.call(rbind, lapply(c("All", "Train", "Test"), function(split) {
    d <- if (split == "All") df else if (split == "Train") train else test
    data.frame(Split = split, N = nrow(d), Hypertension = sum(d$status == 1),
               Competing_death = sum(d$status == 2), Censored = sum(d$status == 0),
               Median_followup_years = median(d$time)/365,
               Hypertension_by_12y = sum(d$status == 1 & d$time <= tau),
               Competing_death_by_12y = sum(d$status == 2 & d$time <= tau),
               Censored_before_12y = sum(d$status == 0 & d$time < tau),
               Observed_event_free_beyond_12y = sum(d$time > tau))
  }))
  utils::write.csv(outcome_summary, file.path(out, "outcome_counts.csv"), row.names = FALSE)

  # Keep the four original model specifications and the common held-out sample.
  models <- list(
    FG = cmprsk::crr(train$time, train$status, as.matrix(train[, covs]), failcode = 1, cencode = 0),
    FG.Reduced = cmprsk::crr(train$time, train$status, as.matrix(train[, reduced_covs]),
                            failcode = 1, cencode = 0),
    AFT.Hypertension = survival::survreg(survival::Surv(time, status == 1) ~ .,
                                        data = train, x = TRUE, y = TRUE, dist = "weibull"),
    AFT.Death = survival::survreg(survival::Surv(time, status == 2) ~ .,
                                data = train, x = TRUE, y = TRUE, dist = "weibull"),
    AFT.Fixed.Hypertension = survival::survreg(survival::Surv(time, status == 1) ~ .,
                                              data = train, x = TRUE, y = TRUE,
                                              dist = "weibull", scale = 5),
    AFT.Fixed.Death = survival::survreg(survival::Surv(time, status == 2) ~ .,
                                      data = train, x = TRUE, y = TRUE,
                                      dist = "weibull", scale = 5))
  if (!isTRUE(models$FG$converged) || !isTRUE(models$FG.Reduced$converged)) {
    stop("A Fine-Gray model did not converge; do not use its metrics.")
  }
  predictions <- list(
    FG = tm_predict_cif(fg_model = models$FG, newdata = test, covs = covs, tau = tau),
    FG.Reduced = tm_predict_cif(fg_model = models$FG.Reduced, newdata = test,
                               covs = reduced_covs, tau = tau),
    CSH.AFT = tm_predict_cif(model1 = models$AFT.Hypertension, model2 = models$AFT.Death,
                            newdata = test, covs = covs, tau = tau),
    CSH.AFT.FIX = tm_predict_cif(model1 = models$AFT.Fixed.Hypertension,
                                model2 = models$AFT.Fixed.Death,
                                newdata = test, covs = covs, tau = tau))
  for (p in predictions) {
    stopifnot(identical(p$times, test$time), identical(p$status, test$status),
              ncol(p$cif_pred) == nrow(test) + 1L,
              all(is.finite(p$cif_pred)), all(is.finite(p$pred)),
              all(p$pred > 0 & p$pred <= tau),
              all(p$cif_pred[, -1] >= -1e-8 & p$cif_pred[, -1] <= 1 + 1e-8))
  }
  # Explicit metrics retain R/L components and omit the withdrawn R_SH/R_E.
  metrics <- tm_summarize_cr(predictions, metrics = c("r_square", "l_square", "pseudo_r2",
                              "pseudo_r2_point", "c_index", "brier_score", "td_auc"),
                             event_type = 1, tau = tau, t_star = tau, digits = 4,
                             censoring = "cox")
  stopifnot(all(is.finite(as.matrix(metrics[, -1]))))
  utils::write.csv(metrics, file.path(out, "metrics.csv"), row.names = FALSE)
  paper_write_md(metrics, file.path(out, "metrics.md"))

  # Record actual support and the CIF step used at the requested 12-year horizon.
  # The package uses the last CIF knot <= t_star, not nearest-time matching.
  support <- do.call(rbind, lapply(c("Train", "Test"), function(split) {
    d <- if (split == "Train") train else test
    km <- survival::survfit(survival::Surv(time, status == 0) ~ 1, data = d)
    data.frame(Split = split, N = nrow(d), tau_days = tau, t_star_days = tau,
               Maximum_followup_days = max(d$time), At_risk_at_12y = sum(d$time >= tau),
               G_12y = if (tau <= max(d$time))
                 summary(km, times = tau, extend = TRUE)$surv else NA_real_)
  }))
  stopifnot(all(support$Maximum_followup_days >= tau), all(support$At_risk_at_12y > 0),
            all(is.finite(support$G_12y)), all(support$G_12y > 0))
  utils::write.csv(support, file.path(out, "followup_support.csv"), row.names = FALSE)
  # Event-observation R2 weights can depend on censoring AFTER the horizon.
  # G(12y) alone is therefore not a diagnostic for this weighting strategy.
  r2_ipcw <- getFromNamespace("tm_cr_event_ipcw", "TimeMetric")(
    test$time, test$status, censoring = "cox")
  r2_weights <- data.frame(RANDID = cohort$RANDID[-idx],
    original_time = test$time, original_status = test$status,
    censoring_survival_at_original_time = r2_ipcw$probability,
    normalized_event_ipcw = r2_ipcw$weight)
  stopifnot(all(r2_weights$normalized_event_ipcw[test$status == 0] == 0))
  utils::write.csv(r2_weights, file.path(out, "r2_original_followup_weights.csv"), row.names = FALSE)
  horizon_knots <- do.call(rbind, lapply(names(predictions), function(m) {
    grid <- predictions[[m]]$cif_pred[, 1]
    data.frame(Model = m, tau_days = tau, t_star_days = tau,
               CIF_knot_used_at_t_star = max(grid[grid <= tau]),
               Maximum_CIF_grid_time = max(grid))
  }))
  utils::write.csv(horizon_knots, file.path(out, "analysis_settings.csv"), row.names = FALSE)
  predictions_12y <- data.frame(RANDID = cohort$RANDID[-idx], time = test$time, status = test$status)
  for (m in names(predictions)) {
    p <- predictions[[m]]
    k <- max(which(p$cif_pred[, 1] <= tau))
    predictions_12y[[paste0(m, "_CIF_12y")]] <- p$cif_pred[k, -1]
    predictions_12y[[paste0(m, "_restricted_mean_days")]] <- p$pred
  }
  utils::write.csv(predictions_12y, file.path(out, "test_predictions_12y.csv"), row.names = FALSE)

  # Use the same 200 held-out subjects in every panel, in years since examination 3.
  set.seed(seed)
  plot_idx <- sample(seq_len(nrow(test)), min(200L, nrow(test)))
  plot_data <- do.call(rbind, lapply(names(predictions), function(m) {
    p <- predictions[[m]]
    inv <- startsWith(m, "CSH")
    data.frame(Model = m, RANDID = cohort$RANDID[-idx][plot_idx],
               Risk_score = (if (inv) -p$linear.pred else p$linear.pred)[plot_idx],
               Observed_years = pmin(p$times[plot_idx], tau)/365,
               Status = p$status[plot_idx], Predicted_years = p$pred[plot_idx]/365)
  }))
  utils::write.csv(plot_data, file.path(out, "plot_data.csv"), row.names = FALSE)
  plot_data$Model <- factor(plot_data$Model, levels = names(predictions),
                            labels = c("Fine-Gray", "Fine-Gray (Reduced)",
                                       "CSH Weibull AFT", "CSH Weibull AFT (Scale = 5)"))
  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(Risk_score)) +
    ggplot2::geom_line(ggplot2::aes(y = Predicted_years), linewidth = 0.5) +
    ggplot2::geom_point(ggplot2::aes(y = Observed_years, shape = factor(Status)),
                        size = 1.2, alpha = 0.7) +
    ggplot2::scale_shape_manual(name = "Observed status", values = c("0" = 2, "1" = 19, "2" = 1),
                                 breaks = c("0", "1", "2"),
                                 labels = c("Censored", "Hypertension", "Competing death")) +
    ggplot2::geom_hline(yintercept = 10, linetype = "dashed", colour = "grey50") +
    ggplot2::facet_wrap(~ Model, ncol = 2, scales = "free_x") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(x = "Model risk score", y = "Years since third examination",
      title = "Framingham: prediction from the third examination",
      subtitle = "12-year horizon; 200 held-out participants",
      caption = paste("Line: integral of 1 - hypertension CIF from 0 to 12 years.",
                      "Points: observed follow-up capped at 12 years; shapes retain original status.",
                      "Descriptive display, not a calibration curve.", sep = "\n"))
  ggplot2::ggsave(file.path(out, "framingham_visit3_12y.png"), plot, width = 10, height = 8, dpi = 150)
  ggplot2::ggsave(file.path(out, "framingham_visit3_12y.pdf"), plot, width = 10, height = 8)
  result <- list(settings = settings, cohort = membership, flow = prepared$flow,
                 outcome_counts = outcome_summary, models = models, predictions = predictions,
                 metrics = metrics, support = support, r2_weights = r2_weights,
                 horizon_knots = horizon_knots, plot = plot)
  saveRDS(result, file.path(out, "framingham_visit3_12y.rds"))
  writeLines(c(
    "Third-examination illustration; the visit is prespecified, not chosen relative to a future outcome.",
    "Population: baseline-eligible attendees without hypertension at examination 3, complete current covariates, and positive remaining follow-up.",
    "Examination 3 occurs at different actual times (approximately year 12); it is not an exact year-12 landmark.",
    "Follow-up starts at each participant's actual third-examination date. No requirement to remain under observation for 12 years is imposed.",
    "Full observed follow-up is used to fit models; prediction/evaluation use tau = t_star = 3650 days.",
    "Status 0 = censoring, 1 = hypertension, 2 = competing death. Endpoint ties retain the original death-first rule.",
    "Generalized R2 uses event-observation IPCW from original follow-up, an intercept-only Cox censoring model and G floor 1e-4; original censoring rows have zero R2 weight even after 12 years.",
    "R2 weights and censoring survival at original follow-up are in r2_original_followup_weights.csv; these weights can depend on follow-up beyond the 12-year outcome horizon.",
    "All models use the same full-model complete-case cohort and seed-1234567 2:1 split.",
    "The reduced Fine-Gray model uses SEX, TOTCHOL, AGE, DIABP, matching the original code.",
    "Package CIF grids/quadrature are retained, including AFT's evaluation-time grid; actual probability knots are in analysis_settings.csv.",
    "Restricted mean integrates 1 - hypertension CIF through 12 years, with competing deaths contributing to the no-hypertension state.",
    "Outcome counts over full follow-up and within 12 years are distinguished in outcome_counts.csv.",
    "Results are an illustration in third-examination attendees, not a baseline-cohort performance estimate."),
    file.path(out, "README.txt"))
  print(prepared$flow, row.names = FALSE)
  print(outcome_summary, row.names = FALSE)
  print(metrics, row.names = FALSE)
  print(support, row.names = FALSE)
  cat("Framingham third-examination 12-year analysis PASS. Results:", out, "\n")
  invisible(result)
}
