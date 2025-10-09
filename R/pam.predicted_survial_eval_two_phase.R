
# Helper (internal)
#' @keywords internal
#' @noRd
km_surv <- function(t, km_cens) {
  if (is.na(t)) return(NA_real_)
  summary(km_cens, times = t, extend = TRUE)$surv
}

#' Evaluate Survival Model Metrics for Two-Phase Sampling Designs
#'
#' This function constructs a standardized evaluation dataset (`dat1`) from the
#' analysis dataset, predicted survival probabilities, and sampling weights, and
#' computes key performance metrics in a tidy table format.
#' Works for common two-phase sampling designs:
#'   - Case–cohort (unstratified or stratified)
#'   - Nested case–control (matched or unmatched)
#'
#' @param analysis_data A data frame containing the analysis dataset used in the fitted Cox model.
#'   Must include `time` and `status` columns.
#' @param case_weights Numeric vector of case weights (from `weighted_param()`).
#' @param pred_results A numeric vector of predicted survival probabilities with (can be calculate by `pam.coxph_restricted()`) :
#'   - `Prob`: predicted survival probabilities over time.
#'   - `pred`: predicted risk scores.
#'   - `time`: time points for `Prob` rows.
#' @param eval_time Time horizon at which to evaluate performance.
#' @param km_cens Kaplan–Meier estimate of censoring distribution (`survival::survfit` object).
#'
#' @return A tibble with columns:
#'   - `Metric`: name of metric
#'   - `Value`: numeric value
#'
#'

#' @keywords internal
pam.predicted_survial_eval_two_phase <- function(analysis_data, pred_results, eval_time, km_cens_fit, case_weights, metrics = c("Harrell’s C", "Uno’s C", "Brier Score", "Time Dependent Auc", "R_square", "L_square", "Pesudo_R")) {
  dat1_intermediate <- analysis_data %>%
    mutate(
      surv_obj = Surv(time, status),
      # Extract the survival probability at the specific evaluation time
      pred = pred_results$surv_prob[max(which(pred_results$time <= eval_time)), ],
      pred.t = pred_results$pred,
      case_weights = case_weights,
      # Determine the time for IPCW calculation based on event status and time
      weight_time = case_when(
        status == 1 & time <= eval_time ~ pmax(time, 0), # Category 1: Events before eval_time
        time > eval_time              ~ eval_time,    # Category 2: Censored after eval_time
        TRUE                          ~ NA_real_      # Category 3: Events after eval_time (no weight)
      )
    )
  
  # Now, create the nested '.pred' column using the columns from the step above.
  dat1 <- dat1_intermediate %>%
    mutate(
      .pred = purrr::map2(
        pred, weight_time,
        ~ tibble(
          .eval_time = eval_time,
          .pred_survival = .x,
          # Calculate Inverse Probability of Censoring Weight (IPCW)
          .weight_censored = if (!is.na(.y)) 1 / km_surv(.y, km_cens_fit) else 0
        )
      )
    )
    results <- list()
    
    if ("Harrell’s C" %in% metrics) {
      c_index <- survival::concordance(
        dat1$surv_obj ~ dat1$pred.t,
        weights = dat1$case_weights,
      )$concordance
      results <- append(results, list("Harrell’s C" = round(c_index, 4)))
    }
    
    if ("Uno’s C" %in% metrics) {
      c_index <- survival::concordance(
        dat1$surv_obj ~ dat1$pred.t,
        weights = dat1$case_weights,
        timewt = "n/G2",
      )$concordance
      results <- append(results, list("Uno’s C" = round(c_index,4 )))
    }
    
    # --- Brier score
    if ("Brier Score" %in% metrics) {
      brier <- yardstick::brier_survival(
        data         = dat1,
        truth        = surv_obj,
        .pred,
        case_weights = case_weights
      )$.estimate
      results <- append(results, list("Brier Score" = round(brier, 4)))
    }
    
    # --- Time-dependent AUC
    if ("Time Dependent Auc" %in% metrics) {
      auc <- yardstick::roc_auc_survival(
        data         = dat1,
        truth        = surv_obj,
        .pred,
        case_weights = case_weights
      )$.estimate
      results <- append(results, list("Time Dependent AUC" = round(auc, 4)))
    }
    
    if (any(c("R_square", "L_square", "Pesudo_R") %in% metrics)) {
      r_l_list <- pam.r2_metrics(pred.results$pred, 
                                 analysis_data$time, 
                                 analysis_data$status,
                                 tau = eval_time,
                                 case_weight = case_weights)

      if ("R_square" %in% metrics) results <- append(results, list("R_square" = round(as.numeric(r_l_list$R_square), 2)))
      if ("L_square" %in% metrics) results <- append(results, list("L_square" = round(as.numeric(r_l_list$L_square), 2)))
      if ("Pesudo_R" %in% metrics) results <- append(results, list("Pesudo_R" = round(as.numeric(r_l_list$Pseudo_R_squared), 2)))
    }

    result_df <- data.frame(
      Metric = names(results),
      Value = unlist(results, use.names = FALSE),
      stringsAsFactors = FALSE
    )
    
    print(result_df)
    
}


#' @rdname pam.predicted_survial_eval_two_phase
#' @export
pam.predicted_survial_eval_casecohort <- pam.predicted_survial_eval_two_phase

#' @rdname pam.predicted_survial_eval_two_phase
#' @export
pam.predicted_survial_eval_ncc <- pam.predicted_survial_eval_two_phase
