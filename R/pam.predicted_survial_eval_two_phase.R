
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
#'   - `status`: indicating event occurrence (1 for event, 0 for censoring, 2 for competing events).
#' @param t_star (Optional) A numeric value specifying the evaluation time for Brier Score and time-dependent AUC. 
#'   If NULL, it defaults to the median survival time.
#' @param tau (Optional) A numeric value specifying the truncation time for calculating explained variation metrics. 
#'   If NULL, it defaults to the maximum observed survival time.
#' @param km_cens Kaplan–Meier estimate of censoring distribution (`survival::survfit` object).
#'
#' @return A tibble with columns:
#'   - `Metric`: name of metric
#'   - `Value`: numeric value
#'
#'
#' @keywords internal
#' @importFrom dplyr mutate case_when
#' @importFrom tibble tibble
#' @importFrom purrr map2
pam.predicted_survial_eval_two_phase <- function(pred_results, 
                                                 t_star = NULL, tau = 10e10, 
                                                 km_cens_fit, case_weights, 
                                                 metrics = c("Pesudo_R", "Harrell’s C", "Uno’s C", "Brier Score", "Time Dependent Auc")) {
  
  if (is.null(t_star)) t_star <- quantile(pred_results$time, 0.5)
  if (is.null(tau)) tau <- max(pred_results$time)
  
  dat1_intermediate <- tibble(time = pred_results$time, 
                                  status = pred_results$status,
                                  pred = pred_results$surv_prob[, max(which(pred_results$time <= t_star))],
                                  pred.t = pred_results$pred,
                                  case_weights = case_weights)
  dat1_intermediate <- dat1_intermediate %>%
    mutate(
      surv_obj = Surv(pred_results$time, 
                      pred_results$status),
      # Determine the time for IPCW calculation based on event status and time
      weight_time = case_when(
        status == 1 & time <= t_star ~ pmax(time, 0), # Category 1: Events before t_star
        time > t_star              ~ t_star,    # Category 2: Censored after t_star
        TRUE                          ~ NA_real_      # Category 3: Events after t_star (no weight)
      )
    )
  
  # Now, create the nested '.pred' column using the columns from the step above.
  dat1 <- dat1_intermediate %>%
    mutate(
      .pred = purrr::map2(
        pred, weight_time,
        ~ tibble(
          .eval_time = t_star,
          .pred_survival = .x,
          # Calculate Inverse Probability of Censoring Weight (IPCW)
          .weight_censored = if (!is.na(.y)) 1 / km_surv(.y, km_cens_fit) else 0
        )
      )
    )
    results <- list()
    
    if (any(c("R_square", "L_square", "Pesudo_R") %in% metrics)) {
      r_l_list <- pam.r2_metrics(pred_results$pred, 
                                 pred_results$time, 
                                 pred_results$status,
                                 tau = tau,
                                 case_weight = case_weights)
      
      if ("R_square" %in% metrics) results <- append(results, list("R_square" = round(as.numeric(r_l_list$R_square), 2)))
      if ("L_square" %in% metrics) results <- append(results, list("L_square" = round(as.numeric(r_l_list$L_square), 2)))
      if ("Pesudo_R" %in% metrics) results <- append(results, list("Pesudo_R" = round(as.numeric(r_l_list$Pseudo_R_squared), 2)))
    }
    
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
    

    result_df <- data.frame(
      Metric = names(results),
      Value = unlist(results, use.names = FALSE),
      stringsAsFactors = FALSE
    )
    
    result_df
    
}


#' Summarize multiple two-phase (case-cohort) models into a wide table
#'
#' @description
#' For a list of models (each providing `time`, `status`, `surv_prob`, `pred`),
#' calls `pam.predicted_survial_eval_two_phase()` for each and returns a wide
#' comparison table with metrics as rows and model names as columns.
#'
#' @param models A **named list** in which each element corresponds to a fitted model.  
#'   Each model entry must itself be a list containing:
#'   \itemize{
#'     \item \code{times} — numeric vector of observed follow-up times.
#'     \item \code{surv_prob} — an \eqn{n \times K} matrix (or data frame) of
#'           subject-specific predicted survival probabilities on a common time grid.
#'     \item \code{status} — event indicator (1 = event, 0 = censored).
#'     \item \code{pred} — (optional) predicted mean survival time (restricted or unrestricted).
#'     \item \code{new_data} — (optional) dataset used for prediction.
#'     \item \code{covs} — character vector of covariate names used for prediction.
#'     \item \code{model} — (optional) the underlying fitted survival model object.
#'   }
#'   Note: If \code{pred} is provided, it will be used to calculate R2 and concordence measure.
#' @param case_weights Numeric vector of Prentice (or other) sampling weights for all subjects.
#' @param km_cens A Kaplan–Meier fit for censoring (based on training data) used to compute IPCW.
#' @param metrics Character vector of metrics to compute (default shown below).
#' @param t_star Optional numeric scalar, specify the time point to evaluate Brier score and AUC.Default is median of observation time.
#' @param tau Optional numeric scalar, specify the max time horizon for R2 measure and concordence measure (default = 10e10).
#' @param digits Integer number of decimal places to round values (default 2).
#'
#' @return A data.frame: rows = metrics, columns = model names.
#' @export
pam.sample_design <- function(models,
                              case_weights,
                              km_cens,
                              metrics = c("Pesudo_R", "Harrell’s C", "Uno’s C",
                                          "Brier Score", "Time Dependent Auc"),
                              t_star = NULL, tau = NULL,
                              digits = 2) {
  
  if (!is.list(models) || length(models) == 0)
    stop("'models' must be a non-empty named list.")
  if (is.null(names(models)) || any(names(models) == ""))
    names(models) <- paste0("Model_", seq_along(models))
  
  required <- c("times", "status", "surv_prob", "pred")
  
  eval_one <- function(mod, name) {
    if (!all(required %in% names(mod))) {
      stop(sprintf("Model '%s' must contain: %s",
                   name, paste(required, collapse = ", ")))
    }
    
    # Build the single-model input structure expected by your evaluator
    pred_results <- list(
      time      = mod$times,
      status    = mod$status,
      surv_prob = mod$surv_prob,
      pred      = mod$pred
    )
    
    res <- pam.predicted_survial_eval_two_phase(
      pred_results = pred_results,
      t_star       = t_star,    # global (NULL allowed)
      tau          = tau,       # global (NULL allowed)
      km_cens_fit  = km_cens,
      case_weights = case_weights,
      metrics      = metrics
    )
    
    if (!all(c("Metric", "Value") %in% names(res)))
      stop(sprintf("Unexpected result from pam.predicted_survial_eval_two_phase() for '%s'", name))
    
    # Normalize minor label variant
    res$Metric <- gsub("^Time Dependent AUC$", "Time Dependent Auc", res$Metric)
    
    res$Model <- name
    res[, c("Model", "Metric", "Value")]
  }
  
  # Evaluate all models
  pieces <- mapply(eval_one, models, names(models), SIMPLIFY = FALSE)
  long   <- do.call(rbind, pieces)
  rownames(long) <- NULL
  
  # Pivot to wide (rows = Metric, cols = Model)
  wide <- reshape(
    long,
    idvar = "Metric",
    timevar = "Model",
    direction = "wide"
  )
  names(wide)    <- sub("^Value\\.", "", names(wide))
  rownames(wide) <- NULL
  wide           <- wide[, c("Metric", setdiff(names(wide), "Metric")), drop = FALSE]
  
  # Preferred ordering (keep present ones)
  preferred <- c("Pesudo_R", "R_square", "L_square",
                 "Harrell’s C", "Uno’s C",
                 "Brier Score", "Time Dependent Auc")
  present <- intersect(preferred, wide$Metric)
  others  <- setdiff(wide$Metric, preferred)
  wide$Metric <- factor(wide$Metric, levels = c(present, sort(others)))
  wide <- wide[order(wide$Metric), ]
  wide$Metric <- as.character(wide$Metric)
  
  # Round numeric columns
  num_cols <- setdiff(names(wide), "Metric")
  wide[num_cols] <- lapply(wide[num_cols], function(x) if (is.numeric(x)) round(x, digits) else x)
  
  wide
}
