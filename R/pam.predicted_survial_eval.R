#' Compute Survival Model Evaluation Metrics
#'
#' This function calculates various predictive performance metrics for survival models based on 
#' predicted survival probabilities and observed survival data. It supports a range of evaluation 
#' measures, including explained variation, concordance indices, Brier Score, and time-dependent AUC.
#'
#' @param event_time A numeric vector of observed survival times.
#' @param predicted_probability A numeric vector of predicted survival probabilities.
#' @param status A numeric vector indicating event occurrence (1 for event, 0 for censoring).
#' @param metrics A character vector specifying the evaluation metrics to compute. Options include:
#'   \itemize{
#'     \item "Pseudo_R_square" - Pseudo R-squared measure
#'     \item "R_square" - Explained variation R²
#'     \item "L_square" - L-squared measure
#'     \item "Harrells_C" - Harrell’s concordance index
#'     \item "Unos_C" - Uno’s concordance index
#'     \item "R_sh" - Schemper-Henderson explained variation (R_sh)
#'     \item "Brier Score" - Brier score for calibration
#'     \item "Time Dependent Auc" - Time-dependent area under the curve (AUC)
#'   }
#'   Default is "all", which computes all available metrics.
#' @param t_star (Optional) A numeric value specifying the evaluation time for Brier Score and time-dependent AUC. 
#'   If NULL, it defaults to the median event time.
#' @param tau (Optional) A numeric value specifying the truncation time for calculating explained variation metrics. 
#'   If NULL, it defaults to the maximum observed survival time.
#'
#' @return A data frame with two columns:
#'   \item{Metric}{The name of the computed evaluation metric.}
#'   \item{Value}{The corresponding computed value.}
#'
#' @references
#' Schemper, M. and R. Henderson (2000). Predictive accuracy and explained variation in Cox regression.
#' Biometrics 56, 249--255.
#' 
#' Lusa, L., R. Miceli and L. Mariani (2007). Estimation of predictive accuracy in survival analysis
#' using R and S-PLUS. Computer Methods and Programs in Biomedicine 87, 132--137.
#' 
#' Potapov, S., Adler, W., Schmid, M., Bertrand, F. (2024). survAUC: Estimating Time-Dependent AUC for Censored Survival Data. 
#' R package version 1.3-0. DOI: \doi{10.32614/CRAN.package.survAUC}. Available at \url{https://CRAN.R-project.org/package=survAUC}.
#'
#' @examples
#' library(PAmeasure)
#' predicted_data <- c(0.8, 0.6, 0.4, 0.2)
#' survival_time <- c(5, 8, 3, 10)
#' status <- c(1, 0, 1, 1)
#' pam.rsh_metric(predicted_data, survival_time, status)
#' @export



integrate_survival <- function(predicted_probability, event_time, status, tau = NULL) {
  if (nrow(predicted_probability) != length(event_time)) {
    stop("Number of rows in predicted_probability must equal length of event_time")
  }
  
  if (any(predicted_probability < 0 | predicted_probability > 1, na.rm = TRUE)) {
    stop("All predicted probabilities must be between 0 and 1")
  }
  
  if (is.null(tau)) {
    tau <- max(event_time, na.rm = TRUE)
  }
  
  order_idx <- order(event_time)
  event_time <- event_time[order_idx]
  predicted_probability <- predicted_probability[order_idx, , drop = FALSE]
  event_time <- pmin(event_time, tau)
  t1 <- event_time
  t2 <- c(0, head(t1, -1)) 
  delta <- status[order_idx]
  delta <- ifelse(event_time <= tau, delta, 0)
  delta.t <- t1 - t2       
  
  cumulative_prediction <- colSums(delta.t * predicted_probability)
  
  return(cumulative_prediction)
}

#' @export

pam.predicted_survial_eval <- function (event_time, predicted_probability, status,
                                        metrics = NULL,  t_star = NULL, tau = NULL) 
{
  
  if (missing(event_time) || missing(predicted_probability)) {
    stop("Please provide 'predicted_data',  and 'covariates' arguments.")
  }
  
  
  
  metrics_results <- list()
  
  valid_metrics <- c("Pseudo_R_square", "R_square", "L_square", 
                     "Harrell’s C", "Uno’s C",
                      "R_sph","Brier Score", "Time Dependent Auc")
  default_metrics <- c("Pseudo_R_square", "Harrell’s C", "Uno’s C",
                       "R_sph","Brier Score", "Time Dependent Auc")
  if (is.null(metrics)){
    metrics <- default_metrics
  }
  else if ("all" %in% metrics) {
    metrics <- valid_metrics
  } else {
    invalid <- setdiff(metrics, valid_metrics)
    if (length(invalid) > 0) 
      stop("Invalid metrics: ", paste(invalid, collapse = ", "))
  }
  
  
  
  predicted_data <- integrate_survival(
    predicted_probability, event_time, status, tau)
  
  t_star <- quantile(event_time, 0.5)
  t_idx <- which.min(abs(event_time - t_star))
  risk_scores <- 1 - predicted_probability[, t_idx] 
  
  if("Pseudo_R_square" %in% metrics 
     || "R_square" %in% metrics 
     || "L_square" %in% metrics) {
    r_l_list <- pam.r2_metrics(predicted_data, event_time, status, tau)
  }
  
  if ("Pseudo_R_square" %in% metrics) {
    metrics_results$Pseudo_R_square <- round(r_l_list$Pseudo_R_squared, 4)
  } 
  if ("R_square" %in% metrics) {
    metrics_results$R_square <- round(r_l_list$R_squared,4)
  }
  if ("L_square" %in% metrics) {
    metrics_results$L_square <- round(r_l_list$L_square, 4)
  }
  
  if ("Harrell’s C" %in% metrics) {
    metrics_results$"Harrell’s C" <- round(
      concordancefit(y = Surv(event_time, status), 
                     x = risk_scores, reverse = FALSE)$concordance, 4)
  }
  
  if ("Uno’s C" %in% metrics) {
    metrics_results$"Uno’s C" <- round(
      concordancefit(y = Surv(event_time, status), 
                     x = risk_scores, reverse = FALSE, 
                     timewt = "n/G2")$concordance, 4)
  }
  
#  if ("R_sh" %in% metrics) {
#    metrics_results$R_sh <- round(pam.rsh_metric(predicted_data = predicted_data, survival_time = event_time, status = status)$R_sh, 4)
#  }
  
  if ("R_sph" %in% metrics) {
   metrics_results$R_sph <- pam.rsph_metric(
     time = event_time, status = status, risk_score = risk_scores)$r2
  }
  
  if ("Brier Score" %in% metrics) {
    t_eval <- event_time[t_idx]
    X <- risk_scores
    brier_result <- tdROC(
      X = X,  
      Y = event_time,      
      delta = status,
      tau = t_eval,       
      method = "both", 
      output = "both"   
    )
    metrics_results$"Brier Score" <- round(
      as.numeric(brier_result$calibration_res[1]), 4)
  }
  
  if ("Time Dependent Auc" %in% metrics) {
    t_eval <- event_time[t_idx]
    X <- risk_scores
    AUC_result <- tdROC(
      X = X,  
      Y = event_time,      
      delta = status,
      tau = t_eval,       
      method = "both", 
      output = "both"   
    )
    metrics_results$"Time Dependent Auc" <- round(
      AUC_result$main_res$AUC.integral, 4)
    #metrics_results$"Time Dependent Auc Empirical" <- round(AUC_result$main_res$AUC.empirical, 4)
  }
  result_df <- data.frame(
    Metric = names(metrics_results),
    Value = unlist(metrics_results, use.names = FALSE),
    stringsAsFactors = FALSE
  )
  
  return(result_df)
}

