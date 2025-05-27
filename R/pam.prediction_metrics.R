#' @title Summary of Performance Metrics for Predicted Survival Data
#'
#' @description 
#' This function computes and summarizes various performance metrics for survival data
#' using predicted values and observed survival data. Users can specify the desired metrics
#' or compute all available metrics by default.
#'
#' @param predicted_data A numeric vector of predicted survival probabilities or scores.
#' @param survival_time A numeric vector of observed survival times.
#' @param metrics A character string or vector specifying the metrics to compute. Default is "all" to compute all available metrics. Options include:
#'   \itemize{
#'     \item "R_square": R-squared metric.
#'     \item "L_square": L-squared metric.
#'     \item "Pesudo_R": Pseudo-R-squared metric.
#'     \item "Harrells_C": Harrell's Concordance Index.
#'     \item "Unos_C": Uno's Concordance Index.
#'     \item "R_sph": Explained variation (R_sph).
#'     \item "R_sh": Explained variation (R_sh).
#'     \item "Brier_Score": Brier Score.
#'     \item "Time_Dependent_Auc": Time-dependent AUC.
#'   }
#'
#' @param status A numeric or logical vector indicating event status (1 for event, 0 for censoring).
#' @param tau An optional numeric value for restricted time horizon. Default is NULL.
#' @param t_star An optional numeric value specifying the time point for certain metrics. Default is NULL.
#'
#' @return A data frame summarizing the requested performance metrics.
#'
#' @examples
#' # Example data
#' predicted_data <- c(0.8, 0.6, 0.4, 0.2)
#' survival_time <- c(5, 8, 3, 10)
#' status <- c(1, 0, 1, 1)
#'
#' # Compute a single metric
#' result <- pam.prediction_metrics(
#' predicted_data = predicted_data,
#' survival_time = survival_time,
#' metrics = "R_square",
#' status = status
#' )
#' print(result)
#'
#' # Compute all metrics
#' result_all <- pam.prediction_metrics(
#' predicted_data = predicted_data,
#' survival_time = survival_time,
#' metrics = "all",
#' status = status
#' )
#' print(result_all)
#'
#' @export

pam.prediction_metrics <- function(predicted_data, survival_time, metrics, 
                                             status = NULL, tau = NULL, t_star = NULL, start_time = NULL) {
  
  if (missing(predicted_data) || missing(survival_time) || missing(metrics)) {
    stop("Please provide 'predicted_data','survival_time' and 'metrics'.")
  }
  
  metrics_results <- list()
  
  
  
  
  valid_metrics <- c("R_square", "L_square", "Pesudo_R", "Harrells_C", "Unos_C",
                     "R_sph", "R_sh", "Brier Score", "Time Dependent Auc")
  
  if ("all" %in% metrics) {
    metrics <- valid_metrics
  } else {
    invalid <- setdiff(metrics, valid_metrics)
    if (length(invalid) > 0) stop("Invalid metrics: ", paste(invalid, collapse = ", "))
  }
  
  if (any(c("R_square", "L_square", "Pesudo_R") %in% metrics)) {
    r_l_list <- pam.r2_metrics(predicted_data, survival_time, status, tau)
    if ("R_square" %in% metrics) metrics_results$R_square <- r_l_list$R_square
    if ("L_square" %in% metrics) metrics_results$L_square <- r_l_list$L_square
    if ("Pesudo_R" %in% metrics) metrics_results$Pesudo_R <- r_l_list$Pseudo_R_squared
  }
  if ("Harrells_C" %in% metrics) {
    metrics_results$Harrells_C <- round(pam.concordance_metric(predicted_data, survival_time, status, weight = "H"), 2)
  }
  
  if ("Unos_C" %in% metrics) {
    metrics_results$Unos_C <- round(pam.concordance_metric(predicted_data, survival_time, status, weight = "U", input_tau = tau), 2)
  }
  if ("R_sph" %in% metrics) {
    metrics_results$R_sph <- round(pam.rsph_metric(predicted_data, survival_time, status, start_time)$Re, 2)
  } 
  if ("R_sh" %in% metrics) {
    value_list <-  pam.rsh_metric(predicted_data, survival_time, status)
    metrics_results$R_sh <- round(value_list$Dx, 2)
  } 
  if ("Brier Score" %in% metrics) {
    metrics_results$"Brier Score" <- round(pam.Brier_metric(predicted_data, survival_time, t_star), 2)
  } 
  if ("Time Dependent Auc" %in% metrics) {
    if(!is.null(t_star)){
      pred_time <- t_star
    } else {
      pred_time <- quantile(survival_time$time, 0.5)
    }
    
    auc <- pam.survivalROC(Stime = survival_time$time, status = survival_time$status, 
                           marker = predicted_data, predict.time = pred_time, 
                           method = "KM")$AUC
    metrics_results$"Time Dependent Auc" <- round(max(auc, 1 - auc), 2)
  }
  
  result_df <- data.frame(
    Metric = names(metrics_results),
    Value = unlist(metrics_results, use.names = FALSE),
    stringsAsFactors = FALSE
  )
  
  return(result_df)
}