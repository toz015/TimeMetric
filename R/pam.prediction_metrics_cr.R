#' @title Summary of Performance Metrics for Competing Risks Survival Predictions
#'
#' @description
#' This function computes performance metrics for survival models with competing risks,
#' using predicted cumulative incidence functions (CIFs) and observed event data. It handles
#' competing events by truncating observation times at a specified horizon (tau) and treating
#' non-target events as censored. Metrics are computed using processed predictions and adjusted
#' event statuses.
#'
#' @param ftime Numeric vector of observed event/censoring times.
#' @param fstatus Numeric vector indicating event status:
#' - 0: censored
#' - 1: primary event of interest
#' - Other values: competing events
#' @param metrics Character vector specifying metrics to compute. Available options:
#' \itemize{
#' \item "R_square": Explained variation metric
#' \item "L_square": Loss-based metric
#' \item "Pesudo_R": Pseudo R-squared
#' \item "Harrells_C": Harrell's concordance index
#' \item "Unos_C": Uno's concordance index
#' \item "R_sph": Spherical explained variation
#' \item "R_sh": Schoeffel's explained variation
#' \item "Brier_Score": Brier score at t_star
#' \item "Time_Dependent_Auc": Time-dependent AUC curve
#' \item "all": Compute all available metrics
#' }
#' @param pred.cif Matrix of predicted cumulative incidence functions (CIFs) for the primary event.
#' Rows represent subjects, columns correspond to time points in time.cif.
#' @param time.cif Numeric vector of time points associated with columns in pred.cif.
#' @param event.type Numeric value specifying which event in fstatus is the primary event
#' (default = 1).
#' @param tau truncation time for evaluation. Times beyond tau are truncated, and
#' competing events are treated as censored at tau (default = NULL).
#' @param t_star Optional time point for metrics requiring fixed-time evaluation (e.g., Brier score).
#' @param start_time Optional start time for evaluation window in time-dependent metrics.
#'
#' @return A data frame summarizing computed performance metrics with columns:
#' - Metric: Name of the performance metric
#' - Value: Computed value
#' 
#' @examples
#' left empty
#'
#' @export
#'

pam.prediction_metrics_cr <- function(ftime, fstatus, metrics, pred.cif, time.cif, event.type = 1, tau = NULL, t_star = NULL, start_time = NULL) {
  
  
  if (missing(ftime) || missing(fstatus) || missing(metrics) || missing(pred.cif) || missing(time.cif)) {
    stop("Please provide 'ftime','fstatus',''metrics','pred.cif' and 'time.cif'.")
  }
  
  ftime.new <- ifelse(ftime > tau, tau, ftime)
  ftime.new <- ifelse(ftime.new <= tau & !(fstatus %in% c(event.type, 0)), 
                      tau, ftime.new)
  fstatus.new <- ifelse(ftime.new == tau | fstatus > 0, 1, fstatus)
  
  
  pred <- apply(pred.cif, 2, m_cif, time.cif = time.cif, tau = tau)
  
  out <- prediction_metrics(pred, ftime.new, metrics, status = fstatus.new, tau = tau, t_star = t_star, start_time = start_time)
  
  return(out)
}