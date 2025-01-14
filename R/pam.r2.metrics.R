#' Compute R^2, L^2, and Pseudo-R^2 for Survival Data
#'
#' This function calculates \( R^2 \), \( L^2 \), and pseudo-\( R^2 \) 
#' for survival data using observed and predicted survival times.
#'
#' @param predicted_data A numeric vector of predicted survival times.
#' @param survival_time A numeric vector of observed survival times.
#' @param status A numeric vector indicating the event occurrence (1 for event, 0 for censoring).
#' @param tau An optional numeric value for restricted time horizon. Default is NULL (no restriction).
#'
#' @return A list with \( R^2 \), \( L^2 \), and pseudo-\( R^2 \).
#'
#' @examples
#' 
#' library(PAmeasure)
#' predicted_data <- c(5, 4, 8, 2)
#' survival_time <- c(6, 5, 10, 3)
#' status <- c(1, 1, 0, 1)
#' pam.r2.metrics(predicted_data, survival_time, status)
#'
#' @export
pam.r2.metrics <- function(predicted_data, survival_time, status, tau = NULL) {
  # Apply restriction if tau is provided
  if (!is.null(tau)) {
    restricted <- restricted_data_gen(survival_time, status, tau)
    survival_time <- restricted$time
    status <- restricted$status
  }
  
  
  wls_fit <- tryCatch(
    lm(survival_time ~ predicted_data, weights = status),
    error = function(e) stop("WLS fitting failed: ", e$message)
  )
  
  
  calibrated <- tryCatch(
    predict(wls_fit),
    error = function(e) stop("Calibration failed: ", e$message)
  )
  
  # Calculate R^2
  mean_weighted <- sum(status * survival_time) / sum(status)
  num_r2 <- sum(status * (calibrated - mean_weighted)^2)
  denom_r2 <- sum(status * (survival_time - mean_weighted)^2)
  r2 <- round(num_r2 / denom_r2, 4)
  
  # Calculate L^2
  num_l2 <- sum(status * (survival_time - calibrated)^2)
  denom_l2 <- sum(status * (survival_time - predicted_data)^2)
  l2 <- round(num_l2 / denom_l2, 4)
  
  # Calculate pseudo-R^2
  pseudo_r2 <- round(r2 * l2, 4)
  
  return(list(
    R_squared = format(round(r2, 2), nsmall = 4),
    L_squared = format(round(l2, 2), nsmall = 4),
    Pseudo_R_squared = format(round(pseudo_r2, 2), nsmall = 4)
  ))
}

# Helper function to restrict survival data
restricted_data_gen <- function(time, status, tau) {
  output <- data.frame(
    time = ifelse(time > tau, tau, time),
    status = ifelse(time > tau, 1, status)
  )
  return(output)
}
