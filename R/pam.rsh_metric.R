#' Compute R_sh Metric
#'
#' This function calculates the \( R_{sh} \) metric, a measure of explained variation for survival models,
#' using only predicted survival probabilities and observed survival data.
#'
#' @param predicted_data A numeric vector of predicted survival probabilities.
#' @param survival_time A numeric vector of observed survival times.
#' @param status A numeric vector indicating event occurrence (1 for event, 0 for censoring).
#'
#' @return A list with the following components:
#' \item{D}{The total variation in the survival data.}
#' \item{Dx}{The unexplained variation by the predictions.}
#' \item{V}{The explained variation (\( R_{sh} \)).}
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
#'
#' @export
pam.rsh_metric <- function(predicted_data, survival_time, status) {
  if (length(predicted_data) != length(survival_time) || length(predicted_data) != length(status)) {
    stop("All input vectors must have the same length.")
  }
  
  # Sort data by survival time
  data <- data.frame(predicted_data, survival_time, status)
  data <- data[order(data$survival_time), ]
  
  # Calculate Kaplan-Meier survival probabilities
  km_fit <- survival::survfit(survival::Surv(survival_time, status) ~ 1)
  km_surv <- stats::approx(km_fit$time, km_fit$surv, xout = data$survival_time, 
                           method = "constant", f = 0, yleft = 1, yright = min(km_fit$surv, na.rm = TRUE))$y
  
  # Calculate Mt (observed variation)
  Mt <- (1 - km_surv) * (1 - data$status) + km_surv * data$status
  D <- mean(Mt)
  
  # Calculate Mtx (unexplained variation using predictions)
  Mtx <- (1 - predicted_data) * (1 - data$status) + predicted_data * data$status
  Dx <- mean(Mtx)
  
  # Calculate R_sh
  R_sh <- (D - Dx) / D
  
  return(list(
    D = round(D, 4),
    Dx = round(Dx, 4),
    R_sh = round(R_sh, 4)
  ))
}
