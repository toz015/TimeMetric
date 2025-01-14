#' @title Prediction Accuracy Measures for Predicted Survival Data
#'
#' @description Computes \(R^2\), \(L^2\), and Pseudo \(R^2\) metrics for predicted survival data.. R-squared is an extension of the classical R2 statistic for a linear model, quantifying the amount of variability in the response that is explained by a corrected prediction based on the original prediction function. L-squared is the proportion of the prediction error of the original prediction function that is explained by the corrected prediction function, quantifying the distance between the corrected and uncorrected predictions. When used together, they give a complete summary of the predictive power of a prediction function.
#' @param observed_data A data frame containing the observed survival data. Must include columns:
#'   - \code{time}: Observed survival times.
#'   - \code{status}: Censoring indicator (1 = event occurred, 0 = censored).
#' @param predicted_data A numeric vector of predicted survival times or risk scores. Must have the same length as \code{observed_data}.
#'
#' @return A list containing:
#'   - \code{R.squared}: The \(R^2\) metric.
#'   - \code{L.squared}: The \(L^2\) metric.
#'   - \code{Pseudo.R}: The pseudo \(R^2\) metric, computed as the product of \(R^2\) and \(L^2\).
#'
#' @examples
#' library(survival)
#' library(PAmeasures)
#'
#'# Use Mayo Clinic Primary Biliary Cirrhosis Data
#'data(pbc)
#'
#' head(pbc)
#'



#' @export
#' 
pam.survreg.metrics <- function(observed_data, predicted_data, dist) {
  # Extract observed survival times and censoring status from validation data
  
  if (!all(c("time", "status") %in% names(observed_data))) {
    stop("The observation dataset must contain 'time' and 'status' columns.")
  }
  
  if (length(observed_data$time) != length(predicted_data)) {
    stop("Observed data and predicted data must have the same length.")
  }
  
  if (dist == "exponential") {
    predicted_data <- predicted_data * gamma(2)
  } else if (dist == "weibull") {
    predicted_data <- predicted_data * gamma(1 + fit.survreg$scale)
  } else if (dist == "lognormal") {
    predicted_data <- predicted_data * exp((fit.survreg$scale)^2 / 2)
  } else if (dist == "loglogistic") {
    predicted_data <- predicted_data * gamma(1 + fit.survreg$scale) * gamma(1 - fit.survreg$scale)
  } else if (dist == "coxph"){
    stop("Use pam.coxph to calculate coxph model's R")
  }
  
  # Extract observed survival times and censoring status
  y.unsorted <- observed_data$time
  censor.unsorted <- observed_data$status
  
  # Sort the data by survival times
  order_indices <- order(y.unsorted)
  y <- y.unsorted[order_indices]
  delta <- censor.unsorted[order_indices]
  predicted_data <- predicted_data[order_indices]
  
  # KM estimate for censoring distribution
  fit.km.censoring <- survfit(Surv(y, 1 - delta) ~ 1)
  sum.km.censoring <- summary(fit.km.censoring, times = y, extend = TRUE)
  km.censoring <- sum.km.censoring$surv
  km.censoring.minus <- c(1, km.censoring[-length(km.censoring)])
  ratio.km <- delta / km.censoring.minus
  ratio.km[is.nan(ratio.km)] <- 0
  weight.km <- ratio.km / sum(ratio.km)
  
  # Weighted least squares fit
  wls.fitted <- tryCatch(
    lm(y ~ predicted_data, weights = weight.km),
    error = function(e) {
      return(c(NA, NA))
    }
  )
  calibrate.fitted <- tryCatch(
    predict(wls.fitted),
    error = function(e) {
      return(rep(NA, length(y)))
    }
  )
  
  # Calculate R-squared
  num.rho2 <- sum(weight.km * (calibrate.fitted - sum(weight.km * y))^2)
  denom.rho2 <- sum(weight.km * (y - sum(weight.km * y))^2)
  R2 <- round(num.rho2 / denom.rho2, digits = 4)
  
  # Calculate L-squared
  num.L2 <- sum(weight.km * (y - calibrate.fitted)^2)
  denom.L2 <- sum(weight.km * (y - predicted_data)^2)
  L2 <- round(num.L2 / denom.L2, digits = 4)
  
  # Calculate Pseudo R-squared
  SR <- round(R2 * L2, digits = 4)
  
  # Return results
  return(list(
    R.squared = format(R2, nsmall = 4),
    L.squared = format(L2, nsmall = 4),
    Pseudo.R = format(SR, nsmall = 4)
  ))
}