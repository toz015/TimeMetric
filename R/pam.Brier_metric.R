#' @title Calculate Brier Score for Survival Data
#'
#' @description
#' Computes the Brier score for survival data to evaluate the accuracy of predicted survival probabilities at a specific time point. 
#' The Brier score measures the mean squared error between the predicted probabilities and the actual outcomes, adjusted for censoring using the Kaplan-Meier estimator.
#'
#' @param suvival_time A survival object created using the \code{Surv} function. It must include survival times and event status indicators.
#' @param predicted_data A numeric vector of predicted survival probabilities at the specified time point t_star for each individual in the dataset. The length must match the number of rows in \code{suvival_time}.
#' @param t_star A positive numeric value specifying the time point at which the Brier score is calculated.
#'
#' @return 
#' A single numeric value representing the Brier score, rounded to six decimal places. The score reflects the accuracy of the predicted survival probabilities at the specified time point.
#'
#' @details
#' The Brier score is calculated by partitioning the observations into two groups:
#' \itemize{
#'   \item Observations with event times less than \code{t_star}.
#'   \item Observations with event times greater than or equal to \code{t_star}.
#' }
#' The score is adjusted for censoring using the Kaplan-Meier estimate of the censoring distribution. Observations are weighted accordingly to account for the censoring bias.
#' @references
#' Graf, E., Schmoor, C., Sauerbrei, W., & Schumacher, M. (1999). 
#' Assessment and comparison of prognostic classification schemes for survival data. \emph{Statistics in Medicine}, 18(17-18), 2529-2545.
#'
#' Brier, G. W. (1950). Verification of forecasts expressed in terms of probability. \emph{Monthly Weather Review}, 78, 1-3.
#'
#' Gneiting, T., & Raftery, A. E. (2007). Strictly Proper Scoring Rules, Prediction, and Estimation. \emph{Journal of the American Statistical Association}, 102(477), 359-378.
#'
#' Zhou, H., Cheng, X., Wang, S., Zou, Y., & Wang, H. (2022). 
#' \emph{SurvMetrics: Predictive Evaluation Metrics in Survival Analysis} (R package version 0.5.0). 
#' Available at \url{https://github.com/skyee1/SurvMetrics}.
#'
#' @importFrom survival Surv
#' @importFrom stats median
#' @examples
#' library(survival)
#' # Example dataset
#' data(lung)
#' lung$SurvObj <- with(lung, Surv(time, status == 2))
#' 
#' # Simulated predicted probabilities
#' set.seed(123)
#' predicted_probs <- runif(nrow(lung), 0, 1)
#' 
#' # Choose a specific time point (t_star) for calculating the Brier Score 
#' # For simplicity, we use 200 as t_star.
#' brier_score <- pam.Brier_metric(
#'   predicted_data = predicted_probs,
#'   suvival_time = lung$SurvObj,
#'   t_star = 200
#' )
#' print(brier_score)
#' @keywords internal
#' @noRd

pam.Brier_metric <- function(predicted_data, suvival_time, t_star = -1) {
  if (!inherits(suvival_time, "Surv")) {
    stop("suvival_time must be a survival object created using Surv().")
  }
  
  if (length(suvival_time[, 1]) != length(predicted_data)) {
    stop("Length of suvival_time and predicted_data must match.")
  }
  
  if (any(is.na(predicted_data))) {
    stop("The input probability vector be calculate by predicted_data cannot have NA")
  }

  
  time <- suvival_time[, 1]
  status <- suvival_time[, 2]
  t_order <- order(time)
  time <- time[t_order]
  status <- status[t_order]
  predicted_data <- predicted_data[t_order]
  
  if (t_star < 0){
    t_star <- median(time)
  }
  
  G_t_star <- Gt(suvival_time, t_star)
  
  sum_before_t <- 0
  sum_after_t <- 0
  
  n <- length(time)
  
  for (i in c(1:n)) {
    
    if (time[i] < t_star && status[i] == 1) {
      Gti <- Gt(Surv(time, status), time[i])
      if (is.na(Gti)) {
        next
      }
      sum_before_t <- sum_before_t + 1/Gti * (predicted_data[i])^2
      next
    }
    if (time[i] >= t_star) {
      if (is.na(G_t_star)) {
        next
      }
      sum_after_t <- sum_after_t + 1/Gt(Surv(time, status), 
                                        t_star) * (1 - predicted_data[i])^2
    }
  }
  BSvalue <- (sum_before_t + sum_after_t)/length(time)
  names(BSvalue) <- "Brier Score"
  return(round(BSvalue, 6))
}
