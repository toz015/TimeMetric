#' @title RE Measure of Explained Variation Using Predicted and Observed Data
#'
#' @description This function computes the RE (explained variation) metric using 
#' predicted risk scores and observed survival times. It evaluates the proportion 
#' of variability in survival times explained by the covariates using a simplified 
#' approach that requires only the predicted values and observed survival data.
#'
#' @param predicted_data A numeric vector of predicted risk scores or linear predictors from a model.
#' @param survival_time A numeric vector of observed survival times.
#' @param status A numeric or logical vector indicating event status (1 for event, 0 for censoring).
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item \code{Re}: The RE measure, representing the proportion of explained variation.
#'   \item \code{numerator}: The numerator of the RE calculation, representing explained variation.
#'   \item \code{denominator}: The denominator of the RE calculation, representing total variation.
#' }
#'
#' @details This function assumes the input data is complete and properly formatted.
#' It adjusts for censoring using Kaplan-Meier weights and computes ranks for the predicted values
#' to evaluate the explained variation in survival times.
#' @references
#' Schemper, M., & Henderson, R. (2000). Predictive accuracy and explained variation in Cox regression. 
#' Biometrics, 56, 249–255.
#' 
#' Lusa, L., Miceli, R., & Mariani, L. (2007). Estimation of predictive accuracy in survival analysis using R and S-PLUS. 
#' Computer Methods and Programs in Biomedicine, 87, 132–137.
#' 
#' @examples
#' library(survival)
#' 
#' data("lung")
#' 
#' predicted_data <- lung$ph.ecog       
#' survival_time <- lung$time          
#' status <- lung$status - 1           
#' 
#' complete_cases <- complete.cases(predicted_data, survival_time, status)
#' predicted_data <- predicted_data[complete_cases]
#' survival_time <- survival_time[complete_cases]
#' status <- status[complete_cases]
#' 
#' # Compute the RE measure
#' result <- pam.rsph_metricc(predicted_data, survival_time, status)
#' 
#' cat("RE Measure:", result$Re, "\n")
#' @keywords internal
#' @noRd
pam.rsph_metric <- function(time, status, risk_score) {
    # Input checks
    stopifnot(length(time) == length(status), 
              length(time) == length(risk_score))
    
    # Remove missing data
    complete <- complete.cases(time, status, risk_score)
    time <- time[complete]
    status <- status[complete]
    risk_score <- risk_score[complete]
    
    # Order data by time (ascending) and status (events first)
    ord <- order(time, -status)
    time <- time[ord]
    status <- status[ord]
    risk_score <- risk_score[ord]
    
    # Unique event times (where status = 1)
    event_times <- unique(time[status == 1])
    n_times <- length(event_times)
    n_subjects <- length(time)
    
    # --- Censoring Adjustment (Gmat) ---
    # Fit Kaplan-Meier for censoring distribution (status=0 means censored)
    km_censoring <- survfit(Surv(time, 1 - status) ~ 1)
    # Extract P(uncensored) at event times
    G <- summary(km_censoring, times = event_times)$surv  
    # Weight matrix (subjects x event times)
    Gmat <- matrix(G, nrow = n_subjects, ncol = n_times, byrow = TRUE)
    
    mean_rank <- perfect_rank <- observed_rank <- numeric(n_times)
    
    for (i in 1:n_times) {
      t_i <- event_times[i]
      at_risk <- (time >= t_i)  # Subjects still under observation at t_i
      n_events <- sum(time == t_i & status == 1)  # Number of events at t_i
      
      # --- Rank Calculations ---
      # Higher risk_score = higher predicted risk → lower rank
      ranks <- rank(-risk_score[at_risk])
      
      # (A) Expected rank under null model (average rank)
      r0 <- mean(ranks)
      mean_rank[i] <- r0 * sum(1 / Gmat[at_risk, i][1:n_events])
      
      # (B) Perfect prediction rank (events have highest risk)
      rP <- 0.5 + 0.5 * sum((time[at_risk] == t_i) / Gmat[at_risk, i])
      perfect_rank[i] <- rP * sum(1 / Gmat[at_risk, i][1:n_events])
      
      # (C) Observed rank sum for events
      observed_rank[i] <- sum(ranks[1:n_events] / Gmat[at_risk, i][1:n_events])
    }
    
    # --- Compute R² ---
    numerator <- sum(mean_rank - observed_rank)   # Explained variation
    denominator <- sum(mean_rank - perfect_rank)  # Total possible variation
    r2 <- numerator / denominator
    
    return(list(r2 = r2, numerator = numerator, denominator = denominator))
  }



