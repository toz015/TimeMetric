#' @title Performance Metrics for Survival Analysis Models
#'
#' @description This function computes a comprehensive set of performance metrics for survival analysis models. It provides metrics such as R_square, L_square, Pseudo_R, Harrell’s C, Uno’s C, R_sph (distance-based estimator for survival predictive accuracy), R_sh, Brier Score, and Time-dependent AUC. Users can specify particular metrics and model types, enabling tailored performance evaluation for various survival models.
#'
#' @param data A data frame containing the survival data.
#' @param time_var The name of the time variable in `data` indicating survival time.
#' @param status_var The name of the status variable in `data` indicating event occurrence.
#' @param covariates A character vector of covariate names to include in the model.
#' @param model A character string or vector specifying the model types to fit (e.g., "coxph", "exp", "lognormal", "weibull"). Default is "coxph" to fit all models.
#' @param metrics A character string or vector specifying the metrics to compute. Default is "all" to compute all available metrics.
#' @param newdata (Optional) A data frame containing validation data. If `NULL`, the function 
#' uses the same data as `data` for model evaluation.
#' @param t_star (Optional) A positive numeric value specifying the time point at which the Brier score and AUC score is calculated.
#' @param tau (Optional) A time point for truncating the survival time. If provided, the function evaluates predictions up to this time point.
#'
#' @return A data frame containing the selected model's performance metrics.
#'
#' @examples
#' library(PAmeasures)
#' library(survival)
#' library(rms)
#' library(dplyr)
#'
#' Use Mayo Clinic Primary Biliary Cirrhosis Data
#' data(pbc)
#' pbc <- pbc %>% 
#'   filter(is.na(trt)==F) %>% 
#' mutate(log_albumin = log(albumin),
#'        log_bili = log(bili),
#'        log_protime = log(protime),
#'        status = ifelse(status==2, 1, 0))
#' time_var <- "time"
#' status_var <- "status"
#' covariates <- c("age", "log_albumin", "log_bili", "log_protime", "edema")
#' Call the function with all metrics and all models
#' results <- pam.performance_metrics(data = pbc, 
#'                                    time_var = time_var, 
#'                                    status_var = status_var, 
#'                                    covariates = covariates
#'                                    )
#' results2 <- pam.performance_metrics(data = pbc, 
#'                                     time_var = time_var,  
#'                                     status_var = status_var, 
#'                                     covariates = covariates, 
#'                                     model = c("lognormal", "weibull"),  
#'                                     metrics = c("R_square", "L_square", "Brier Score")
#'                                     )
#' @export

pam.performance_metrics <- function (data, time_var, status_var, covariates, model = "coxph", 
                                 metrics = "all", newdata = NULL, t_star = NULL, tau = NULL) {
  # Validate inputs
  if (missing(data) || missing(time_var) || missing(status_var) || missing(covariates)) {
    stop("Please provide 'data', 'time_var', 'status_var', and 'covariates' arguments.")
  }
  if (!is.null(newdata)) {
      test_data <- newdata
      } else {
      test_data <- data
      }
  
  # Fit models based on user input
  fits <- list()
  
  # Ensure the formula is a single string
  formula_text <- paste(
    "Surv(", time_var, ", ", status_var, ") ~ ", 
    paste(covariates, collapse = " + "), 
    sep = ""
  )
  # Convert to formula
  formula <- as.formula(formula_text)
  
  model_types <- if (("all" %in% model)) c("coxph", "exp", "lognormal", "weibull") else model
  metrics <- if (("all" %in% metrics))c("R_square", "L_square", "Pesudo_R", "Harrell’s C", "Uno’s C", "R_sph", "R_sh", "Brier Score", "Time Dependent Auc") else metrics
  # Define a list to hold metrics
  metrics_results <- list()
  
  if ("coxph" %in% model_types) {
    fits$coxph <- survival::coxph(formula, data = data, x = TRUE, y = TRUE)
  }
  if ("exp" %in% model_types) {
    fits$exp <- survival::survreg(formula, data = data, dist = "exponential", x = TRUE, y = TRUE)
  }
  if ("lognormal" %in% model_types) {
    fits$lognormal <- survival::survreg(formula, data = data, dist = "lognormal", x = TRUE, y = TRUE)
  }
  if ("weibull" %in% model_types) {
    fits$weibull <- survival::survreg(formula, data = data, dist = "weibull", x = TRUE, y = TRUE)
  }
  for (fit_name in names(fits)) {
    metrics_results[[fit_name]] <- list()
    if (is.null(tau)) {
      event_times <- data[[time_var]]
      if (length(event_times) == 0) {
        stop("No observed events to determine default tau.")
      }
      tau <- max(event_times)
    }
    if (fit_name == "coxph") {
      r_l_list <- pam.coxph_restricted(fits[[fit_name]], covariates = covariates, time_var = time_var, status_var = status_var, tau = tau, newdata = test_data) %>% Reduce("c", .) %>% as.numeric()
    } 
      else {
      r_l_list <- pam.surverg_restricted(fits[[fit_name]], covariates = covariates, time_var = time_var, status_var = status_var, tau = tau, newdata = test_data) %>% Reduce("c", .) %>% as.numeric()
    }
    # Extract metrics if requested
    if ("R_square" %in% metrics) {
      metrics_results[[fit_name]]$R_square <- round(r_l_list[1], 2)
    }
    if ("L_square" %in% metrics) {
      metrics_results[[fit_name]]$L_square <- round(r_l_list[2], 2)
    }
    
    if ( "Pesudo_R" %in% metrics ){
      metrics_results[[fit_name]]$Pesudo_R <- round(r_l_list[1] * r_l_list[2], 2)
    }
    
    if ("Harrell’s C" %in% metrics) {
      metrics_results[[fit_name]]$"Harrell’s C" <- pam.concordance(fits[[fit_name]], newdata = test_data)$concordance
    }
    
    if ("Uno’s C" %in% metrics) {
      metrics_results[[fit_name]]$"Uno’s C" <- pam.concordance(fits[[fit_name]], newdata = test_data, timewt="n/G2")$concordance
    }
    
    if ("R_sph" %in% metrics) {
      metrics_results[[fit_name]]$R_sph <- pam.rsph(fits[[fit_name]])$Re
    }
    
    
    if ("R_sh" %in% metrics) {
      if (fit_name == "coxph" ) {
        rms_coxph <- rms::cph(formula, data = data, x = TRUE, y = TRUE)
        check_factors <- function(data) {
          factors <- sapply(data, is.factor)
          if (any(factors)) {
            factor_cols <- names(data)[factors]
            warning("The following columns are factors: ", 
                    paste(factor_cols, collapse = ", "),
                    ". This function to calculate R_sph does not support factor variables. Returning NA.")
            return(TRUE)
          }
          return(FALSE)
        }
        
        # Notify users about factors in both datasets and return NA if any are found
        if (check_factors(data) || check_factors(test_data)) {
          metrics_results[[fit_name]]$R_sh <- NA
        } else {
          R_sh_coxph <- pam.schemper(rms_coxph, traindata = data, 
                                     newdata = test_data)$Dx
          metrics_results[[fit_name]]$R_sh <- R_sh_coxph 
        }
      } else {
        metrics_results[[fit_name]]$R_sh <- NA
      }
    }
    
    if ("Brier Score" %in% metrics) {
      if (is.null(t_star)) {
        metrics_results[[fit_name]]$"Brier Score" <- pam.Brier(fits[[fit_name]], test_data)
      }
      else{
        metrics_results[[fit_name]]$"Brier Score" <- pam.Brier(fits[[fit_name]], test_data, t_star)
      }
      
    }
    
    if ("Time Dependent Auc" %in% metrics) {
      if(!is.null(t_star)){
        pred_time <- t_star
      } else {
        pred_time <- quantile(test_data[[time_var]], 0.5, na.rm = TRUE)
      }
      auc <- pam.survivalROC(Stime = test_data[[time_var]], status = test_data[[status_var]], marker = predict(fits[[fit_name]], newdata = test_data, type = "lp"), predict.time = pred_time, method = "KM")$AUC
      auc <- max(auc, 1- auc)
      metrics_results[[fit_name]]$"Time Dependent Auc" <- auc
    }
  }
  
  # Format the result as a data frame
  metrics_df <- do.call(rbind, lapply(names(metrics_results), function(model) {
    model_metrics <- metrics_results[[model]]
    
    row_data <- c(Model = model, unlist(model_metrics))
    
    as.data.frame(t(row_data), stringsAsFactors = FALSE)
  }))
  metrics_df[-1] <- lapply(metrics_df[-1], function(x) round(as.numeric(x), 2))
  colnames(metrics_df) <- c("Model", metrics)
  return(metrics_df)
}
