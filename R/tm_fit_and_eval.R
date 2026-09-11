#' @title Performance Metrics for Survival Analysis Models
#'
#' @description This function computes a comprehensive set of performance metrics for survival analysis models. It provides metrics such as R_square, L_square, Pseudo_R, Harrell's C, Uno's C, R_sph (distance-based estimator for survival predictive accuracy), R_sh, Brier Score, and Time-dependent AUC. Users can specify particular metrics and model types, enabling tailored performance evaluation for various survival models.
#'
#' @param train_data A data frame containing the survival data.
#' @param covariates A character vector of covariate names to include in the model.
#' @param models A character string or vector specifying the model types to fit (e.g., "coxph", "exp", "lognormal", "weibull"). Default is "coxph" to fit all models.
#' @param metrics A character string or vector specifying the metrics to compute. Default is "all" to compute all available metrics. Options include:
#'   \itemize{
#'     \item "r_square": R-squared metric.
#'     \item "l_square": L-squared metric.
#'     \item "pseudo_r2": Pseudo-R-squared metric.
#'     \item "harrell_c": Harrell's Concordance Index.
#'     \item "uno_c": Uno's Concordance Index.
#'     \item "r_e": Explained variation (R_sph).
#'     \item "r_sh": Explained variation (R_sh).
#'     \item "brier_score": Brier Score.
#'     \item "td_auc": Time-dependent AUC.
#'   }
#'
#' @param predicted_data (Optional) A data frame containing validation data. If `NULL`, the function 
#' uses the same data as `data` for model evaluation.
#' @param t_star (Optional) A positive numeric value specifying the time point at which the Brier score and AUC score is calculated.
#' @param tau (Optional) A time point for truncating the survival time. If provided, the function evaluates predictions up to this time point.
#'
#' @return A data frame containing the selected model's performance metrics.
#'
#' @examples
#' data(pbc, package = "survival")
#' pbc <- pbc[!is.na(pbc$trt), ]
#' pbc$log_albumin <- log(pbc$albumin)
#' pbc$log_bili    <- log(pbc$bili)
#' pbc$log_protime <- log(pbc$protime)
#' pbc$status <- ifelse(pbc$status == 2, 1, 0)
#' covariates <- c("age", "log_albumin", "log_bili", "log_protime", "edema")
#' dat <- pbc[, c("time", "status", covariates)]
#' dat <- dat[stats::complete.cases(dat), ]
#'
#' # All available models and metrics
#' results <- tm_fit_and_eval(train_data = dat, covariates = covariates)
#' results
#'
#' # A specific subset of models and metrics
#' results2 <- tm_fit_and_eval(
#'   train_data = dat,
#'   covariates = covariates,
#'   models  = c("lognormal", "weibull"),
#'   metrics = c("r_square", "l_square", "brier_score")
#' )
#' results2
#'
#' @export

tm_fit_and_eval <- function (train_data, covariates, models = "coxph", 
                                 metrics = "all", predicted_data = NULL, t_star = NULL, tau = NULL) {
  time_var <- "time"
  status_var <- "status"
  # Validate inputs
  if (missing(train_data) || missing(covariates)) {
    stop("Please provide 'train_data', 'time_var', 'status_var', and 'covariates' arguments.")
  }
  if (!is.null(predicted_data)) {
      test_data <- predicted_data
      } else {
      test_data <- train_data
      }
  
  # Fit models based on user input
  fits <- list()
  
  # Ensure the formula is a single string
  formula_text <- paste(
    "Surv(", time_var, ", ", status_var, ") ~ ", 
    paste(covariates, collapse = " + "), 
    sep = ""
  )
  formula <- as.formula(formula_text)
  
  model_types <- if (("all" %in% models)) c("coxph", "exp", "lognormal", "weibull") else models
  metrics <- tm_normalize_metrics(metrics)
  metrics <- if (("all" %in% metrics))c("pseudo_r2", "r_square", "l_square", "harrell_c", "uno_c", "r_e", "r_sh", "brier_score", "td_auc") else metrics
  # Define a list to hold metrics
  metrics_results <- list()
  
  if ("coxph" %in% model_types) {
    fits$coxph <- survival::coxph(formula, data = train_data, x = TRUE, y = TRUE)
  }
  if ("exp" %in% model_types) {
    fits$exp <- survival::survreg(formula, data = train_data, dist = "exponential", x = TRUE, y = TRUE)
  }
  if ("lognormal" %in% model_types) {
    fits$lognormal <- survival::survreg(formula, data = train_data, dist = "lognormal", x = TRUE, y = TRUE)
  }
  if ("weibull" %in% model_types) {
    fits$weibull <- survival::survreg(formula, data = train_data, dist = "weibull", x = TRUE, y = TRUE)
  }
  for (fit_name in names(fits)) {
    metrics_results[[fit_name]] <- list()
    if (is.null(tau)) {
      event_times <- train_data[[time_var]]
      if (length(event_times) == 0) {
        stop("No observed events to determine default tau.")
      }
      tau <- max(event_times)
    }
    if (fit_name == "coxph") {
      r_l_list <- tm_predict_coxph(fits[[fit_name]], covs = covariates,
                                       tau = tau, new_data = test_data,
                                       predict = FALSE) %>%
        Reduce("c", .) %>% as.numeric()
    } 
      else {
      r_l_list <- tm_predict_survreg(fits[[fit_name]], covs = covariates,
                                         tau = tau, new_data = test_data,
                                         predict = FALSE) %>%
        Reduce("c", .) %>% as.numeric()
    }
    # Extract metrics if requested
    if ( "pseudo_r2" %in% metrics ){
      metrics_results[[fit_name]]$Pesudo_R <- round(r_l_list[1] * r_l_list[2], 2)
    }
    if ("r_square" %in% metrics) {
      metrics_results[[fit_name]]$r_square <- round(r_l_list[1], 2)
    }
    if ("l_square" %in% metrics) {
      metrics_results[[fit_name]]$l_square <- round(r_l_list[2], 2)
    }
    
    if ("harrell_c" %in% metrics) {
      metrics_results[[fit_name]]$"harrell_c" <- pam.concordance(fits[[fit_name]], newdata = test_data)$concordance
    }
    
    if ("uno_c" %in% metrics) {
      metrics_results[[fit_name]]$"uno_c" <- pam.concordance(fits[[fit_name]], newdata = test_data, timewt="n/G2")$concordance
    }
    
    if ("r_e" %in% metrics) {
      metrics_results[[fit_name]]$r_e <- pam.rsph(fits[[fit_name]])$Re
    }
    
    
    if ("r_sh" %in% metrics) {
      if (fit_name == "coxph" ) {
        sh_coxph <- survival::coxph(formula, data = train_data,
                                    x = TRUE, y = TRUE)
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
        if (check_factors(train_data) || check_factors(test_data)) {
          metrics_results[[fit_name]]$r_sh <- NA
        } else {
          R_sh_coxph <- pam.schemper(sh_coxph, traindata = train_data, 
                                     newdata = test_data)$Dx
          metrics_results[[fit_name]]$r_sh <- R_sh_coxph 
        }
      } else {
        metrics_results[[fit_name]]$r_sh <- NA
      }
    }
    
    if ("brier_score" %in% metrics) {
      if (is.null(t_star)) {
        metrics_results[[fit_name]]$"brier_score" <- pam.Brier(fits[[fit_name]], test_data)
      }
      else{
        metrics_results[[fit_name]]$"brier_score" <- pam.Brier(fits[[fit_name]], test_data, t_star)
      }
      
    }
    
    if ("td_auc" %in% metrics) {
      if(!is.null(t_star)){
        pred_time <- t_star
      } else {
        pred_time <- quantile(test_data[[time_var]], 0.5, na.rm = TRUE)
      }
      auc <- pam.survivalROC(Stime = test_data[[time_var]], status = test_data[[status_var]], marker = predict(fits[[fit_name]], newdata = test_data, type = "lp"), predict.time = pred_time, method = "KM")$AUC
      auc <- max(auc, 1- auc)
      metrics_results[[fit_name]]$"td_auc" <- auc
    }
  }
  
  # Format the result as a data frame
  metrics_df <- do.call(rbind, lapply(names(metrics_results), function(models) {
    model_metrics <- metrics_results[[models]]
    
    row_data <- c(Model = models, unlist(model_metrics))
    
    as.data.frame(t(row_data), stringsAsFactors = FALSE)
  }))
  metrics_df[-1] <- lapply(metrics_df[-1], function(x) round(as.numeric(x), 2))
  colnames(metrics_df) <- c("Model", metrics)
  return(metrics_df)
}
