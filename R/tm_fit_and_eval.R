#' @title Performance Metrics for Survival Analysis Models
#'
#' @description This function computes a comprehensive set of performance metrics for survival analysis models. It provides metrics such as R_square, L_square, Pseudo_R, Harrell's C, Uno's C, Brier Score, and Time-dependent AUC. Users can specify particular metrics and model types, enabling tailored performance evaluation for various survival models.
#'
#' @param train_data A data frame containing the survival data.
#' @param covariates A character vector of covariate names to include in the model.
#' @param models Character vector of model types: "coxph", "exp", "lognormal",
#'   or "weibull". The default fits only "coxph"; use "all" to fit all four.
#' @param metrics A character string or vector specifying the metrics to compute. Default is "all" to compute all available metrics. Options include:
#'   \itemize{
#'     \item "r_square": R-squared metric.
#'     \item "l_square": L-squared metric.
#'     \item "pseudo_r2": Pseudo-R-squared metric.
#'     \item "harrell_c": Harrell's Concordance Index.
#'     \item "uno_c": Uno's Concordance Index.
#'     \item "brier_score": Brier Score.
#'     \item "td_auc": Time-dependent AUC.
#'   }
#'
#' @param predicted_data (Optional) A data frame containing validation data. If `NULL`, the function 
#' uses the same data as `data` for model evaluation.
#' @param t_star (Optional) A positive numeric value specifying the time point at which the Brier score and AUC score is calculated.
#' @param tau (Optional) A time point for truncating the survival time. If provided, the function evaluates predictions up to this time point.
#'
#' @details Rows with missing values in \code{time}, \code{status}, or a
#'   requested covariate are excluded from training and evaluation with a
#'   warning. All requested metrics use the same complete evaluation rows.
#'   The warning can be silenced with \code{suppressWarnings()}.
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
#' # Default Cox model and all its available metrics
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
  required <- c(time_var, status_var, covariates)
  train_data <- tm_complete_survival_rows(train_data, required, "train_data")
  test_data <- if (is.null(predicted_data)) train_data else
    tm_complete_survival_rows(predicted_data, required, "predicted_data")
  
  # Fit models based on user input
  fits <- list()
  
  # Ensure the formula is a single string
  formula_text <- paste(
    "Surv(", time_var, ", ", status_var, ") ~ ", 
    paste(covariates, collapse = " + "), 
    sep = ""
  )
  formula <- as.formula(formula_text)
  
  valid_models <- c("coxph", "exp", "lognormal", "weibull")
  invalid_models <- setdiff(models, c("all", valid_models))
  if (length(invalid_models) > 0 || length(models) == 0) {
    stop("Invalid models: ", paste(invalid_models, collapse = ", "))
  }
  model_types <- if ("all" %in% models) valid_models else models
  metrics <- tm_normalize_metrics(metrics)
  valid_metrics <- c("pseudo_r2", "r_square", "l_square", "harrell_c",
                     "uno_c", "brier_score", "td_auc")
  invalid_metrics <- setdiff(metrics, c("all", valid_metrics))
  if (length(invalid_metrics) > 0 || length(metrics) == 0) {
    stop("Invalid metrics: ", paste(invalid_metrics, collapse = ", "))
  }
  metrics <- if ("all" %in% metrics) valid_metrics else unique(metrics)
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
      metrics_results[[fit_name]]$pseudo_r2 <- round(r_l_list[1] * r_l_list[2], 2)
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
      marker <- predict(fits[[fit_name]], newdata = test_data, type = "lp")
      # Cox scores increase with risk; AFT location scores increase with time.
      # Do not choose orientation using validation outcomes: an AUC below 0.5
      # is evidence of poor discrimination and must not be folded above 0.5.
      if (fit_name != "coxph") marker <- -marker
      auc <- survivalROC::survivalROC(
        Stime = test_data[[time_var]], status = test_data[[status_var]],
        marker = marker, predict.time = pred_time, method = "KM"
      )$AUC
      metrics_results[[fit_name]]$"td_auc" <- auc
    }
  }
  
  # Format the result as a data frame
  metrics_df <- do.call(rbind, lapply(names(metrics_results), function(models) {
    model_metrics <- metrics_results[[models]]
    
    row_data <- c(Model = models, unlist(model_metrics[metrics]))
    
    as.data.frame(t(row_data), stringsAsFactors = FALSE)
  }))
  metrics_df[-1] <- lapply(metrics_df[-1], function(x) round(as.numeric(x), 2))
  colnames(metrics_df) <- c("Model", metrics)
  return(metrics_df)
}

# Keep model fitting and every metric on the same rows. In particular, do not
# let survivalROC drop missing rows itself: it reports that through cat().
tm_complete_survival_rows <- function(data, columns, label) {
  if (!is.data.frame(data) || !all(columns %in% names(data))) {
    stop("'", label, "' must be a data frame containing: ",
         paste(columns, collapse = ", "), call. = FALSE)
  }
  keep <- stats::complete.cases(data[, columns, drop = FALSE])
  if (!any(keep)) {
    stop("No complete rows remain in '", label, "'.", call. = FALSE)
  }
  if (any(!keep)) {
    warning("Dropped ", sum(!keep), " incomplete row(s) from '", label,
            "'.", call. = FALSE)
  }
  data[keep, , drop = FALSE]
}
