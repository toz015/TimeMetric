#' Compute Survival Model Evaluation Metrics
#'
#' This function calculates various predictive performance metrics for survival models based on 
#' predicted survival probabilities and observed survival data. It supports a range of evaluation 
#' measures, including explained variation, concordance indices, Brier Score, and time-dependent AUC.
#' @param model A fitted survival model object used by certain metrics:
#'   \itemize{
#'     \item For \code{"R_sh"} (Schemper–Henderson), a Cox model fitted with \code{x=TRUE, y=TRUE}
#'           (e.g., \code{survival::coxph}) or an \code{rms::cph} model.
#'     \item For \code{"R_E"} (rank-based \(R^2\)), a Cox model compatible with \code{pam.rsph()}.
#'   }
#'
#' @param event_time A numeric vector of observed survival times.
#' @param pred_mean_survival A numeric vector of predicted mean survival. If input this value, \code{predicted_probability} value would be ignored.
#' @param predicted_probability A numeric vector of predicted survival probabilities (Note: Only models with discrete estimated survival probabilities can use this input.). with
#'   \itemize{
#'     \item rows = subjects.
#'     \item columns = observed survival times
#'   }
#' @param status A numeric vector indicating event occurrence (1 for event, 0 for censoring).
#' @param covariates A character vector specifying the names of the covariates used in the model.
#' @param new_data Optional data frame used by metrics that require refitting or
#'   prediction from \code{model} (e.g., Schemper–Henderson \code{"R_sh"} and
#'   rank-based \code{"R_E"}). If supplied, it should contain the variables
#'   needed by those procedures and columns \code{time} and \code{status}
#'   coded as above.
#'   
#' @param metrics A character vector specifying the evaluation metrics to compute. Options include:
#'   \itemize{
#'     \item "Pseudo_R_square" - Pseudo R-squared measure
#'     \item "R_square" - Explained variation R²
#'     \item "L_square" - L-squared measure
#'     \item "Harrell’s C" - Harrell’s concordance index
#'     \item "Uno's C" - Uno’s concordance index
#'     \item "R_sh" - Schemper-Henderson explained variation (R_sh)
#'     \item "R_E" - Rank-based R²
#'     \item "Brier Score" - Brier score for calibration
#'     \item "Time Dependent Auc" - Time-dependent area under the curve (AUC)
#'   }
#'   Default is "all", which computes all available metrics.
#'   
#'   
#' @param t_star (Optional) A numeric value specifying the evaluation time for Brier Score and time-dependent AUC. 
#'   If NULL, it defaults to the median survival time.
#' @param tau (Optional) A numeric value specifying the truncation time for calculating explained variation metrics. 
#'   If NULL, it defaults to the maximum observed survival time.
#'
#' @return A data frame with two columns:
#'   \item{Metric}{The name of the computed evaluation metric.}
#'   \item{Value}{The corresponding computed value.}
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
#' @export

pam.predicted_survial_eval <- function (model, event_time, predicted_probability, 
                                        pred_mean_survival = NULL,
                                        status, covariates, new_data = NULL,
                                        metrics = NULL,  t_star = NULL, tau = NULL) 
{
  
  if (missing(covariates) || missing(new_data)) {
    stop("Please provide 'new_data',  and 'covariates' arguments.")
  }
  
  
  
  metrics_results <- list()
  
  valid_metrics <- c("Pseudo_R_square", "R_square", "L_square", 
                     "Harrell’s C", "Uno’s C",
                      "R_E","R_sh", "Brier Score", "Time Dependent Auc")
  default_metrics <- c("Pseudo_R_square", "Harrell’s C", "Uno’s C",
                       "R_E","R_sh", "Brier Score", "Time Dependent Auc")
  if (is.null(metrics)){
    metrics <- default_metrics
  }
  else if ("all" %in% metrics) {
    metrics <- valid_metrics
  } else {
    invalid <- setdiff(metrics, valid_metrics)
    if (length(invalid) > 0) 
      stop("Invalid metrics: ", paste(invalid, collapse = ", "))
  }
  
  
  if(is.null(pred_mean_survival)){
    predicted_data <- integrate_survival(
      predicted_probability, event_time, status, tau)
  }else{
    predicted_data <- pred_mean_survival
  }
  
  #print(predicted_data)
  if (is.null(t_star)) t_star <- quantile(event_time, 0.5)
  t_idx <- which.min(abs(event_time - t_star))
  risk_scores <- 1 - predicted_probability[, t_idx] 

  if("Pseudo_R_square" %in% metrics 
     || "R_square" %in% metrics 
     || "L_square" %in% metrics) {
    r_l_list <- pam.r2_metrics(predicted_data, event_time, status, tau)
  }
  
  if ("Pseudo_R_square" %in% metrics) {
    metrics_results$Pseudo_R_square <- round(r_l_list$Pseudo_R_squared, 4)
  } 
  if ("R_square" %in% metrics) {
    metrics_results$R_square <- round(r_l_list$R_squared,4)
  }
  if ("L_square" %in% metrics) {
    metrics_results$L_square <- round(r_l_list$L_square, 4)
  }
  
  if ("Harrell’s C" %in% metrics) {
    metrics_results$"Harrell’s C" <- round(
      concordancefit(y = Surv(event_time, status), 
                     x = predicted_data, reverse = FALSE)$concordance, 4)
  }
  
  if ("Uno’s C" %in% metrics) {
    metrics_results$"Uno’s C" <- round(
      concordancefit(y = Surv(event_time, status), 
                     x = predicted_data, reverse = FALSE, 
                     timewt = "n/G2")$concordance, 4)
  }
  
  if ("R_sh" %in% metrics) {
    if (!inherits(model, "coxph")) {
      message("R_sh is only defined for Cox PH models (coxph). Returning NA.")
      metrics_results$R_sh <- NA
    } else {
      message("Note: 'R_sh' metric is only valid for Cox model.")
      
      check_factors <- function(data) {
        factors <- sapply(data, is.factor)
        if (any(factors)) {
          factor_cols <- names(data)[factors]
          message(
            "The following columns are factors: ",
            paste(factor_cols, collapse = ", "),
            ". This function to calculate R_sh does not support factor variables. Returning NA."
          )
          return(TRUE)
        }
        return(FALSE)
      }
      if (check_factors(new_data)) {
        metrics_results$"R_sh" <- NA
      } else {
        train_data <- data.frame(
          time   = model$y[, 1],
          status = model$y[, 2],
          as.data.frame(model$x, check.names = FALSE)
        )
        
        formula_text <- paste(
          "Surv(", "time", ", ", "status", ") ~ ", 
          paste(covariates, collapse = " + "), 
          sep = ""
        )
    
        formula <- as.formula(formula_text)
        rms_coxph <- rms::cph(formula, data = train_data, x = TRUE, y = TRUE)
        metrics_results$"R_sh" <- pam.schemper(
          train.fit = rms_coxph,
          traindata = train_data,
          newdata   = new_data
        )$Dx
      }
    }
  }
  if ("R_E" %in% metrics) {
   metrics_results$"R_E" <- pam.rsph(model, test_data = new_data)$Re
  }
  
  if ("Brier Score" %in% metrics) {
    t_eval <- event_time[t_idx]
    X <- risk_scores
    brier_result <- suppressMessages(tdROC::tdROC(
      X = X,  
      Y = event_time,      
      delta = status,
      tau = t_eval,       
      method = "both", 
      output = "both"   
    ))
    metrics_results$"Brier Score" <- round(
      as.numeric(brier_result$calibration_res[1]), 4)
  }
  
  if ("Time Dependent Auc" %in% metrics) {
    t_eval <- event_time[t_idx]
    X <- risk_scores
    AUC_result <- suppressMessages(tdROC::tdROC(
      X = X,  
      Y = event_time,      
      delta = status,
      tau = t_eval,       
      method = "both", 
      output = "both"   
    ))
    metrics_results$"Time Dependent Auc" <- round(AUC_result$main_res$AUC.empirical, 4)
  }
  result_df <- data.frame(
    Metric = names(metrics_results),
    Value = unlist(metrics_results, use.names = FALSE),
    stringsAsFactors = FALSE
  )
  
  return(result_df)
}

#' Summarize multiple survival models into a wide comparison table
#'
#' @description
#' Calls \code{pam.predicted_survial_eval()} for each model in a named list and
#' outputs a wide table: each row corresponds to one metric and each column
#' corresponds to a model.
#'
#' @param models A **named list** where each element is a list with:
#' \itemize{
#'   \item \code{times} — numeric vector of observed times.
#'   \item \code{surv_prob} — \eqn{n \times K} matrix (or data frame) of
#'         subject-specific survival probabilities on a common time grid.
#'   \item \code{status} — 0/1 event indicator (1 = event, 0 = censored).
#' }
#' @param metrics Optional character vector of metrics to compute (passed through).
#' @param t_star Optional numeric scalar, passed to \code{pam.predicted_survial_eval()}.
#' @param tau Optional numeric scalar, passed to \code{pam.predicted_survial_eval()}.
#' @param digits Number of decimal places to round the metric values (default = 2).
#'
#' @return A data frame with:
#'   \itemize{
#'     \item Each row = metric
#'     \item Each column = model name
#'     \item Cell values = corresponding rounded metric values
#'   }
#'
#' @export
pam.summary <- function(models,
                        metrics = NULL,
                        t_star = NULL,
                        tau = NULL,
                        digits = 2) {
  if (!is.list(models) || length(models) == 0) {
    stop("'models' must be a non-empty named list.")
  }
  if (is.null(names(models)) || any(names(models) == "")) {
    names(models) <- paste0("Model_", seq_along(models))
  }
  
  required <- c("times", "surv_prob", "status")
  
  eval_one <- function(mod, name) {
    if (!all(required %in% names(mod))) {
      stop(sprintf("Model '%s' must contain: %s",
                   name, paste(required, collapse = ", ")))
    }
    
    res <- pam.predicted_survial_eval(
      model = mod$model,
      event_time = mod$times,
      predicted_probability = mod$surv_prob,
      pred_mean_survival = mod$pred,
      status = mod$status,
      metrics = metrics,
      new_data = mod$new_data,
      covariates = mod$covs,
      t_star = t_star,
      tau = tau
    )
    
    # Ensure expected columns exist and add Model label
    if (!all(c("Metric", "Value") %in% names(res))) {
      stop(sprintf("Unexpected result structure from pam.predicted_survial_eval() for model '%s'.", name))
    }
    res$Model <- name
    res[, c("Model", "Metric", "Value")]
  }
  
  # run single-model evaluation for all
  res_list <- mapply(eval_one, models, names(models), SIMPLIFY = FALSE)
  res_long <- do.call(rbind, res_list)
  rownames(res_long) <- NULL
  
  # pivot wider: rows = Metric, cols = Model, cells = Value
  res_wide <- reshape(
    res_long,
    idvar = "Metric",
    timevar = "Model",
    direction = "wide"
  )
  
  # clean column names like "Value.Cox" → "Cox"
  names(res_wide) <- sub("^Value\\.", "", names(res_wide))
  rownames(res_wide) <- NULL
  res_wide <- res_wide[, c("Metric", setdiff(names(res_wide), "Metric"))]
  
  # enforce preferred metric display order
  preferred_order <- c(
    "Pseudo_R_square",
    "R_square",
    "L_square",
    "Harrell’s C",
    "Uno’s C",
    "R_sh",
    "R_E",
    "Brier Score",
    "Time Dependent Auc"
  )
  
  res_wide$Metric <- factor(res_wide$Metric, levels = preferred_order)
  res_wide <- res_wide[order(res_wide$Metric), ]
  res_wide$Metric <- as.character(res_wide$Metric)
  
  # round numeric columns
  numeric_cols <- setdiff(names(res_wide), "Metric")
  res_wide[numeric_cols] <- lapply(res_wide[numeric_cols], function(x) round(x, digits))
  
  return(res_wide)
}

#' Integrate Predicted Survival Probabilities Over Time
#'
#' @description
#' Computes an integrated survival (or risk) measure by summing predicted probabilities 
#' across observed follow-up times, weighted by the time increments between events.
#' This function is intended for internal use within the package.
#'
#' @param predicted_probability A numeric matrix or data frame of predicted survival
#' probabilities at observed times, with one row per subject (matching the order of \code{event_time})
#' and one or more columns for different prediction models or time points.
#' @param event_time A numeric vector of observed event or censoring times.
#' @param status A binary vector indicating event occurrence 
#' (1 = event, 0 = censored) for each subject.
#' @param tau Optional numeric value specifying the maximum truncation time. 
#' If \code{NULL} (default), the maximum observed event time is used.
#'
#' @return
#' A numeric vector of cumulative integrated predictions, with one element per
#' column in \code{predicted_probability}.
#'
#' @keywords internal

integrate_survival <- function(predicted_probability, event_time, status, tau = NULL) {
  if (nrow(predicted_probability) != length(event_time)) {
    stop("Number of rows in predicted_probability must equal length of event_time")
  }
  
  if (any(predicted_probability < 0 | predicted_probability > 1, na.rm = TRUE)) {
    stop("All predicted probabilities must be between 0 and 1")
  }
  
  if (is.null(tau)) {
    tau <- max(event_time, na.rm = TRUE)
  }
  
  order_idx <- order(event_time)
  event_time <- event_time[order_idx]
  predicted_probability <- predicted_probability[order_idx, , drop = FALSE]
  event_time <- pmin(event_time, tau)
  t1 <- event_time
  t2 <- c(0, head(t1, -1)) 
  delta <- status[order_idx]
  delta <- ifelse(event_time <= tau, delta, 0)
  delta.t <- t1 - t2       
  
  cumulative_prediction <- colSums(delta.t %*% t(predicted_probability))
  #print(order_idx)
  return(cumulative_prediction)
}