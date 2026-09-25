#' Compute Survival Model Evaluation Metrics
#'
#' This function calculates various predictive performance metrics for survival models based on 
#' predicted survival probabilities and observed survival data. It supports a range of evaluation 
#' measures, including explained variation, concordance indices, Brier Score, and time-dependent AUC.
#' @param model A fitted survival model object. Optional for the metrics
#'   computed directly from \code{predicted_probability}.
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
#'   prediction from \code{model}. If supplied, it should contain the variables
#'   needed by those procedures and columns \code{time} and \code{status}
#'   coded as above.
#'   
#' @param metrics A character vector specifying the evaluation metrics to compute. Options include:
#'   \itemize{
#'     \item "pseudo_r2" - Pseudo R-squared measure
#'     \item "r_square" - Explained variation R^2
#'     \item "l_square" - L-squared measure
#'     \item "Harrell's C" - Harrell's concordance index
#'     \item "Uno's C" - Uno's concordance index
#'     \item "brier_score" - Brier score for calibration
#'     \item "td_auc" - Time-dependent area under the curve (AUC)
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

#' @examples
#' # A Cox model evaluated on right-censored data. The simulator is seeded,
#' # so the reported metric values are reproducible.
#' d <- tm_sim_cox_weibull(n = 150, pi_c = 0.3, v = 2,
#'                         beta = c(0.5, -0.5), seed = 2025)
#' d <- d[, c("time", "status", "x1", "x2")]
#' fit <- survival::coxph(survival::Surv(time, status) ~ x1 + x2,
#'                        data = d, x = TRUE, y = TRUE)
#' pred <- tm_predict_coxph(model = fit, covs = c("x1", "x2"), new_data = d)
#'
#' # the default metric set
#' tm_survival_eval(
#'   model = fit, event_time = pred$times,
#'   predicted_probability = pred$surv_prob, status = pred$status,
#'   covariates = c("x1", "x2"), new_data = d
#' )
#'
#' # a chosen subset; see tm_metric_names() for the accepted identifiers
#' tm_survival_eval(
#'   model = fit, event_time = pred$times,
#'   predicted_probability = pred$surv_prob, status = pred$status,
#'   covariates = c("x1", "x2"), new_data = d,
#'   metrics = c("harrell_c", "brier_score")
#' )
tm_survival_eval <- function (model, event_time, 
                                        predicted_probability, 
                                        pred_mean_survival = NULL,
                                        status, covariates, new_data = NULL,
                                        metrics = NULL,  t_star = NULL, tau = NULL) 
{
  
  #if (missing(covariates) || missing(new_data)) {
  #  stop("Please provide 'new_data',  and 'covariates' arguments.")
  #}
  
  
  
  metrics_results <- list()
  
  valid_metrics <- c("pseudo_r2", "r_square", "l_square",
                     "pseudo_r2_point", "r2_point", "l2_point",
                     "harrell_c", "uno_c", "brier_score", "td_auc")
  default_metrics <- c("pseudo_r2", "pseudo_r2_point",
                       "harrell_c", "uno_c", "brier_score", "td_auc")
  
  metrics <- tm_normalize_metrics(metrics)
  
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
  
  y.order <- order(event_time)
  event_time <- event_time[y.order]
  predicted_probability <- predicted_probability[y.order, ]
  status <- status[y.order]
  
  if(is.null(pred_mean_survival)){
    predicted_data <- integrate_survival(
      predicted_probability, event_time, status, tau)
  }else{
    predicted_data <- pred_mean_survival[y.order]
  }
  
  #lp.pred <- predict(model, newdata = new_data, type="lp")
  #if(inherits(model, "coxph")){
  #  lp.pred <- -predict(model, newdata = new_data, type="lp")
  #}
  
  
  #print(predicted_data)
  if (is.null(t_star)) t_star <- quantile(event_time, 0.5)
  t_idx <- which.min(abs(event_time - t_star))
  risk_scores <- 1 - predicted_probability[, t_idx] 

  if("pseudo_r2" %in% metrics 
     || "r_square" %in% metrics 
     || "l_square" %in% metrics) {
    r_l_list <- pam.r2_metrics(predicted_data, event_time, status, tau)
  }
  
  if ("pseudo_r2" %in% metrics) {
    metrics_results$pseudo_r2 <- round(r_l_list$Pseudo_R_squared, 4)
  } 
  if ("r_square" %in% metrics) {
    metrics_results$r_square <- round(r_l_list$R_squared,4)
  }
  if ("l_square" %in% metrics) {
    metrics_results$l_square <- round(r_l_list$L_square, 4)
  }
  
  if("pseudo_r2_point" %in% metrics ||
     "r2_point" %in% metrics || 
     "l2_point" %in% metrics) {
    i.obs <- ifelse(event_time < t_star & status == 1, 1, 0)
    i.predict <- predicted_probability[, t_idx]
    restricted <- restricted_data_gen(event_time, status, t_star)
    event_time.R2.point <- restricted$time
    status.R2.point <- restricted$status
    r_l_p <- pam.censor.point(i.obs = i.obs, i.predict = i.predict,
                              y = event_time.R2.point, delta = status.R2.point)
  }

  if ("pseudo_r2_point" %in% metrics) {
    metrics_results$pseudo_r2_point <- round(as.numeric(r_l_p$R.squared) * 
                                               as.numeric(r_l_p$L.square), 4)
  } 
  if ("r2_point" %in% metrics) {
    metrics_results$r2_point <- round(as.numeric(r_l_p$R.squared),4)
  }
  if ("l2_point" %in% metrics) {
    metrics_results$l2_point <- round(as.numeric(r_l_p$L.square), 4)
  }
  
  if ("harrell_c" %in% metrics) {
    metrics_results$"harrell_c" <- round(
      concordancefit(y = Surv(event_time, status), 
                     x = predicted_data, 
                     reverse = FALSE)$concordance, 4)
  }
  
  if ("uno_c" %in% metrics) {
    metrics_results$"uno_c" <- round(
      concordancefit(y = Surv(event_time, status), 
                     x = predicted_data, ymax = tau,
                     reverse = FALSE, 
                     timewt = "n/G2")$concordance, 4)
  }
  
  if ("brier_score" %in% metrics) {
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
    metrics_results$"brier_score" <- round(
      as.numeric(brier_result$calibration_res[1]), 4)
  }
  
  if ("td_auc" %in% metrics) {
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
    metrics_results$"td_auc" <- round(AUC_result$main_res$AUC.empirical, 4)
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
#' Calls \code{tm_survival_eval()} for each model in a named list and
#' outputs a wide table: each row corresponds to one metric and each column
#' corresponds to a model.
#'
#' @param models A **named list** in which each element corresponds to a fitted model.  
#'   Each model entry must itself be a list containing:
#'   \itemize{
#'     \item \code{times} -- numeric vector of observed follow-up times.
#'     \item \code{surv_prob} -- an \eqn{n \times K} matrix (or data frame) of
#'           subject-specific predicted survival probabilities on a common time grid.
#'     \item \code{status} -- event indicator (1 = event, 0 = censored).
#'     \item \code{pred} -- (optional) predicted mean survival time (restricted or unrestricted).
#'     \item \code{new_data} -- (optional) dataset used for prediction.
#'     \item \code{covs} -- character vector of covariate names used for prediction.
#'     \item \code{model} -- (optional) the underlying fitted survival model object.
#'   }
#'   Note: If \code{pred} is provided, it will be used to calculate R2 and concordence measure.
#' @param metrics Optional character vector of metrics to compute (passed through).
#' @param t_star Optional numeric scalar, specify the time point to evaluate Brier score and AUC.Default is median of observation time.
#' @param tau Optional numeric scalar, specify the max time horizon for R2 measure and concordence measure (default = 10e10).
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
#' @examples
#' d <- tm_sim_cox_weibull(n = 150, pi_c = 0.3, v = 2,
#'                         beta = c(0.5, -0.5), seed = 2025)
#' d <- d[, c("time", "status", "x1", "x2")]
#' fit <- survival::coxph(survival::Surv(time, status) ~ x1 + x2,
#'                        data = d, x = TRUE, y = TRUE)
#' pred <- tm_predict_coxph(model = fit, covs = c("x1", "x2"), new_data = d)
#'
#' # one column per model, metrics down the rows
#' tm_summarize(list(cox = pred))
tm_summarize <- function(models,
                        metrics = NULL,
                        t_star = NULL,
                        tau = 10e10,
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
    
    res <- tm_survival_eval(
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
      stop(sprintf("Unexpected result structure from tm_survival_eval() for model '%s'.", name))
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
  
  # clean column names like "Value.Cox" ? "Cox"
  names(res_wide) <- sub("^Value\\.", "", names(res_wide))
  rownames(res_wide) <- NULL
  res_wide <- res_wide[, c("Metric", setdiff(names(res_wide), "Metric"))]
  
  # enforce preferred metric display order
  preferred_order <- c(
    "pseudo_r2",
    "r_square",
    "l_square",
    "pseudo_r2_point", 
    "r2_point", 
    "l2_point",
    "harrell_c",
    "uno_c",
    "brier_score",
    "td_auc"
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

#' @keywords internal
#' @noRd
integrate_survival <- function(predicted_probability, event_time, status, tau = NULL) {
  if (nrow(predicted_probability) != length(event_time)) {
    stop("Number of rows in predicted_probability must equal length of event_time")
  }
  
  if (any(predicted_probability < 0 | predicted_probability > 1, na.rm = TRUE)) {
    stop("All predicted probabilities must be between 0 and 1")
  }
  if (!is.null(tau)) {
    restricted <- restricted_data_gen(event_time, status, tau)
    event_time <- restricted$time
    status <- restricted$status
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
