#' @title Summary of Performance Metrics for Survival Model
#' 
#' @description This function calculates various performance metrics for survival models, such as Cox proportional hazards models (`coxph`) 
#' and parametric survival models (`survreg`). It supports metrics like pseudo R-squared, Harrell's C, Uno's C, Brier Score, 
#' and time-dependent AUC. The function is designed to work with both training and test datasets, allowing users to evaluate 
#' model performance on new data.
#'
#' @param object A fitted survival model object. Supported classes include `coxph` (from the `survival` package) 
#'               and `survreg` (from the `survival` package).
#' @param train_data A data frame containing the training data used to fit the model. This is required for certain metrics 
#'                   like `R_sh`.
#' @param predicted_data A data frame containing the test data for which predictions are to be evaluated. This data must 
#'                       include the time-to-event and event status columns.
#' @param covariates A character vector specifying the names of the covariates used in the model.
#' @param metrics A character string or vector specifying the metrics to compute. Default is "all" to compute all available metrics. Options include:
#'   \itemize{
#'     \item "R_square": R-squared metric.
#'     \item "L_square": L-squared metric.
#'     \item "Pesudo_R": Pseudo-R-squared metric.
#'     \item "Harrells_C": Harrell's Concordance Index.
#'     \item "Unos_C": Uno's Concordance Index.
#'     \item "R_sph": Explained variation (R_sph).
#'     \item "R_sh": Explained variation (R_sh).
#'     \item "Brier_Score": Brier Score.
#'     \item "Time_Dependent_Auc": Time-dependent AUC.
#'   }
#'
#' @param t_star A numeric value specifying the time point at which to evaluate time-dependent metrics (e.g., Brier Score 
#'               and AUC). If `NULL`, the median survival time of the test data is used.
#' @param tau A numeric value specifying the time horizon for restricted mean survival time (RMST) calculations. If `NULL`, 
#'            the maximum observed time in the test data is used.
#'
#' @return A data frame with two columns:
#' \itemize{
#'   \item `Metric`: The name of the performance metric.
#'   \item `Value`: The calculated value of the metric, rounded to two decimal places.
#' }
#'
#' @references
#' 
#' Li, G., & Wang, X. (2016). Prediction Accuracy Measures for a Nonlinear Model and for Right-Censored Time-to-Event Data. arXiv preprint arXiv:1611.03063. Available at https://arxiv.org/abs/1611.03063
#' F Harrell, R Califf, D Pryor, K Lee and R Rosati, Evaluating the yield of medical tests, J Am Medical Assoc, 1982.
#' 
#' R Peto and J Peto, Asymptotically efficient rank invariant test procedures (with discussion), J Royal Stat Soc A, 1972.
#' 
#' M Schemper, Cox analysis of survival data with non-proportional hazard functions, The Statistician, 1992.
#' 
#' H Uno, T Cai, M Pencina, R D'Agnostino and Lj Wei, On the C-statistics for evaluating overall adequacy of risk prediction procedures with censored survival data, Statistics in Medicine, 2011.
#' 
#' Therneau, T. M., Lumley, T., Atkinson, E., Crowson, C. (2024). survival: Survival Analysis. 
#' R package version 3.7-0. DOI: \doi{10.32614/CRAN.package.survival}. Available at \url{https://CRAN.R-project.org/package=survival}.
#' 
#' Schemper, M. and R. Henderson (2000). Predictive accuracy and explained variation in Cox regression. Biometrics 56, 249--255.
#' 
#' LUSA, L., MICELI, R. and MARIANI, L. (2007). Estimation of predictive accuracy in survival analysis using R and S-PLUS. Computer methods and programs in biomedicine 87 132–137.
#' Graf, E., Schmoor, C., Sauerbrei, W., & et al. (1999). Assessment and comparison of prognostic classification schemes for survival data. *Statistical Medicine*, 18(17-18), 2529-2545.
#' Schemper, M. and R. Henderson (2000). Predictive accuracy and explained variation in Cox regression.
#' Biometrics 56, 249--255.
#' 
#' Lusa, L., R. Miceli and L. Mariani (2007). Estimation of predictive accuracy in survival analysis
#' using R and S-PLUS. Computer Methods and Programs in Biomedicine 87, 132--137.
#' 
#' Potapov, S., Adler, W., Schmid, M., Bertrand, F. (2024). survAUC: Estimating Time-Dependent AUC for Censored Survival Data. 
#' R package version 1.3-0. DOI: \doi{10.32614/CRAN.package.survAUC}. Available at \url{https://CRAN.R-project.org/package=survAUC}.
#'
#' Brier, G. W. (1950). Verification of forecasts expressed in terms of probability. *Monthly Weather Review*, 78.
#'
#' Gneiting, T., & Raftery, A. E. (2007). Strictly Proper Scoring Rules, Prediction, and Estimation.
#' 
#' Zhou, H., Cheng, X., Wang, S., Zou, Y., & Wang, H. (2022). SurvMetrics: Predictive Evaluation Metrics in Survival Analysis. 
#' R package version 0.5.0. Available at \url{https://github.com/skyee1/SurvMetrics}.
#' 
#' Heagerty, P. J., Lumley, T., & Pepe, M. S. (2000). Time-dependent ROC curves for censored
#' survival data and a diagnostic marker. \emph{Biometrics}, 56(2), 337-344.
#'
#' Heagerty, P. J., Saha-Chaudhuri, P. (2022). survivalROC: Time-Dependent ROC Curve Estimation from Censored Survival Data. 
#' R package version 1.0.3.1. DOI: \doi{10.32614/CRAN.package.survivalROC}. Available at \url{https://CRAN.R-project.org/package=survivalROC}.
#' 
#' @examples
#' library(PAmeasures)
#' library(survival)
#'
#' # Use Mayo Clinic Primary Biliary Cirrhosis Data
#' data(pbc)
#' pbc <- pbc %>% 
#'   filter(is.na(trt)==F) %>% 
#' mutate(log_albumin = log(albumin),
#'        log_bili = log(bili),
#'        log_protime = log(protime),
#'        status = ifelse(status==2, 1, 0))
#' #Schemper and Henderson's estimator of the absolute deviation between survival functions
#' schemper(train.fit.full, pbc, pbc)$Dx       
#' @keywords internal
#' @noRd
pam.prediction_survial_eval <- function (object, train_data, predicted_data, covariates,
                               metrics = "all",  t_star = NULL, tau = NULL) 
{
  time_var <- "time"
  status_var <- "status"
  if (missing(predicted_data) || missing(covariates)) {
    stop("Please provide 'predicted_data',  and 'covariates' arguments.")
  }
  
  test_data <- predicted_data
  metrics_results <- list()
  
  valid_metrics <- c("Pseudo_R_square", "R_square", "L_square", "Harrell’s C", "Uno’s C",
                     "R_sph", "R_sh", "Brier Score", "Time Dependent Auc")
  
  if ("all" %in% metrics) {
    metrics <- valid_metrics
  } else {
    invalid <- setdiff(metrics, valid_metrics)
    if (length(invalid) > 0) stop("Invalid metrics: ", paste(invalid, collapse = ", "))
  }
  
  if (inherits(object, "coxph")){
    r_l_list <- pam.coxph_restricted(object, 
                                     covariates = covariates, tau = tau, predicted_data = test_data) %>% 
      Reduce("c", .) %>% as.numeric()
  }
  else if (inherits(object, c("survreg"))) {
    r_l_list <- pam.surverg_restricted(object, 
                                       covariates = covariates, tau = tau, predicted_data = test_data) %>% 
      Reduce("c", .) %>% as.numeric()
  }
  
  if ("Pseudo_R_square" %in% metrics) {
    metrics_results$Pesudo_R_square <- round(r_l_list[1] * 
                                        r_l_list[2], 2)
  } 
  if ("R_square" %in% metrics) {
    metrics_results$R_square <- round(r_l_list[1], 
                                      2)
  }
  if ("L_square" %in% metrics) {
    metrics_results$L_square <- round(r_l_list[2], 
                                      2)
  }

  if ("Harrell’s C" %in% metrics) {
    metrics_results$"Harrell’s C" <- pam.concordance(object, 
                                                     newdata = test_data)$concordance
  }
  
  if ("Uno’s C" %in% metrics) {
    metrics_results$"Uno’s C" <- pam.concordance(object, 
                                                 newdata = test_data, timewt = "n/G2")$concordance
  }
  
  if ("R_sph" %in% metrics) {
    metrics_results$R_sph <- pam.rsph(object)$Re
  }
  
  if ("R_sh" %in% metrics) {
    if (inherits(object, "coxph")){
      check_factors <- function(data) {
        factors <- sapply(data, is.factor)
        if (any(factors)) {
          factor_cols <- names(data)[factors]
          warning("The following columns are factors: ", 
                  paste(factor_cols, collapse = ", "), ". This function to calculate R_sph does not support factor variables. Returning NA.")
          return(TRUE)
        }
        return(FALSE)
      }
      if (check_factors(predicted_data)) {
        metrics_results$R_sh <- NA
      }
      else {
        formula_text <- paste(
          "Surv(", time_var, ", ", status_var, ") ~ ", 
          paste(covariates, collapse = " + "), 
          sep = ""
        )
        # Convert to formula
        formula <- as.formula(formula_text)
        rms_coxph <- rms::cph(formula, data = train_data, x = TRUE, y = TRUE)
        R_sh_coxph <- pam.schemper(rms_coxph, traindata = train_data, 
                                   newdata = test_data)$Dx
        metrics_results$R_sh <- R_sh_coxph
      }
    }
    else {
      metrics_results$R_sh <- NA
    }
  }
  if ("Brier Score" %in% metrics) {
    if (is.null(t_star)) {
      metrics_results$"Brier Score" <- pam.Brier(object, 
                                                 test_data)
    }
    else {
      metrics_results$"Brier Score" <- pam.Brier(object, 
                                                 test_data, t_star)
    }
  }
  if ("Time Dependent Auc" %in% metrics) {
    if (!is.null(t_star)) {
      pred_time <- t_star
    }
    else {
      pred_time <- quantile(test_data[[time_var]], 
                            0.5, na.rm = TRUE)
    }
    auc <- pam.survivalROC(Stime = test_data[[time_var]], 
                           status = test_data[[status_var]], marker = predict(object, 
                                                                              newdata = test_data, type = "lp"), predict.time = pred_time, 
                           method = "KM")$AUC
    auc <- max(auc, 1 - auc)
    metrics_results$"Time Dependent Auc" <- auc
  }
  
  result_df <- data.frame(
    Metric = names(metrics_results),
    Value = unlist(metrics_results, use.names = FALSE),
    stringsAsFactors = FALSE
  )
  result_df$Value <- round(result_df$Value, 2)
  
  return(result_df)
  
}

