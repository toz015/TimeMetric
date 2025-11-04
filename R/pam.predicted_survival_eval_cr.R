#' Evaluate Survival Prediction Performance with Competing Risks
#'
#' This function evaluates the predictive performance of a survival model in the presence of competing risks.
#' It calculates various performance metrics including Pseudo R-squared, R-squared, L-squared, C-index, 
#' Brier Score, and Time-Dependent AUC.
#'
#' @param pred_cif A matrix of predicted cumulative incidence functions (CIFs) for different subjects.
#' @param event_time A numeric vector of observed event or censoring times for subjects.
#' @param time.cif A numeric vector of time points corresponding to the predicted CIFs.
#' @param status A numeric vector indicating the event type for each subject: by default
#'   \itemize{
#'     \item `0` for censored cases.
#'     \item `1` for the primary event of interest.
#'     \item `2` (or other values) for competing events.
#'   }
#' @param metrics A character vector specifying the performance metrics to calculate. 
#'   The default is `"all"`, which computes all available metrics. 
#'   Valid options include:
#'   \itemize{
#'     \item `"Pseudo_R_square"`
#'     \item `"R_square"`
#'     \item `"L_square"`
#'     \item `"Pseudo_R2_point"`
#'     \item `"R2_point"`
#'     \item `"L2_point"`
#'     \item `"C_index"`
#'     \item `"Brier Score"`
#'     \item `"Time Dependent Auc"`
#'   }
#' @param t_star Optional numeric value specifying a reference time for evaluating point-version Pseudo_R_square, 
#'   the Brier Score and Time-Dependent AUC. Defaults to the median of `event_time`.
#' @param tau Optional numeric value specifying a restricted time horizon for model evaluation. 
#'   Defaults to the maximum observed event time plus one.
#' @param event_type An integer indicating the primary event type of interest (default = `1`).
#'
#' @return A data frame with two columns:
#'   \itemize{
#'     \item `"Metric"`: Name of the evaluated metric.
#'     \item `"Value"`: Computed value of the corresponding metric.
#'   }
#'
#' 
#' 
#' @export
pam.predicted_survial_eval_cr <- function (pred_cif, event_time, time.cif, status,
                                           metrics = NULL,  t_star = NULL, 
                                           tau = NULL, event_type = 1) 
  
{
  
  if (missing(event_time) || missing(pred_cif)) {
    stop("Please provide 'event_time',  and 'pred_cif' arguments.")
  }
  
  if(is.null(tau)) tau <- max(event_time)
  if(is.null(t_star)) t_star <- median(event_time)
  time_idx <- max(which(time.cif <= t_star))
  
  metrics_results <- list()
  
  valid_metrics <- c("Pseudo_R_square", "R_square", "L_square", 
                     "Pseudo_R2_point", "R2_point", "L2_point",
                     "C_index", "Brier Score", "Time Dependent Auc")
  default_metrics <- c("Pseudo_R_square", "Pseudo_R2_point", "C_index",
                       "Brier Score", "Time Dependent Auc")
  if ("all" %in% metrics) {
    metrics <- valid_metrics
  } else {
    invalid <- setdiff(metrics, valid_metrics)
    if (length(invalid) > 0) 
      stop("Invalid metrics: ", paste(invalid, collapse = ", "))
  }
  


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
  
  
  if("Pseudo_R_square" %in% metrics ||
     "R_square" %in% metrics || 
     "L_square" %in% metrics) {
    r_l_list <- pam.censor.cr(event_time, status, tau, 
                              pred_cif, time.cif, event_type)
  }
  
  if ("Pseudo_R_square" %in% metrics) {
    metrics_results$Pseudo_R_square <- r_l_list$Pseudo.R
  } 
  if ("R_square" %in% metrics) {
    metrics_results$R_square <- r_l_list$R.squared
  }
  if ("L_square" %in% metrics) {
    metrics_results$L_square <- r_l_list$L.squared
  }
  
  if("Pseudo_R2_point" %in% metrics ||
     "R2_point" %in% metrics || 
     "L2_point" %in% metrics) {
    
    r_p_list <- pam.censor.cr.point(ftime = event_time, 
                                    fstatus = status,
                                    tau = t_star, pred.cif = pred_cif[time_idx, ],
                                    event.type = event_type)
  }
  
  if ("Pseudo_R2_point" %in% metrics) {
    metrics_results$Pseudo_R2_point <- r_p_list$Pseudo.R
  } 
  if ("R2_point" %in% metrics) {
    metrics_results$R2_point <- r_p_list$R.squared
  }
  if ("L2_point" %in% metrics) {
    metrics_results$L2_point <- r_p_list$L.squared
  }
  
  pred_risk <- apply(pred_cif, 2, m_cif, time.cif=time.cif, tau=tau)
  
  C_index <- C_cr(
    time = event_time,
    status = status,
    predicted = pred_risk,
    tau = tau,
    Cause_int = event_type
  )
  if ("C_index" %in% metrics) {
    metrics_results$"C_index" <- round(C_index, 4)
  }
  
  status.recode <- ifelse(status == event_type, -1, status)
  status.recode <- ifelse(status.recode > 0, 2, status.recode)
  status.recode <- ifelse(status.recode == -1, 1, status.recode)
  
  if ("Brier Score" %in% metrics) {
    X <- pred_cif[time_idx, ]
    
    brier_result <- suppressMessages(tdROC::tdROC.cr(
      X = X,  
      Y = event_time,      
      delta = status.recode,
      tau = t_star, 
      nboot = 0))
    
    metrics_results$"Brier Score" <- round(
      as.numeric(brier_result$calibration_res[1]), 4)
  }
  
  
  if ("Time Dependent Auc" %in% metrics) {
    AUC_result <- suppressMessages(tdROC::tdROC.cr(
      X = X,  
      Y = event_time,      
      delta = status.recode,
      tau = t_star,       
      method = "both", 
      output = "AUC"   
    ))
    metrics_results$"Time Dependent Auc" <- round(AUC_result$main_res$AUC.A.integral, 4)
    #metrics_results$"Time Dependent Auc Empirical" <- round(AUC_result$main_res$AUC.B.integral, 4)
  }
  
  result_df <- data.frame(
    Metric = names(metrics_results),
    Value = unlist(metrics_results),
    row.names = NULL,  
    stringsAsFactors = FALSE
  )
  
  return(result_df)
}


##################################################
## calculate integral for int_0^tau t dF(t|x) 
##################################################
### input:
# x: input estimated CIF for subject i
# time.cif: time point of observation event 1
# tau: restricted time
#' @noRd
m_cif <- function(x, time.cif, tau){
  index <- which(time.cif <= tau)
  t(diff(c(0, time.cif[index], tau))) %*% c(1, (1-x)[index])
}




##################################################
### calculate PA measure for competing risk model
##################################################
### input:
# ftime: vector of failure/censoring times
# fstatus: vector with a unique code for each failure type (1,2,3...)
#          and censored observations (0)
# tau: restricted time
# time.cif: time point of observation event of interest
# pred.cif: estimated CIF for all subject
#' @noRd
pam.censor.cr <- function(ftime, fstatus, tau, 
                          pred.cif, time.cif,
                          event.type = 1){
  

  ftime.new <- ifelse(ftime > tau, tau, ftime)
  ftime.new <- ifelse(ftime.new <= tau & !(fstatus %in% c(event.type, 0)), 
                      tau, ftime.new)
  fstatus.new <- ifelse(ftime.new == tau | fstatus > 0, 1, fstatus)
  

  pred <- apply(pred.cif, 2, m_cif, time.cif = time.cif, tau = tau)
  

  out <- pam.censor(ftime.new, pred, fstatus.new)
  
  out <- c(out, 
           Pseudo.R = format(round(as.numeric(out$R.squared) * 
                                     as.numeric(out$L.squared), 
                                   digits = 4), nsmall = 4))
  return(out)
}


####################################################################
### calculate competing risk PA measure for a given time point
####################################################################
#' @noRd
pam.censor.cr.point <- function(ftime, fstatus, tau, 
                                pred.cif,
                                event.type = 1){
  
  # create new variable T(1,tau)
  ftime.new <- ifelse(ftime > tau, tau, ftime)
  ftime.new <- ifelse(ftime.new <= tau & !(fstatus %in% c(event.type, 0)), 
                      tau, ftime.new)
  fstatus.new <- ifelse(ftime.new == tau | fstatus > 0, 1, fstatus)
  
  
  # obtain PA measure
  i.obs <- ifelse(ftime < tau & fstatus == event.type, 1, 0)
  out <- pam.censor.point(i.obs = i.obs, i.predict = pred.cif,
                          y =  ftime.new, delta = fstatus.new)
  out <- c(out, 
           Pseudo.R = format(round(as.numeric(out$R.squared) * 
                                     as.numeric(out$L.squared), 
                                   digits = 4), nsmall = 4))
  return(out)
}


####################################################################
### calculate PA measure for a given time point
####################################################################
#' @noRd
pam.censor.point <- function(i.obs, i.predict, y, delta){
  
  y.sorted<-sort(y)
  delta.sorted<-delta[order(y)]
  
  #KM estimate for censoring distribution
  fit.km.censoring <- survfit(Surv(y.sorted, 1-delta.sorted) ~ 1  )
  #sum.km.censoring<-summary(fit.km.censoring,  censored=TRUE)
  sum.km.censoring<-summary(fit.km.censoring,times=y.sorted,extend=TRUE)
  km.censoring<-sum.km.censoring$surv
  
  km.censoring.minus<-c(1,km.censoring[-length(km.censoring)])
  ratio.km<- delta.sorted/km.censoring.minus
  
  ratio.km[is.nan(ratio.km)]<-0
  
  weight.km<-ratio.km/(sum(ratio.km))
  
  ### now change y (observe event time) to indicator version
  
  y.sorted<-i.obs[order(y)]
  y.predict.sorted<-i.predict[order(y)]
  
  wls.fitted<-tryCatch(lm(y.sorted~y.predict.sorted,weights=weight.km),error=function(e){return(c(NA,NA))})
  calibrate.fitted<-tryCatch(predict(wls.fitted),error=function(e){return(c(NA,NA))})
  
  num.rho2<-sum(weight.km*(calibrate.fitted-sum(weight.km*y.sorted))^2)
  denom.rho2<-sum(weight.km*(y.sorted-sum(weight.km*y.sorted))^2)
  
  R2 <-format(round(num.rho2/denom.rho2,digits = 4) ,nsmall=4)
  
  
  num.L2<- sum(weight.km*(y.sorted-calibrate.fitted)^2)
  denom.L2<- sum(weight.km*(y.sorted-y.predict.sorted)^2)
  L2 <-format(round(num.L2/denom.L2,digits = 4),nsmall=4)
  
  return(list(R.squared=R2,L.squared=L2))
}


#' @noRd
pam.censor<-function(y,y.predict,delta){
  
  
  
  y.sorted<-sort(y)
  delta.sorted<-delta[order(y)]
  y.predict.sorted<-y.predict[order(y)]
  
  fit.km.censoring <- survfit(Surv(y.sorted, 1-delta.sorted) ~ 1  )
  sum.km.censoring<-summary(fit.km.censoring,times=y.sorted,extend=TRUE)
  km.censoring<-sum.km.censoring$surv
  
  km.censoring.minus<-c(1,km.censoring[-length(km.censoring)])
  ratio.km<- delta.sorted/km.censoring.minus
  
  ratio.km[is.nan(ratio.km)]<-0
  
  weight.km<-ratio.km/(sum(ratio.km))
  
  
  wls.fitted<-tryCatch(lm(y.sorted~y.predict.sorted,weights=weight.km),error=function(e){return(c(NA,NA))})
  calibrate.fitted<-tryCatch(predict(wls.fitted),error=function(e){return(c(NA,NA))})
  
  
  
  num.rho2<-sum(weight.km*(calibrate.fitted-sum(weight.km*y.sorted))^2)
  denom.rho2<-sum(weight.km*(y.sorted-sum(weight.km*y.sorted))^2)
  
  R2 <-format(round(num.rho2/denom.rho2,digits = 4) ,nsmall=4)
  
  num.L2<- sum(weight.km*(y.sorted-calibrate.fitted)^2)
  denom.L2<- sum(weight.km*(y.sorted-y.predict.sorted)^2)
  L2 <-format(round(num.L2/denom.L2,digits = 4),nsmall=4)
  
  return(list(R.squared=R2,L.squared=L2))
}

#' @noRd
C_cr <- function(time, status, predicted, tau = NULL, Cause_int = 1,
                 print.value = FALSE) {
  if (any(is.na(time))) {
    stop("The input vector cannot have NA")
  }
  if (any(is.na(status))) {
    stop("The input vector cannot have NA")
  }
  if (any(!(status %in% c(0, 1, 2)))) {
    stop("The status must be 0 or 1 or 2")
  }
  if (any(is.na(predicted))) {
    stop("The input vector cannot have NA")
  }
  if (!(Cause_int %in% status)) {
    stop("Invalid input of Cause_int")
  }
  if (min(time) <= 0) {
    stop("Survival time must be positive")
  }
  
  KM <- summary(survfit(Surv(time, status==0)~1))
  df_KM <- data.frame(time = c(0, KM$time), surv = c(1, KM$surv))
  
  Time_survival <- time
  Censoring <- ifelse(status == 0, 0, 1)
  Cause <- ifelse(status == 2, 2, 1)
  Prediction <- -log(predicted)
  if(is.null(tau)) tau <- max(Time_survival) + 1
  n <- length(Prediction)
  A <- matrix(0, nrow = n, ncol = n)
  B <- matrix(0, nrow = n, ncol = n)
  Q <- matrix(0, nrow = n, ncol = n)
  N_t <- matrix(0, nrow = n, ncol = n)
  Num_mat <- matrix(0, nrow = n, ncol = n)
  Den_mat <- matrix(0, nrow = n, ncol = n)
  Num <- 0
  Den <- 0
  
  w1 <- rep(0, n)
  w2 <- rep(0, n)
  for(i in 1:n){
    w1[i] <- min(df_KM$surv[df_KM$time < Time_survival[i]]) 
    w2[i] <- min(df_KM$surv[df_KM$time <= Time_survival[i]])
  }
  w3 <- w1 %*% t(w1)
  w1 <- w1 * w2
  
  for (i in 1:n) {
    A[i, which(Time_survival[i] < Time_survival)] <- 1
    A[i, ] <- A[i, ] * w1[i]
    B[i, intersect(intersect(which((Time_survival[i] >= 
                                      Time_survival)), which(Cause != Cause_int)), 
                   which(Censoring == 1))] <- 1
    B[i, ] <- B[i, ] * w3[i, ]
    Q[i, which(Prediction[i] > Prediction)] <- 1
  }
  for (i in 1:n) {
    if (Time_survival[i] <= tau && Cause[i] == Cause_int && 
        Censoring[i] == 1) {
      N_t[i, ] <- 1
    }
  }
  Num_mat <- (A + B) * Q * N_t
  Den_mat <- (A + B) * N_t
  Num <- sum(Num_mat)
  Den <- sum(Den_mat)
  if(print.value==TRUE) print(c(Num, Den))
  return(Num/Den)
}

#' Summarize multiple competing-risks models into a wide comparison table
#'
#' @description
#' Calls \code{pam.predicted_survial_eval_cr()} for each model in a named list and
#' returns a wide table where rows are metrics and columns are model names.
#'
#' @param models A **named list**; each element is a list with components:
#'   \itemize{
#'     \item \code{pred_cif}   — numeric matrix of predicted CIF values (rows: time grid; cols: subjects)
#'     \item \code{time.cif}   — numeric vector of the time grid corresponding to \code{pred_cif} rows
#'     \item \code{event_time} — numeric vector of observed times
#'     \item \code{status}     — integer vector of event codes (0=censoring; \code{event_type}=target; others=competing)
#'   }
#' @param metrics Character vector of metrics to compute (passed through to
#'   \code{pam.predicted_survial_eval_cr()}); use \code{"all"} for all supported.
#' @param t_star Optional numeric evaluation time (passed through).
#' @param tau Optional numeric truncation time (passed through).
#' @param event_type Integer code for the cause of interest (passed through).
#'
#' @return A data.frame with:
#' \itemize{
#'   \item each row = metric
#'   \item each column = model name
#'   \item cell values = metric values
#' }
#' @export
pam.summary_cr <- function(models,
                           metrics = NULL,
                           t_star = NULL,
                           tau = NULL,
                           event_type = 1,
                           digits = 2) {
  if (!is.list(models) || length(models) == 0)
    stop("'models' must be a non-empty named list.")
  if (is.null(names(models)) || any(names(models) == ""))
    names(models) <- paste0("Model_", seq_along(models))
  
  required <- c("cif_pred", "times", "status")
  
  eval_one <- function(mod, name) {
    if (!all(required %in% names(mod))) {
      stop(sprintf("Model '%s' must contain: %s",
                   name, paste(required, collapse = ", ")))
    }
    
    res <- pam.predicted_survial_eval_cr(
      pred_cif   = mod$cif_pred[, -1],
      event_time = mod$times,
      time.cif   = mod$cif_pred[, 1],
      status     = mod$status,
      metrics    = metrics,
      t_star     = t_star,
      tau        = tau,
      event_type = event_type
    )
    
    if (!all(c("Metric", "Value") %in% names(res)))
      stop(sprintf("Unexpected result from pam.predicted_survial_eval_cr() for '%s'", name))
    
    res$Model <- name
    res[, c("Model", "Metric", "Value")]
  }
  
  # evaluate each model
  res_list <- mapply(eval_one, models, names(models), SIMPLIFY = FALSE)
  res_long <- do.call(rbind, res_list)
  rownames(res_long) <- NULL
  
  # pivot to wide: rows = Metric, cols = Model, cells = Value
  res_wide <- reshape(
    res_long,
    idvar = "Metric",
    timevar = "Model",
    direction = "wide"
  )
  
  # clean column names "Value.Model" -> "Model"
  names(res_wide) <- sub("^Value\\.", "", names(res_wide))
  rownames(res_wide) <- NULL
  res_wide <- res_wide[, c("Metric", setdiff(names(res_wide), "Metric"))]
  
  # optional: enforce a pleasant metric order if present
  preferred_order <- c("Pseudo_R_square", "R_square", "L_square",
                       "Pseudo_R2_point", "R2_point", "L2_point",
                       "C_index", "Brier Score", "Time Dependent Auc")
  present <- intersect(preferred_order, res_wide$Metric)
  others  <- setdiff(res_wide$Metric, preferred_order)
  res_wide$Metric <- factor(res_wide$Metric, levels = c(present, sort(others)))
  res_wide <- res_wide[order(res_wide$Metric), ]
  res_wide$Metric <- as.character(res_wide$Metric)
  
  # round numeric columns
  numeric_cols <- setdiff(names(res_wide), "Metric")
  res_wide[numeric_cols] <- lapply(res_wide[numeric_cols], function(x) round(as.numeric(x), digits))
  
  res_wide
}
