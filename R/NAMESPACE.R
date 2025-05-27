#' Evaluate Survival Prediction Performance with Competing Risks
#'
#' This function evaluates the predictive performance of a survival model in the presence of competing risks.
#' It calculates various performance metrics including Pseudo R-squared, R-squared, L-squared, C-index, 
#' Brier Score, and Time-Dependent AUC.
#'
#' @param pred_cif A matrix of predicted cumulative incidence functions (CIFs) for different subjects.
#' @param event_time A numeric vector of observed event or censoring times for subjects.
#' @param time.cif A numeric vector of time points corresponding to the predicted CIFs.
#' @param status A numeric vector indicating the event type for each subject:
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
#'     \item `"C_index"`
#'     \item `"Brier Score"`
#'     \item `"Time Dependent Auc"`
#'   }
#' @param t_star Optional numeric value specifying a reference time for evaluating 
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
#' @examples
#' \dontrun{
#' result <- pam.predicted_survial_eval_cr(
#'     pred_cif = pred_matrix,
#'     event_time = event_times,
#'     time.cif = pred_time_points,
#'     status = event_status,
#'     metrics = c("C_index", "Brier Score"),
#'     tau = 5
#' )
#' print(result)
#' }
#'
#' 
#' 
#' 
#' @export
pam.predicted_survial_eval_cr <- function (pred_cif, event_time, time.cif, status,
                                           metrics = "all",  t_star = NULL, tau = NULL, event_type = 1) 
  
{
  
  if (missing(event_time) || missing(pred_cif)) {
    stop("Please provide 'event_time',  and 'pred_cif' arguments.")
  }
  
  
  
  metrics_results <- list()
  
  valid_metrics <- c("Pseudo_R_square", "R_square", "L_square", "C_index", "Brier Score", "Time Dependent Auc")
  
  if (metrics == "all") {
    metrics <- valid_metrics
  } else {
    invalid <- setdiff(metrics, valid_metrics)
    if (length(invalid) > 0) stop("Invalid metrics: ", paste(invalid, collapse = ", "))
  }
  
  
  
  if("Pseudo_R_square" %in% metrics || "R_square" %in% metrics || "L_square" %in% metrics) {
    r_l_list <- pam.censor.cr(event_time, status, tau, pred_cif, time.cif, event_type)
  }
  
  if ("Pseudo_R_square" %in% metrics) {
    metrics_results$Pseudo_R_square <- r_l_list$Pseudo_R_squared
  } 
  if ("R_square" %in% metrics) {
    metrics_results$R_square <- r_l_list$R_squared
  }
  if ("L_square" %in% metrics) {
    metrics_results$L_square <- r_l_list$L_squared
  }
  
  pred_risk <- apply(pred_cif, 2, m_cif, time.cif=time.cif, tau=tau)
  
  C_index <- C_cr(
    time = event_time,
    status = status,
    predicted = pred_risk,
    tau = tau,
    Cause_int = 1
  )
  if ("C_index" %in% metrics) {
    metrics_results$"C_index" <- round(C_index, 4)
  }
  
  if ("Brier Score" %in% metrics) {
    if (is.null(t_star)){
      t_star = median(event_time)
    }
    time_idx <- max(which(time.cif <= tau))
    X <- pred.cif[time_idx, ]
    
    brier_result <- tdROC.cr(
      X = X,  
      Y = event_time,      
      delta = status,
      tau = tau, 
      nboot = 0)
    
    metrics_results$"Brier Score" <- round(as.numeric(brier_result$calibration_res[1]), 4)
  }
  
  
  if ("Time Dependent Auc" %in% metrics) {
    if (is.null(t_star)){
      t_star = median(event_time)
    }
    AUC_result <- tdROC.cr(
      X = X,  
      Y = event_time,      
      delta = status,
      tau = tau,       
      method = "both", 
      output = "AUC"   
    )
    metrics_results$"Time Dependent Auc Integral" <- round(AUC_result$main_res$AUC.integral, 4)
    metrics_results$"Time Dependent Auc Empirical" <- round(AUC_result$main_res$AUC.empirical, 4)
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


