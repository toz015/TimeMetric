
#' Predict Subject-Specific Cumulative Incidence Function (CIF)
#'
#' This function calculates cumulative incidence function (CIF) predictions for individual subjects 
#' across multiple types of competing risks models. It supports cause-specific Cox models, 
#' accelerated failure time (AFT) models (Weibull, Log-normal, Log-logistic), Fine-Gray models, 
#' and competing risks random survival forests.
#'
#' @param model1 The fitted model for the target cause (either a `coxph` or `survreg` object).
#' @param model2 The fitted model for competing causes (either a `coxph` or `survreg` object).
#' @param final_data A data frame containing observed event times (`obs.times`), event indicators (`obs.event`), and covariates.
#' @param event.type The event type of interest for which CIF is to be predicted.
#' @param tau The maximum time horizon for prediction.
#' @param num.grid Number of evaluation time points (i.e., grid points) for computing CIF. 
#'        This parameter only affects calculations used in **pseudo R² metrics**. For standard use cases, keep `num.grid = 1` (the default).
#' @param fg_model Optional Fine-Gray model object (`crr` from `cmprsk` package).
#' @param cr_model Optional competing risks random forest model object.
#' @return A list with two components:
#' \describe{
#'   \item{tau_list}{Vector of time points used for evaluation.}
#'   \item{m_pred}{A matrix where the first column is time, and subsequent columns are predicted CIFs for each subject.}
#' }
#'


handle_error <- function(e) {
  print(e)
  return()
}


#######################################################################
# Internal helper function: calculate CIF based on Cox cause-specific hazard models
#######################################################################


get_CIF <- function(cox.focus, cox.other, X){
  
  fit <- tryCatch({survfit(cox.focus, newdata = X, se.fit = FALSE)}, 
                  error = handle_error)
  Ht <- fit$cumhaz[fit$n.event>=1,]
  if(is.null(dim(Ht))) return()
  hi <- apply(Ht, 2, function(x) diff(c(0, x)))
  
  for(cox.model in cox.other){
    fit_other <- tryCatch({survfit(cox.model, newdata = X, se.fit = FALSE)}, 
                          error = handle_error)
    Ht_other <- fit_other$cumhaz[fit$n.event>=1,]
    Ht <- Ht + Ht_other
    if(is.null(dim(Ht))) return()
  }
  
  CIF <- apply(exp(-Ht) * hi, 2, cumsum)
  return(list(y = fit$time[fit$n.event>=1], CIF = CIF))
}


#######################################################################
# Internal helper function: calculate CIF based on AFT models
#######################################################################


get_CIF_aft <- function(ftime, fstatus, aft.focus, aft.other, newX = NULL){
  if(is.null(newX)) newX <- aft.focus$x[,-1]
  lambda <- exp(-aft.focus$coefficients[1]/aft.focus$scale)
  alpha <- 1/aft.focus$scale
  rs1 <- (-newX %*% aft.focus$coefficients[-1]) %>% exp
  if(aft.focus$dist %in% c("weibull", "exponential")){
    hi <- sapply(sort(ftime) %>% unique,
                 function(x) rs1 * lambda * alpha * (rs1 * x)^(alpha-1)) %>% t
    Ht <- sapply(sort(ftime) %>% unique,
                 function(x) lambda * (rs1 * x)^(alpha)) %>% t
  }else if(aft.focus$dist == "lognormal"){
    St <- sapply(sort(ftime) %>% unique,
                 function(x){1 - pnorm((log(x) -
                                          as.matrix(cbind(1, newX)) %*% 
                                          aft.focus$coefficients)/
                                         aft.focus$scale)}) %>% t
    fx <- sapply(sort(ftime) %>% unique,
                 function(x){exp(-1/2*((log(x)-
                                          as.matrix(cbind(1, newX)) %*% 
                                          aft.focus$coefficients)/
                                         aft.focus$scale)^2)/
                     (x*sqrt(2*pi)*aft.focus$scale)}) %>% t
    hi <- fx / St 
    Ht <- - log(St)
  }else if(aft.focus$dist == "loglogistic"){
    hi <- sapply(sort(ftime) %>% unique,
                 function(x){(rs1 * lambda * alpha * (rs1 * x)^(alpha-1))/
                     (1 + lambda *(rs1 * x)^(alpha))}) %>% t
    Ht <- sapply(sort(ftime) %>% unique,
                 function(x) -log(1/(1 + lambda * (rs1 * x)^(alpha)))) %>% t
  }
  
  
  for(aft.model in aft.other){
    lambda <- exp(-aft.model$coefficients[1]/aft.model$scale)
    alpha <- 1/aft.model$scale
    rs2 <- (-newX %*% aft.model$coefficient[-1]) %>% exp
    if(aft.model$dist %in% c("weibull", "exponential")){
      Ht_other <- sapply(sort(ftime) %>% unique,
                         function(x) lambda * (rs2 * x)^(alpha)) %>% t
      Ht <- Ht + Ht_other
    }else if(aft.model$dist == "lognormal"){
      St <- sapply(sort(ftime) %>% unique,
                   function(x){1 - pnorm((log(x) -
                                            as.matrix(cbind(1, newX)) %*% 
                                            aft.model$coefficient)/
                                           aft.model$scale)}) %>% t
      Ht_other <- - log(St)
      Ht <- Ht + Ht_other
    }else if(aft.model$dist == "loglogistic"){
      Ht_other <- sapply(sort(ftime) %>% unique,
                         function(x) -log(1/(1 + lambda * (rs2 * x)^(alpha)))) %>% t
      Ht <- Ht + Ht_other
    }
    
    if(is.null(dim(Ht))) return()
  }
  
  temp.matrix <- rbind(rep(0, dim(Ht)[2]), exp(-Ht) * hi)
  new_matrix <- (temp.matrix[-nrow(temp.matrix), ] + temp.matrix[-1, ]) / 2
  new_matrix <- diff(c(0, sort(ftime) %>% unique)) * new_matrix
  new_matrix <- apply(new_matrix, 2, cumsum)
  return(list(y = sort(ftime) %>% unique, CIF = new_matrix))
}




#' @export


pam.predict_subject_cif <- function(model1, model2, final_data, event.type, tau = NULL, num.grid = 1, fg_model = NULL, cr_model = NULL) {
  
  quantile_list <- seq(0, 1, length.out = num.grid)
  LB_tau <- min(final_data$obs.times[final_data$obs.event == event.type])
  UB_tau <- max(final_data$obs.times[final_data$obs.event == event.type])
  tau_list <- qunif(quantile_list, min = LB_tau, max = min(UB_tau, tau))
  
  
  ###  Handle Fine-Gray model (crr model fg_model)
  if (!is.null(fg_model)) {
    m_pred <- predict(fg_model, as.matrix(final_data[, -c(1,2)]))
    if (is.null(dim(m_pred[, -1]))) return(NULL)
    return(list(tau_list = tau_list, m_pred = m_pred))
  }
  
  if (!is.null(cr_model)) {
    cr_pred <- predict(cr_model, newdata = final_data)
    cif_event <- cr_pred$cif[,,event.type]
    m_pred <- cbind(cr_pred$time.interest, t(cif_event))
    return(list(tau_list = tau_list, m_pred = m_pred))
  }
  
  if (inherits(model1, "coxph")) {
    model.info1 <- get_CIF(model1, list(model2), final_data)
    m_pred <- cbind(model.info1$y, model.info1$CIF)
    return(list(tau_list = tau_list, m_pred = m_pred))
  }
  
  if (inherits(model1, "survreg")) {
    model.info1 <- get_CIF_aft(final_data$obs.times, final_data$obs.event, model1, list(model2),
                               newX = as.matrix(final_data[,-c(1,2)]))
    m_pred <- cbind(model.info1$y, model.info1$CIF)
    
    return(list(tau_list = tau_list, m_pred = m_pred))
  }
  
  stop("Unknown model type")
}
