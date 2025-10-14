#' @title Prediction Accuracy Measures for Parametric Survival Regression Models
#'
#' @description This function computes Calculates a set of measures, $R^2$ and $L^2$, $Pseudo-R^2$ for evaluating the predictive accuracy of parametric survival regression models. 
#' R-squared quantifies the proportion of variability in the response that is explained by a corrected prediction function derived from the original model. L-squared measures the proportion of the prediction error from the original model that is explained by the corrected prediction function, effectively capturing the improvement in prediction accuracy. Together, these metrics provide a comprehensive assessment of the predictive performance of a survival regression model.
#'
#' @param model An object of class survreg representing a fitted parametric survival regression model. The survreg call must include x = TRUE and y = TRUE to ensure the design matrix and response vector are stored in the model object.

#' @param covs A vector of covariate names to be used for prediction when newdata is provided.
#' @param tau (Optional) A time point for truncating the survival time. If provided, the function evaluates predictions up to this time point.
#' @param new_data (Optional) A new dataset for evaluating the model. If NULL, the function uses the training data stored in model.
#' @param predict (Optional) A logical value indicating whether to return individual predictions. Default is FALSE. 
#' @return A list containing three components:
#' - R.squared: The R-squared measure, quantifying the explained variability in the response.
#' - L.squared: The L-squared measure, quantifying the reduction in prediction error.
#' - Psuedo.R: A pseudo R-squared measure, calculated as the product of R-squared and L-squared.
#'
#' @examples
#' library(survival)
#' library(PAmeasures)
#'
#' # Use Mayo Clinic Primary Biliary Cirrhosis Data
#' data(pbc)
#'
#' # Fit an exponential model with bilirubin
#' fit1 <- survreg(Surv(time, status == 2) ~ bili, data = pbc, dist = "exponential", x = TRUE, y = TRUE)
#' pam.survreg(fit1)
#'
#' # Fit a lognormal model with standardised blood clotting time
#' fit2 <- survreg(Surv(time, status == 2) ~ protime, data = pbc, dist = "lognormal", x = TRUE, y = TRUE)
#' pam.survreg(fit2)
#'
#' # Fit a Weibull model with bilirubin and standardised blood clotting time
#' fit3 <- survreg(Surv(time, status == 2) ~ bili + protime, data = pbc, dist = "weibull", x = TRUE, y = TRUE)
#' pam.survreg(fit3)
#'
#' @export

pam.surverg_restricted <- function(model, covs, tau = NULL,  new_data = NULL, predict = T) 
{
  
  # Check inputs
  if (!inherits(model, "survreg")) {
    stop("model must be an object of class 'survreg'.")
  }
  if (is.null(model$x) || is.null(model$y)) {
    stop("model must include x = TRUE and y = TRUE.")
  }
  
  if(is.null(new_data)){
    x.matrix.unsorted <- model$x
    y.unsorted <- model$y[, 1]
    censor.unsorted <- model$y[, 2]
    y.order.new <- NULL
    
  } else {
    if (!all(covs %in% colnames(new_data))) {
      stop("All covariates must be present in new_data.")
    }
    time_var = "time"
    status_var = "status"
    if (any(!c(time_var, status_var) %in% colnames(new_data))) {
      stop("new_data require time and status columns")
    }
    y.unsorted <- model$y[, 1]
    censor.unsorted <- model$y[, 2]
    x.matrix.unsorted <- new_data[, covs, drop = FALSE]
    y.unsorted.new <- new_data[[time_var]]
    censor.unsorted.new <- new_data[[status_var]]
    y.order.new <- order(y.unsorted.new)
    
  }
  
  if (is.null(tau)) {
    tau <- max(new_data[[time_var]]) + 1  
  }
  nsize <- length(y.unsorted)
  y.order <- order(y.unsorted)
  y <- y.unsorted[y.order]
  delta <- censor.unsorted[y.order]
  p <- dim(as.matrix(x.matrix.unsorted))[2]
  
  if(is.null(y.order.new)==F) {
    y.order <- y.order.new
  }
  
  if (p == 1) {
    x.matrix <- as.matrix(x.matrix.unsorted[y.order])
  } else {
    x.matrix <- x.matrix.unsorted[y.order, ]
  }
  
  
  nsize <- length(y)
  fit.km.censoring <- survfit(Surv(y, 1 - delta) ~ 1)
  
  if(is.null(y.order.new)==F){
    y <- y.unsorted.new[y.order]
    delta <- censor.unsorted.new[y.order]
  } 
  
  sum.km.censoring <- summary(fit.km.censoring, times = y, extend = TRUE)
  km.censoring <- sum.km.censoring$surv
  km.censoring.minus <- c(1, km.censoring[-length(km.censoring)])
  
  ratio.km <- delta/km.censoring.minus
  ratio.km[is.nan(ratio.km)] <- 0
  weight.km <- ratio.km/(sum(ratio.km))
  
  if (model$dist %in% c("exponential", "weibull")) {
    if(is.null(new_data)){
      temp.pred <- predict(model, type = "response")[y.order]
    }else{
      temp.pred <- predict(model, newdata = as.data.frame(x.matrix), type = "response")
    }
    
    gtau <- (tau/temp.pred)^(1/model$scale)
    t.predicted <- temp.pred * gamma(1 + model$scale) *
      (1-expint::gammainc((1 + model$scale), gtau))
  } else if (model$dist == "lognormal") {
    if(is.null(new_data)){
      temp.pred <- predict(model, type = "response")[y.order]
    }else{
      temp.pred <- predict(model, newdata = as.data.frame(x.matrix), type = "response")
    }
    gtau <- (log(tau)-log(temp.pred))/model$scale
    t.predicted <- temp.pred * exp((model$scale)^2/2) * pnorm(gtau)
  }
  
  if (!is.null(tau)) {
    delta <- ifelse(y <= tau, delta, 0)
    y <- pmin(y, tau)
    y.input <- ifelse(tau > y, y, tau)
  }
  
  if(predict == T){
    time_points <- y
    time_points <- time_points[time_points <= tau] 
    
    pre_sp <- predictSurvProb2survreg(model, as.data.frame(x.matrix), time_points)
    t.predicted2 <- integrate_survival(pre_sp, y, delta, tau)
    return(list(pred = t.predicted2, #t.predicted,
                times = y,
                status = delta,
                surv_prob = pre_sp,
                linear.pred = predict(model, 
                                      newdata = as.data.frame(x.matrix), 
                                      type = "lp")))
  }
  
  wls.fitted <- tryCatch(lm(y ~ t.predicted, weights = weight.km), 
                         error = function(e) {
                           cat("Error in WLS fitting:", e$message, "\n")
                           return(c(NA, NA))
                         })
  calibrate.fitted <- tryCatch(predict(wls.fitted), error = function(e) {
    return(c(NA, NA))
  })
  
  num.rho2 <- sum(weight.km * (calibrate.fitted - sum(weight.km * y))^2)
  denom.rho2 <- sum(weight.km * (y - sum(weight.km * y))^2)
  R2 <- round(num.rho2/denom.rho2, digits = 4)
  num.L2 <- sum(weight.km * (y - calibrate.fitted)^2)
  denom.L2 <- sum(weight.km * (y - t.predicted)^2)
  L2 <- round(num.L2/denom.L2, digits = 4)
  SR <- round(R2 * L2, digits = 4)
  
  return(list(R.squared = format(R2, nsmall = 4), 
              L.squared = format(L2, nsmall = 4), 
              Psuedo.R = format(SR, nsmall = 4)))
}