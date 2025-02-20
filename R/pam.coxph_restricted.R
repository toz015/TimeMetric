#' @title Prediction Accuracy Measures for Cox Proportional Hazards Model
#'
#' @description This function calculates a pair of measures, R-squared and L-squared, for a Cox proportional hazards model. 
#' R-squared quantifies the proportion of variability in the observed survival times explained by the model's predictions,
#' while L-squared measures the proportion of prediction error explained by a corrected prediction function. Together, these metrics provide a comprehensive evaluation of the predictive power of the Cox model.
#'
#' @param fit.cox An object of class `coxph` representing a fitted Cox proportional hazards regression model. The model must be fitted with `x = TRUE` and `y = TRUE` to include the design matrix and response vector.
#' @param covariates A character vector specifying the names of the covariates used in the model.
#' @param time_var A string specifying the name of the time variable in the data.
#' @param status_var A string specifying the name of the event indicator variable in the data (1 for event, 0 for censored).
#' @param tau (Optional) A numeric value specifying a time horizon for truncating the observed survival times. If provided, the function calculates metrics up to time `tau`.
#' @param newdata (Optional) A new dataset for validation. If provided, the function calculates metrics on this new dataset instead of the training data.
#'
#' @return A list containing the following components:
#' \item{R.squared}{The R-squared measure, quantifying the proportion of variability explained by the model.}
#' \item{L.squared}{The L-squared measure, quantifying the proportion of prediction error explained by the corrected prediction.}
#' \item{Psuedo.R}{A pseudo-R measure, calculated as the product of R-squared and L-squared.}
#'
#' @references
#' Li, G., & Wang, X. (2016). Prediction Accuracy Measures for a Nonlinear Model and for Right-Censored Time-to-Event Data.
#' arXiv preprint arXiv:1611.03063. Available at \url{https://arxiv.org/abs/1611.03063}.
#'
#' @examples
#' library(survival)
#' library(PAmeasure)
#'
#' # Use Mayo Clinic Primary Biliary Cirrhosis Data
#' data(pbc)
#'
#' # Fit a univariate Cox PH model with bilirubin
#' fit1 <- coxph(Surv(time, status == 2) ~ bili, data = pbc, x = TRUE, y = TRUE)
#'
#' # Calculate prediction accuracy measures
#' pam.cox_metrics(fit1, covariates = "bili", time_var = "time", status_var = "status")
#'
#' # Fit a multivariate Cox PH model with bilirubin and albumin
#' fit2 <- coxph(Surv(time, status == 2) ~ bili + albumin, data = pbc, x = TRUE, y = TRUE)
#'
#' # Calculate prediction accuracy measures
#' pam.cox_metrics(fit2, covariates = c("bili", "albumin"), time_var = "time", status_var = "status")
#'
#' @export

pam.coxph_restricted <- function(fit.cox, covariates, time_var, status_var, tau = NULL, newdata = NULL) 
{
  if(is.null(newdata)){
    x.matrix.unsorted <- fit.cox$x
    y.unsorted <- fit.cox$y[, 1]
    censor.unsorted <- fit.cox$y[, 2]
    y.order.new <- NULL
  }else{
    y.unsorted <- fit.cox$y[, 1]
    censor.unsorted <- fit.cox$y[, 2]
    x.matrix.unsorted <- newdata[, covariates, drop = FALSE]
    y.unsorted.new <- newdata[[time_var]]
    censor.unsorted.new <- newdata[[status_var]]
    y.order.new <- order(y.unsorted.new)
  }
  my.beta <- fit.cox$coeff
  y.order <- order(y.unsorted)
  nsize <- length(y.unsorted)
  y <- y.unsorted[y.order]
  delta <- censor.unsorted[y.order]
  p <- dim(as.matrix(x.matrix.unsorted))[2]
  if(is.null(y.order.new)==F) y.order <- y.order.new
  if (p == 1) {
    x.matrix <- as.matrix(x.matrix.unsorted[y.order])
  } else {
    x.matrix <- x.matrix.unsorted[y.order, ]
  }
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
  y.length <- length(y)
  yi.matrix <- matrix(rep(y, each = y.length), nrow = y.length)
  yj.matrix <- t(yi.matrix)
  R <- ((yj.matrix >= yi.matrix) * 1)
  if(is.null(newdata)){
    temp_rs <- x.matrix %*% my.beta
  }else{
    temp_rs <- predict(fit.cox, newdata[, covariates, drop = FALSE], type = "lp")[y.order]
  }
  my.Lambda <- R %*% ((delta)/t(t(exp(temp_rs)) %*% R))
  rm(R)
  my.power <- matrix(rep(t(exp(temp_rs)), each = y.length), 
                     nrow = y.length)
  my.factor <- (max(my.Lambda)/100)
  my.Lambda2 <- t(matrix(rep(exp(-my.Lambda/my.factor), each = y.length), 
                         nrow = y.length))
  S.hat.x <- my.Lambda2^(my.factor * my.power)
  rm(my.Lambda)
  rm(my.power)
  rm(my.Lambda2)
  
  # Truncate data if tau is provided
  if (!is.null(tau)) {
    delta <- ifelse(y <= tau, delta, 0)
    y <- pmin(y, tau)
    y.input <- ifelse(tau > y, y, tau)
  }
  else {
    y.input <- NULL
  }
  
  
  if(is.null(y.input)){ 
    t1 <- y
  }else{
    t1 <- y.input[y.order]
  }
  
  t2 <- c(0, t1[1:length(t1) - 1])
  delta.t <- t1 - t2
  t.predicted <- colSums(delta.t * S.hat.x)
  
  wls.fitted <- tryCatch(lm(y ~ t.predicted, weights = weight.km), 
                         error = function(e) {
                           return(c(NA, NA))
                         })
  calibrate.fitted <- tryCatch(predict(wls.fitted), error = function(e) {
    return(c(NA, NA))
  })
  num.rho2 <- sum(weight.km * (calibrate.fitted - sum(weight.km * 
                                                        y))^2)
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

