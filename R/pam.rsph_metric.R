#' @title RE Measure of Explained Variation Using Predicted and Observed Data
#'
#' @description This function computes the RE (explained variation) metric using 
#' predicted risk scores and observed survival times. It evaluates the proportion 
#' of variability in survival times explained by the covariates using a simplified 
#' approach that requires only the predicted values and observed survival data.
#'
#' @param predicted_data A numeric vector of predicted risk scores or linear predictors from a model.
#' @param survival_time A numeric vector of observed survival times.
#' @param status A numeric or logical vector indicating event status (1 for event, 0 for censoring).
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item \code{Re}: The RE measure, representing the proportion of explained variation.
#'   \item \code{numerator}: The numerator of the RE calculation, representing explained variation.
#'   \item \code{denominator}: The denominator of the RE calculation, representing total variation.
#' }
#'
#' @details This function assumes the input data is complete and properly formatted.
#' It adjusts for censoring using Kaplan-Meier weights and computes ranks for the predicted values
#' to evaluate the explained variation in survival times.
#' @references
#' Schemper, M., & Henderson, R. (2000). Predictive accuracy and explained variation in Cox regression. 
#' Biometrics, 56, 249–255.
#' 
#' Lusa, L., Miceli, R., & Mariani, L. (2007). Estimation of predictive accuracy in survival analysis using R and S-PLUS. 
#' Computer Methods and Programs in Biomedicine, 87, 132–137.
#' 
#' @examples
#' library(survival)
#' 
#' data("lung")
#' 
#' predicted_data <- lung$ph.ecog       
#' survival_time <- lung$time          
#' status <- lung$status - 1           
#' 
#' complete_cases <- complete.cases(predicted_data, survival_time, status)
#' predicted_data <- predicted_data[complete_cases]
#' survival_time <- survival_time[complete_cases]
#' status <- status[complete_cases]
#' 
#' # Compute the RE measure
#' result <- pam.rsph_metricc(predicted_data, survival_time, status)
#' 
#' cat("RE Measure:", result$Re, "\n")
#' @keywords internal
#' @noRd
pam.rsph_metric <- function(predicted_data, survival_time, status, start_time = NULL) {
if (is.null(start_time)) {
  start_time <- rep(0, length(survival_time))
}
# Input validation
if (length(predicted_data) != length(survival_time) || 
    length(survival_time) != length(status) || 
    length(start_time) != length(status)) {
  stop("Lengths of predicted_data, start_time, survival_time, and status must match.")
}

if (!all(status %in% c(0, 1))) {
  stop("Status vector must only contain 0 (censored) and 1 (event) values.")
}

# Sort data by survival time and event status
sorted_indices <- order(survival_time, -status)
predicted_data <- predicted_data[sorted_indices]
survival_time <- survival_time[sorted_indices]
start_time <- start_time[sorted_indices]
status <- status[sorted_indices]

km_fit <- my.survfit(start = start_time, stop = survival_time, event = status)

Gmat <- matrix(km_fit$surv.i2, nrow = length(survival_time), ncol = length(unique(survival_time)), byrow = TRUE)

numerator <- 0
denominator <- 0

# Extract unique event times
event_times <- unique(survival_time[status == 1])

# Compute individual hazards (pi)
pi <- exp(predicted_data) / sum(exp(predicted_data)) # Hazard probabilities
# Loop through each event time
for (time in event_times) {
  # Identify at-risk individuals
  at_risk <- (survival_time >= time & start_time < time)
  events <- (survival_time == time & status == 1)
  
  # Compute ranks for at-risk individuals
  rangi <- rank(-predicted_data[at_risk]) / Gmat[at_risk, event_times == time]+0.5
  r0 <- mean(rangi)
  rP <- 0.5 + 0.5 * sum((survival_time[at_risk] == time & status[at_risk] == 1) / Gmat[at_risk, event_times == time])
  
  # Adjust ranks using pi and Gmat
  r0_pi <- mean(rangi)
  
  numerator <- numerator + sum((r0_pi - rangi) * pi[at_risk] / Gmat[at_risk, event_times == time])
  denominator <- denominator + sum((r0_pi - rP) * pi[at_risk] / Gmat[at_risk, event_times == time])
}

# Compute Re measure
re <- numerator / denominator


return(list(
  Re = re,
  numerator = numerator,
  denominator = denominator
))
}



"my.survfit" <- function(start,stop,event){
  #internal function for calculation of weights (allows counting process notation)
  
  if(missing(start)) start <- rep(0,length(stop))   #if no count. process notation, set everyone to start at 0
  
  data <- data.frame(start=start,stop=stop,event=event) 
  
  data <- data[order(data$start),]      #ordered according to the starting time
  data$Y <- unlist(lapply(data$stop,function(x)sum(data$start<=x))) #for each time point: the number of indiv. that started before (or at) it
  
  data <- data[order(data$stop,-data$event),]   #order in time, put events before censoring
  
  howmany <- nrow(data)         #total number of rows
  
  lt.t1 <- c(0,data$stop[-nrow(data)])      #aux. vector for searching the stop times 
  lt.t2 <- c(data$stop[-1],data$stop[nrow(data)]+1)
  lt1 <- which(data$stop!=lt.t1)        #first case at a time point
  lt2 <- which(data$stop!=lt.t2)        #last case at a time point
  
  ct.t1 <- c(-1,data$start[-nrow(data)])      #aux. vector for searching the start times
  ct.t2 <- c(data$start[-1],data$start[nrow(data)]+1)
  ct1 <- which(data$start!=ct.t1)       #first case at a time point
  ct2 <- which(data$start!=ct.t2)       #last case at a time point
  n.incom <- rep(0,length(lt1))       #prepare a vector that will contain the no. of incoming indiv. at a time point (length=no. of unique time points)
  inx <- unlist(lapply(data$start[ct1],function(x)min(which(data$stop[lt1]>=x)))) #the index of the first one to fail after a starting point
  incom.sum <- diff(c(0,ct2))       #the number of incoming
  if(length(inx)>1){
    if(inx[1]==inx[2]){         #if indiv. come in at time 0 and first stop time
      incom.sum <- c(sum(incom.sum[1:2]),incom.sum[-(1:2)])    #sum them up
      inx <- unique(inx)          #only one value per each stop time
    }
  }
  n.incom[inx] <- incom.sum       #at each event time, we need the number of patients that got in just before it.
  n.incom2 <- n.incom         
  n.incom2[1] <- n.incom2[1] - sum(data$start==0)   #the number of additional patients to come in at time 0 is set to 0
  
  
  evcum <- cumsum(data$event)     #number of events up to a certain time point
  evcum <- evcum[lt2]       #number of events for t<=t_i
  cencum <- cumsum(data$event==0)     #number of censored
  cencum <- cencum[lt2]       #number of censored t<=t_i
  n.event <- diff(c(0,evcum))     #number of events at each unique event time
  n.cens <- diff(c(-sum(data$start==0),cencum)) #number of censored at each unique event time, the first number is the total of the patients
  n.t.cens <- n.cens - n.incom      #number truly censored: number censored - number of incoming (0 at time 0)
  n.risk <- (data$Y - (0:(howmany-1)))[lt1]-n.incom2 #number at risk at each event time: no.started before or at this time - no. of incoming at this time - 1 per row before this time
  
  kmji <- 1- n.event/n.risk     #proportion of survivors
  
  kmji.i <- 1- n.t.cens/(n.risk-n.event)    #the events are taken out of the risk set before calculating the G(t_i)
  
  km <- (cumprod(kmji))       #km 
  km.i <- cumprod(kmji.i)       #km for censoring     
  
  #correction for the last value (if division by 0) (this happens if the last time is an event time)
  inx <- which(is.na(km.i))       #missing at the end, only needed if not lagged?
  if(length(inx)){
    km.i[inx] <- km.i[min(inx)-1] # carrie last value forward
  }
  
  km.i2 <- c(1,km.i[-length(km.i)])     #this gives the lagged weight (weight at time just before t)
  
  list(surv=km,surv.i=km.i,n.event=n.event,n.cens=n.t.cens,time=data$stop[lt1],n.risk=n.risk,surv.i2=km.i2)
}

