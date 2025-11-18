library(TimeMetric)
library(tidyverse)

set.seed(123)
event_type <- 1 # define event type you want to focus
n <- 300 # number of subjects  
ran.time <- runif(n) # generate a vector of random prob (you can put your predicted CIF here) 
cif_pred <- matrix(data = ran.time, nrow = 1) # put predicted CIF in matrix  
tau <- 50; t_star <- 50; time.cif <- 50 # set time point you evaluated
event_time <- ran.time*100 # input your event time
status <- rbinom(n, size = 2, prob = 0.5) # input your status

pam.predicted_survial_eval_cr(
  pred_cif   = cif_pred,
  event_time = event_time,
  time.cif   = time.cif,
  status     = status,
  metrics    = "Pseudo_R2_point",
  t_star     = t_star,
  tau        = tau,
  event_type = event_type
)
