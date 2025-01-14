#' Compute the Concordance Statistic for Survival Data
#'
#'
#' This function calculates the Concordance Index (C-index), a performance metric 
#' used to evaluate survival prediction models. The C-index measures the ability 
#' of a model to correctly rank predicted survival times with observed survival outcomes.
#'
#' @param predicted_time Numeric vector of predicted survival times.
#' @param survival_time Numeric vector of observed survival times.
#' @param status Numeric vector indicating censoring status (1 = event occurred, 0 = censored).
#' @param weight Character string specifying the weighting method to use. Options are:
#'   - `"H"`: Harrell's C-index.
#'   - `"U"`: Unweighted C-index with KM survival probability.
#'   - `"U_tau"`: Unweighted C-index truncated at \code{input_tau}.
#' @param input_tau Optional numeric value specifying the truncation point for \code{"U_tau"} weighting.
#'
#' @details The function performs pairwise comparisons of predicted survival times and observed outcomes. 
#' For weighted calculations (`"U"` and `"U_tau"`), the Kaplan-Meier survival probability 
#' is used to weight the pairs. The `"U_tau"` method additionally truncates comparisons 
#' at a specified time point (\code{input_tau}).
#'
#' @return A numeric value representing the computed C-index. 
#' A value closer to 1 indicates better model performance, while a value closer to 0.5 
#' indicates no discriminative power.
#'
#' @references 
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
#' @examples
#' library(PAmeasure)
#' predicted_time <- c(2.5, 3.2, 1.8, 4.1)
#' survival_time <- c(3, 4, 2, 5)
#' status <- c(1, 0, 1, 1)
#' pam.concordance.metric(predicted_time, survival_time, status, weight = "H")
#'
#' @export
pam.concordance.metric <- function(predicted_time, survival_time, status, 
                           weight = "H", input_tau = NULL) {
  concordant <- 0
  usable_pairs <- 0
  if(weight == "U"){
    KM <- summary(survfit(Surv(survival_time, status==0)~1))
    df_KM <- data.frame(time = c(0, KM$time), surv = c(1, KM$surv))
    for (i in 1:length(survival_time)) {
      for (j in 1:length(survival_time)) {
        if(status[i] == 1){
          usable_pairs <- usable_pairs + 
            (survival_time[i] < survival_time[j]) * 
            min(df_KM$surv[df_KM$time <= survival_time[i]])^(-2)
          concordant <- concordant +
            (predicted_time[i] < predicted_time[j]) * 
            (survival_time[i] < survival_time[j]) * 
            min(df_KM$surv[df_KM$time <= survival_time[i]])^(-2)
        }
      }
    }
  } else if(weight == "H"){
    for (i in 1:length(survival_time)) {
      for (j in 1:length(survival_time)) {
        if(i != j & status[i] == 1){
          usable_pairs <- usable_pairs + 
            (survival_time[i] < survival_time[j])
          concordant <- concordant +
            (predicted_time[i] < predicted_time[j]) * 
            (survival_time[i] < survival_time[j])
        }
      }
    }
  } else if(weight == "U_tau"){
    KM <- summary(survfit(Surv(survival_time, status==0)~1))
    df_KM <- data.frame(time = c(0, KM$time), surv = c(1, KM$surv))
    for (i in 1:length(survival_time)) {
      for (j in 1:length(survival_time)) {
        if(status[i] == 1){
          usable_pairs <- usable_pairs + 
            (survival_time[i] < survival_time[j]) * 
            (survival_time[i] < input_tau) *
            min(df_KM$surv[df_KM$time <= survival_time[i]])^(-2)
          concordant <- concordant +
            (predicted_time[i] < predicted_time[j]) * 
            (survival_time[i] < survival_time[j]) *
            (survival_time[i] < input_tau) *
            min(df_KM$surv[df_KM$time <= survival_time[i]])^(-2)
        }
      }
    }
  } 
  
  cindex <- concordant / usable_pairs
  return(cindex)
}
