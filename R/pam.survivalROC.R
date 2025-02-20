#' @title Time-dependent ROC Curve from Censored Survival Data
#'
#' @description This function creates a time-dependent ROC curve from censored survival data
#' using the Kaplan-Meier (KM) or Nearest Neighbor Estimation (NNE) method of
#' Heagerty, Lumley, and Pepe (2000).
#'
#' @param Stime Event time or censoring time for subjects.
#' @param status Indicator of status, 1 if death or event, 0 otherwise.
#' @param marker Predictor or marker value.
#' @param entry Entry time for the subjects.
#' @param predict.time Time point of the ROC curve.
#' @param cut.values Marker values to use as a cut-off for calculating sensitivity and specificity.
#' @param method Method for fitting joint distribution of (marker, t), either "KM" or "NNE". Default is "NNE".
#' @param lambda Smoothing parameter for NNE.
#' @param span Span for the NNE. Either \code{lambda} or \code{span} is required for NNE.
#' @param window Window type for NNE, either "symmetric" or "asymmetric".
#' 
#' @return A list containing the following elements:
#' \describe{
#'   \item{cut.values}{Unique marker values for calculating TP and FP.}
#'   \item{TP}{True Positive corresponding to the cut-offs in the marker.}
#'   \item{FP}{False Positive corresponding to the cut-offs in the marker.}
#'   \item{predict.time}{Time point of interest.}
#'   \item{Survival}{Kaplan-Meier survival estimate at \code{predict.time}.}
#'   \item{AUC}{Area Under the ROC Curve at \code{predict.time}.}
#' }
#'
#' @examples
#' library(survival)
#' library(PAmeasures)
#' # Use Mayo Clinic Primary Biliary Cirrhosis Data
#' data(pbc)
#' pbc <- pbc %>% 
#'   filter(is.na(trt)==F) %>% 
#'   mutate(log_albumin = log(albumin),
#'          log_bili = log(bili),
#'          log_protime = log(protime),
#'          status = ifelse(status==2, 1, 0))
#' 
#' time_dep_auc_full <- survivalROC(Stime = pbc$time,
#'             status = pbc$status,
#'             marker = predict(fit.coxph.full, newdata=pbc, type = "lp") ,
#'             predict.time = quantile(pbc$time, 0.5), 
#'             method="KM")$AUC
#' @references
#' Heagerty, P. J., Lumley, T., & Pepe, M. S. (2000). Time-dependent ROC curves for censored
#' survival data and a diagnostic marker. \emph{Biometrics}, 56(2), 337-344.
#'
#' Heagerty, P. J., Saha-Chaudhuri, P. (2022). survivalROC: Time-Dependent ROC Curve Estimation from Censored Survival Data. 
#' R package version 1.0.3.1. DOI: \doi{10.32614/CRAN.package.survivalROC}. Available at \url{https://CRAN.R-project.org/package=survivalROC}.
#' @keywords internal
#' @noRd

pam.survivalROC <- function (Stime, status, marker, entry = NULL, predict.time, 
                             cut.values = NULL, method = "NNE", lambda = NULL, span = NULL, 
                             window = "symmetric") 
{
  times = Stime
  x <- marker
  if (is.null(entry)) 
    entry <- rep(0, length(times))
  bad <- is.na(times) | is.na(status) | is.na(x) | is.na(entry)
  entry <- entry[!bad]
  times <- times[!bad]
  status <- status[!bad]
  x <- x[!bad]
  if (sum(bad) > 0) 
    cat(paste("\n", sum(bad), "records with missing values dropped. \n"))
  if (is.null(cut.values)) 
    cut.values <- unique(x)
  cut.values <- cut.values[order(cut.values)]
  ncuts <- length(cut.values)
  ooo <- order(times)
  times <- times[ooo]
  status <- status[ooo]
  x <- x[ooo]
  s0 <- 1
  unique.t0 <- unique(times)
  unique.t0 <- unique.t0[order(unique.t0)]
  n.times <- sum(unique.t0 <= predict.time)
  for (j in 1:n.times) {
    n <- sum(entry <= unique.t0[j] & times >= unique.t0[j])
    d <- sum((entry <= unique.t0[j]) & (times == unique.t0[j]) & 
               (status == 1))
    if (n > 0) 
      s0 <- s0 * (1 - d/n)
  }
  s.pooled <- s0
  roc.matrix <- matrix(NA, ncuts, 2)
  roc.matrix[ncuts, 1] <- 0
  roc.matrix[ncuts, 2] <- 1
  if (method == "KM") {
    for (c in 1:(ncuts - 1)) {
      s0 <- 1
      subset <- as.logical(x > cut.values[c])
      e0 <- entry[subset]
      t0 <- times[subset]
      c0 <- status[subset]
      if (!is.null(t0)) {
        unique.t0 <- unique(t0)
        unique.t0 <- unique.t0[order(unique.t0)]
        n.times <- sum(unique.t0 <= predict.time)
        if (n.times > 0) {
          for (j in 1:n.times) {
            n <- sum(e0 <= unique.t0[j] & t0 >= unique.t0[j])
            d <- sum((e0 <= unique.t0[j]) & (t0 == unique.t0[j]) & 
                       (c0 == 1))
            if (n > 0) 
              s0 <- s0 * (1 - d/n)
          }
        }
      }
      p0 <- mean(subset)
      roc.matrix[c, 1] <- (1 - s0) * p0/(1 - s.pooled)
      roc.matrix[c, 2] <- 1 - s0 * p0/s.pooled
    }
  }
  if (method == "NNE") {
    if (is.null(lambda) & is.null(span)) {
      cat("method = NNE requires either lambda or span! \n")
      stop(0)
    }
    x.unique <- unique(x)
    x.unique <- x.unique[order(x.unique)]
    S.t.x <- rep(0, length(x.unique))
    t.evaluate <- unique(times[status == 1])
    t.evaluate <- t.evaluate[order(t.evaluate)]
    t.evaluate <- t.evaluate[t.evaluate <= predict.time]
    for (j in 1:length(x.unique)) {
      if (!is.null(span)) {
        if (window == "symmetric") {
          ddd <- (x - x.unique[j])
          n <- length(x)
          ddd <- ddd[order(ddd)]
          index0 <- sum(ddd < 0) + 1
          index1 <- index0 + trunc(n * span + 0.5)
          if (index1 > n) 
            index1 <- n
          lambda <- ddd[index1]
          wt <- as.integer(((x - x.unique[j]) <= lambda) & 
                             ((x - x.unique[j]) >= 0))
          index0 <- sum(ddd <= 0)
          index2 <- index0 - trunc(n * span/2)
          if (index2 < 1) 
            index2 <- 1
          lambda <- abs(ddd[index1])
          set.index <- ((x - x.unique[j]) >= -lambda) & 
            ((x - x.unique[j]) <= 0)
          wt[set.index] <- 1
        }
        if (window == "asymmetric") {
          ddd <- (x - x.unique[j])
          n <- length(x)
          ddd <- ddd[order(ddd)]
          index0 <- sum(ddd < 0) + 1
          index <- index0 + trunc(n * span)
          if (index > n) 
            index <- n
          lambda <- ddd[index]
          wt <- as.integer(((x - x.unique[j]) <= lambda) & 
                             ((x - x.unique[j]) >= 0))
        }
      }
      else {
        wt <- exp(-(x - x.unique[j])^2/lambda^2)
      }
      s0 <- 1
      for (k in 1:length(t.evaluate)) {
        n <- sum(wt * (entry <= t.evaluate[k]) & (times >= 
                                                    t.evaluate[k]))
        d <- sum(wt * (entry <= t.evaluate[k]) & (times == 
                                                    t.evaluate[k]) * (status == 1))
        if (n > 0) 
          s0 <- s0 * (1 - d/n)
      }
      S.t.x[j] <- s0
    }
    S.all.x <- S.t.x[match(x, x.unique)]
    n <- length(times)
    S.marginal <- sum(S.all.x)/n
    for (c in 1:(ncuts - 1)) {
      p1 <- sum(x > cut.values[c])/n
      Sx <- sum(S.all.x[x > cut.values[c]])/n
      roc.matrix[c, 1] <- (p1 - Sx)/(1 - S.marginal)
      roc.matrix[c, 2] <- 1 - Sx/S.marginal
    }
  }
  sensitivity = roc.matrix[, 1]
  specificity = roc.matrix[, 2]
  x <- 1 - c(0, specificity)
  y <- c(1, sensitivity)
  n <- length(x)
  dx <- x[-n] - x[-1]
  mid.y <- (y[-n] + y[-1])/2
  area <- sum(dx * mid.y)
  list(cut.values = c(-Inf, cut.values), TP = y, FP = x, predict.time = predict.time, 
       Survival = s.pooled, AUC = area)
}
