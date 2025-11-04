#' Simulate Cox survival data with targeted (or fixed) censoring
#'
#' @description
#' Generates survival data under a Cox proportional hazards model with a Weibull
#' baseline cumulative hazard \eqn{H_0(t) = (0.5 t)^v}. Event times are sampled via
#' inverse-transform:
#' \deqn{Y = H_0^{-1}(-\log U \cdot \exp(-\beta^T X)), \quad U \sim \mathrm{Unif}(0,1),}
#' where \eqn{H_0^{-1}(s) = 2 s^{1/v}} and covariates \eqn{X_j = 10 \cdot \mathrm{Bernoulli}(0.5)}.
#'
#' Right censoring times are generated as \eqn{C_i = \exp(\mu + r_i)} with
#' \eqn{r_i \sim \mathcal{N}(0, \mathrm{sd}^2)}. By default, \code{sd} is taken to be
#' \eqn{\mathrm{sd}(\log Y)} and \code{mu} is calibrated so that the expected censoring
#' proportion equals \code{pi_c}. If \code{mu} and/or \code{sd} are supplied, calibration
#' is skipped appropriately:
#' \itemize{
#'   \item \strong{mu and sd provided:} use them as-is (no solving).
#'   \item \strong{sd provided only:} solve for \code{mu} using this \code{sd}.
#'   \item \strong{mu provided only:} estimate \eqn{\mathrm{sd}(\log Y)} and use the given \code{mu}.
#'   \item \strong{neither provided:} estimate \eqn{\mathrm{sd}(\log Y)} and solve for \code{mu}.
#' }
#'
#' @param n Integer. Number of individuals.
#' @param pi_c Numeric in \eqn{[0, 0.99]}. Target censoring proportion. Ignored if both
#'   \code{mu} and \code{sd} are provided (since no calibration is done).
#' @param v Positive numeric. Weibull shape in \eqn{H_0(t) = (0.5 t)^v}.
#' @param beta Numeric vector of regression coefficients; length defines \eqn{p}.
#'   Each covariate is \eqn{X_j = 10 \cdot \mathrm{Bernoulli}(0.5)}.
#' @param mu Optional numeric. Mean of \eqn{\log C}. If supplied with \code{sd}, no
#'   calibration is performed.
#' @param sd Optional positive numeric. Standard deviation used for \eqn{r_i} in
#'   \eqn{\log C = \mu + r_i}. If \code{NULL}, defaults to \eqn{\mathrm{sd}(\log Y)}.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A \code{data.frame} with columns:
#' \itemize{
#'   \item \code{time}: observed time \eqn{T_i = \min(Y_i, C_i)}.
#'   \item \code{status}: event indicator \eqn{1\{Y_i \le C_i\}}.
#'   \item \code{x1, x2, ...}: generated covariates.
#'   \item \code{y_true}: latent event time \eqn{Y_i}.
#'   \item \code{cens_time}: censoring time \eqn{C_i}.
#' }
#' Attributes include:
#' \itemize{
#'   \item \code{mu}: censoring shift used.
#'   \item \code{sd_log}: sd used for \eqn{r_i} in \eqn{\log C}.
#'   \item \code{sd_logY}: \eqn{\mathrm{sd}(\log Y)} from the generated events.
#'   \item \code{pi_c_target}: requested censoring rate.
#'   \item \code{pi_c_observed}: realized censoring rate in the sample.
#'   \item \code{v}, \code{beta}: model parameters.
#' }
#'
#' @examples
#' # 1) Default: estimate sd=sd(logY), solve for mu to hit pi_c
#' set.seed(1)
#' d1 <- sim_cox_weibull_censored(n = 500, pi_c = 0.30, v = 1.2, beta = c(0.4, -0.2))
#' attr(d1, "mu"); attr(d1, "sd_log"); attr(d1, "pi_c_observed")
#'
#' # 2) Provide sd only: solve for mu using this sd
#' d2 <- sim_cox_weibull_censored(n = 500, pi_c = 0.30, v = 1.2, beta = 0.5, sd = 0.8)
#'
#' # 3) Provide mu only: use mu, estimate sd = sd(logY)
#' d3 <- sim_cox_weibull_censored(n = 500, pi_c = 0.30, v = 1.2, beta = 0.5, mu = 1.5)
#'
#' # 4) Provide both mu and sd: skip calibration entirely
#' d4 <- sim_cox_weibull_censored(n = 500, pi_c = 0.30, v = 1.2, beta = 0.5, mu = 1.5, sd = 0.7)
#'
#' @export
sim_cox_weibull_censored <- function(n, pi_c, v, beta, mu = NULL, sd = NULL, seed = NULL) {
  stopifnot(length(n) == 1L, n > 0, is.finite(n))
  stopifnot(length(v) == 1L, v > 0, is.finite(v))
  stopifnot(length(pi_c) == 1L, is.finite(pi_c), pi_c >= 0, pi_c < 1)
  stopifnot(is.numeric(beta), length(beta) >= 1L, all(is.finite(beta)))
  if (!is.null(mu)) stopifnot(length(mu) == 1L, is.finite(mu))
  if (!is.null(sd)) stopifnot(length(sd) == 1L, is.finite(sd), sd > 0)
  
  if (!is.null(seed)) set.seed(seed)
  
  p <- length(beta)
  
  # Generate covariates: X_j = 10 * Bernoulli(0.5)
  # X <- matrix(rbinom(n * p, size = 1, prob = 0.5), nrow = n, ncol = p)
  # Use normal covariates:
  X <- matrix(rnorm(n * p, mean = 0, sd = 1), nrow = n, ncol = p)
  colnames(X) <- paste0("x", seq_len(p))
  
  # Event times: Y = H0^{-1}(-log U * exp(-eta)), H0^{-1}(s) = 2 * s^{1/v}
  U <- stats::runif(n)
  eta <- as.numeric(X %*% beta)
  s <- -log(U) * exp(-eta)
  Y <- 2 * (s)^(1 / v)
  
  logY <- log(Y)
  sd_logY <- stats::sd(logY)
  if (!is.finite(sd_logY) || sd_logY <= .Machine$double.eps) sd_logY <- 1e-6
  
  if(pi_c == 0){
    out <- data.frame(
      time = Y,
      status = 1,
      X
    )
    return(out)
  }
  # Determine sd to use for log-censoring noise
  sd_use <- if (is.null(sd)) sd_logY else sd
  
  # Determine mu to use
  mu_use <- mu
  if (is.null(mu_use)) {
    # Solve for mu so that mean(Phi((logY - mu)/sd_use)) = 1 - pi_c
    target_obs <- 1 - pi_c
    f_mu <- function(m) {
      mean(stats::pnorm((m - logY) / sd_use)) - target_obs
    }
    lo <- min(logY) - 10 * sd_use - 5
    hi <- max(logY) + 10 * sd_use + 5
    f_lo <- f_mu(lo); f_hi <- f_mu(hi)
    expand <- 0
    while (f_lo * f_hi > 0 && expand < 5) {
      lo <- lo - 10; hi <- hi + 10
      f_lo <- f_mu(lo); f_hi <- f_mu(hi)
      expand <- expand + 1
    }
    mu_use <- tryCatch(
      stats::uniroot(f_mu, lower = lo, upper = hi)$root,
      error = function(e) stats::quantile(logY, probs = target_obs, names = FALSE)
    )
  }
  
  # Generate censoring and observed
  r <- stats::rnorm(n, mean = 0, sd = sd_use)
  C <- exp(mu_use + r)
  
  time <- pmin(Y, C)
  status <- as.integer(Y <= C)
  pi_c_obs <- mean(status == 0)
  
  out <- data.frame(
    time = time,
    status = status,
    X,
    y_true = Y,
    cens_time = C,
    check.names = FALSE
  )
  
  attr(out, "mu") <- as.numeric(mu_use)
  attr(out, "sd_log") <- as.numeric(sd_use)
  attr(out, "sd_logY") <- as.numeric(sd_logY)
  attr(out, "pi_c_target") <- as.numeric(pi_c)
  attr(out, "pi_c_observed") <- as.numeric(pi_c_obs)
  attr(out, "v") <- as.numeric(v)
  attr(out, "beta") <- as.numeric(beta)
  
  out
}

