# Deterministic fixtures shared by all characterization tests.
# Seeds are fixed; never change them without regenerating every snapshot.

# TimeMetric calls survival::concordancefit unqualified without importing it
# (see findings.md #11), so the package only works when survival is attached --
# which is what paper.code.Rmd and every real user session does. Attaching it
# here characterizes the package as users actually experience it. The missing
# import is asserted separately in test-eval-survival.R.
library(survival)

fx_covs <- function() c("x1", "x2")

fx_surv <- function() {
  d <- sim_cox_weibull_censored(
    n = 200, pi_c = 0.3, v = 2, beta = c(0.5, -0.5), seed = 1001
  )
  # pi_c > 0 also returns y_true and cens_time; drop them so model
  # formulas built with "." cannot pick them up as covariates.
  d[, c("time", "status", "x1", "x2")]
}

fx_surv_uncensored <- function() {
  sim_cox_weibull_censored(
    n = 200, pi_c = 0, v = 2, beta = c(0.5, -0.5), seed = 1002
  )[, c("time", "status", "x1", "x2")]
}

fx_cr <- function() {
  simulateTwoCauseFineGrayModel(
    n = 200, v = 2, beta1 = c(0.5, -0.5), beta2 = c(-0.3, 0.3),
    censor = 0.3, seed = 1003
  )
}

fx_cox <- function() {
  survival::coxph(
    survival::Surv(time, status) ~ x1 + x2,
    data = fx_surv(), x = TRUE, y = TRUE
  )
}

fx_survreg <- function() {
  # x = TRUE, y = TRUE is required: pam.surverg_restricted errors without it
  survival::survreg(
    survival::Surv(time, status) ~ x1 + x2,
    data = fx_surv(), dist = "weibull", x = TRUE, y = TRUE
  )
}
