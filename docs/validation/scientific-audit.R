# Run from the repository root: Rscript docs/validation/scientific-audit.R
# This is a diagnostic audit, not a snapshot of acceptable behaviour.
# It exits nonzero while independently checked scientific contracts fail.
# The ordinary testthat suite remains separate; do not mark these as skips.
pkgload::load_all(".", quiet = TRUE)
dir.create("docs/validation/results", recursive = TRUE, showWarnings = FALSE)
checks <- list()
record <- function(id, actual, expected, tolerance = 1e-8) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = id, actual = as.numeric(actual), expected = as.numeric(expected),
    pass = isTRUE(all.equal(as.numeric(actual), as.numeric(expected),
                            tolerance = tolerance)))
}

train <- tm_sim_cox_weibull(160, pi_c = 0.3, v = 2,
                            beta = c(0.5, -0.5), seed = 5101)
test <- tm_sim_cox_weibull(80, pi_c = 0.3, v = 2,
                           beta = c(0.5, -0.5), seed = 5102)
train <- train[, c("time", "status", "x1", "x2")]
test <- test[, c("time", "status", "x1", "x2")]
covs <- c("x1", "x2")
cox <- survival::coxph(survival::Surv(time, status) ~ x1 + x2,
                       train, x = TRUE, y = TRUE)
p <- tm_predict_coxph(cox, covs, new_data = test)
changed <- test
changed$status <- 1 - changed$status
p_changed <- tm_predict_coxph(cox, covs, new_data = changed)
record("Cox predictions invariant to validation event labels",
       max(abs(p$surv_prob - p_changed$surv_prob)), 0)
sf <- survival::survfit(cox, newdata = test[order(test$time), ], se.fit = FALSE)
reference <- t(summary(sf, times = sort(test$time), extend = TRUE)$surv)
record("Cox probabilities agree with fitted-model survfit predictions",
       max(abs(p$surv_prob - reference)), 0)

tau <- unname(stats::median(test$time))
for (distribution in c("weibull", "lognormal")) {
  fit <- survival::survreg(survival::Surv(time, status) ~ x1 + x2,
                           train, x = TRUE, y = TRUE, dist = distribution)
  aft <- tm_predict_survreg(fit, covs, tau = tau, new_data = test)
  lp <- stats::predict(fit, test[order(test$time), ], type = "lp")
  ref <- vapply(lp, function(eta) {
    stats::integrate(function(t) {
      if (distribution == "weibull") {
        exp(-(t / exp(eta))^(1 / fit$scale))
      } else {
        stats::pnorm((log(t) - eta) / fit$scale, lower.tail = FALSE)
      }
    }, lower = 0, upper = tau)$value
  }, numeric(1))
  record(paste(distribution, "restricted means equal integral of survival"),
         max(abs(aft$pred - ref)), 0, tolerance = 1e-6)
  record(paste(distribution, "prediction preserves evaluation event labels"),
         sum(aft$status != test$status[order(test$time)]), 0)
}

# One case at time 1, two eligible controls: with one control per case,
# each control has inclusion probability 1/2 and inverse weight 2.
w <- tm_nested_case_control_weights(c(1, 2, 3), c(1, 0, 0), m = 1)
record("NCC one control per case: control inverse probability weight", w[2], 2)
w2 <- tm_nested_case_control_weights(c(1, 2, 3), c(1, 0, 0), m = 2)
record("NCC sampling all eligible controls: control weight", w2[2], 1)

# Hand-computable uncensored explained variation must not depend on row order.
y <- c(1, 2, 4, 8, 10)
predicted <- c(1.2, 2.3, 3.8, 7.5, 9.7)
idx <- c(5, 3, 1, 4, 2)
a <- TimeMetric:::pam.r2_metrics(predicted, y, rep(1, length(y)))
b <- TimeMetric:::pam.r2_metrics(predicted[idx], y[idx], rep(1, length(y)))
record("Explained variation is invariant to joint subject permutation",
       b$Pseudo_R_squared, a$Pseudo_R_squared)

d <- tm_sim_cox_weibull(120, 0, 2, c(1, -1), seed = 9191)
d <- d[, c("time", "status", "x1", "x2")]
fit <- survival::coxph(survival::Surv(time, status) ~ x1 + x2,
                       d, x = TRUE, y = TRUE)
pr <- tm_predict_coxph(fit, covs, new_data = d)
t_star <- pr$times[60]
event <- as.numeric(pr$times <= t_star)
risk <- 1 - pr$surv_prob[, 60]
calibration <- stats::lm(event ~ risk)
r2 <- round(summary(calibration)$r.squared, 4)
l2 <- round(sum(stats::residuals(calibration)^2) / sum((event - risk)^2), 4)
point <- tm_survival_eval(fit, pr$times, pr$surv_prob, status = pr$status,
  covariates = covs, t_star = t_star, metrics = "pseudo_r2_point")
record("Point pseudo R-squared agrees with event-probability regression",
       point$Value, round(r2 * l2, 4))

validation <- d
validation$x1 <- -validation$x1
validation$x2 <- -validation$x2
t <- stats::median(validation$time)
score <- stats::predict(fit, validation, type = "lp")
cases <- validation$time <= t
comparisons <- outer(score[cases], score[!cases], "-")
reference_auc <- mean((comparisons > 0) + 0.5 * (comparisons == 0))
auc <- tm_fit_and_eval(d, covs, predicted_data = validation,
                       metrics = "td_auc", t_star = t)
record("Reversed held-out AUC agrees with empirical case/control ranking",
       auc$td_auc, round(reference_auc, 2))

out <- do.call(rbind, checks)
utils::write.csv(out, "docs/validation/results/scientific-audit.csv", row.names = FALSE)
print(out, row.names = FALSE)
writeLines(capture.output(sessionInfo()), "docs/validation/results/session-info.txt")
if (any(!out$pass)) quit(status = 1L)
