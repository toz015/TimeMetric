# Validation of r_sh (Schemper-Henderson) against an INDEPENDENT implementation
# written from the published definition, not against a second TimeMetric code
# path. See findings.md #35 for the defect this guards against: the baseline
# survival curve was read from an rms::cph object fitted without surv = TRUE,
# so it collapsed to a two-valued step function instead of the fitted hazard.
#
# Schemper, M. & Henderson, R. (2000). Predictive accuracy and explained
# variation in Cox regression. Biometrics 56, 249-255.

# With no censoring the inverse-probability-of-censoring machinery drops out and
# the estimator reduces to the mean absolute deviation between the 0/1 survival
# indicator and the predicted survival probability, averaged over event times:
#
#   M(t)  = (1/n) sum_i | I(T_i > t) - S_i(t) |
#         = (1/n) sum_i [ I(T_i > t)(1 - S_i(t)) + I(T_i <= t) S_i(t) ]
#   D     = mean_j M(t_j) using the marginal Kaplan-Meier for every subject
#   D_x   = mean_j M(t_j) using each subject's covariate-specific survival
#   r_sh  = (D - D_x) / D
sh_reference <- function(d) {
  fit <- survival::coxph(survival::Surv(time, status) ~ .,
                         data = d, x = TRUE, y = TRUE)
  tj <- sort(unique(d$time[d$status == 1]))
  km <- survival::survfit(survival::Surv(time, status) ~ 1, data = d)

  s_null <- stats::approx(km$time, km$surv, xout = tj, method = "constant",
                          f = 0, yleft = 1, rule = 2)$y
  s_mod <- summary(survival::survfit(fit, newdata = d),
                   times = tj, extend = TRUE)$surv

  alive <- outer(tj, d$time, function(a, b) as.numeric(b > a))
  M <- function(S) rowMeans(alive * (1 - S) + (1 - alive) * S)

  D  <- mean(M(matrix(s_null, nrow = length(tj), ncol = nrow(d))))
  Dx <- mean(M(s_mod))
  (D - Dx) / D
}

# r_sh as the package computes it
tm_rsh <- function(d) {
  m <- survival::coxph(survival::Surv(time, status) ~ x1 + x2,
                       data = d, x = TRUE, y = TRUE)
  p <- tm_predict_coxph(model = m, covs = c("x1", "x2"),
                        new_data = d, tau = 10e10)
  r <- tm_survival_eval(model = m, event_time = p$times,
                        predicted_probability = p$surv_prob, status = p$status,
                        covariates = c("x1", "x2"), new_data = d, tau = 10e10)
  r$Value[r$Metric == "r_sh"]
}

sim <- function(beta, pi_c, seed, n = 200) {
  tm_sim_cox_weibull(n = n, pi_c = pi_c, v = 2, beta = beta,
                     seed = seed)[, c("time", "status", "x1", "x2")]
}

test_that("r_sh matches an independent implementation of the SH definition", {
  # Uncensored, where the definition is directly computable without IPCW.
  for (sc in list(list(b = c(1.5, -1.5), seed = 7001, label = "strong predictor"),
                  list(b = c(0.3, -0.3), seed = 7002, label = "weak predictor"),
                  list(b = c(0, 0),      seed = 7003, label = "null model"))) {
    d <- sim(sc$b, pi_c = 0, seed = sc$seed, n = 150)
    expect_equal(tm_rsh(d), sh_reference(d),
                 tolerance = 1e-8, info = sc$label)
  }
})

test_that("r_sh is approximately zero for a null model", {
  # Definitional: with no covariate effect the model survival equals the
  # marginal Kaplan-Meier, so D_x = D and the explained variation vanishes.
  for (pi_c in c(0, 0.2, 0.4)) {
    v <- tm_rsh(sim(c(0, 0), pi_c, seed = 8001))
    expect_lt(abs(v), 0.05)
    expect_gt(v, -0.05)
  }
})

test_that("r_sh increases strictly with predictor strength", {
  vals <- vapply(c(0, 0.25, 0.5, 1.0, 2.0),
                 function(b) tm_rsh(sim(c(b, -b), 0.3, seed = 8004)),
                 numeric(1))

  expect_true(all(diff(vals) > 0))
  expect_lt(vals[1], 0.05)   # null end
  expect_gt(vals[5], 0.4)    # strong end
})

test_that("r_sh stays within [0, 1] across censoring levels", {
  # The pre-correction implementation produced negative values (-0.026 on the
  # standard fixture), which is impossible for an explained-variation measure.
  for (pi_c in c(0, 0.2, 0.4, 0.6)) {
    v <- tm_rsh(sim(c(1, -1), pi_c, seed = 9001))
    expect_gte(v, -1e-6)
    expect_lte(v, 1)
  }
})

test_that("r_sh is computable with heavily tied event times", {
  d <- sim(c(1, -1), 0.3, seed = 8003)
  d$time <- round(d$time, 1)   # force many ties

  expect_lt(length(unique(d$time)), nrow(d) * 0.5)   # ties really present
  v <- suppressWarnings(tm_rsh(d))
  expect_true(is.finite(v))
  expect_gte(v, -1e-6)
  expect_lte(v, 1)
})

test_that("r_sh no longer depends on rms", {
  # rms was used only to fit a model for pam.schemper and to extract linear
  # predictors. Both are now survival:: calls, so no rms reference remains.
  root <- skip_without_source_tree()
  src <- unlist(lapply(list.files(file.path(root, "R"), full.names = TRUE),
                       readLines, warn = FALSE))
  code <- grep("^\\s*#", src, value = TRUE, invert = TRUE)

  expect_false(any(grepl("rms::", code, fixed = TRUE)))
})
