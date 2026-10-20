test_that('competing-risk AUC and Brier work without attaching survival', {
  # Namespace-only calls must work; library(survival) would hide the upstream
  # tdROC formula lookup defect with survival < 3.8.
  expect_false('package:survival' %in% search())
  withr::local_seed(436)
  n <- 120L
  event_time <- rexp(n)
  status <- sample(0:2, n, replace = TRUE)
  times <- c(0.25, 0.5, 1)
  risks <- runif(n, 0.05, 0.7)
  cif <- outer(times, risks, function(t, p) 1 - (1 - p)^t)
  result <- tm_survival_eval_cr(cif, event_time, times, status,
    metrics = c('td_auc', 'brier_score'), t_star = 0.5)
  expect_setequal(result$Metric, c('td_auc', 'brier_score'))
  expect_true(all(is.finite(result$Value)))
  expect_true(all(result$Value >= 0 & result$Value <= 1))
  expect_false('package:survival' %in% search())
})
