# Run after installing survival 3.8-3 and TimeMetric; never attach survival.
stopifnot(packageVersion('survival') == '3.8.3')
library(TimeMetric)
stopifnot(!'package:survival' %in% search())
set.seed(436)
n <- 120L
event_time <- rexp(n)
status <- sample(0:2, n, replace = TRUE)
times <- c(0.25, 0.5, 1)
risks <- runif(n, 0.05, 0.7)
cif <- outer(times, risks, function(t, p) 1 - (1 - p)^t)
result <- tm_survival_eval_cr(cif, event_time, times, status,
  metrics = c('td_auc', 'brier_score'), t_star = 0.5)
stopifnot(setequal(result$Metric, c('td_auc', 'brier_score')),
          all(is.finite(result$Value)),
          all(result$Value >= 0 & result$Value <= 1),
          !'package:survival' %in% search())
print(result)
cat('PASS: competing-risk AUC and Brier with survival 3.8-3, namespace only.\n')
