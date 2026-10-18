test_that("held-out outcomes cannot change fixed-grid Cox probabilities", {
  fit <- fx_cox()
  d <- fx_surv_uncensored()
  changed <- d
  changed$status[seq(1, nrow(d), 2)] <- 0
  a <- tm_predict_coxph(fit, fx_covs(), new_data = d, tau = 2)
  b <- tm_predict_coxph(fit, fx_covs(), new_data = changed, tau = 2)
  expect_identical(a$surv_prob, b$surv_prob)
  expect_identical(a$pred, b$pred)
  ref <- survival::survfit(fit, newdata = d[order(d$time), ], se.fit = FALSE)
  expected <- t(summary(ref, times = sort(d$time), extend = TRUE)$surv)
  expect_equal(unname(a$surv_prob), unname(expected), tolerance = 1e-12)
  # Changing observed times only moves the output grid, not the fitted curve.
  changed$time[which.max(changed$time)] <- max(d$time) + 10
  c <- tm_predict_coxph(fit, fx_covs(), new_data = changed, tau = 2)
  common <- intersect(a$times, c$times)
  expect_equal(a$surv_prob[, match(common, a$times)],
               c$surv_prob[, match(common, c$times)], tolerance = 1e-12)
})

test_that("training weights, ties, factors, and strata match survfit", {
  d <- fx_surv()
  d$time <- ceiling(d$time * 3)
  d$g <- factor(rep(c("A", "B"), length.out = nrow(d)))
  d$w <- rep(c(1, 2, 4), length.out = nrow(d))
  for (formula in list(survival::Surv(time, status) ~ x1 + g,
                      survival::Surv(time, status) ~ x1 + survival::strata(g))) {
    fit <- survival::coxph(formula, data = d, weights = w, ties = "efron",
                            x = TRUE, y = TRUE, model = TRUE)
    new <- d[c(41, 6, 78, 3, 22, 5), ]
    a <- tm_predict_coxph(fit, c("x1", "g"), new_data = new)
    ordered <- new[order(new$time), ]
    expected <- t(vapply(seq_len(nrow(ordered)), function(i) {
      ref <- survival::survfit(fit, newdata = ordered[i, , drop = FALSE], se.fit = FALSE)
      as.numeric(summary(ref, times = sort(new$time), extend = TRUE)$surv)
    }, numeric(nrow(ordered))))
    expect_equal(unname(a$surv_prob), unname(expected), tolerance = 1e-12)
    changed <- new; changed$status <- 1 - changed$status
    b <- tm_predict_coxph(fit, c("x1", "g"), new_data = changed)
    expect_identical(a$surv_prob, b$surv_prob)
    in_sample <- tm_predict_coxph(fit, c("x1", "g"))
    explicit <- tm_predict_coxph(fit, c("x1", "g"), new_data = d)
    expect_equal(in_sample$surv_prob, explicit$surv_prob, tolerance = 1e-12)
  }
})

test_that("one-covariate and zero-event evaluation sets remain predictable", {
  d <- fx_surv()
  fit <- survival::coxph(survival::Surv(time, status) ~ x1, d, x = TRUE, y = TRUE)
  d$status <- 0
  p <- tm_predict_coxph(fit, "x1", new_data = d)
  expect_true(all(is.finite(p$surv_prob)))
  expect_equal(dim(p$surv_prob), c(nrow(d), nrow(d)))
})

test_that("matched NCC excludes other-stratum sampling occasions", {
  time <- c(1, 4, 5, 2, 4, 5)
  status <- c(1, 0, 0, 1, 0, 0)
  strata <- c("A", "A", "A", "B", "B", "B")
  w <- tm_nested_case_control_weights(time, status, strata, m = 2)
  expect_equal(w, c(1, 2, 2, 1, 2, 2))
  o <- c(6, 2, 4, 1, 5, 3)
  expect_equal(tm_nested_case_control_weights(time[o], status[o], strata[o], 2), w[o])
  # Equal-time events in separate strata are separate valid occasions.
  time[c(1, 4)] <- 1
  expect_equal(tm_nested_case_control_weights(time, status, strata, 2), w)
})

test_that("NCC accounts for repeated occasions and small risk sets", {
  # Two cases at t=1,2. A non-case survives both draws: p=1-(1-1/3)(1-1/2).
  expect_equal(tm_nested_case_control_weights(c(1, 2, 4, 5), c(1, 1, 0, 0), m = 2),
               c(1, 1, 1.5, 1.5))
  expect_equal(tm_nested_case_control_weights(c(1, 4, 5), c(1, 0, 0), m = 9), c(1, 1, 1))
  expect_equal(tm_nested_case_control_weights(c(1, 1, 4, 5), c(1, 1, 0, 0),
                                               strata = rep("A", 4), m = 2), c(1, 1, 1.8, 1.8))
  expect_error(tm_nested_case_control_weights(c(1, 2), c(1, 0), m = 1), "integer >= 2")
  expect_error(tm_nested_case_control_weights(c(1, 2), c(1, 0), strata = "A", m = 2),
                "one non-missing")
})
