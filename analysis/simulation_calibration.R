# Evaluate simulation calibration at the requested outcome horizon, using the
# original fitted model rather than a grid defined by validation follow-up.
paper_simulation_risk_at <- function(model, new_data, t_star) {
  if (!is.numeric(t_star) || length(t_star) != 1L || !is.finite(t_star) ||
      t_star <= 0) stop("Calibration time must be one positive finite value.")
  if (inherits(model, "coxph")) {
    lp <- as.numeric(stats::predict(model, newdata = new_data,
      type = "lp", reference = "sample"))
    survival_at <- getFromNamespace("tm_cox_survival", "TimeMetric")
    survival_probability <- survival_at(model, new_data, lp, t_star)
  } else if (inherits(model, "survreg")) {
    survival_at <- getFromNamespace("predictSurvProb2survreg", "TimeMetric")
    survival_probability <- survival_at(model, new_data, t_star)
  } else {
    survival_probability <- pec::predictSurvProb(model, newdata = new_data,
                                                times = t_star)
  }
  p <- 1 - as.numeric(survival_probability)
  if (length(p) != nrow(new_data) || any(!is.finite(p)) ||
      any(p < 0 | p > 1)) stop("Invalid exact-time simulation probabilities.")
  p
}

# Keep its unrounded estimate for replication means, errors and MCSEs. Model
# adapters sort Cox/AFT subjects by observed time; RF retains input row order.
paper_simulation_calibration <- function(pred, t_star, calibration_args,
                                          model, new_data) {
  grid_time <- NA_real_
  method <- NA_character_
  error <- ""
  ici <- tryCatch({
    if (inherits(model, c("coxph", "survreg")))
      new_data <- new_data[order(new_data$time), , drop = FALSE]
    if (!identical(as.numeric(new_data$time), as.numeric(pred$times)) ||
        !identical(as.numeric(new_data$status), as.numeric(pred$status)))
      stop("Calibration subjects and outcomes are not aligned.")
    p <- paper_simulation_risk_at(model, new_data, t_star)
    # Retain the existing diagnostic column name, now always equal to t_star.
    grid_time <- t_star
    cal <- do.call(tm_calibration, c(list(p = p, event_time = pred$times,
      status = pred$status, t_star = t_star), calibration_args))
    method <- cal$method
    unname(cal$metrics[["ici"]])
  }, error = function(e) {
    error <<- conditionMessage(e)
    warning("Simulation ICI unavailable: ", error, call. = FALSE)
    NA_real_
  })
  list(ici = ici, diagnostics = data.frame(
    t_star = t_star, prediction_grid_time = grid_time,
    requested_knots = calibration_args$nk, method = method,
    subjects = length(pred$times), observed_events = sum(pred$status == 1),
    events_by_horizon = sum(pred$status == 1 & pred$times <= t_star),
    at_risk = sum(pred$times >= t_star), ici = ici, error = error))
}

# Missing ICI rows must never be accepted as a completed old checkpoint.
paper_simulation_record_complete <- function(record, metrics, finite = TRUE) {
  if (!isTRUE(record$ok) || !is.data.frame(record$metrics) ||
      !all(c("Model", "Metric", "Value") %in% names(record$metrics))) return(FALSE)
  models <- switch(record$job$part,
    sim1 = "Cox", sim3 = "Cox",
    sim2_linear = c("Cox", "Cox.reduce", "Ln", "Ln.fix"),
    sim2_nonlinear = c("Cox", "Cox.reduce", "Ln", "Ln.fix", "RF"),
    character())
  d <- record$metrics
  if (!length(models) || anyNA(d[c("Model", "Metric")]) ||
      !setequal(d$Model, models)) return(FALSE)
  complete <- all(vapply(models, function(model) {
    x <- d$Metric[d$Model == model]
    length(x) == length(metrics) && setequal(x, metrics)
  }, logical(1)))
  complete && (!finite || all(is.finite(d$Value)))
}
