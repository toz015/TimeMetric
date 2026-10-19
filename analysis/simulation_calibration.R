# Simulation ICI uses the same grid convention as tm_survival_eval().
# Keep its unrounded estimate for replication means, errors and MCSEs.
paper_simulation_calibration <- function(pred, t_star, calibration_args) {
  grid_time <- NA_real_
  method <- NA_character_
  error <- ""
  ici <- tryCatch({
    grid <- sort(pred$times)
    grid_risk <- getFromNamespace("tm_calibration_grid_risk", "TimeMetric")
    p <- grid_risk(pred$surv_prob, grid, t_star, length(pred$times),
                   survival = TRUE)
    at <- findInterval(t_star, grid)
    grid_time <- if (at == 0L) 0 else grid[at]
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
