#!/usr/bin/env Rscript
# Default: print the full job plan. No simulation is run without --mode=smoke/full.
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
script <- sub("^--file=", "", file_arg[1])
if (!file.exists(script)) script <- gsub("~+~", " ", script, fixed = TRUE)
source(file.path(dirname(normalizePath(script)), "common.R"))
root <- paper_root(script)
source(file.path(root, "analysis", "simulation_calibration.R"))

options_from <- function(args) {
  o <- list(mode = "plan", part = "all", repeats = NULL, output = NULL)
  for (a in args) {
    if (!grepl("^--[a-z]+=.\u002a$", a)) stop("Use --name=value arguments.")
    key <- sub("^--([^=]+)=.*$", "\\1", a)
    if (!key %in% names(o)) stop("Unknown option: ", key)
    o[[key]] <- sub("^--[^=]+=", "", a)
  }
  if (!o$mode %in% c("plan", "smoke", "full")) stop("mode: plan, smoke, or full")
  parts <- c("sim1", "sim2_linear", "sim2_nonlinear", "sim3")
  o$part <- if (o$part == "all") parts else strsplit(o$part, ",", fixed = TRUE)[[1]]
  if (any(!o$part %in% parts)) stop("Unknown simulation part.")
  if (is.null(o$repeats)) o$repeats <- if (o$mode == "smoke") "1" else "1:100"
  if (!grepl("^[0-9]+(:[0-9]+)?(,[0-9]+(:[0-9]+)?)*$", o$repeats)) stop("Invalid repeats.")
  o$repeats <- sort(unique(unlist(lapply(strsplit(o$repeats, ",")[[1]], function(x) {
    v <- as.integer(strsplit(x, ":", fixed = TRUE)[[1]])
    if (length(v) == 1) v else seq.int(v[1], v[2])
  }))))
  if (any(!o$repeats %in% 1:100)) stop("repeats must be in 1:100")
  if (o$mode == "smoke" && length(o$repeats) > 3) stop("Smoke mode is limited to 3 repetitions.")
  o
}
opt <- options_from(commandArgs(trailingOnly = TRUE))
smoke <- opt$mode == "smoke"
settings <- list(mode = opt$mode, train_n = if (smoke) 240L else 10000L,
                 test_n = if (smoke) 120L else 3000L,
                 sample_n = if (smoke) c(80L, 120L) else c(300L, 3000L),
                 ntree = if (smoke) 30L else 1000L,
                 rf_seed = 1234L, legacy_integration = TRUE,
                 metrics = paper_simulation_metrics,
                 calibration_args = list(nk = 5L),
                 calibration_grid = "last at or before t_star",
                 calibration_precision = "unrounded")
scenarios <- expand.grid(v = 1/c(0.5, 1, 1.5, 2), beta = c(0.5, 0.7, 1, 1.5))
scenarios$scenario <- seq_len(nrow(scenarios))
sim1 <- expand.grid(scenario = seq_len(nrow(scenarios)), rep = 1:100)
sim1$part <- "sim1"; sim1$seed <- seq_len(nrow(sim1))
sim1$v <- scenarios$v[sim1$scenario]; sim1$beta <- scenarios$beta[sim1$scenario]
sim1$n <- settings$test_n; sim1$censor <- 0; sim1$quantile <- 0
make_sim2_jobs <- function(name, v) {
  data.frame(scenario = 1L, rep = 1:100, part = name, seed = 1:100,
             v = v, beta = 1, n = settings$test_n, censor = 0, quantile = 0)
}
sample <- expand.grid(n = settings$sample_n, censor = c(0, 0.3, 0.7), rep = 1:100)
sample$seed <- seq_len(nrow(sample))
sim3 <- do.call(rbind, lapply(c(50L, 90L), function(q) {
  data.frame(scenario = 1L, rep = sample$rep, part = "sim3", seed = sample$seed,
             v = 1, beta = 1, n = sample$n, censor = sample$censor, quantile = q)
}))
jobs <- rbind(sim1, make_sim2_jobs("sim2_linear", 1),
               make_sim2_jobs("sim2_nonlinear", 3), sim3)
jobs$id <- with(jobs, sprintf("%s_s%02d_r%03d_n%d_c%02d_q%02d", part, scenario, rep, n,
                              round(100*censor), quantile))
stopifnot(!anyDuplicated(jobs$id))
chosen <- jobs[jobs$part %in% opt$part & jobs$rep %in% opt$repeats, ]
print(aggregate(id ~ part, chosen, length), row.names = FALSE)
if (opt$mode == "plan") {
  cat("PLAN ONLY. Full settings: train=10000; test=3000; repetitions=100.\n",
       "The complete design contains", nrow(jobs), "prediction/evaluation jobs.\n",
       "Use --mode=smoke for a small execution; --mode=full explicitly starts the full workload.\n")
  quit(save = "no", status = 0)
}
paper_load(root)
out <- if (is.null(opt$output)) file.path(root, "results", "simulations", opt$mode) else opt$output
out <- normalizePath(out, mustWork = FALSE)
paper_manifest(root, out, settings)
for (d in c("checkpoints", "models", "tables", "figures")) dir.create(file.path(out, d), showWarnings = FALSE)
utils::write.csv(jobs, file.path(out, "full_job_plan.csv"), row.names = FALSE)

# The old design fits once and regenerates test data. Keep exactly that structure.
clean_data <- function(d) d[, c("time", "status", "x1", "x2")]
simulate <- function(n, censor, v, beta, seed, nonlinear = FALSE) {
  clean_data(tm_sim_cox_weibull(n = n, pi_c = censor, v = v,
                                beta = beta, seed = seed, nonlinear = nonlinear))
}
compact_model <- function(fit) {
  attr(fit$terms, ".Environment") <- baseenv()
  if (!is.null(fit$formula)) environment(fit$formula) <- baseenv()
  fit
}
model_cache <- new.env(parent = emptyenv())
get_models <- function(job) {
  # Simulation 3 uses Simulation 1's beta=(1,-1), v=1 training fit (scenario 10).
  key <- if (job$part == "sim3") "sim1_10" else
    if (job$part == "sim1") paste0("sim1_", job$scenario) else job$part
  if (exists(key, model_cache, inherits = FALSE)) return(model_cache[[key]])
  path <- file.path(out, "models", paste0(key, ".rds"))
  if (file.exists(path)) {
    models <- readRDS(path)
  } else {
    if (startsWith(key, "sim1_")) {
      k <- as.integer(sub("sim1_", "", key))
      train <- simulate(settings$train_n, 0, scenarios$v[k],
                         scenarios$beta[k] * c(1, -1), k)
      models <- list(Cox = compact_model(survival::coxph(
        survival::Surv(time, status) ~ x1 + x2, train, x = TRUE, y = TRUE)))
    } else {
      nonlinear <- key == "sim2_nonlinear"
      train <- simulate(settings$train_n, 0, if (nonlinear) 3 else 1,
                         if (nonlinear) c(1,-1,-0.5,0.5) else c(1,-1),
                         1234, nonlinear)
      models <- list(
        Cox = compact_model(survival::coxph(survival::Surv(time, status) ~ x1 + x2,
                                            train, x = TRUE, y = TRUE)),
        Cox.reduce = compact_model(survival::coxph(survival::Surv(time, status) ~ x2,
                                                   train, x = TRUE, y = TRUE)),
        Ln = compact_model(survival::survreg(survival::Surv(time, status) ~ x1 + x2,
                                             train, dist = "lognormal", x = TRUE, y = TRUE)),
        Ln.fix = compact_model(survival::survreg(survival::Surv(time, status) ~ x1 + x2,
                                                 train, dist = "lognormal", scale = if (nonlinear) 5 else 2,
                                                 x = TRUE, y = TRUE)))
      if (nonlinear) {
        if (!requireNamespace("randomForestSRC", quietly = TRUE)) stop("Install randomForestSRC.")
        set.seed(settings$rf_seed)
        Surv <- survival::Surv
        # RF has only x1,x2 inputs; no interaction or quadratic input columns.
        models$RF <- randomForestSRC::rfsrc(Surv(time, status) ~ x1 + x2,
                                             data = train, ntree = settings$ntree,
                                             seed = -settings$rf_seed)
      }
    }
    saveRDS(models, path)
  }
  model_cache[[key]] <- models
  models
}
ref <- simulate(settings$train_n, 0, 1, c(1,-1), 1234)
horizons <- stats::quantile(ref$time, c(0.5, 0.9), names = FALSE)
utils::write.csv(data.frame(quantile = c(50, 90), tau = horizons),
                 file.path(out, "horizons.csv"), row.names = FALSE)

evaluate_job <- function(job) {
  warnings <- character()
  result <- tryCatch(withCallingHandlers({
    models <- get_models(job)
    nonlinear <- job$part == "sim2_nonlinear"
    test <- simulate(job$n, job$censor, job$v,
                      if (nonlinear) c(1,-1,-0.5,0.5) else job$beta * c(1,-1),
                      job$seed, nonlinear)
    tau <- if (job$part == "sim3") horizons[match(job$quantile, c(50,90))] else 10e10
    # The linear/nonlinear comparison uses the training median, as in the Rmd.
    t_star <- if (job$part == "sim3") tau else if (startsWith(job$part, "sim2")) {
      train_ref <- simulate(settings$train_n, 0, job$v,
                             if (nonlinear) c(1,-1,-0.5,0.5) else c(1,-1), 1234, nonlinear)
      median(train_ref$time)
    } else median(test$time)
    rows <- list(); support <- list(); calibration <- list()
    for (m in names(models)) {
      fit <- models[[m]]
      covs <- if (m == "Cox.reduce") "x2" else c("x1", "x2")
      if (inherits(fit, "coxph")) {
        pred <- tm_predict_coxph(fit, covs, tau = tau, new_data = test)
      } else if (inherits(fit, "survreg")) {
        pred <- tm_predict_survreg(fit, covs, tau = tau, new_data = test)
      } else {
        grid <- sort(test$time)
        sp <- pec::predictSurvProb(fit, newdata = test, times = grid)
        pred <- list(model = NULL, times = test$time, status = test$status,
                     surv_prob = sp, pred = NULL, new_data = test, covs = covs)
      }
      # Preserve the existing metrics and their historical rounding.
      table <- tm_summarize(setNames(list(pred), m), metrics = setdiff(settings$metrics, "ici"),
                            tau = tau, t_star = t_star, digits = 4)
      cal <- paper_simulation_calibration(pred, t_star, settings$calibration_args)
      rows[[m]] <- rbind(
        data.frame(Model = m, Metric = table$Metric, Value = table[[m]]),
        data.frame(Model = m, Metric = "ici", Value = cal$ici))
      calibration[[m]] <- cbind(Model = m, cal$diagnostics)
      if (inherits(fit, "coxph")) support[[m]] <- cbind(Model = m, paper_support(fit, test, tau))
    }
    list(ok = TRUE, job = job, metrics = do.call(rbind, rows), tau = tau,
         t_star = t_star, support = do.call(rbind, support),
         calibration = do.call(rbind, calibration),
         observed_censoring = mean(test$status == 0))
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) list(ok = FALSE, job = job, error = conditionMessage(e)))
  result$warnings <- unique(warnings)
  result
}

for (i in seq_len(nrow(chosen))) {
  j <- chosen[i, ]
  path <- file.path(out, "checkpoints", paste0(j$id, ".rds"))
  if (file.exists(path) &&
      paper_simulation_record_complete(readRDS(path), settings$metrics)) {
    message("Reuse ", j$id); next
  }
  message("Run ", j$id)
  record <- evaluate_job(j)
  temporary <- paste0(path, ".tmp")
  saveRDS(record, temporary)
  if (!file.rename(temporary, path)) stop("Could not commit checkpoint: ", path)
  if (!record$ok) message("FAILED: ", record$error)
  invisible(gc())
}

records <- lapply(list.files(file.path(out, "checkpoints"), pattern = "\\.rds$", full.names = TRUE), readRDS)
status <- do.call(rbind, lapply(records, function(r) data.frame(
  id = r$job$id, part = r$job$part, rep = r$job$rep, ok = r$ok,
  metrics_complete = paper_simulation_record_complete(r, settings$metrics, finite = FALSE),
  finite = paper_simulation_record_complete(r, settings$metrics),
  error = if (r$ok) "" else r$error, warnings = paste(r$warnings, collapse = " | "))))
utils::write.csv(status, file.path(out, "run_status.csv"), row.names = FALSE)
successful <- Filter(function(r)
  paper_simulation_record_complete(r, settings$metrics, finite = FALSE), records)
if (!length(successful)) stop("No job succeeded. See run_status.csv.")
raw <- do.call(rbind, lapply(successful, function(r) {
  meta <- r$job[rep(1, nrow(r$metrics)), ]
  cbind(meta, tau = r$tau, t_star = r$t_star,
         observed_censoring = r$observed_censoring, r$metrics)
}))
utils::write.csv(raw, file.path(out, "tables", "metrics_raw.csv"), row.names = FALSE)
support <- do.call(rbind, lapply(successful, function(r) {
  if (!is.null(r$support)) cbind(id = r$job$id, r$support)
}))
utils::write.csv(support, file.path(out, "tables", "support_diagnostics.csv"), row.names = FALSE)
calibration <- do.call(rbind, lapply(successful, function(r) {
  if (!is.null(r$calibration)) cbind(id = r$job$id, part = r$job$part,
    rep = r$job$rep, r$calibration)
}))
utils::write.csv(calibration, file.path(out, "tables", "calibration_diagnostics.csv"), row.names = FALSE)
group <- c("part", "scenario", "v", "beta", "n", "censor", "quantile", "Model", "Metric")
keys <- interaction(raw[group], drop = TRUE, lex.order = TRUE)
summary <- do.call(rbind, lapply(split(seq_len(nrow(raw)), keys), function(ix) {
  d <- raw[ix, ]; value <- d$Value[is.finite(d$Value)]
  B <- length(value)
  cbind(d[1, group], completed_jobs = nrow(d), finite_B = B,
         mean = if (B) mean(value) else NA_real_,
         empirical_SD = if (B > 1) stats::sd(value) else NA_real_,
         MCSE = if (B > 1) stats::sd(value)/sqrt(B) else NA_real_)
}))
utils::write.csv(summary, file.path(out, "tables", "metrics_summary.csv"), row.names = FALSE)

# Recompute the same finite-sample uncensored references used by the manuscript.
# When only some repetitions exist, report the reference as partial explicitly.
s3 <- raw[raw$part == "sim3", ]
errors <- data.frame()
if (nrow(s3)) {
  ref_summary <- summary[summary$part == "sim3" & summary$n == max(settings$sample_n) & summary$censor == 0, ]
  refs <- ref_summary[, c("quantile", "Model", "Metric", "mean", "MCSE", "finite_B")]
  names(refs)[4:6] <- c("reference", "reference_MCSE", "reference_B")
  utils::write.csv(refs, file.path(out, "tables", "sim3_references.csv"), row.names = FALSE)
  errors <- merge(s3, refs, by = c("quantile", "Model", "Metric"), all.x = TRUE, sort = FALSE)
  errors$error <- errors$Value - errors$reference
  utils::write.csv(errors, file.path(out, "tables", "sim3_errors.csv"), row.names = FALSE)
  error_group <- interaction(errors[c("n", "censor", "quantile", "Metric")], drop = TRUE)
  accuracy <- do.call(rbind, lapply(split(errors, error_group), function(d) {
    z <- d$error[is.finite(d$error)]; B <- length(z)
    cbind(d[1, c("n", "censor", "quantile", "Metric", "reference_B", "reference_MCSE")],
           finite_B = B, bias = if (B) mean(z) else NA_real_,
           RMSE = if (B) sqrt(mean(z^2)) else NA_real_,
           bias_MCSE = if (B > 1) stats::sd(z)/sqrt(B) else NA_real_)
  }))
  utils::write.csv(accuracy, file.path(out, "tables", "sim3_accuracy.csv"), row.names = FALSE)
}

plotted <- c("pseudo_r2", "brier_score", "harrell_c", "td_auc", "ici")
save_plot <- function(p, name, width, height) {
  p <- p + ggplot2::labs(caption = if (smoke) "SMOKE TEST: reduced sample sizes and repetitions; not manuscript results." else
                          "See configuration.rds and run_status.csv for versions and completion.")
  ggplot2::ggsave(file.path(out, "figures", paste0(name, ".png")), p, width = width, height = height, dpi = 120)
  ggplot2::ggsave(file.path(out, "figures", paste0(name, ".pdf")), p, width = width, height = height)
}
s1 <- summary[summary$part == "sim1" & summary$Metric %in% plotted, ]
if (nrow(s1)) {
  p <- ggplot2::ggplot(s1, ggplot2::aes(1/v, mean, group = factor(beta), colour = factor(beta))) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - 1.96*MCSE, ymax = mean + 1.96*MCSE), width = 0.03, na.rm = TRUE) +
    ggplot2::facet_wrap(~Metric, scales = "free_y") + ggplot2::theme_bw() +
    ggplot2::labs(x = "1 / Weibull shape", y = "Replication mean", colour = "Beta magnitude")
  save_plot(p, "paper.fig4.updated", 10, 7)
}
for (part in c("sim2_linear", "sim2_nonlinear")) {
  d <- summary[summary$part == part, ]
  if (!nrow(d)) next
  bars <- d[d$Metric %in% plotted, ]
  if (part == "sim2_nonlinear") bars <- bars[bars$Model %in% c("RF", "Cox", "Ln"), ]
  components <- d[d$Metric %in% c("r_square", "l_square") & d$Model %in% bars$Model, ]
  components$Component <- components$Metric; components$Metric <- "pseudo_r2"
  p <- ggplot2::ggplot(bars, ggplot2::aes(Model, mean)) + ggplot2::geom_col(fill = "grey65") +
    ggplot2::geom_point(data = components, ggplot2::aes(shape = Component), size = 2.5) +
    ggplot2::scale_shape_manual(values = c(r_square = 16, l_square = 1)) +
    ggplot2::facet_wrap(~Metric, scales = "free_y") + ggplot2::theme_bw() +
    ggplot2::labs(y = "Replication mean", x = NULL)
  save_plot(p, if (part == "sim2_linear") "paper.fig5.updated" else "paper.fig5.2.updated", 10, 7)
}
if (nrow(errors)) {
  d <- errors[errors$Metric %in% c(plotted, "uno_c"), ]
  p <- ggplot2::ggplot(d, ggplot2::aes(factor(n), error, fill = factor(quantile))) +
    ggplot2::geom_boxplot(outlier.size = 0.4) +
    ggplot2::facet_grid(censor ~ Metric, scales = "free_y") + ggplot2::theme_bw() +
    ggplot2::labs(x = "Test sample size", y = "Estimate minus uncensored reference", fill = "Horizon percentile")
  save_plot(p, "paper.fig6.updated", 15.5, 7)
}
completion <- aggregate(id ~ part, jobs, length); names(completion)[2] <- "planned_jobs"
done <- aggregate(id ~ part, status, length); names(done)[2] <- "attempted_jobs"
completion <- merge(completion, done, all.x = TRUE)
completion$successful_finite_jobs <- vapply(completion$part, function(p) {
  sum(status$part == p & status$ok & status$finite)
}, integer(1))
completion$failed_jobs <- vapply(completion$part, function(p) sum(status$part == p & !status$ok), integer(1))
completion$nonfinite_jobs <- vapply(completion$part, function(p) {
  sum(status$part == p & status$ok & !status$finite)
}, integer(1))
completion$failure_rate <- (completion$failed_jobs + completion$nonfinite_jobs)/completion$attempted_jobs
completion$complete <- !is.na(completion$attempted_jobs) &
  completion$successful_finite_jobs == completion$planned_jobs
completion$mode <- opt$mode
utils::write.csv(completion, file.path(out, "completion.csv"), row.names = FALSE)
cat("Outputs:", out, "\n")
cat("Mode:", opt$mode, "; attempted:", nrow(status), "; failed:", sum(!status$ok),
     "; non-finite:", sum(status$ok & !status$finite), "\n")
if (any(!status$ok | !status$finite)) quit(save = "no", status = 1)
