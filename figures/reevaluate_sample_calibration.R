#!/usr/bin/env Rscript
# Validate and add weighted ICI to saved GBSG predictions, without model refits.
script <- sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE)[1])
if (!file.exists(script)) script <- gsub('~+~', ' ', script, fixed = TRUE)
root <- normalizePath(file.path(dirname(script), '..'))
source(file.path(root, 'analysis', 'common.R'))
source(file.path(root, 'figures', 'original_format.R'))
source(file.path(root, 'figures', 'paper_outputs.R'))
input <- NULL
output <- NULL
for (arg in commandArgs(trailingOnly = TRUE)) {
  if (startsWith(arg, '--input=')) input <- sub('^--input=', '', arg) else
  if (startsWith(arg, '--output=')) output <- sub('^--output=', '', arg) else
    stop('Use --input=/path/to/source_run and --output=/path/to/new_run.')
}
if (is.null(input)) input <- paper_current_real_run(root)
input <- normalizePath(input, mustWork = TRUE)
if (is.null(output)) output <- file.path(root, 'results', 'real_data', 'runs',
  paste0(format(Sys.time(), '%Y-%m-%d_%H%M%S', tz = 'America/Los_Angeles'), '_sample_ici'))
if (file.exists(output)) stop('Choose a new output directory; existing runs are preserved.')
paper_load(root)
original <- readRDS(file.path(input, 'gbsg', 'gbsg2.rds'))
result <- original

# Recover the actual inclusion weights independently from phase-one membership.
full <- result$test
selected <- full$subco | full$status == 1
sample <- full[selected, ]
index <- match(result$cc_test$id, sample$id)
stopifnot(!anyNA(index), !anyDuplicated(result$cc_test$id),
          nrow(sample) == nrow(result$cc_test),
          identical(sample$time[index], result$cc_test$time),
          identical(sample$status[index], result$cc_test$status))
expected_weights <- ifelse(result$cc_test$status == 1, 1,
  sum(full$status == 0)/sum(sample$status == 0))
stopifnot(isTRUE(all.equal(result$weights, expected_weights)),
          isTRUE(all.equal(result$cc_test$weight, expected_weights)),
          isTRUE(all.equal(sum(result$weights), nrow(full))))
weight_file <- utils::read.csv(file.path(input, 'gbsg', 'evaluation_subjects_and_weights.csv'))
stopifnot(identical(weight_file$id, result$cc_test$id),
          isTRUE(all.equal(weight_file$weight, expected_weights)))

calibration <- list()
diagnostics <- list()
subjects <- list()
sensitivity <- list()
t <- unname(stats::median(result$cc_test$time))
for (name in names(result$predictions)) {
  p <- result$predictions[[name]]
  stopifnot(identical(p$times, result$cc_test$time),
            identical(as.numeric(p$status), as.numeric(result$cc_test$status)))
  k <- findInterval(t, sort(p$times))
  risk <- 1-p$surv_prob[, k]
  cal <- tm_calibration(risk, p$times, p$status, t, case_weights = result$weights)
  error <- abs(risk-cal$curve$calibrated)
  hand_ici <- sum(expected_weights*error)/sum(expected_weights)
  stopifnot(nrow(cal$curve) == nrow(result$cc_test),
            isTRUE(all.equal(unname(cal$metrics['ici']), hand_ici)))

  # Independent natural-cubic basis spanning the same spline space, with
  # survival curves predicted by survfit rather than predict(type='expected').
  stopifnot(cal$method == 'weighted_cox_rcs')
  z <- log(-log1p(-risk))
  knots <- cal$knots
  basis <- splines::ns(z, knots = knots[-c(1, length(knots))],
                       Boundary.knots = knots[c(1, length(knots))])
  colnames(basis) <- paste0('B', seq_len(ncol(basis)))
  frame <- data.frame(time = p$times, event = p$status, basis)
  weights <- expected_weights
  fit <- survival::coxph(survival::Surv(time, event) ~ ., data = frame,
                          weights = weights, x = TRUE, model = TRUE)
  curves <- survival::survfit(fit, newdata = as.data.frame(basis), se.fit = FALSE)
  reference <- 1-as.numeric(summary(curves, times = t, extend = TRUE)$surv)
  reference_error <- max(abs(reference-cal$curve$calibrated))
  reference_ici <- sum(weights*abs(risk-reference))/sum(weights)
  stopifnot(reference_error < 1e-7, abs(reference_ici-hand_ici) < 1e-8)

  # Check weight scaling, subject order and reduction to the ordinary method.
  scaled <- tm_calibration(risk, p$times, p$status, t,
                             case_weights = 10*weights)
  o <- rev(seq_along(risk))
  permuted <- tm_calibration(risk[o], p$times[o], p$status[o], t,
                               case_weights = weights[o])
  unit <- tm_calibration(risk, p$times, p$status, t,
                           case_weights = rep(1, length(weights)))
  ordinary <- tm_calibration(risk, p$times, p$status, t)
  scale_error <- max(abs(scaled$metrics-cal$metrics))
  order_error <- max(abs(permuted$curve$calibrated-rev(cal$curve$calibrated)))
  unit_error <- max(abs(unit$curve$calibrated-ordinary$curve$calibrated))
  stopifnot(scale_error < 1e-9, order_error < 1e-7, unit_error < 1e-7)
  calibration[[name]] <- cal
  diagnostics[[name]] <- data.frame(Model = name, Evaluation_time_days = t,
    Prediction_grid_time_days = sort(p$times)[k], N = length(risk),
    Phase_one_test_n = nrow(full), Weight_sum = sum(weights),
    Effective_weighted_n = sum(weights)^2/sum(weights^2),
    Method = cal$method, ICI = hand_ici, E50 = unname(cal$metrics['e50']),
    E90 = unname(cal$metrics['e90']), Independent_basis_ICI = reference_ici,
    Independent_curve_max_error = reference_error,
    Scaling_metric_max_error = scale_error, Permutation_curve_max_error = order_error,
    Unit_weight_curve_max_error = unit_error)
  subjects[[name]] <- data.frame(Model = name, id = result$cc_test$id,
    Time = p$times, Status = p$status, Sampling_weight = weights,
    Predicted_risk = risk, Calibrated_risk = cal$curve$calibrated,
    Absolute_error = error, ICI_contribution = weights*error/sum(weights))
  for (nk in 3:5) {
    alternative <- tm_calibration(risk, p$times, p$status, t,
                                   nk = nk, case_weights = weights)
    sensitivity[[length(sensitivity)+1L]] <- data.frame(Model = name,
      Requested_knots = nk, Actual_knots = length(alternative$knots),
      ICI = unname(alternative$metrics['ici']))
  }
}
ici <- tm_sample_design(result$predictions, case_weights = result$weights,
                         metrics = 'ici', t_star = t, digits = 4)
previous <- original$metrics[original$metrics$Metric != 'ici', , drop = FALSE]
table <- rbind(previous, ici[, names(previous), drop = FALSE])
for (name in names(result$predictions))
  stopifnot(table[[name]][table$Metric == 'ici'] == round(calibration[[name]]$metrics['ici'], 4))
retained <- table[match(previous$Metric, table$Metric), names(previous), drop = FALSE]
rownames(retained) <- rownames(previous) <- NULL
stopifnot(identical(retained, previous), sum(table$Metric == 'ici') == 1L)

# Write a fresh complete run so the paper manifest has a single source root.
dir.create(output, recursive = TRUE)
output <- normalizePath(output)
stopifnot(all(file.copy(list.files(input, full.names = TRUE, all.files = TRUE,
  no.. = TRUE), output, recursive = TRUE)))
provenance <- list(source_run = input, output_run = output, models_refitted = FALSE,
  evaluated_at = format(Sys.time(), tz = 'America/Los_Angeles', usetz = TRUE),
  source_signature = paper_signature(root, list(calibration = 'sampling-weighted ICI')),
  script_md5 = unname(tools::md5sum(script)))
result$metrics <- table
result$calibration <- calibration
result$calibration_reevaluation <- provenance
for (field in setdiff(names(original), c('metrics', 'calibration', 'calibration_reevaluation')))
  stopifnot(identical(result[[field]], original[[field]]))
saveRDS(result, file.path(output, 'gbsg', 'gbsg2.rds'))
write_metrics <- function(table, stem) {
  utils::write.csv(table, paste0(stem, '.csv'), row.names = FALSE)
  paper_write_md(table, paste0(stem, '.md'))
  tex <- c(paste0('\\begin{tabular}{l', strrep('r', ncol(table)-1), '}'),
    '\\toprule', paste0(paste(gsub('_', '\\_', names(table), fixed = TRUE),
      collapse = ' & '), ' \\\\'), '\\midrule')
  for (i in seq_len(nrow(table))) tex <- c(tex, paste0(paste(c(
    gsub('_', '\\_', table$Metric[i], fixed = TRUE),
    formatC(as.numeric(table[i, -1]), format = 'f', digits = 4)), collapse = ' & '), ' \\\\'))
  writeLines(c(tex, '\\bottomrule', '\\end{tabular}'), paste0(stem, '.tex'))
}
write_metrics(table, file.path(output, 'gbsg', 'metrics'))
write_metrics(table, file.path(output, 'tables', 'gbsg_metrics'))
settings_file <- file.path(output, 'gbsg', 'analysis_settings.csv')
settings <- utils::read.csv(settings_file)
settings$Calibration_time <- t
settings$Calibration_method <- 'sampling-weighted Cox restricted cubic spline; five knots'
utils::write.csv(settings, settings_file, row.names = FALSE)
utils::write.csv(do.call(rbind, diagnostics),
  file.path(output, 'gbsg', 'calibration_diagnostics.csv'), row.names = FALSE)
utils::write.csv(do.call(rbind, subjects),
  file.path(output, 'gbsg', 'calibration_subjects.csv'), row.names = FALSE)
utils::write.csv(do.call(rbind, sensitivity),
  file.path(output, 'gbsg', 'calibration_knot_sensitivity.csv'), row.names = FALSE)
saveRDS(provenance, file.path(output, 'gbsg', 'calibration_reevaluation.rds'))
capture.output(sessionInfo(), file = file.path(output, 'gbsg', 'calibration_sessionInfo.txt'))

# Relocate the FHS active-horizon pointer to its byte-identical copied result.
active_file <- file.path(output, 'framingham', 'active_horizon.csv')
active <- utils::read.csv(active_file, stringsAsFactors = FALSE)
stopifnot(startsWith(active$Source_predictions, paste0(input, '/')))
active$Source_predictions <- file.path(output,
  substring(active$Source_predictions, nchar(input)+2L))
active$Fit_source_md5 <- unname(tools::md5sum(file.path(output, 'framingham',
  'framingham_visit3_10y.rds')))
utils::write.csv(active, active_file, row.names = FALSE)
unchanged <- c(list.files(file.path(input, 'pbc'), recursive = TRUE, full.names = TRUE),
  list.files(file.path(input, 'framingham'), recursive = TRUE, full.names = TRUE),
  list.files(input, pattern = 'paper\\.fig[123].*', recursive = TRUE, full.names = TRUE),
  file.path(input, 'tables', paste0(rep(c('pbc', 'framingham'), each = 3),
    '_metrics.', rep(c('csv', 'md', 'tex'), 2))))
unchanged <- unchanged[unchanged != file.path(input, 'framingham', 'active_horizon.csv')]
copied <- file.path(output, substring(unchanged, nchar(input)+2L))
stopifnot(identical(unname(tools::md5sum(unchanged)), unname(tools::md5sum(copied))))
historical_note <- file.path(output, 'CALIBRATION_UPDATE.md')
if (file.exists(historical_note)) writeLines(c(
  'Historical cohort-update note copied from the source run. See SAMPLE_CALIBRATION_UPDATE.md for the subsequent GBSG update.',
  '', readLines(historical_note)), historical_note)
writeLines(c('# GBSG weighted calibration update', '', paste('Source run:', input),
  'Saved models, splits, sampled IDs, sampling weights, predictions and figures are unchanged.',
  'Only weighted ICI was added to the existing GBSG metric table; old metric values are retained exactly.',
  'The evaluation horizon and prediction grid time are both 637 days, on 115 case-cohort test subjects.',
  'Weights reconstructed from the 229-person phase-one test cohort match saved weights and sum to 229.',
  'Reference calculations use a separate natural-cubic basis and survfit predictions at the same knots.',
  'Diagnostics also check unit weights, weight scaling, subject order and 3/4/5-knot sensitivity.',
  'The five-knot estimator was prespecified; sensitivity results do not select the reported estimator.',
  'All PBC/FHS results and all figure data remain byte-identical; only the FHS active pointer is relocated.',
  'Historical configuration files describe original model fits; gbsg/calibration_reevaluation.rds records this update.'),
  file.path(output, 'SAMPLE_CALIBRATION_UPDATE.md'))
paper_publish_real_outputs(root, output)
print(do.call(rbind, diagnostics), row.names = FALSE)
print(do.call(rbind, sensitivity), row.names = FALSE)
cat('GBSG weighted ICI validation PASS. Output:', output, '\n')
