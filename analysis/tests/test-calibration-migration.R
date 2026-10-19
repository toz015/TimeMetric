# Run from the checkout: Rscript -e 'testthat::test_file("analysis/tests/test-calibration-migration.R")'
# This integration test is outside the CRAN tarball because figures/ is excluded.
testthat::test_that("pre-ICI saved runs migrate completely before publication", {
  root <- normalizePath(file.path("..", ".."))
  scratch <- tempfile("calibration-migration-")
  dir.create(scratch)
  withr::defer(unlink(scratch, recursive = TRUE))
  # Run the actual entry points in an isolated checkout, never publish fixtures
  # into the user's manuscript folder. No clinical data or network required.
  files <- c("R", "man", "data", "analysis", "figures", "DESCRIPTION", "NAMESPACE", "LICENSE")
  stopifnot(all(file.copy(file.path(root, files), scratch, recursive = TRUE)))
  old <- file.path(scratch, "old")
  for (d in c("pbc", "framingham", "gbsg", "tables"))
    dir.create(file.path(old, d), recursive = TRUE)
  old <- normalizePath(old)
  withr::local_seed(994)
  n <- 180L
  x <- rnorm(n)
  rate <- exp(.4*x)/5000
  event <- rexp(n, rate)
  censor <- rexp(n, 1/10000)
  time <- pmin(event, censor)
  status <- as.integer(event <= censor)
  survival <- list(times = time, status = status,
    surv_prob = outer(rate, sort(time), function(r,t) exp(-r*t)))
  make_models <- function(pred, n) setNames(rep(list(pred), n), paste0("Model", seq_len(n)))
  make_table <- function(n) {
    tab <- data.frame(Metric = c("r_square", "l_square", "pseudo_r2", "brier_score", "td_auc"))
    for (i in seq_len(n)) tab[[paste0("Model", i)]] <- seq(.1, .5, length.out = 5)
    tab
  }
  write_table <- function(tab, stem) {
    write.csv(tab, paste0(stem, ".csv"), row.names = FALSE)
    for (ext in c("md", "tex")) writeLines("synthetic old table", paste0(stem, ".", ext))
  }
  pbc <- make_models(survival, 3)
  saveRDS(list(before = pbc, after = pbc, models = "saved synthetic fits"),
    file.path(old, "pbc", "pbc_comparison.rds"))
  for (version in c("before", "after")) write_table(make_table(3),
    file.path(old, "pbc", paste0("metrics_", version)))
  write_table(make_table(3), file.path(old, "tables", "pbc_metrics"))
  write.csv(data.frame(setting = "synthetic"), file.path(old, "pbc", "analysis_settings.csv"), row.names = FALSE)

  competing <- rexp(n, 1/8000)
  cr_time <- pmin(event, competing, censor)
  cr_status <- ifelse(censor < pmin(event, competing), 0L, ifelse(event <= competing, 1L, 2L))
  grid <- c(0, 3650, 4380)
  cif <- outer(grid, rate, function(t,r) r/(r+1/8000)*(1-exp(-(r+1/8000)*t)))
  pred <- list(times = cr_time, status = cr_status, cif_pred = cbind(time = grid, cif))
  for (year in c(10, 12)) saveRDS(list(settings = list(tau = year*365),
    predictions = make_models(pred, 4), metrics = make_table(4), models = "saved fits", cohort = "synthetic"),
    file.path(old, "framingham", paste0("framingham_visit3_", year, "y.rds")))
  write.csv(data.frame(Horizon_years = 12, tau_days = 4380,
    Source_predictions = file.path(old, "framingham", "framingham_visit3_12y.rds")),
    file.path(old, "framingham", "active_horizon.csv"), row.names = FALSE)
  write_table(make_table(4), file.path(old, "tables", "framingham_metrics"))

  full <- data.frame(id = seq_len(n), time = time, status = status, subco = rep(c(TRUE, FALSE), n/2))
  sampled <- full[full$subco | full$status == 1, ]
  weights <- ifelse(sampled$status == 1, 1, sum(full$status == 0)/sum(sampled$status == 0))
  sampled$weight <- weights
  sample_pred <- list(times = sampled$time, status = sampled$status,
    surv_prob = outer(rate[sampled$id], sort(sampled$time), function(r,t) exp(-r*t)))
  saveRDS(list(test = full, cc_test = sampled, weights = weights,
    predictions = make_models(sample_pred, 2), metrics = make_table(2), models = "saved fits"),
    file.path(old, "gbsg", "gbsg2.rds"))
  write.csv(sampled[c("id", "weight")], file.path(old, "gbsg", "evaluation_subjects_and_weights.csv"), row.names = FALSE)
  write.csv(data.frame(setting = "synthetic"), file.path(old, "gbsg", "analysis_settings.csv"), row.names = FALSE)
  write_table(make_table(2), file.path(old, "tables", "gbsg_metrics"))
  for (i in 1:3) for (ext in c("png", "pdf")) writeLines("synthetic figure placeholder",
    file.path(old, c("pbc", "framingham", "gbsg")[i], paste0("paper.fig", i, ".", ext)))
  originals <- list.files(old, full.names = TRUE, recursive = TRUE)
  hashes <- tools::md5sum(originals)
  cohort <- file.path(scratch, "cohort")
  final <- file.path(scratch, "final")
  published <- file.path(scratch, "results", "real_data", "paper")
  # An existing manifest must not advance until all three datasets are ready.
  dir.create(published, recursive = TRUE)
  write.csv(data.frame(Dataset = c("pbc", "framingham", "gbsg"), Source_run = old),
    file.path(published, "source_manifest.csv"), row.names = FALSE)
  before_manifest <- tools::md5sum(file.path(published, "source_manifest.csv"))
  run <- function(script, input, output) {
    log <- system2(file.path(R.home("bin"), "Rscript"),
      c("--vanilla", shQuote(file.path(scratch, "figures", script)),
        shQuote(paste0("--input=", input)), shQuote(paste0("--output=", output))),
      stdout = TRUE, stderr = TRUE)
    testthat::expect_null(attr(log, "status"), info = paste(log, collapse = "\n"))
    if (!is.null(attr(log, "status"))) stop(paste(log, collapse = "\n"))
  }
  run("reevaluate_real_calibration.R", old, cohort)
  testthat::expect_equal(tools::md5sum(file.path(published, "source_manifest.csv")), before_manifest)
  testthat::expect_equal(read.csv(file.path(cohort, "tables", "gbsg_metrics.csv")), make_table(2))
  # Explicit --input also works if no published run exists (or its pointer is stale).
  unlink(file.path(published, "source_manifest.csv"))
  run("reevaluate_sample_calibration.R", cohort, final)
  for (d in c("pbc", "framingham", "gbsg")) {
    previous <- read.csv(file.path(old, "tables", paste0(d, "_metrics.csv")))
    updated <- read.csv(file.path(final, "tables", paste0(d, "_metrics.csv")))
    retained <- updated[updated$Metric != "ici", ]; rownames(retained) <- NULL
    testthat::expect_equal(retained, previous)
    testthat::expect_equal(sum(updated$Metric == "ici"), 1)
    testthat::expect_true(all(is.finite(as.numeric(updated[updated$Metric == "ici", -1]))))
  }
  manifest <- read.csv(file.path(published, "source_manifest.csv"))
  testthat::expect_setequal(manifest$Dataset, c("pbc", "framingham", "gbsg"))
  testthat::expect_identical(unique(manifest$Source_run), normalizePath(final))
  testthat::expect_true(all(file.exists(file.path(published, manifest$Table))))
  testthat::expect_equal(tools::md5sum(originals), hashes)
})
