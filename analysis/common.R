# Shared helpers. Always load this checkout, never the globally installed 0.1.0.
paper_root <- function(script) normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)

paper_load <- function(root) {
  if (!requireNamespace("pkgload", quietly = TRUE)) stop("Install pkgload first.")
  suppressPackageStartupMessages(pkgload::load_all(root, quiet = TRUE))
  actual <- normalizePath(getNamespaceInfo(asNamespace("TimeMetric"), "path"))
  stopifnot(identical(actual, normalizePath(root)))
  invisible(actual)
}

paper_metrics <- c("r_square", "l_square", "pseudo_r2", "harrell_c",
                   "uno_c", "brier_score", "td_auc")

paper_git <- function(root, ...) {
  ans <- system2("git", c("-C", shQuote(root), ...), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(ans, "status"))) stop(paste(ans, collapse = "\n"))
  ans
}

paper_signature <- function(root, settings) {
  paths <- c(list.files(file.path(root, "R"), full.names = TRUE),
             list.files(file.path(root, "analysis"), pattern = "\\.R$", full.names = TRUE),
             file.path(root, "DESCRIPTION"))
  packages <- c("TimeMetric", "survival", "tdROC", "yardstick", "randomForestSRC",
                "cmprsk", "riskCommunicator", "expint", "pec", "dplyr", "pkgload")
  versions <- setNames(vapply(packages, function(p) {
    if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_
  }, character(1)), packages)
  list(upstream = "426272c2bdcbf94079c035166fc7e3804ad2b632",
       source_md5 = tools::md5sum(paths), settings = settings,
       R = R.version.string, RNG = RNGkind(), packages = versions)
}

paper_manifest <- function(root, out, settings) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  file <- file.path(out, "configuration.rds")
  sig <- paper_signature(root, settings)
  if (file.exists(file) && !identical(readRDS(file), sig)) {
    stop("Source, settings or versions differ from this output directory. Use a new --output directory.")
  }
  saveRDS(sig, file)
  capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
  invisible(sig)
}

paper_legacy_cox <- function(root) {
  code <- paper_git(root, "show", "426272c2bdcbf94079c035166fc7e3804ad2b632:R/tm_predict_coxph.R")
  e <- new.env(parent = asNamespace("TimeMetric"))
  eval(parse(text = code), envir = e)
  e$tm_predict_coxph
}

paper_write_md <- function(d, file) {
  rows <- apply(d, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
  writeLines(c(paste0("| ", paste(names(d), collapse = " | "), " |"),
               paste0("| ", paste(rep("---", ncol(d)), collapse = " | "), " |"), rows), file)
}

paper_support <- function(model, test, horizon) {
  km <- survival::survfit(survival::Surv(time, 1 - status) ~ 1, data = test)
  data.frame(train_n = nrow(model$y), test_n = nrow(test),
             train_last_time = max(model$y[, 1]), test_last_time = max(test$time),
             requested_horizon = horizon, integration_end = min(horizon, max(test$time)),
             at_risk = sum(test$time >= horizon),
             G_horizon = if (horizon <= max(test$time))
               summary(km, times = horizon, extend = TRUE)$surv else NA_real_,
             evaluation_beyond_training = any(test$time > max(model$y[, 1])))
}
