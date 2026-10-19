# Run in a fresh R session after installing the source package and hard deps.
# Restrict library discovery to the declared Depends/Imports closure,
# so an incidental suggested package in a runner cache cannot hide a defect.
# LinkingTo supplies build-time headers, which binary installations need not
# retain; those packages are not required to run an already installed package.
db <- installed.packages()
db <- db[!duplicated(db[, "Package"]), , drop = FALSE]
required <- unique(c("TimeMetric", unlist(tools::package_dependencies(
  "TimeMetric", db = db, which = c("Depends", "Imports"),
  recursive = TRUE))))
required <- setdiff(required, "R")
missing <- setdiff(required, db[, "Package"])
if (length(missing)) stop("Missing runtime dependencies: ", paste(missing, collapse = ", "))
stopifnot("tidyr" %in% required)
lib <- tempfile("timemetric-hard-deps-")
dir.create(lib)
check <- function() {
  on.exit(unlink(lib, recursive = TRUE), add = TRUE)
  # Symlinks avoid copying compiled libraries; this job runs on Linux/macOS.
  paths <- file.path(db[required, "LibPath"], required)
  stopifnot(all(file.symlink(paths, file.path(lib, required))))
  .libPaths(c(lib, .Library), include.site = FALSE)
  library(TimeMetric)
  d <- tm_sim_cox_weibull(120, pi_c = .3, v = 2,
                          beta = c(.5, -.5), seed = 1051)
  fit <- survival::coxph(survival::Surv(time, status) ~ x1 + x2,
                         d, x = TRUE, y = TRUE)
  pred <- tm_predict_coxph(fit, c("x1", "x2"), new_data = d)
  km <- survival::survfit(survival::Surv(time, 1 - status) ~ 1, d)
  result <- tm_evaluate_two_phase(pred, km_cens_fit = km,
                                  case_weights = rep(1, nrow(d)))
  stopifnot(nrow(result) == 5L, all(is.finite(result$Value)),
            "td_auc" %in% result$Metric)
  print(result)
  cat("Default two-phase evaluation passed with required dependencies only.\n")
}
check()
