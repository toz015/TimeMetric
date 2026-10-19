#!/usr/bin/env Rscript
# Run the three updated real-data illustrations in separate clean R sessions.
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
script <- sub("^--file=", "", file_arg[1])
if (!file.exists(script)) script <- gsub("~+~", " ", script, fixed = TRUE)
source(file.path(dirname(normalizePath(script)), "common.R"))
root <- paper_root(script)
out <- file.path(root, "results", "real_data", "full")
datasets <- c("pbc", "framingham", "gbsg")
selected <- datasets
for (a in commandArgs(trailingOnly = TRUE)) {
  if (startsWith(a, "--output=")) out <- sub("^--output=", "", a) else
    if (startsWith(a, "--dataset=")) selected <- strsplit(sub("^--dataset=", "", a), ",", fixed = TRUE)[[1]] else
      stop("Use --dataset=pbc,framingham,gbsg or --output=/path.")
}
if (!nzchar(out) || !length(selected) || any(!selected %in% datasets)) stop("Invalid output or dataset.")
paper_load(root)
required <- unique(c(
  if ("framingham" %in% selected) c("cmprsk", "riskCommunicator"),
  if ("gbsg" %in% selected) "pec"))
for (p in required) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Install ", p, " first.")
}
paper_manifest(root, out, list(analyses = c("PBC seed 1234567", "Framingham examination 3, 12 years",
                                           "GBSG2 seed 2025, 10% case-cohort, aligned weights")))
out <- normalizePath(out)
dir.create(file.path(out, "logs"), showWarnings = FALSE)
dir.create(file.path(out, "tables"), showWarnings = FALSE)
status_file <- file.path(out, "run_status.csv")
status <- if (file.exists(status_file)) utils::read.csv(status_file) else
  data.frame(Dataset = character(), Status = character(), Exit_code = integer(),
              Output = character(), Log = character())
table_files <- c(pbc = "metrics_after.csv", framingham = "metrics.csv", gbsg = "metrics.csv")
for (d in unique(selected)) {
  message("Running ", d, " with this checkout...")
  target <- file.path(out, d)
  log <- file.path(out, "logs", paste0(d, ".log"))
  code <- system2(file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(file.path(root, "analysis", paste0("run_", d, ".R"))),
       shQuote(paste0("--output=", target))), stdout = log, stderr = log)
  metric_file <- file.path(target, table_files[[d]])
  okay <- identical(as.integer(code), 0L) && file.exists(metric_file)
  status <- status[status$Dataset != d, ]
  status <- rbind(status, data.frame(Dataset = d, Status = if (okay) "PASS" else "FAILED",
                                     Exit_code = as.integer(code), Output = target, Log = log))
  utils::write.csv(status, status_file, row.names = FALSE)
  if (okay) {
    tab <- utils::read.csv(metric_file, check.names = FALSE)
    stopifnot(all(is.finite(as.matrix(tab[, -1]))))
    utils::write.csv(tab, file.path(out, "tables", paste0(d, "_metrics.csv")), row.names = FALSE)
    paper_write_md(tab, file.path(out, "tables", paste0(d, "_metrics.md")))
    # A tabular fragment that can be included in the manuscript.
    header <- gsub("_", "\\_", names(tab), fixed = TRUE)
    tex <- c(paste0("\\begin{tabular}{l", paste(rep("r", ncol(tab)-1), collapse = ""), "}"),
              "\\toprule", paste0(paste(header, collapse = " & "), " \\\\"), "\\midrule")
    for (i in seq_len(nrow(tab))) {
      row <- c(gsub("_", "\\_", tab$Metric[i], fixed = TRUE),
                formatC(as.numeric(tab[i, -1]), format = "f", digits = 4))
      tex <- c(tex, paste0(paste(row, collapse = " & "), " \\\\"))
    }
    writeLines(c(tex, "\\bottomrule", "\\end{tabular}"),
                file.path(out, "tables", paste0(d, "_metrics.tex")))
    message(d, " PASS")
  } else {
    message(d, " FAILED. See ", log)
  }
}
print(status, row.names = FALSE)
cat("Real-data outputs:", out, "\n")
if (any(status$Status[status$Dataset %in% selected] != "PASS")) quit(save = "no", status = 1)
