#!/usr/bin/env Rscript
# Preferred real-data entry: recompute all requested analyses, then draw the
# original manuscript figures using the package's own plotting functions.
f <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
if (!file.exists(f)) f <- gsub("~+~", " ", f, fixed = TRUE)
root <- normalizePath(file.path(dirname(f), ".."))
args <- commandArgs(trailingOnly = TRUE)
# Computation runs stay separate; final figures/tables always share one folder.
paper_out <- file.path(root, "results", "real_data", "paper")
for (a in args) if (startsWith(a, "--paper-output="))
  paper_out <- sub("^--paper-output=", "", a)
args <- args[!startsWith(args, "--paper-output=")]
if (!any(startsWith(args, "--output="))) {
  run_base <- file.path(root, "results", "real_data", "runs",
    format(Sys.time(), "%Y-%m-%d_%H%M%S", tz = "America/Los_Angeles"))
  run_dir <- run_base
  index <- 1L
  while (dir.exists(run_dir)) {
    run_dir <- paste0(run_base, "_", index)
    index <- index + 1L
  }
  args <- c(args, paste0("--output=", run_dir))
}
status <- system2(file.path(R.home("bin"), "Rscript"),
  c("--vanilla", shQuote(file.path(root,"analysis","run_real_data.R")), vapply(args, shQuote, character(1))))
if (status != 0) quit(save = "no", status = status)
source(file.path(root,"analysis","common.R"))
source(file.path(root,"figures","original_format.R"))
source(file.path(root,"figures","paper_outputs.R"))
paper_load(root)
out <- NULL
selected <- c("pbc","framingham","gbsg")
for (a in args) {
  if (startsWith(a,"--output=")) out <- sub("^--output=","",a)
  if (startsWith(a,"--dataset=")) selected <- strsplit(sub("^--dataset=","",a),",",fixed=TRUE)[[1]]
}
if (is.null(out)) stop("The computation output directory was not selected.")
if ("framingham" %in% selected) {
  status <- system2(file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(file.path(root, "figures", "run_fhs_horizon.R")), "--years=12",
      shQuote(paste0("--input=", file.path(out, "framingham", "framingham_visit3_10y.rds")))))
  if (status != 0) quit(save = "no", status = status)
}
print(paper_original_real(root,out,selected), row.names=FALSE)
print(paper_publish_real_outputs(root,out,paper_out,selected), row.names=FALSE)
cat("Analysis complete. All manuscript figures and tables:", normalizePath(paper_out), "\n")
