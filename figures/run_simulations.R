#!/usr/bin/env Rscript
# Preferred simulation entry: same computation/resume driver, original graphics.
f <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
if (!file.exists(f)) f <- gsub("~+~", " ", f, fixed = TRUE)
root <- normalizePath(file.path(dirname(f), ".."))
args <- commandArgs(trailingOnly = TRUE)
status <- system2(file.path(R.home("bin"), "Rscript"),
  c("--vanilla", shQuote(file.path(root,"analysis","run_simulations.R")), vapply(args, shQuote, character(1))))
mode <- "plan"; out <- NULL
for (a in args) {
  if (startsWith(a,"--mode=")) mode <- sub("^--mode=","",a)
  if (startsWith(a,"--output=")) out <- sub("^--output=","",a)
}
if (mode == "plan") quit(save="no",status=0)
source(file.path(root,"analysis","common.R"))
source(file.path(root,"figures","original_format.R"))
if (is.null(out)) out <- file.path(root,"results","simulations",mode)
paper_preserve_simulation_diagnostics(out)
paper_original_simulations(root,out,preview = mode == "smoke")
cat("Original-format simulation figures complete when required replications are available.\n")
if (status != 0L) {
  message("Some simulation results are incomplete. Completed parts have manuscript figures; incomplete parts have labeled previews. See completion.csv and run_status.csv.")
  quit(save = "no", status = status)
}
