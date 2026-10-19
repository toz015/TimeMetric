#!/usr/bin/env Rscript
# Preferred simulation entry: same computation/resume driver, original graphics.
f <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
if (!file.exists(f)) f <- gsub("~+~", " ", f, fixed = TRUE)
root <- normalizePath(file.path(dirname(f), ".."))
args <- commandArgs(trailingOnly = TRUE)
status <- system2(file.path(R.home("bin"), "Rscript"),
  c("--vanilla", shQuote(file.path(root,"analysis","run_simulations.R")), vapply(args, shQuote, character(1))))
if (status != 0) quit(save = "no", status = status)
mode <- "plan"; out <- NULL
for (a in args) {
  if (startsWith(a,"--mode=")) mode <- sub("^--mode=","",a)
  if (startsWith(a,"--output=")) out <- sub("^--output=","",a)
}
if (mode == "plan") quit(save="no",status=0)
source(file.path(root,"analysis","common.R"))
source(file.path(root,"figures","original_format.R"))
paper_load(root)
if (is.null(out)) out <- file.path(root,"results","simulations",mode)
paper_original_simulations(root,out,preview = mode == "smoke")
cat("Original-format simulation figures complete when required replications are available.\n")
