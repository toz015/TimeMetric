#!/usr/bin/env Rscript
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
script <- sub("^--file=", "", file_arg[1])
if (!file.exists(script)) script <- gsub("~+~", " ", script, fixed = TRUE)
source(file.path(dirname(normalizePath(script)), "common.R"))
source(file.path(dirname(normalizePath(script)), "framingham.R"))
root <- paper_root(script)
out <- file.path(root, "results", "framingham_visit3_12y")
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  if (length(args) != 1 || !startsWith(args, "--output=")) stop("Use --output=/path.")
  out <- sub("^--output=", "", args)
  if (!nzchar(out)) stop("The output path cannot be empty.")
}
paper_load(root)
paper_run_framingham(root, out)
