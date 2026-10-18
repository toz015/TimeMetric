#!/usr/bin/env Rscript
# Redraw ONLY; never refit or alter statistical outputs/checkpoints.
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
script <- sub("^--file=", "", file_arg[1])
if (!file.exists(script)) script <- gsub("~+~", " ", script, fixed = TRUE)
root <- normalizePath(file.path(dirname(script), ".."))
source(file.path(root, "analysis", "common.R"))
source(file.path(root, "figures", "original_format.R"))
source(file.path(root, "figures", "paper_outputs.R"))
section <- "all"; real_out <- NULL
paper_out <- file.path(root,"results","real_data","paper")
sim_out <- file.path(root,"results","simulations","full"); preview <- FALSE
for (arg in commandArgs(trailingOnly = TRUE)) {
  if (startsWith(arg,"--section=")) section <- sub("^--section=","",arg) else
    if (startsWith(arg,"--paper-output=")) paper_out <- sub("^--paper-output=","",arg) else
    if (startsWith(arg,"--real-output=")) real_out <- sub("^--real-output=","",arg) else
      if (startsWith(arg,"--sim-output=")) sim_out <- sub("^--sim-output=","",arg) else
        if (arg == "--preview") preview <- TRUE else stop("Unknown argument: ",arg)
}
if (!section %in% c("all","real","simulations")) stop("section: all, real, simulations")
if (section %in% c("all", "real") && !any(startsWith(commandArgs(trailingOnly = TRUE), "--real-output="))) {
  real_out <- paper_current_real_run(root, paper_out)
}
paper_load(root)
if (section %in% c("all","real")) {
  print(paper_original_real(root,real_out), row.names=FALSE)
  print(paper_publish_real_outputs(root,real_out,paper_out), row.names=FALSE)
}
if (section %in% c("all","simulations")) paper_original_simulations(root,sim_out,preview)
cat("Original-format redraw finished. Statistical results were not recomputed.\n")
