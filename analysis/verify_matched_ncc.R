#!/usr/bin/env Rscript
arg <- grep("^--file=", commandArgs(), value = TRUE)
script <- sub("^--file=", "", arg[1])
if (!file.exists(script)) script <- gsub("~+~", " ", script, fixed = TRUE)
source(file.path(dirname(normalizePath(script)), "common.R"))
root <- paper_root(script)
paper_load(root)
out <- file.path(root, "results", "matched_ncc")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
code <- paper_git(root, "show", "426272c2bdcbf94079c035166fc7e3804ad2b632:R/pam.weighted_param.R")
e <- new.env(parent = asNamespace("TimeMetric")); eval(parse(text = code), e)
d <- data.frame(time = c(1,4,5,2,4,5), status = c(1,0,0,1,0,0),
                 strata = c("A","A","A","B","B","B"))
d$before <- e$weighted_param(d$time, d$status, design = "matched_ncc", strata = d$strata, m = 2)
# The historical implementation used m - 1 despite documenting controls.
# Compare the same one-control design under the corrected public API.
d$after <- tm_nested_case_control_weights(d$time, d$status, d$strata, m = 1)
d$expected <- c(1,2,2,1,2,2)
stopifnot(isTRUE(all.equal(d$after, d$expected)))
utils::write.csv(d, file.path(out, "hand_calculated_check.csv"), row.names = FALSE)
paper_write_md(d, file.path(out, "hand_calculated_check.md"))
print(d)
cat("MATCHED NCC PASS: m=1 means one sampled control per case.\n")
