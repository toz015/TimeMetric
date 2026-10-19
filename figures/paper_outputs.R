# Publish manuscript-ready real-data figures and tables into one flat folder.
# Dataset-specific model fits and diagnostics remain in the recorded run folder.
paper_publish_real_outputs <- function(root, run_dir,
    destination = file.path(root, "results", "real_data", "paper"),
    datasets = c("pbc", "framingham", "gbsg")) {
  datasets <- unique(datasets)
  if (!length(datasets) || any(!datasets %in% c("pbc", "framingham", "gbsg")))
    stop("Unknown real-data dataset to publish.")
  run_dir <- normalizePath(run_dir, mustWork = TRUE)
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  destination <- normalizePath(destination, mustWork = TRUE)
  stage <- tempfile("paper_export_")
  dir.create(stage)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  numbers <- c(pbc = 2L, framingham = 3L, gbsg = 4L)
  figures <- c(pbc = "paper.fig1", framingham = "paper.fig2", gbsg = "paper.fig3")
  metric_order <- c("r_square", "l_square", "pseudo_r2", "harrell_c", "uno_c",
                    "c_index", "brier_score", "td_auc", "ici")
  labels <- c(r_square = "$R^2$", l_square = "$L^2$", pseudo_r2 = "Pseudo-$R^2$",
    harrell_c = "Harrell's C-index", uno_c = "Uno's C-index", c_index = "C-index",
    brier_score = "Brier Score", td_auc = "Time-dependent AUC", ici = "ICI")
  headers <- list(pbc = c("Cox", "Cox (Reduced)", "Weibull AFT"),
    framingham = c("Fine--Gray", "Fine--Gray (Reduced)", "CSH-AFT",
      "CSH-AFT ($\\text{scale} = 5$)"),
    gbsg = c("Cox (Case--cohort)", "Cox (Random sample)"))
  records <- list()
  for (d in datasets) {
    source_table <- file.path(run_dir, "tables", paste0(d, "_metrics.csv"))
    horizon <- NA_real_
    if (d == "framingham") {
      active <- utils::read.csv(file.path(run_dir, d, "active_horizon.csv"))
      if (nrow(active) != 1L || active$Horizon_years != 12)
        stop("The manuscript export requires the selected 12-year FHS result.")
      horizon <- active$tau_days
      result <- readRDS(active$Source_predictions)
      table <- result$metrics
      consolidated <- utils::read.csv(source_table, check.names = FALSE)
      if (!isTRUE(all.equal(table, consolidated, check.attributes = FALSE)))
        stop("The FHS active-horizon result and consolidated table disagree.")
    } else {
      table <- utils::read.csv(source_table, check.names = FALSE)
    }
    wanted <- metric_order[metric_order %in% table$Metric]
    table <- table[match(wanted, table$Metric), , drop = FALSE]
    if (!all(c("r_square", "l_square", "pseudo_r2", "brier_score", "td_auc") %in% wanted) ||
        ncol(table) - 1L != length(headers[[d]]) ||
        !all(is.finite(as.matrix(table[, -1])))) stop("Incomplete manuscript table for ", d)
    if (!"ici" %in% wanted)
      stop("The manuscript table requires ICI for ", d)
    # The manuscript displays interval R/L/pseudo-R2; point-version diagnostics
    # remain in the analysis folder, rather than adding an unrequested table row.
    stem <- paste0("paper.table", numbers[[d]], ".", d)
    utils::write.csv(table, file.path(stage, paste0(stem, ".csv")), row.names = FALSE)
    readable <- table
    readable$Metric <- unname(labels[readable$Metric])
    names(readable)[-1] <- headers[[d]]
    readable[-1] <- lapply(readable[-1], function(x)
      formatC(round(x, digits = 2), format = "f", digits = 2))
    paper_write_md(readable, file.path(stage, paste0(stem, ".md")))
    tex <- c(paste0("\\begin{tabular}{l", strrep("c", ncol(table) - 1L), "}"),
      "\\toprule", paste0("\\textbf{Metric} & ",
        paste(paste0("\\textbf{", headers[[d]], "}"), collapse = " & "), " \\\\"),
      "\\midrule")
    for (i in seq_len(nrow(table))) {
      label <- labels[[table$Metric[i]]]
      if (table$Metric[i] %in% c("r_square", "l_square", "pseudo_r2", "ici"))
        label <- paste0("\\rev{", label, "}")
      values <- paste0("\\rev{", formatC(round(as.numeric(table[i, -1]), digits = 2),
        format = "f", digits = 2), "}")
      tex <- c(tex, paste0(paste(c(label, values), collapse = " & "), " \\\\"))
    }
    writeLines(c(tex, "\\bottomrule", "\\end{tabular}"),
      file.path(stage, paste0(stem, ".tex")))
    for (ext in c("pdf", "png")) {
      from <- file.path(run_dir, d, paste0(figures[[d]], ".", ext))
      if (!file.exists(from) || file.info(from)$size <= 0) stop("Missing figure: ", from)
      if (!file.copy(from, file.path(stage, basename(from)))) stop("Cannot stage figure: ", from)
    }
    records[[d]] <- data.frame(Dataset = d, Figure = paste0(figures[[d]], ".pdf"),
      Table = paste0(stem, ".tex"), Source_run = run_dir,
      Source_table = normalizePath(source_table), Source_table_md5 = unname(tools::md5sum(source_table)),
      Horizon_days = horizon, Exported_at = format(Sys.time(), tz = "America/Los_Angeles", usetz = TRUE))
  }
  # Validate every requested dataset before replacing any managed output file.
  files <- list.files(stage, full.names = TRUE)
  if (!all(file.copy(files, file.path(destination, basename(files)), overwrite = TRUE)))
    stop("Could not publish all manuscript files.")
  manifest <- do.call(rbind, records)
  manifest_file <- file.path(destination, "source_manifest.csv")
  if (file.exists(manifest_file)) {
    previous <- utils::read.csv(manifest_file, stringsAsFactors = FALSE)
    manifest <- rbind(previous[!previous$Dataset %in% datasets, ], manifest)
  }
  utils::write.csv(manifest, manifest_file, row.names = FALSE)
  managed <- unlist(lapply(seq_len(nrow(manifest)), function(i) {
    c(paste0(sub("\\.pdf$", "", manifest$Figure[i]), c(".pdf", ".png")),
      paste0(sub("\\.tex$", "", manifest$Table[i]), c(".csv", ".md", ".tex")))
  }))
  paths <- file.path(destination, managed)
  if (any(!file.exists(paths))) stop("An existing manuscript export is incomplete.")
  utils::write.csv(data.frame(File = managed, Bytes = file.info(paths)$size,
    MD5 = unname(tools::md5sum(paths))), file.path(destination, "file_checksums.csv"), row.names = FALSE)
  writeLines(c("# Real-data manuscript figures and tables", "",
    "All final real-data paper figures and tables are in this folder.", "",
    "- Figure 1 / Table 2: PBC (paper.fig1 and paper.table2.pbc).",
    "- Figure 2 / Table 3: Framingham, examination 3, 12 years (paper.fig2 and paper.table3.framingham).",
    "- Figure 3 / Table 4: GBSG2 case-cohort (paper.fig3 and paper.table4.gbsg).", "",
    "Figures: PDF for the manuscript, PNG for preview.",
    "Tables: CSV values, Markdown preview, and LaTeX tabular fragments.",
    "LaTeX fragments use the manuscript's existing rev macro for blue revision text.",
    "Point-version metrics, model files and diagnostics remain in the source run folders.",
    "See source_manifest.csv for the source run of each dataset and file_checksums.csv for integrity.", "",
    "Recompute and refresh: Rscript --vanilla figures/run_real_data.R (from TimeMetric-paper)."),
    file.path(destination, "README.md"))
  invisible(manifest)
}
