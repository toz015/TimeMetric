# Presentation only. Kept outside analysis/ so changing figure styling cannot
# invalidate or interrupt the existing simulation computation/checkpoints.
# Resolve the recorded source run rather than an obsolete result directory.
paper_current_real_run <- function(root,
    paper_out = file.path(root, "results", "real_data", "paper")) {
  manifest_file <- file.path(paper_out, "source_manifest.csv")
  if (!file.exists(manifest_file))
    stop("No published real-data run. Run figures/run_real_data.R first, or supply an explicit input/output path.")
  manifest <- utils::read.csv(manifest_file, stringsAsFactors = FALSE)
  current_runs <- unique(manifest$Source_run)
  if (length(current_runs) != 1L || is.na(current_runs) || !dir.exists(current_runs))
    stop("Published datasets need one existing source run; supply an explicit input/output path.")
  normalizePath(current_runs, mustWork = TRUE)
}

paper_original_save <- function(plot, dir, stem, width, height) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path(dir, paste0(stem, ".pdf")), plot, width = width, height = height,
                   encoding = "WinAnsi.enc")
  ggplot2::ggsave(file.path(dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = 200)
  invisible(plot)
}

paper_original_real_plot <- function(dataset, predictions, tau = 3650) {
  if (dataset == "pbc") {
    return(TimeMetric::tm_plot_summary(predictions, ncol = 3,
      titles = c("Cox PH", "Cox PH (Reduced)", "Weibull AFT"),
      invert_linear = c(FALSE, FALSE, TRUE)))
  }
  if (dataset == "framingham") {
    return(TimeMetric::tm_plot_summary(predictions,
      titles = c("Fine\u2013Gray", "Fine\u2013Gray (Reduced)",
                 "CSH Weibull AFT", "CSH Weibull AFT (Scale = 5)"),
      invert_linear = c(FALSE, FALSE, TRUE, TRUE),
      label_name = c("Censored", "Event of interest", "Competing event"),
      ncol = 2, sample_size = 200, restrict_time = tau, seed = 1234567))
  }
  if (dataset == "gbsg") {
    return(TimeMetric::tm_plot_summary(predictions,
      titles = c("Cox (Case-cohort)", "Cox (Random sample)"), ncol = 2,
      tag_levels = "a", invert_linear = c(FALSE, FALSE)))
  }
  stop("Unknown real-data figure: ", dataset)
}

paper_original_real <- function(root, out, datasets = c("pbc", "framingham", "gbsg")) {
  files <- c(pbc = "pbc/pbc_comparison.rds",
             framingham = "framingham/framingham_visit3_10y.rds", gbsg = "gbsg/gbsg2.rds")
  stems <- c(pbc = "paper.fig1", framingham = "paper.fig2", gbsg = "paper.fig3")
  sizes <- list(pbc = c(9, 3.5), framingham = c(9, 7), gbsg = c(9, 4))
  aliases <- c(pbc = "pbc", framingham = "framingham_visit3_10y", gbsg = "gbsg2")
  active_file <- file.path(out, "framingham", "active_horizon.csv")
  if ("framingham" %in% datasets && file.exists(active_file)) {
    active <- utils::read.csv(active_file, stringsAsFactors = FALSE)
    fit <- file.path(out, files[["framingham"]])
    if (nrow(active) != 1L || !file.exists(fit) ||
        !identical(active$Fit_source_md5, unname(tools::md5sum(fit)))) {
      stop("FHS fit changed: rerun figures/run_fhs_horizon.R --years=12 with this --input.")
    }
    files[["framingham"]] <- active$Source_predictions
    aliases[["framingham"]] <- paste0("framingham_visit3_", active$Horizon_years, "y")
  }
  records <- list()
  for (d in datasets) {
    path <- if (startsWith(files[[d]], "/")) files[[d]] else file.path(out, files[[d]])
    if (!file.exists(path)) stop("Run the updated real-data analysis first: missing ", path)
    result <- readRDS(path)
    predictions <- if (d == "pbc") result$after else result$predictions
    tau <- if (d == "framingham") result$settings$tau else NULL
    plot <- paper_original_real_plot(d, predictions, tau)
    dest <- file.path(out, d)
    paper_original_save(plot, dest, stems[[d]], sizes[[d]][1], sizes[[d]][2])
    # Keep the previous stable output links, now using the requested native style.
    for (ext in c("pdf", "png")) {
      file.copy(file.path(dest, paste0(stems[[d]], ".", ext)),
                 file.path(dest, paste0(aliases[[d]], ".", ext)), overwrite = TRUE)
    }
    saveRDS(plot, file.path(dest, paste0(stems[[d]], ".plot.rds")))
    plot_data <- do.call(rbind, lapply(seq_along(predictions), function(i) {
      cbind(Model = names(predictions)[i], plot[[i]]$data)
    }))
    utils::write.csv(plot_data, file.path(dest, paste0(stems[[d]], ".plot_data.csv")), row.names = FALSE)
    # The earlier PBC before/after display is a diagnostic, not the paper figure.
    if (d == "pbc") {
      for (ext in c("pdf", "png")) {
        old <- file.path(dest, paste0("pbc_before_after.", ext))
        if (file.exists(old)) {
          dir.create(file.path(dest, "diagnostics"), showWarnings = FALSE)
          moved <- file.path(dest, "diagnostics", basename(old))
          if (!file.copy(old, moved, overwrite = TRUE)) stop("Could not preserve diagnostic.")
          unlink(old)
        }
      }
    }
    records[[d]] <- data.frame(Dataset = d, Figure = paste0(stems[[d]], ".pdf"),
      Function = "TimeMetric::tm_plot_summary", Width_inches = sizes[[d]][1],
      Height_inches = sizes[[d]][2], Units = "Days", Source_predictions = path,
      Plot_format = "Original manuscript", stringsAsFactors = FALSE)
  }
  records <- do.call(rbind, records)
  record_file <- file.path(out, "original_figure_format.csv")
  if (file.exists(record_file)) {
    previous <- utils::read.csv(record_file, stringsAsFactors = FALSE)
    records <- rbind(previous[!previous$Dataset %in% datasets, ], records)
  }
  utils::write.csv(records, record_file, row.names = FALSE)
  invisible(records)
}

paper_original_metric_label <- function(metric) switch(metric,
  pseudo_r2 = expression("Pseudo-" * R^2), r_square = expression(R^2),
  l_square = expression(L^2), brier_score = "Brier score", harrell_c = "Harrell's C",
  uno_c = "Uno's C", td_auc = "AUC", ici = "ICI", metric)

# The following layouts follow the existing manuscript redraw script: bars,
# paired R/L dots in Simulation 2, and three labeled rows of boxplots in Sim 3.
paper_original_simulations <- function(root, out, preview = FALSE,
    parts = c("sim1", "sim2_linear", "sim2_nonlinear", "sim3")) {
  raw_file <- file.path(out, "tables", "metrics_raw.csv")
  if (!file.exists(raw_file)) {
    message("Simulation still computing: no aggregated results yet. Plotting skipped.")
    return(invisible(NULL))
  }
  raw <- utils::read.csv(raw_file, check.names = FALSE)
  scalar <- c("pseudo_r2", "brier_score", "harrell_c", "td_auc", "ici")
  destination <- if (preview) file.path(out, "figures", "original_format_preview") else file.path(out, "figures")
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  counts <- stats::aggregate(is.finite(raw$Value),
    raw[c("part", "scenario", "n", "censor", "quantile", "Model", "Metric")], sum)
  names(counts)[ncol(counts)] <- "finite_repetitions"
  counts <- counts[counts$part %in% parts, ]
  # Simulation 3 ICI is not estimable when no test subject is followed to
  # t_star. Only these recorded cases may be missing from a final figure; their
  # finite counts are written beside Figure 6 for the manuscript caption.
  diagnostics_file <- file.path(out, "tables", "calibration_diagnostics.csv")
  unestimable <- if (file.exists(diagnostics_file)) {
    x <- utils::read.csv(diagnostics_file)
    x <- x[x$part == "sim3" & x$at_risk == 0 &
           grepl("within observed follow-up", x$error, fixed = TRUE), ]
    paste(x$id, x$Model)
  } else character()
  completion_file <- file.path(out, "completion.csv")
  smoke <- file.exists(completion_file) &&
    any(utils::read.csv(completion_file)$mode == "smoke")
  save_plot <- function(p, name, width, height) {
    if (preview) {
      incomplete <- counts[counts$finite_repetitions < 100L, ]
      caption <- if (smoke) {
        "SMOKE PREVIEW: reduced sample sizes and repetitions; not manuscript results."
      } else if (nrow(incomplete)) {
        details <- vapply(seq_len(nrow(incomplete)), function(i) {
          x <- incomplete[i, ]
          sprintf("%s: %s, n=%d, censoring=%g%%, horizon=%g%%: %d/100 finite",
            x$part, x$Metric, x$n, 100*x$censor, x$quantile,
            x$finite_repetitions)
        }, character(1))
        paste("INCOMPLETE PREVIEW -", paste(details, collapse = "; "))
      } else "PREVIEW: design or replication checks are incomplete; not manuscript results."
      p <- p + patchwork::plot_annotation(caption = caption,
        theme = ggplot2::theme(plot.caption = ggplot2::element_text(
          hjust = 0, size = 9, margin = ggplot2::margin(t = 8))))
      height <- height + 0.3
    }
    paper_original_save(p, destination, name, width, height)
    # Keep both stable manuscript links on the selected .updated.pdf design.
    for (ext in c("pdf", "png")) {
      if (!file.copy(file.path(destination, paste0(name, ".", ext)),
                     file.path(destination, paste0(sub(".updated", "", name, fixed = TRUE), ".", ext)),
                     overwrite = TRUE)) stop("Could not refresh figure alias: ", name)
    }
  }
  ready <- function(part) {
    if (!part %in% parts) return(FALSE)
    d <- raw[raw$part == part, ]
    if (!nrow(d)) return(FALSE)
    required <- c(scalar, "r_square", "l_square", "uno_c")
    designs <- split(d$Metric, interaction(
      d[c("scenario", "n", "censor", "quantile", "Model")], drop = TRUE))
    if (!all(vapply(designs, function(x) setequal(x, required), logical(1)))) {
      message(part, ": required metric rows including ICI are missing. Figure skipped.")
      return(FALSE)
    }
    if (preview) return(TRUE)
    keys <- c("scenario", "n", "censor", "quantile", "Model", "Metric")
    d$usable <- is.finite(d$Value) |
      (d$Metric == "ici" & paste(d$id, d$Model) %in% unestimable)
    cells <- split(d, interaction(d[keys], drop = TRUE))
    if (!all(vapply(cells, function(x) nrow(x) == 100L && setequal(x$rep, 1:100) &&
                   all(x$usable), logical(1)))) {
      message(part, ": incomplete results; manuscript figure skipped; drawing a labeled preview instead.")
      paper_original_simulations(root, out, preview = TRUE, parts = part)
      return(FALSE)
    }
    if (part == "sim1" && length(unique(d$scenario)) != 16L) return(FALSE)
    if (part == "sim3" && nrow(unique(d[c("n", "censor", "quantile")])) != 12L) return(FALSE)
    TRUE
  }
  row_header <- function(expr_label, hjust) {
    ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::annotate(
      "text", x = 0, y = 1, label = expr_label, parse = TRUE,
      hjust = hjust, vjust = 1, fontface = "bold", size = 4.2) +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
      ggplot2::theme(plot.margin = ggplot2::margin(t = 2, r = 0, b = 4, l = 20))
  }
  if (ready("sim1")) {
    d <- raw[raw$part == "sim1" & raw$Metric %in% scalar, ]
    means <- stats::aggregate(Value ~ Metric + v + beta, d, mean)
    historical <- means$Metric != "ici"
    means$Value[historical] <- round(means$Value[historical], 3)
    groups <- lapply(c("beta", "inv_v"), function(variable) {
      plots <- lapply(scalar, function(metric) {
        d <- means[means$Metric == metric, ]
        d <- if (variable == "beta") d[d$v == 1, ] else d[d$beta == 1, ]
        d$x <- if (variable == "beta") d$beta else 1/d$v
        d$x <- factor(d$x, levels = sort(unique(d$x)))
        limits <- range(means$Value[means$Metric == metric])
        ggplot2::ggplot(d, ggplot2::aes(x, Value)) + ggplot2::geom_col(width = 0.7) +
          ggplot2::theme_bw() +
          ggplot2::labs(x = if (variable == "beta") expression(beta) else "1/v",
            y = if (metric == "harrell_c") "C-index" else paper_original_metric_label(metric)) +
          ggplot2::coord_cartesian(ylim = limits) +
          ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.1, 0.17))) +
          ggplot2::scale_x_discrete(drop = FALSE) + ggplot2::theme(
            axis.title.x = ggplot2::element_text(size = 9), axis.title.y = ggplot2::element_text(size = 9))
      })
      header <- if (variable == "beta")
        row_header("bold('(a) Vary the effect size through ') * beta", 0.3) else
          row_header("bold('(b) Vary the variance of event times through ') * 1/v", 0.203)
      # Keep the manuscript's four existing panels in their 2x2 positions;
      # append ICI in the right-hand column without moving C-index or AUC.
      patchwork::wrap_plots(header,
        patchwork::wrap_plots(plots, design = "ABE\nCD#"),
        ncol = 1, heights = c(0.08, 1))
    })
    save_plot(patchwork::wrap_plots(groups, ncol = 1), "paper.fig4.updated", 10.5, 7)
    utils::write.csv(means, file.path(destination, "paper.fig4.plot_data.csv"), row.names = FALSE)
  }
  for (part in c("sim2_linear", "sim2_nonlinear")) {
    if (!ready(part)) next
    linear <- part == "sim2_linear"
    order <- if (linear) c("Cox", "Cox.reduce", "Ln", "Ln.fix") else c("RF", "Cox", "Ln")
    labels <- if (linear) c("Cox", "Cox (Reduced)", "Ln", "Ln (Scale = 2)") else order
    d <- raw[raw$part == part & raw$Model %in% order & raw$Metric %in% c(scalar, "r_square", "l_square"), ]
    means <- stats::aggregate(Value ~ Metric + Model, d, mean)
    means$Model <- factor(means$Model, levels = order, labels = labels)
    dots <- means[means$Metric %in% c("r_square", "l_square"), ]
    dots$Component <- factor(dots$Metric, levels = c("r_square", "l_square"))
    original_theme <- ggplot2::theme_bw() + ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 11),
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, size = 9),
      axis.title.y = ggplot2::element_text(size = 9), legend.position = "none")
    bw <- if (linear) 0.7 else 0.5
    headroom <- if (linear) 0.15 else 0.2
    plots <- list(ggplot2::ggplot(means[means$Metric == "pseudo_r2", ], ggplot2::aes(Model, Value)) +
      ggplot2::geom_col(width = bw) +
      ggplot2::geom_point(data = dots, ggplot2::aes(shape = Component), size = 2.5,
                         colour = "black", show.legend = FALSE) +
      ggplot2::scale_shape_manual(values = c(16, 1), guide = "none") +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, headroom))) +
      ggplot2::labs(x = NULL, y = "Mean Value", title = "Pseudo-R\u00b2") + original_theme)
    for (metric in scalar[-1]) {
      plots[[length(plots)+1L]] <- ggplot2::ggplot(means[means$Metric == metric, ], ggplot2::aes(Model, Value)) +
        ggplot2::geom_col(width = bw) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, headroom))) +
        ggplot2::labs(x = NULL, y = "Mean Value",
          title = if (metric == "harrell_c") "C-index" else paper_original_metric_label(metric)) + original_theme
    }
    p <- patchwork::wrap_plots(plots, design = "ABE\nCD#") + patchwork::plot_annotation(
      theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold")))
    stem <- if (linear) "paper.fig5.updated" else "paper.fig5.2.updated"
    # Preserve each existing panel's physical dimensions as ICI is appended.
    width <- if (linear) 15 else 12
    height <- if (linear) 10 * 3.5/5.5 else 8 * 3.5/5.5
    save_plot(p, stem, width, height)
    utils::write.csv(means, file.path(destination, paste0(stem, ".plot_data.csv")), row.names = FALSE)
  }
  if (ready("sim3")) {
    error_file <- file.path(out, "tables", "sim3_errors.csv")
    d <- utils::read.csv(error_file)
    box_metrics <- c("pseudo_r2", "brier_score", "harrell_c", "uno_c", "td_auc", "ici")
    if (!all(box_metrics %in% d$Metric))
      stop("Simulation 3 error table is missing metrics including ICI. Recompute its summaries.")
    d <- d[d$Metric %in% box_metrics, ]
    # Use the NEW error calculations. Styling never recalculates their reference.
    d$quantile <- factor(d$quantile, levels = c(90, 50))
    censoring <- sort(unique(d$censor))
    finite <- stats::aggregate(list(B = is.finite(d$error)),
      d[c("Metric", "n", "censor", "quantile")], sum)
    partial <- finite[finite$B < 100L, ]
    rows <- lapply(seq_along(censoring), function(i) {
      label <- paste0("(", letters[i], ") censoring ", round(100*censoring[i]), "%")
      label_plot <- ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 1, label = label,
        hjust = -0.2, vjust = -1.5, fontface = "bold", size = 4.8) +
        ggplot2::xlim(0,1) + ggplot2::ylim(0,1) + ggplot2::coord_cartesian(clip = "off") +
        ggplot2::theme_void(base_size = 12) +
        ggplot2::theme(plot.margin = ggplot2::margin(t = 22, r = 16, b = 4, l = 14))
      plots <- lapply(box_metrics, function(metric) {
        sub <- d[d$Metric == metric & d$censor == censoring[i], ]
        limits <- range(d$error[d$Metric == metric], finite = TRUE)
        ggplot2::ggplot(sub, ggplot2::aes(factor(n), error, fill = quantile)) +
          ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
          ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6, colour = "black", na.rm = TRUE,
            position = ggplot2::position_dodge(width = 0.75)) +
          ggplot2::scale_fill_manual(values = c("90" = "grey70", "50" = "white"),
            labels = c("90%", "50%"), name = "Quantile") +
          ggplot2::scale_y_continuous(limits = limits, expand = ggplot2::expansion(mult = c(0.05,0.05))) +
          ggplot2::labs(x = NULL, y = paper_original_metric_label(metric)) +
          ggplot2::theme_bw(base_size = 12) + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(6,6,6,6))
      })
      patchwork::wrap_plots(c(list(label_plot), plots), ncol = length(box_metrics)+1,
                            widths = c(0.1, rep(1,length(box_metrics))))
    })
    big <- (Reduce(`/`, rows) + patchwork::plot_layout(heights = rep(1,length(rows)))) &
      ggplot2::theme(legend.position = "none")
    xlab <- ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Sample size", size = 5) +
      ggplot2::theme_void()
    save_plot(patchwork::wrap_plots(big, xlab, ncol = 1, heights = c(1,0.05)), "paper.fig6.updated", 15.5, 6)
    # Boxes from fewer than 100 finite repetitions are reported for the
    # manuscript caption rather than drawn on the figure.
    utils::write.csv(partial, file.path(destination, "paper.fig6.partial_repetitions.csv"),
                     row.names = FALSE)
  }
  if (preview) {
    utils::write.csv(counts, file.path(destination, "replication_counts.csv"), row.names = FALSE)
    writeLines(if (smoke)
      "Reduced smoke-test results; NOT manuscript simulation results." else
      "Incomplete simulation results. See replication_counts.csv and the figure caption; NOT final manuscript results.",
      file.path(destination, "PREVIEW_ONLY.txt"))
  }
  invisible(destination)
}

# The calculation driver creates generic diagnostics before it reports any
# failed jobs. Preserve those separately so they cannot be mistaken for paper
# figures when a later manuscript-format figure is incomplete.
paper_preserve_simulation_diagnostics <- function(out) {
  figures <- file.path(out, "figures")
  stems <- c("paper.fig4", "paper.fig5", "paper.fig5.2", "paper.fig6")
  for (stem in stems) {
    for (ext in c("pdf", "png")) {
      source <- file.path(figures, paste0(stem, ".updated.", ext))
      alias <- file.path(figures, paste0(stem, ".", ext))
      # The native renderer always saves both aliases. The calculation driver
      # saves only .updated; unequal files also identify a fresh diagnostic.
      if (!file.exists(source)) next
      if (file.exists(alias) && identical(unname(tools::md5sum(source)),
                                         unname(tools::md5sum(alias)))) next
      destination <- file.path(figures, "diagnostics", "analysis_format")
      dir.create(destination, recursive = TRUE, showWarnings = FALSE)
      target <- file.path(destination, basename(source))
      if (!file.copy(source, target, overwrite = TRUE))
        stop("Could not preserve simulation diagnostic: ", source)
      if (file.exists(alias)) {
        previous <- file.path(figures, "diagnostics", "previous_manuscript")
        dir.create(previous, recursive = TRUE, showWarnings = FALSE)
        if (!file.copy(alias, file.path(previous, basename(alias)), overwrite = TRUE))
          stop("Could not preserve previous manuscript figure: ", alias)
      }
      unlink(source)
      if (file.exists(alias)) unlink(alias)
    }
  }
  invisible(NULL)
}
