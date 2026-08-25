#!/usr/bin/env Rscript

# Plot a BUSCO completeness summary from an explicit tabular source file.

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])
helper_path <- file.path("analysis", "figures", "export_utils.R")
if (!file.exists(helper_path)) helper_path <- file.path(dirname(script_path), "..", "figures", "export_utils.R")
source(helper_path)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The ggplot2 package is required: install.packages('ggplot2')")
}

read_busco_summary <- function(path) {
  data <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("sample", "single_copy", "duplicated", "fragmented", "missing")
  if (!all(required %in% names(data))) {
    stop(paste("BUSCO summary must contain", paste(required, collapse = ", ")))
  }
  for (name in required[-1L]) data[[name]] <- as.numeric(data[[name]])
  if (any(!is.finite(as.matrix(data[, required[-1L]]))) || any(data[, required[-1L]] < 0)) {
    stop("BUSCO counts must be finite non-negative numbers")
  }
  if (anyDuplicated(data$sample)) stop("BUSCO sample names must be unique")
  data
}

build_busco_plot <- function(data) {
  count_columns <- c("single_copy", "duplicated", "fragmented", "missing")
  totals <- rowSums(data[, count_columns, drop = FALSE])
  if (any(totals <= 0)) stop("Every BUSCO sample must contain at least one assessed ortholog")
  categories <- c("Single-copy", "Duplicated", "Fragmented", "Missing")
  long <- data.frame(
    sample = rep(data$sample, each = length(categories)),
    category = factor(rep(categories, times = nrow(data)), levels = rev(categories)),
    count = as.vector(t(as.matrix(data[, count_columns, drop = FALSE]))),
    total = rep(totals, each = length(categories)),
    stringsAsFactors = FALSE
  )
  long$percentage <- 100 * long$count / long$total
  annotation <- data.frame(
    sample = data$sample,
    label = sprintf(
      "C:%d [S:%d, D:%d], F:%d, M:%d, n:%d",
      data$single_copy + data$duplicated,
      data$single_copy,
      data$duplicated,
      data$fragmented,
      data$missing,
      totals
    ),
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot(long, ggplot2::aes(x = sample, y = percentage, fill = category)) +
    ggplot2::geom_col(width = 0.68, position = "stack") +
    ggplot2::geom_text(
      data = annotation,
      ggplot2::aes(x = sample, y = 2, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      size = 2.7,
      family = "sans"
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 20),
      labels = function(value) paste0(value, "%"),
      expand = c(0, 0)
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "Single-copy" = "#56B4E9",
        "Duplicated" = "#3492C7",
        "Fragmented" = "#F0E442",
        "Missing" = "#F04442"
      ),
      breaks = categories
    ) +
    ggplot2::labs(title = "BUSCO assessment results", x = NULL, y = "% BUSCOs", fill = NULL) +
    ggplot2::theme_classic(base_size = 7, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 8),
      legend.position = "top",
      legend.text = ggplot2::element_text(size = 6),
      axis.text = ggplot2::element_text(color = "black"),
      axis.line = ggplot2::element_line(linewidth = 0.35)
    )
}

main <- function() {
  args <- parse_cli(
    commandArgs(trailingOnly = TRUE),
    list(width_mm = "120", height_mm = "75", dpi = "600", formats = "pdf,svg,png,tiff")
  )
  require_cli(args, c("summary", "output_prefix"))
  plot <- build_busco_plot(read_busco_summary(args$summary))
  export_figure(
    function() print(plot),
    args$output_prefix,
    args$width_mm,
    args$height_mm,
    args$dpi,
    args$formats
  )
}

if (sys.nframe() == 0L) main()
