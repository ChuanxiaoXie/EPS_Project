#!/usr/bin/env Rscript

# Plot a chromosome-partitioned whole-genome Hi-C contact heatmap from X/Y/Z data.

parse_args <- function(values) {
  result <- list(
    width = 10,
    height = 10,
    dpi = 300,
    label_offset = 25,
    separator_size = 0.2,
    max_log_contact = 4
  )
  index <- 1
  while (index <= length(values)) {
    key <- sub("^--", "", values[[index]])
    if (index == length(values)) stop(paste("Missing value for", values[[index]]))
    result[[gsub("-", "_", key)]] <- values[[index + 1]]
    index <- index + 2
  }
  required <- c("matrix", "breaks", "labels", "output_pdf", "output_png")
  missing <- required[!required %in% names(result)]
  if (length(missing) > 0) stop(paste("Missing arguments:", paste(missing, collapse = ", ")))
  for (name in c("width", "height", "dpi", "label_offset", "separator_size", "max_log_contact")) {
    result[[name]] <- as.numeric(result[[name]])
  }
  result
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
suppressPackageStartupMessages({
  library(ggplot2)
  library(RColorBrewer)
  library(grid)
})

contacts <- read.table(args$matrix, header = TRUE, sep = "", check.names = FALSE)
if (!all(c("X", "Y", "Z") %in% names(contacts))) {
  stop("The contact table must contain X, Y and Z columns")
}
chromosome_breaks <- scan(args$breaks, what = numeric(), quiet = TRUE)
chromosome_labels <- scan(args$labels, what = character(), quiet = TRUE)
if (length(chromosome_breaks) != length(chromosome_labels) + 1) {
  stop("The breaks file must contain one more value than the labels file")
}

contacts$plot_value <- contacts$Z
legend_label <- "Contact value"
nonnegative_contacts <- min(contacts$Z, na.rm = TRUE) >= 0
if (nonnegative_contacts) {
  contacts$plot_value <- log10(contacts$Z + 1)
  legend_label <- "log10(KR + 1)"
}

palette_seed <- colorRampPalette(c(brewer.pal(9, "YlOrRd"), "black"), bias = 1)(9)
palette_values <- colorRampPalette(palette_seed)(100)
label_positions <- (chromosome_breaks[-1] + chromosome_breaks[-length(chromosome_breaks)]) / 2
label_data <- data.frame(label = chromosome_labels, position = label_positions)
maximum_coordinate <- max(c(contacts$X, contacts$Y), na.rm = TRUE)
maximum_axis <- max(maximum_coordinate, chromosome_breaks)
effective_offset <- min(abs(args$label_offset), maximum_coordinate * 0.08)
minimum_axis <- -effective_offset

fill_scale <- if (nonnegative_contacts) {
  scale_fill_gradientn(
    colours = palette_values,
    name = legend_label,
    limits = c(0, args$max_log_contact),
    oob = scales::squish
  )
} else {
  scale_fill_gradientn(colours = palette_values, name = legend_label)
}

plot <- ggplot(contacts, aes(x = X, y = Y, fill = plot_value)) +
  geom_tile() +
  geom_vline(xintercept = chromosome_breaks, colour = "black", linetype = "dashed", linewidth = args$separator_size) +
  geom_hline(yintercept = chromosome_breaks, colour = "black", linetype = "dashed", linewidth = args$separator_size) +
  geom_text(
    data = label_data,
    aes(x = position, y = minimum_axis, label = label),
    inherit.aes = FALSE,
    size = 2.8,
    fontface = "bold",
    angle = 45,
    hjust = 1
  ) +
  geom_text(
    data = label_data,
    aes(x = minimum_axis, y = position, label = label),
    inherit.aes = FALSE,
    size = 2.8,
    fontface = "bold",
    hjust = 1
  ) +
  fill_scale +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(minimum_axis, maximum_axis), ylim = c(minimum_axis, maximum_axis), clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.7, "cm"),
    plot.margin = margin(10, 10, 60, 60)
  )

dir.create(dirname(args$output_pdf), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(args$output_png), recursive = TRUE, showWarnings = FALSE)
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
ggsave(args$output_pdf, plot = plot, width = args$width, height = args$height, device = pdf_device)
grDevices::png(args$output_png, width = args$width, height = args$height,
               units = "in", res = args$dpi)
print(plot)
grDevices::dev.off()
