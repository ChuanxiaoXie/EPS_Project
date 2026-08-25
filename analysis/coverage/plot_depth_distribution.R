#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3 || length(args) > 4) {
  stop("Usage: plot_depth_distribution.R INPUT_TABLE OUTPUT_PDF OUTPUT_PNG [TAIL_FRACTION]")
}

input_table <- args[[1]]
output_pdf <- args[[2]]
output_png <- args[[3]]
tail_fraction <- if (length(args) == 4) as.numeric(args[[4]]) else 0.005

depth_table <- read.table(input_table, header = FALSE, sep = "", stringsAsFactors = FALSE)
if (ncol(depth_table) < 2) {
  stop("Input must contain depth and base-fraction columns")
}
depth <- as.numeric(depth_table[[1]])
fraction <- as.numeric(depth_table[[2]])
if (any(!is.finite(depth)) || any(!is.finite(fraction)) || any(fraction < 0)) {
  stop("Depth and fraction columns must contain finite non-negative values")
}
ordering <- order(depth)
depth <- depth[ordering]
fraction <- fraction[ordering]

cutoff <- 100
for (candidate in seq(40, 100, by = 10)) {
  if (sum(fraction[depth > candidate]) <= tail_fraction) {
    cutoff <- candidate
    break
  }
}
inside <- depth <= cutoff
plot_depth <- c(depth[inside], cutoff)
plot_fraction <- c(fraction[inside], sum(fraction[!inside])) * 100
plot_cumulative <- cumsum(plot_fraction)

draw_depth_plot <- function() {
  par(mar = c(4.5, 4.5, 2.5, 4.5))
  plot(
    plot_depth,
    plot_fraction,
    col = "#2166ac",
    type = "l",
    lwd = 3,
    xlab = "Sequencing depth",
    ylab = "Fraction of bases (%)",
    bty = "l",
    xlim = c(0, cutoff)
  )
  par(new = TRUE)
  plot(
    plot_depth,
    plot_cumulative,
    col = "#b2182b",
    type = "l",
    lwd = 3,
    axes = FALSE,
    xlab = "",
    ylab = "",
    ylim = c(0, 100)
  )
  axis(side = 4)
  mtext("Cumulative fraction of bases (%)", side = 4, line = 3)
  legend(
    "right",
    legend = c("Depth distribution", "Cumulative distribution"),
    lty = 1,
    lwd = 3,
    col = c("#2166ac", "#b2182b"),
    bty = "n"
  )
}

dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_png), recursive = TRUE, showWarnings = FALSE)
pdf(output_pdf, width = 8, height = 6)
draw_depth_plot()
dev.off()
png(output_png, width = 1600, height = 1200, res = 200)
draw_depth_plot()
dev.off()
