#!/usr/bin/env Rscript

# Draw the SNP calling, joint-genotyping and hotspot-filtering workflow.

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])
helper_path <- file.path("analysis", "figures", "export_utils.R")
if (!file.exists(helper_path)) helper_path <- file.path(dirname(script_path), "..", "figures", "export_utils.R")
source(helper_path)

box <- function(x, y, width, height, label, fill = "#DCEAF6", border = "#3182BD", cex = 0.5) {
  graphics::rect(x - width / 2, y - height / 2, x + width / 2, y + height / 2,
                 col = fill, border = border, lwd = 0.8)
  graphics::text(x, y, label, cex = cex, family = "sans")
}

diamond <- function(x, y, width, height, label) {
  graphics::polygon(c(x, x + width / 2, x, x - width / 2),
                    c(y + height / 2, y, y - height / 2, y),
                    col = "#FFF0E6", border = "#F39C6B", lwd = 0.8)
  graphics::text(x, y, label, cex = 0.46)
}

arrow <- function(x0, y0, x1, y1) {
  graphics::arrows(x0, y0, x1, y1, length = 0.045, angle = 22, lwd = 0.75)
}

draw_workflow <- function(window_kb, step_kb, threshold) {
  graphics::par(mar = rep(0.35, 4), family = "sans", xpd = NA)
  graphics::plot.new(); graphics::plot.window(c(0, 100), c(0, 100))
  graphics::abline(v = 50, col = "#333333", lwd = 0.7)

  box(24, 95, 38, 5, "NGS paired-end reads", "#D8E8F5", "#2E9BD6", 0.62)
  box(11, 87, 20, 7, "QC: N >10%; bases Q <20\n>50%; adapter contamination", "#BCD3EA", "#91B2D1", 0.43)
  diamond(33, 87, 14, 8, "Quality\ncontrol")
  arrow(24, 92.5, 31, 91); arrow(21, 87, 26, 87)
  box(24, 79, 38, 5, "Clean reads", "#E7F1F8", "#2E9BD6", 0.57); arrow(33, 83, 28, 81.5)

  box(11, 70, 20, 6, "Reference genome\nCA7301-WT", "#BCD3EA", "#91B2D1", 0.47)
  diamond(33, 70, 15, 8, "Sentieon BWA mem\n-k 32 -M")
  arrow(24, 76.5, 31, 74); arrow(21, 70, 25.5, 70)
  box(24, 61, 38, 5, "Coordinate-sorted BAM and duplicate removal", "#E7F1F8", "#2E9BD6", 0.5)
  arrow(33, 66, 28, 63.5)

  box(11, 51, 21, 8,
      "samtools view: paired, mapped,\nprimary, MAPQ >=20; exclude\nsupplementary; NM <=1; CIGAR ops <=2",
      "#BCD3EA", "#91B2D1", 0.39)
  diamond(33, 51, 15, 8, "Post-dedup\nBAM filter")
  arrow(24, 58.5, 31, 55); arrow(21.5, 51, 25.5, 51)
  box(24, 42, 38, 5, "Filtered BAM", "#E7F1F8", "#2E9BD6", 0.57); arrow(33, 47, 28, 44.5)

  box(11, 33, 21, 7, "Sentieon Haplotyper\nploidy 2; emit/call confidence 30;\ngVCF mode", "#BCD3EA", "#91B2D1", 0.41)
  diamond(33, 33, 15, 8, "Per-sample\ngVCF")
  arrow(24, 39.5, 31, 37); arrow(21.5, 33, 25.5, 33)
  box(24, 24, 38, 5, "Sentieon GVCFtyper: merged population genotypes", "#E7F1F8", "#2E9BD6", 0.5)
  arrow(33, 29, 28, 26.5)
  box(24, 15, 30, 5, "Reference sample genotype = 0/0", "#BCD3EA", "#91B2D1", 0.5)
  arrow(24, 21.5, 24, 17.5)

  workflow_title <- sprintf("Density windows: %d kb; step: %d kb", window_kb, step_kb)
  box(74, 88, 38, 6, workflow_title, "#BCD3EA", "#91B2D1", 0.53)
  graphics::segments(39, 15, 48, 15, lwd = 0.75)
  graphics::segments(48, 15, 48, 88, lwd = 0.75)
  arrow(48, 88, 55, 88); arrow(74, 85, 74, 80)
  box(64, 76, 18, 6, sprintf("Windows with\n> %d SNPs", threshold), "#DCEAF6", "#91B2D1", 0.5)
  box(84, 76, 18, 6, sprintf("Windows with\n<= %d SNPs", threshold), "#DCEAF6", "#91B2D1", 0.5)
  arrow(74, 80, 64, 79); arrow(74, 80, 84, 79)
  box(64, 66, 18, 5, "Genomic hotspot", "#E7F1F8", "#2E9BD6", 0.53)
  box(84, 66, 19, 5, "Non-hotspot window", "#E7F1F8", "#2E9BD6", 0.5)
  arrow(64, 73, 64, 68.5); arrow(84, 73, 84, 68.5)

  filters <- paste("QD <2.0; MQ <40.0; FS >40.0;",
                   "HaplotypeScore >13.0; MQRankSum <-12.5;",
                   "ReadPosRankSum <-4.0", sep = "\n")
  box(64, 54, 31, 8, filters, "#BCD3EA", "#91B2D1", 0.43)
  arrow(64, 63.5, 64, 58); arrow(84, 63.5, 84, 46)
  box(74, 43, 38, 6, "Merge filtered hotspot and non-hotspot SNPs", "#E7F1F8", "#2E9BD6", 0.51)
  arrow(64, 50, 70, 46); arrow(84, 46, 80, 46)
  box(74, 32, 34, 6, "Publication SNP set and five-track Circos inputs", "#D8E8F5", "#2E9BD6", 0.5)
  arrow(74, 40, 74, 35)

  graphics::text(2, 98, "a", font = 2, cex = 0.9, adj = c(0, 1))
  graphics::text(50, 5,
                 "Defaults reproduce the final analysis: 200-kb windows, 20-kb steps and a strict count >100 hotspot rule.",
                 cex = 0.45, col = "#444444")
}

main <- function() {
  args <- parse_cli(commandArgs(trailingOnly = TRUE),
                    list(width_mm = "183", height_mm = "110", dpi = "600",
                         formats = "pdf,svg,png,tiff", window_kb = "200",
                         step_kb = "20", hot_count_threshold = "100"))
  require_cli(args, c("output_prefix"))
  window_kb <- as.integer(args$window_kb); step_kb <- as.integer(args$step_kb)
  threshold <- as.integer(args$hot_count_threshold)
  if (any(!is.finite(c(window_kb, step_kb, threshold))) || any(c(window_kb, step_kb, threshold) <= 0L)) {
    stop("Window, step and hotspot threshold values must be positive integers")
  }
  export_figure(function() draw_workflow(window_kb, step_kb, threshold),
                args$output_prefix, args$width_mm, args$height_mm, args$dpi, args$formats)
}

if (sys.nframe() == 0L) main()
