#!/usr/bin/env Rscript

# Draw the genome-annotation workflow used for the manuscript supplement.

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])
helper_path <- file.path("analysis", "figures", "export_utils.R")
if (!file.exists(helper_path)) helper_path <- file.path(dirname(script_path), "..", "figures", "export_utils.R")
source(helper_path)

node <- function(x, y, width, height, label, fill = "#DCEAF6", border = "#2878B5", cex = 0.55) {
  graphics::rect(x - width / 2, y - height / 2, x + width / 2, y + height / 2,
                 col = fill, border = border, lwd = 0.8)
  graphics::text(x, y, label, cex = cex, family = "sans")
}

edge <- function(x0, y0, x1, y1, lty = 1, color = "#222222") {
  graphics::arrows(x0, y0, x1, y1, length = 0.055, angle = 22, lwd = 0.75,
                   lty = lty, col = color)
}

draw_annotation_workflow <- function() {
  graphics::par(mar = rep(0.35, 4), family = "sans", xpd = NA)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 100), ylim = c(0, 100), asp = 1)

  node(50, 94, 20, 5.5, "Genome assembly", "#B9D2E8", "#245B84", 0.72)

  node(17, 84, 16, 5, "Repeat annotation", "#D8EEF7", "#1B708F")
  node(17, 75, 15, 6, "De novo and\nhomology-based prediction", "#E5F4FA", "#1B708F", 0.40)
  node(17, 65, 14, 5, "Repeat-masked\ngenome", "#D8EEF7", "#1B708F")
  edge(45, 92, 17, 86.5); edge(17, 81.5, 17, 78); edge(17, 72, 17, 67.5)

  node(45, 84, 17, 5, "Structural annotation", "#B9D2E8", "#245B84")
  node(35, 74, 15, 6, "Homology-based\ngene prediction", "#DCEAF6", "#245B84", 0.48)
  node(50, 74, 14, 6, "De novo gene\nprediction", "#DCEAF6", "#245B84", 0.48)
  node(61, 74, 7, 6, "Other\nevidence", "#DCEAF6", "#245B84", 0.40)
  edge(48, 91.3, 45, 86.5); edge(45, 81.5, 35, 77); edge(45, 81.5, 50, 77); edge(45, 81.5, 61, 77)
  edge(24, 65, 43, 71.5); edge(35, 71, 47, 67); edge(61, 71, 53, 67)

  node(76, 84, 14, 5, "ncRNA annotation", "#DCEAF6", "#245B84")
  for (i in seq_along(c("tRNA", "rRNA", "miRNA", "snRNA"))) {
    x <- c(67, 74, 81, 88)[[i]]
    node(x, 74, 5.5, 4.5, c("tRNA", "rRNA", "miRNA", "snRNA")[[i]], "#E8F1F8", "#245B84", 0.43)
    edge(76, 81.5, x, 76.25)
  }
  edge(55, 92, 76, 86.5)

  node(83, 95, 12, 5, "PacBio Iso-Seq", "#E2EEF8", "#78A6CC", 0.5)
  node(95, 95, 9, 5, "RNA-seq", "#E2EEF8", "#78A6CC", 0.52)
  node(83, 88, 11, 4.5, "Transcripts", "#E2EEF8", "#78A6CC")
  node(95, 88, 8, 4.5, "Trinity", "#E2EEF8", "#78A6CC")
  node(91, 81, 12, 5, "PASA models", "#E2EEF8", "#78A6CC")
  node(91, 66, 12, 5, "Training set", "#E2EEF8", "#78A6CC")
  edge(83, 92.5, 83, 90.25); edge(95, 92.5, 95, 90.25)
  edge(83, 85.75, 88, 83.5); edge(95, 85.75, 94, 83.5); edge(91, 78.5, 91, 68.5)
  edge(85, 66, 55, 65, lty = 2)

  node(50, 64, 15, 5.5, "EvidenceModeler", "#B9D2E8", "#245B84", 0.62)
  node(50, 55, 10, 4.5, "PASA refinement", "#B9D2E8", "#245B84", 0.52)
  node(50, 47, 11, 4.5, "Quality control", "#B9D2E8", "#245B84", 0.52)
  node(50, 39, 11, 4.5, "Final gene set", "#B9D2E8", "#245B84", 0.54)
  edge(50, 61.25, 50, 57.25); edge(50, 52.75, 50, 49.25); edge(50, 44.75, 50, 41.25)

  node(69, 48, 14, 5, "Functional annotation", "#E3F2D8", "#4C8A2B", 0.55)
  edge(56, 39, 67, 45.5)
  databases <- c("InterPro", "KEGG", "Swiss-Prot", "TrEMBL")
  for (i in seq_along(databases)) {
    x <- c(57, 66, 75, 84)[[i]]
    node(x, 36, 8, 4.5, databases[[i]], "#EAF6E2", "#4C8A2B", 0.46)
    edge(69, 45.5, x, 38.25)
  }

  graphics::text(3, 98, "c", font = 2, cex = 0.9, adj = c(0, 1))
  graphics::text(50, 25, "Genome annotation workflow", font = 2, cex = 0.72)
  graphics::text(50, 20,
                 "Repeat, structural, non-coding RNA, transcript-assisted and functional annotation evidence are integrated.",
                 cex = 0.46, col = "#444444")
}

main <- function() {
  args <- parse_cli(commandArgs(trailingOnly = TRUE),
                    list(width_mm = "183", height_mm = "145", dpi = "600",
                         formats = "pdf,svg,png,tiff"))
  require_cli(args, c("output_prefix"))
  export_figure(draw_annotation_workflow, args$output_prefix, args$width_mm,
                args$height_mm, args$dpi, args$formats)
}

if (sys.nframe() == 0L) main()
