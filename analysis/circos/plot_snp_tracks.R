#!/usr/bin/env Rscript

# Draw the manuscript five-track SNP Circos figure from prepared interval files.

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])
helper_path <- file.path("analysis", "figures", "export_utils.R")
if (!file.exists(helper_path)) helper_path <- file.path(dirname(script_path), "..", "figures", "export_utils.R")
source(helper_path)

if (!requireNamespace("circlize", quietly = TRUE)) {
  stop("The circlize package is required: install.packages('circlize')")
}

read_chromosome_sizes <- function(path) {
  data <- utils::read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  if (ncol(data) < 2L) stop("Chromosome sizes require two columns")
  data <- data[, 1:2]
  names(data) <- c("chromosome", "length")
  data$length <- as.numeric(data$length)
  if (anyDuplicated(data$chromosome) || any(!is.finite(data$length)) || any(data$length <= 0)) {
    stop("Chromosome names must be unique and lengths must be positive")
  }
  data
}

validate_track <- function(track, chromosome_sizes, path) {
  unknown <- setdiff(unique(track$sector), chromosome_sizes$chromosome)
  if (length(unknown) > 0L) stop(paste("Unknown chromosome in", path, ":", paste(unknown, collapse = ", ")))
  limits <- chromosome_sizes$length[match(track$sector, chromosome_sizes$chromosome)]
  if (any(track$start < 0) || any(track$end > limits)) stop(paste("Out-of-range interval in", path))
  track
}

plot_point_track <- function(track, color, point_cex = 0.22) {
  track$value <- normalize_values(track$value)
  circlize::circos.genomicTrackPlotRegion(
    track,
    ylim = c(0, 1),
    track.height = 0.105,
    bg.border = NA,
    panel.fun = function(region, value, ...) {
      circlize::circos.genomicPoints(region, value, pch = 16, cex = point_cex, col = color)
    }
  )
}

plot_bar_track <- function(track, color, height = 0.09) {
  track$value <- normalize_values(track$value)
  circlize::circos.genomicTrackPlotRegion(
    track,
    ylim = c(0, 1),
    track.height = height,
    bg.border = NA,
    panel.fun = function(region, value, ...) {
      circlize::circos.genomicRect(
        region,
        ybottom = rep(0, nrow(region)),
        ytop = value[, 1L],
        col = color,
        border = NA
      )
    }
  )
}

draw_snp_circos <- function(chromosome_sizes, tracks) {
  circlize::circos.clear()
  chromosome_count <- nrow(chromosome_sizes)
  gaps <- rep(1.8, chromosome_count)
  gaps[[chromosome_count]] <- 13
  circlize::circos.par(
    start.degree = 82,
    gap.after = gaps,
    track.margin = c(0.003, 0.003),
    cell.padding = c(0, 0, 0, 0),
    canvas.xlim = c(-1.28, 1.28),
    canvas.ylim = c(-1.28, 1.28),
    points.overflow.warning = FALSE
  )
  circlize::circos.initialize(
    factors = chromosome_sizes$chromosome,
    xlim = cbind(rep(0, chromosome_count), chromosome_sizes$length)
  )
  circlize::circos.trackPlotRegion(
    ylim = c(0, 1),
    track.height = 0.075,
    bg.col = "#A66F00",
    bg.border = "#A66F00",
    panel.fun = function(x, y) {
      sector <- circlize::get.cell.meta.data("sector.index")
      xlim <- circlize::get.cell.meta.data("xlim")
      sector_length <- diff(xlim)
      major_step <- if (sector_length >= 100000000) 50000000 else max(100000, round(sector_length / 2 / 100000) * 100000)
      major_at <- seq(0, sector_length, by = major_step)
      circlize::circos.axis(
        h = "top",
        major.at = major_at,
        labels = format(round(major_at / 1000000, 1), trim = TRUE),
        labels.cex = 0.42,
        major.tick.length = 0.035,
        minor.ticks = 0,
        labels.facing = "clockwise",
        labels.niceFacing = TRUE
      )
      circlize::circos.text(
        mean(xlim),
        1.85,
        sub("^Chr", "Chr ", sector),
        facing = "bending.outside",
        niceFacing = TRUE,
        cex = 0.62
      )
    }
  )
  plot_point_track(tracks$genome_snp, "#148DB7", 0.23)
  plot_point_track(tracks$genic_snp, "#42B7C9", 0.21)
  plot_bar_track(tracks$gene_density, grDevices::adjustcolor("#98D48D", alpha.f = 0.55), 0.09)
  plot_bar_track(tracks$gc_density, "#9DD8CF", 0.095)
  plot_bar_track(tracks$gc_to_at, "#61B5E0", 0.095)

  graphics::text(0, 0.70, "a", cex = 0.75)
  graphics::text(0, 0.59, "b", cex = 0.75)
  graphics::text(0, 0.48, "c", cex = 0.75)
  graphics::text(0, 0.37, "d", cex = 0.75)
  graphics::text(0, 0.26, "e", cex = 0.75)
  legend_labels <- c(
    "a  Genome SNP density",
    "b  Genic SNP density",
    "c  Gene density",
    "d  Genome GC density",
    "e  GC>AT SNP density"
  )
  graphics::text(-0.43, 0.02, paste(legend_labels, collapse = "\n"), adj = c(0, 0.5), cex = 0.55)
  circlize::circos.clear()
}

main <- function() {
  args <- parse_cli(
    commandArgs(trailingOnly = TRUE),
    list(width_mm = "140", height_mm = "140", dpi = "600", formats = "pdf,svg,png,tiff")
  )
  require_cli(args, c("track_dir", "chrom_sizes", "output_prefix"))
  chromosome_sizes <- read_chromosome_sizes(args$chrom_sizes)
  track_paths <- list(
    genome_snp = file.path(args$track_dir, "1_wai.txt"),
    genic_snp = file.path(args$track_dir, "1_nei.txt"),
    gene_density = file.path(args$track_dir, "2_genedensity.txt"),
    gc_density = file.path(args$track_dir, "3_GC_wai.txt"),
    gc_to_at = file.path(args$track_dir, "3_GC2AT_nei.txt")
  )
  tracks <- Map(
    function(path) validate_track(read_interval_track(path), chromosome_sizes, path),
    track_paths
  )
  export_figure(
    function() draw_snp_circos(chromosome_sizes, tracks),
    args$output_prefix,
    args$width_mm,
    args$height_mm,
    args$dpi,
    args$formats
  )
}

if (sys.nframe() == 0L) main()
