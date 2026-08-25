#!/usr/bin/env Rscript

# Draw paired CA7301/B73 chromosome sectors with TE, gene, LAI and synteny tracks.

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])
helper_path <- file.path("analysis", "figures", "export_utils.R")
if (!file.exists(helper_path)) helper_path <- file.path(dirname(script_path), "..", "figures", "export_utils.R")
source(helper_path)

if (!requireNamespace("circlize", quietly = TRUE)) {
  stop("The circlize package is required: install.packages('circlize')")
}

read_sectors <- function(path) {
  data <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("sector", "label", "genome", "chromosome", "length", "color")
  if (!all(required %in% names(data))) stop(paste("Sector table must contain", paste(required, collapse = ", ")))
  data$length <- as.numeric(data$length)
  if (anyDuplicated(data$sector) || any(data$length <= 0)) stop("Sector identifiers must be unique and lengths positive")
  data
}

read_links <- function(path, sector_names) {
  data <- utils::read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  if (ncol(data) < 6L) stop("Link table requires at least six columns")
  if (ncol(data) == 6L) data[[7L]] <- 1
  data <- data[, 1:7, drop = FALSE]
  names(data) <- c("sector1", "start1", "end1", "sector2", "start2", "end2", "score")
  if (any(!c(data$sector1, data$sector2) %in% sector_names)) stop("Link table contains unknown sectors")
  data
}

track_colors <- function(track, sectors, low_alpha = 0.18, high_alpha = 0.90) {
  normalized <- normalize_values(track$value)
  base_colors <- sectors$color[match(track$sector, sectors$sector)]
  mapply(
    function(color, alpha) grDevices::adjustcolor(color, alpha.f = alpha),
    base_colors,
    low_alpha + normalized * (high_alpha - low_alpha),
    USE.NAMES = FALSE
  )
}

draw_heat_track <- function(track, sectors, height) {
  colors <- track_colors(track, sectors)
  track$normalized <- normalize_values(track$value)
  circlize::circos.genomicTrackPlotRegion(
    track[, c("sector", "start", "end", "normalized")],
    ylim = c(0, 1),
    track.height = height,
    bg.border = NA,
    panel.fun = function(region, value, ...) {
      sector <- circlize::get.cell.meta.data("sector.index")
      indices <- which(track$sector == sector)
      circlize::circos.genomicRect(region, 0, 1, col = colors[indices], border = NA)
    }
  )
}

draw_line_track <- function(track, color, height) {
  track$normalized <- normalize_values(track$value)
  circlize::circos.genomicTrackPlotRegion(
    track[, c("sector", "start", "end", "normalized")],
    ylim = c(0, 1),
    track.height = height,
    bg.border = NA,
    panel.fun = function(region, value, ...) {
      centers <- rowMeans(region)
      circlize::circos.lines(centers, value[, 1L], col = color, lwd = 0.7, type = "l")
    }
  )
}

draw_comparative_circos <- function(sectors, te, gene, lai, links) {
  circlize::circos.clear()
  sector_count <- nrow(sectors)
  gaps <- rep(1.0, sector_count)
  chromosome_ends <- seq(2L, sector_count, by = 2L)
  gaps[chromosome_ends] <- 3.0
  gaps[[sector_count]] <- 13
  circlize::circos.par(
    start.degree = 83,
    gap.after = gaps,
    track.margin = c(0.003, 0.003),
    cell.padding = c(0, 0, 0, 0),
    canvas.xlim = c(-1.28, 1.28),
    canvas.ylim = c(-1.28, 1.28),
    points.overflow.warning = FALSE
  )
  circlize::circos.initialize(sectors$sector, xlim = cbind(0, sectors$length))
  circlize::circos.trackPlotRegion(
    ylim = c(0, 1),
    track.height = 0.075,
    bg.col = NA,
    bg.border = NA,
    panel.fun = function(x, y) {
      sector <- circlize::get.cell.meta.data("sector.index")
      index <- match(sector, sectors$sector)
      xlim <- circlize::get.cell.meta.data("xlim")
      circlize::circos.rect(xlim[[1L]], 0, xlim[[2L]], 1,
                            col = sectors$color[[index]], border = sectors$color[[index]])
      circlize::circos.text(
        mean(xlim),
        1.78,
        sectors$label[[index]],
        facing = "bending.outside",
        niceFacing = TRUE,
        cex = 0.37,
        col = sectors$color[[index]]
      )
    }
  )
  draw_heat_track(te, sectors, 0.12)
  draw_line_track(gene, "#2B6EA6", 0.095)
  draw_line_track(lai, "#70AD47", 0.095)
  for (row_index in seq_len(nrow(links))) {
    link <- links[row_index, ]
    circlize::circos.link(
      link$sector1,
      c(link$start1, link$end1),
      link$sector2,
      c(link$start2, link$end2),
      col = grDevices::adjustcolor("#70AD47", alpha.f = 0.32),
      border = NA
    )
  }
  graphics::text(0, 0.26, "a  TE density\nb  Gene density\nc  LAI\nd  Gene collinearity", cex = 0.48)
  circlize::circos.clear()
}

main <- function() {
  args <- parse_cli(
    commandArgs(trailingOnly = TRUE),
    list(width_mm = "140", height_mm = "140", dpi = "600", formats = "pdf,svg,png,tiff")
  )
  require_cli(args, c("sectors", "te", "gene", "lai", "links", "output_prefix"))
  sectors <- read_sectors(args$sectors)
  tracks <- lapply(c(args$te, args$gene, args$lai), read_interval_track)
  known <- sectors$sector
  if (any(!unlist(lapply(tracks, function(track) unique(track$sector))) %in% known)) stop("A comparative track contains an unknown sector")
  links <- read_links(args$links, known)
  export_figure(
    function() draw_comparative_circos(sectors, tracks[[1L]], tracks[[2L]], tracks[[3L]], links),
    args$output_prefix,
    args$width_mm,
    args$height_mm,
    args$dpi,
    args$formats
  )
}

if (sys.nframe() == 0L) main()
