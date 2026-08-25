# Shared command-line and publication-export helpers for R figures.

parse_cli <- function(values, defaults = list()) {
  result <- defaults
  index <- 1L
  while (index <= length(values)) {
    key <- values[[index]]
    if (!startsWith(key, "--")) stop(paste("Unexpected positional argument:", key))
    if (index == length(values)) stop(paste("Missing value for", key))
    key <- gsub("-", "_", substring(key, 3L), fixed = TRUE)
    result[[key]] <- values[[index + 1L]]
    index <- index + 2L
  }
  result
}

require_cli <- function(args, names) {
  missing <- names[!names %in% names(args) | vapply(args[names], function(x) !nzchar(x), logical(1))]
  if (length(missing) > 0L) stop(paste("Missing arguments:", paste(missing, collapse = ", ")))
}

script_directory <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) != 1L) stop("Unable to determine the current script directory")
  dirname(sub("^--file=", "", file_arg[[1L]]))
}

open_figure_device <- function(path, width_in, height_in, dpi) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  extension <- tolower(tools::file_ext(path))
  if (extension == "pdf") {
    if (capabilities("cairo")) {
      grDevices::cairo_pdf(path, width = width_in, height = height_in, family = "sans")
    } else {
      grDevices::pdf(path, width = width_in, height = height_in, family = "sans")
    }
  } else if (extension == "svg") {
    if (!requireNamespace("svglite", quietly = TRUE)) {
      stop("Editable SVG export requires the svglite package: install.packages('svglite')")
    }
    svglite::svglite(path, width = width_in, height = height_in, bg = "white")
  } else if (extension == "png") {
    grDevices::png(path, width = width_in, height = height_in, units = "in", res = dpi)
  } else if (extension %in% c("tif", "tiff")) {
    grDevices::tiff(
      path,
      width = width_in,
      height = height_in,
      units = "in",
      res = dpi,
      compression = "lzw"
    )
  } else {
    stop(paste("Unsupported figure extension:", extension))
  }
}

export_figure <- function(draw_function, output_prefix, width_mm, height_mm, dpi, formats) {
  width_in <- as.numeric(width_mm) / 25.4
  height_in <- as.numeric(height_mm) / 25.4
  dpi <- as.numeric(dpi)
  requested <- trimws(strsplit(formats, ",", fixed = TRUE)[[1L]])
  allowed <- c("pdf", "svg", "png", "tiff")
  if (length(requested) == 0L || any(!requested %in% allowed)) {
    stop(paste("formats must be a comma-separated subset of", paste(allowed, collapse = ", ")))
  }
  for (format in requested) {
    extension <- if (format == "tiff") "tiff" else format
    path <- paste0(output_prefix, ".", extension)
    open_figure_device(path, width_in, height_in, dpi)
    tryCatch(
      draw_function(),
      finally = grDevices::dev.off()
    )
    if (!file.exists(path) || file.info(path)$size <= 0L) stop(paste("Empty figure output:", path))
    message(path)
  }
}

read_interval_track <- function(path) {
  data <- utils::read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  if (ncol(data) < 4L) stop(paste("Interval track requires four columns:", path))
  data <- data[, 1:4]
  names(data) <- c("sector", "start", "end", "value")
  data$start <- as.numeric(data$start)
  data$end <- as.numeric(data$end)
  data$value <- as.numeric(data$value)
  if (any(!is.finite(data$value)) || any(data$end <= data$start)) {
    stop(paste("Invalid interval values in", path))
  }
  data
}

normalize_values <- function(values) {
  observed <- range(values, finite = TRUE)
  if (diff(observed) == 0) return(rep(0.5, length(values)))
  (values - observed[[1L]]) / diff(observed)
}
