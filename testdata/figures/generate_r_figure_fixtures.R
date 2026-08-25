#!/usr/bin/env Rscript

# Generate deterministic ten-chromosome inputs for executable R figure tests.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: generate_r_figure_fixtures.R OUTPUT_DIR")
output_dir <- args[[1L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

chromosomes <- paste0("Chr", seq_len(10L))
chromosome_lengths <- seq(1200000, 2100000, length.out = 10L)
chromosome_lengths <- as.integer(round(chromosome_lengths / 100000) * 100000)

snp_dir <- file.path(output_dir, "snp_circos")
dir.create(snp_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.table(
  data.frame(chromosome = chromosomes, length = chromosome_lengths),
  file.path(snp_dir, "chrom.sizes.tsv"),
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

make_track <- function(scale, phase, minimum = 0) {
  rows <- list()
  row_index <- 1L
  for (chromosome_index in seq_along(chromosomes)) {
    chromosome <- chromosomes[[chromosome_index]]
    chromosome_length <- chromosome_lengths[[chromosome_index]]
    starts <- seq(0, chromosome_length - 1L, length.out = 24L)
    starts <- as.integer(floor(starts / 1000) * 1000)
    ends <- c(starts[-1L], chromosome_length)
    values <- minimum + scale * (
      0.15 + 0.85 * abs(sin(seq_along(starts) / 3 + chromosome_index / 2 + phase))
    )
    rows[[row_index]] <- data.frame(chromosome, starts, ends, round(values, 3))
    row_index <- row_index + 1L
  }
  do.call(rbind, rows)
}

tracks <- list(
  "1_wai.txt" = make_track(65, 0.0),
  "1_nei.txt" = make_track(25, 0.7),
  "2_genedensity.txt" = make_track(15, 1.2),
  "3_GC_wai.txt" = make_track(70000, 1.8, 45000),
  "3_GC2AT_nei.txt" = make_track(49, 2.4)
)
for (name in names(tracks)) {
  utils::write.table(
    tracks[[name]],
    file.path(snp_dir, name),
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
}

comparative_dir <- file.path(output_dir, "comparative_circos")
dir.create(comparative_dir, recursive = TRUE, showWarnings = FALSE)
sector_rows <- list()
track_rows <- list(te = list(), gene = list(), lai = list())
links <- list()
sector_index <- 1L
link_index <- 1L
for (chromosome_index in seq_along(chromosomes)) {
  for (genome in c("CA7301", "B73")) {
    chromosome <- chromosomes[[chromosome_index]]
    length_multiplier <- if (genome == "CA7301") 1 else 0.94
    sector_length <- as.integer(chromosome_lengths[[chromosome_index]] * length_multiplier)
    sector <- paste(genome, chromosome, sep = "_")
    color <- if (genome == "CA7301") "#2B6EA6" else "#79B75B"
    sector_rows[[sector_index]] <- data.frame(
      sector = sector,
      label = paste(chromosome, genome),
      genome = genome,
      chromosome = chromosome,
      length = sector_length,
      color = color
    )
    starts <- as.integer(seq(0, sector_length - 1L, length.out = 18L))
    ends <- c(starts[-1L], sector_length)
    phase <- chromosome_index + if (genome == "CA7301") 0 else 0.6
    track_rows$te[[sector_index]] <- data.frame(sector, starts, ends, round(15 + 75 * abs(sin(seq_along(starts) / 4 + phase)), 3))
    track_rows$gene[[sector_index]] <- data.frame(sector, starts, ends, round(5 + 35 * abs(cos(seq_along(starts) / 3 + phase)), 3))
    track_rows$lai[[sector_index]] <- data.frame(sector, starts, ends, round(8 + 18 * abs(sin(seq_along(starts) / 5 + phase)), 3))
    sector_index <- sector_index + 1L
  }
  ca_sector <- paste("CA7301", chromosomes[[chromosome_index]], sep = "_")
  b73_sector <- paste("B73", chromosomes[[chromosome_index]], sep = "_")
  for (block in seq_len(3L)) {
    start1 <- as.integer((block - 0.75) * chromosome_lengths[[chromosome_index]] / 3)
    end1 <- start1 + as.integer(chromosome_lengths[[chromosome_index]] * 0.08)
    start2 <- as.integer(start1 * 0.94)
    end2 <- as.integer(end1 * 0.94)
    links[[link_index]] <- data.frame(ca_sector, start1, end1, b73_sector, start2, end2, 1)
    link_index <- link_index + 1L
  }
}
utils::write.table(do.call(rbind, sector_rows), file.path(comparative_dir, "sectors.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
for (track_name in names(track_rows)) {
  utils::write.table(do.call(rbind, track_rows[[track_name]]), file.path(comparative_dir, paste0(track_name, ".tsv")), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}
utils::write.table(do.call(rbind, links), file.path(comparative_dir, "links.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

hic_dir <- file.path(output_dir, "hic")
dir.create(hic_dir, recursive = TRUE, showWarnings = FALSE)
bins_per_chromosome <- 4L
bin_count <- length(chromosomes) * bins_per_chromosome
contacts <- expand.grid(X = 0:(bin_count - 1L), Y = 0:(bin_count - 1L))
same_chromosome <- (contacts$X %/% bins_per_chromosome) == (contacts$Y %/% bins_per_chromosome)
distance <- abs(contacts$X - contacts$Y)
log_contact <- pmin(4, 0.25 + 2.5 * exp(-distance / 2) + 1.2 * same_chromosome)
contacts$Z <- round(10^log_contact - 1)
utils::write.table(contacts, file.path(hic_dir, "contact_matrix.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
writeLines(as.character(seq(0, bin_count, by = bins_per_chromosome)), file.path(hic_dir, "chromosome_breaks.txt"))
writeLines(chromosomes, file.path(hic_dir, "chromosome_labels.txt"))

assembly_dir <- file.path(output_dir, "assembly")
dir.create(assembly_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.table(
  data.frame(
    sample = "SyntheticAssembly",
    single_copy = 80,
    duplicated = 10,
    fragmented = 5,
    missing = 5,
    mode = "genome",
    lineage = "synthetic_lineage"
  ),
  file.path(assembly_dir, "busco_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message(output_dir)
