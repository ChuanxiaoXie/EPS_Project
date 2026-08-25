#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: go_enrichment.R GENE2GO_RDS_OR_TSV GENE_LIST GENE_LENGTHS OUTPUT_DIR PREFIX")
}

gene2go_file <- args[[1]]
gene_list_file <- args[[2]]
gene_length_file <- args[[3]]
output_dir <- args[[4]]
prefix <- args[[5]]

suppressPackageStartupMessages({
  library(goseq)
  library(topGO)
})

read_gene2go <- function(path) {
  if (grepl("[.]rds$", path, ignore.case = TRUE)) {
    return(readRDS(path))
  }
  mapping <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if (!all(c("gene", "go_id") %in% names(mapping))) {
    stop("GENE2GO_TSV must contain gene and go_id columns")
  }
  split(as.character(mapping$go_id), as.character(mapping$gene))
}

gene2GO <- read_gene2go(gene2go_file)
if (!is.list(gene2GO) || is.null(names(gene2GO))) {
  stop("GENE2GO input must define a named mapping from genes to GO identifiers")
}
selected_genes <- scan(gene_list_file, what = "character", quiet = TRUE)
length_table <- read.table(gene_length_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
if (ncol(length_table) < 2) {
  stop("GENE_LENGTHS must contain gene identifier and length columns")
}
gene_lengths <- setNames(as.numeric(length_table[[2]]), as.character(length_table[[1]]))
universe <- intersect(names(gene2GO), names(gene_lengths))
if (!length(universe)) {
  stop("No genes are shared by the GO mapping and length table")
}
gene_vector <- as.integer(universe %in% selected_genes)
names(gene_vector) <- universe

probability_weighting <- nullp(gene_vector, bias.data = gene_lengths[universe], plot.fit = FALSE)
result <- goseq(probability_weighting, gene2cat = gene2GO[universe], method = "Wallenius")
positive <- result$over_represented_pvalue[result$over_represented_pvalue > 0]
if (length(positive)) {
  result$over_represented_pvalue[result$over_represented_pvalue == 0] <- min(positive) / 10000
}
result$adjusted_pvalue <- p.adjust(result$over_represented_pvalue, method = "BH")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
table_path <- file.path(output_dir, paste0(prefix, ".go_enrichment.tsv"))
write.table(result, table_path, row.names = FALSE, quote = FALSE, sep = "\t")

ontology_map <- c(BP = "biological_process", CC = "cellular_component", MF = "molecular_function")
for (ontology in names(ontology_map)) {
  selected <- result[result$ontology == ontology_map[[ontology]] & result$adjusted_pvalue < 0.05, , drop = FALSE]
  if (!nrow(selected)) {
    next
  }
  topgo_data <- new(
    "topGOdata",
    description = ontology,
    ontology = ontology,
    allGenes = factor(gene_vector),
    annot = annFUN.gene2GO,
    gene2GO = gene2GO[universe],
    nodeSize = 10
  )
  graph_nodes <- nodes(graph(topgo_data))
  scores <- setNames(selected$adjusted_pvalue, selected$category)[graph_nodes]
  significant_nodes <- min(10, sum(!is.na(scores)))
  if (significant_nodes < 1) {
    next
  }
  pdf(file.path(output_dir, paste0(prefix, ".", ontology, ".dag.pdf")))
  showSigOfNodes(topgo_data, useInfo = "all", scores, firstSigNodes = significant_nodes)
  dev.off()
}
