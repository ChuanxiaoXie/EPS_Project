# Load required packages
library(readxl)
library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)

# Set input file paths
chrom_size_path <- ""
snp_pos_path <- ""

# Set output PDF path
output_pdf_path <- ""

cat(output_pdf_path, "\n")

# Read chromosome size file
chrom_size <- read_excel(
  chrom_size_path,
  sheet = "Sheet2",
  skip = 1,
  col_names = c("CHROM", "size")
)

# Retrieve all worksheet names from the SNP file
sheet_names <- excel_sheets(snp_pos_path)
cat(paste(sheet_names, collapse = ", "), "\n")

# Define the desired chromosome order
chrom_order <- c(
  "Chr1", "Chr2", "Chr3", "Chr4", "Chr5",
  "Chr6", "Chr7", "Chr8", "Chr9", "Chr10"
)

# Sort chromosome size data and calculate cumulative genomic coordinates
chrom_size <- chrom_size %>%
  mutate(CHROM = factor(CHROM, levels = chrom_order)) %>%
  arrange(CHROM) %>%
  mutate(
    cum_size = cumsum(size),
    start_pos = lag(cum_size, default = 0)
  )

# Define a function to generate an SNP density plot for an individual worksheet
create_density_plot <- function(sheet_name, snp_data) {
  
  # Define sliding-window parameters
  window_size <- 3000000
  step_size <- 200000
  
  results <- data.frame()
  
  # Process SNP density separately for each chromosome
  for (i in 1:nrow(chrom_size)) {
    chrom <- as.character(chrom_size$CHROM[i])
    chr_size <- chrom_size$size[i]
    
    # Generate sliding-window start positions
    window_starts <- seq(0, chr_size - window_size, by = step_size)
    
    # Generate at least one window for chromosomes shorter than the window size
    if (length(window_starts) == 0) {
      window_starts <- 0
    }
    
    # Extract SNP positions for the current chromosome
    snp_chr <- snp_data[snp_data$CHROM == chrom, "POS", drop = TRUE]
    
    # Count SNPs within each sliding window
    for (start_pos in window_starts) {
      end_pos <- start_pos + window_size
      snp_count <- sum(snp_chr >= start_pos & snp_chr < end_pos)
      
      results <- rbind(
        results,
        data.frame(
          CHROM = chrom,
          Window_Start = start_pos,
          Window_End = end_pos,
          SNP_Count = snp_count
        )
      )
    }
  }
  
  # Map chromosome-specific window positions to cumulative genomic coordinates
  plot_data <- results %>%
    left_join(select(chrom_size, CHROM, start_pos), by = "CHROM") %>%
    mutate(
      global_start = start_pos + Window_Start,
      global_end = start_pos + Window_End
    )
  
  # Define SNP count intervals for heatmap color classification
  plot_data$Color_Group <- cut(
    plot_data$SNP_Count,
    breaks = c(-Inf, 0, 10, 20, 30, 40, 50, 60, Inf),
    labels = c(
      "0", "1-10", "11-20", "21-30",
      "31-40", "41-50", "51-60", "60+"
    ),
    right = TRUE
  )
  
  # Define the heatmap color palette
  color_values <- c(
    "0" = "#FFFFFF",
    "1-10" = "#E1F0F7",
    "11-20" = "#A6D2E6",
    "21-30" = "#6BB4D6",
    "31-40" = "#60A2C1",
    "41-50" = "#5690AB",
    "51-60" = "#4B7E96",
    "60+" = "#406C80"
  )
  
  # Calculate chromosome label positions
  chrom_label_pos <- chrom_size %>%
    mutate(mid = start_pos + size / 2)
  
  # Generate chromosome boundary line coordinates
  chrom_borders <- chrom_size %>%
    mutate(end_pos = start_pos + size) %>%
    slice(rep(1:n(), each = 4)) %>%
    mutate(
      type = rep(
        c("top", "bottom", "left", "right"),
        times = nrow(chrom_size)
      ),
      x = case_when(
        type %in% c("top", "bottom") ~ start_pos,
        type == "left" ~ start_pos,
        type == "right" ~ end_pos
      ),
      y = case_when(
        type == "top" ~ 1,
        type == "bottom" ~ 0,
        type %in% c("left", "right") ~ 0
      ),
      xend = case_when(
        type %in% c("top", "bottom") ~ end_pos,
        type == "left" ~ start_pos,
        type == "right" ~ end_pos
      ),
      yend = case_when(
        type == "top" ~ 1,
        type == "bottom" ~ 0,
        type %in% c("left", "right") ~ 1
      )
    )
  
  # Generate SNP density heatmap
  p <- ggplot(plot_data) +
    geom_segment(
      data = chrom_borders,
      aes(x = x, y = y, xend = xend, yend = yend),
      color = "black",
      linewidth = 0.8
    ) +
    geom_rect(
      aes(
        xmin = global_start,
        xmax = global_end,
        ymin = 0.05,
        ymax = 0.95,
        fill = Color_Group
      ),
      color = NA
    ) +
    scale_fill_manual(values = color_values) +
    scale_x_continuous(
      breaks = chrom_label_pos$mid,
      labels = chrom_label_pos$CHROM,
      expand = c(0, 0)
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      x = "Chromosome",
      y = NULL,
      title = paste0(sheet_name)
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5,
        size = 10,
        face = "bold"
      ),
      axis.text.x = element_text(size = 8),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  return(p)
}

# Define a function to create the SNP count legend
create_legend_plot <- function() {
  
  # Create legend data
  legend_data <- data.frame(
    Group = factor(
      c(
        "0", "1-10", "11-20", "21-30",
        "31-40", "41-50", "51-60", "60+"
      ),
      levels = c(
        "0", "1-10", "11-20", "21-30",
        "31-40", "41-50", "51-60", "60+"
      )
    ),
    Value = 1
  )
  
  # Define the legend color palette
  color_values <- c(
    "0" = "#FFFFFF",
    "1-10" = "#E1F0F7",
    "11-20" = "#A6D2E6",
    "21-30" = "#6BB4D6",
    "31-40" = "#60A2C1",
    "41-50" = "#5690AB",
    "51-60" = "#4B7E96",
    "60+" = "#406C80"
  )
  
  # Generate the legend plot
  legend_plot <- ggplot(
    legend_data,
    aes(x = Group, y = Value, fill = Group)
  ) +
    geom_tile(
      color = "black",
      linewidth = 0.2,
      width = 0.5,
      height = 0.4
    ) +
    geom_text(
      aes(label = Group),
      vjust = 3,
      size = 1.8,
      fontface = "bold"
    ) +
    scale_fill_manual(values = color_values) +
    scale_y_continuous(limits = c(0, 1.5)) +
    labs(title = "SNP Count per 3Mb Window\n(200Kb Sliding Step)") +
    theme_void() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 8,
        face = "bold",
        margin = margin(b = 10)
      ),
      legend.position = "none",
      plot.margin = margin(10, 10, 10, 10)
    )
  
  return(legend_plot)
}

# Initialize a list to store valid worksheet plots
plot_list <- list()

# Process each worksheet in the SNP file
for (sheet_name in sheet_names) {
  
  # Exclude non-data worksheets
  if (sheet_name %in% c("Sheet1", "Summary", "Metadata")) next
  
  cat(sheet_name, "\n")
  
  tryCatch({
    
    # Read SNP position data from the current worksheet
    snp_data <- read_excel(
      snp_pos_path,
      sheet = sheet_name,
      skip = 1,
      col_names = c("CHROM", "POS", "Genotype")
    )
    
    # Skip empty worksheets
    if (nrow(snp_data) == 0) {
      cat(sheet_name, "\n")
      next
    }
    
    # Generate and store the SNP density plot
    density_plot <- create_density_plot(sheet_name, snp_data)
    plot_list[[sheet_name]] <- density_plot
    
    cat(sheet_name, "\n")
    
  }, error = function(e) {
    cat(sheet_name, e$message, "\n")
  })
}

# Determine the total number of valid sample plots
n_plots <- length(plot_list)

if (n_plots == 0) {
  stop("ERROR")
}

# Generate the legend plot
legend_plot <- create_legend_plot()

# Append the legend plot to the plot list
plot_list[["Legend"]] <- legend_plot
n_plots <- n_plots + 1

# Determine the multi-panel layout, with the legend placed in the final row
n_cols <- ifelse(max(n_plots - 1, 1) <= 4, 2, 3)
n_rows <- ceiling((n_plots - 1) / n_cols) + 1

# Combine all sample plots and the legend into a multi-panel figure
combined_plot <- wrap_plots(
  plot_list,
  ncol = n_cols,
  nrow = n_rows,
  heights = c(rep(1, n_rows - 1), 0.3)
) +
  plot_annotation(
    title = "SNP Distribution Density Across All Samples\n(3Mb Window with 200Kb Sliding Step)",
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 16,
        face = "bold"
      )
    )
  )

# Export the combined figure to a PDF file
pdf(output_pdf_path, width = 16, height = 4 + n_rows * 1)
print(combined_plot)
dev.off()

cat(output_pdf_path, "\n")