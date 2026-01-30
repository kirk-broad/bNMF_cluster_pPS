
#!/usr/bin/env Rscript
# pPS_summary_plots.R
# Final version: Includes robust column mapping + Summary Statistics output

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggplot2)
  library(gridExtra)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript pPS_summary_plots.R <score_file> <weights_file> <output_dir>")
}

prs_scores_file <- args[1]
weights_file <- args[2]
output_dir <- args[3]

cat("Reading PRS scores from:", prs_scores_file, "\n")
cat("Reading weights file from:", weights_file, "\n")
cat("Saving outputs to:", output_dir, "\n")

# 1. Load Data
prs_scores <- fread(prs_scores_file, header = TRUE)
cat("PRS scores dimensions:", dim(prs_scores), "\n")

# --- FIX 1: Handle PLINK's #IID column name ---
names(prs_scores) <- gsub("^#", "", names(prs_scores))

# 2. Identify Cluster Columns
weights <- fread(weights_file, header = TRUE)
cat("Weights file dimensions:", dim(weights), "\n")

# Extract cluster names (assuming columns 3 onwards are clusters in the weights file)
weight_cols <- colnames(weights)
cluster_cols <- weight_cols[sapply(weights, is.numeric)]
cluster_names <- setdiff(cluster_cols, c("CHR", "POS", "pos", "chr"))

cat("Identified cluster names:", paste(cluster_names, collapse=", "), "\n")

# Map .sscore columns (SCORE1_SUM, etc.) to Cluster Names
score_cols <- grep("_SUM", colnames(prs_scores), value = TRUE)
cat("Score columns found:", paste(score_cols, collapse=", "), "\n")

if(length(score_cols) != length(cluster_names)) {
  warning("Mismatch between number of score columns in .sscore file and clusters in weights file!")
  min_len <- min(length(score_cols), length(cluster_names))
  score_cols <- score_cols[1:min_len]
  cluster_names <- cluster_names[1:min_len]
}

# --- FIX 2: Correct Rename Mapping (NewName = OldName) ---
name_map <- setNames(score_cols, cluster_names) 

prs_plot_data <- prs_scores %>%
  rename(any_of(name_map)) %>%
  dplyr::select(IID, all_of(cluster_names))

# 3. Generate Summary Statistics (NEW SECTION)
cat("Calculating summary statistics...\n")

summary_stats <- prs_plot_data %>%
  pivot_longer(cols = all_of(cluster_names), names_to = "Cluster", values_to = "Score") %>%
  group_by(Cluster) %>%
  summarise(
    N = n(),
    Min = min(Score, na.rm = TRUE),
    Q1 = quantile(Score, 0.25, na.rm = TRUE),
    Median = median(Score, na.rm = TRUE),
    Mean = mean(Score, na.rm = TRUE),
    Q3 = quantile(Score, 0.75, na.rm = TRUE),
    Max = max(Score, na.rm = TRUE),
    SD = sd(Score, na.rm = TRUE)
  ) %>%
  mutate(across(where(is.numeric), ~ round(., 4))) # Round to 4 decimal places

# Save Summary Stats Table
stats_file <- file.path(output_dir, "pPS_summary_stats.txt")
write_tsv(summary_stats, stats_file)
cat("Saved summary statistics to", stats_file, "\n")


# 4. Generate Plots

# A. Histograms
cat("Creating histogram plots...\n")
plot_list <- list()

for (cluster in cluster_names) {
  p <- ggplot(prs_plot_data, aes(x = .data[[cluster]])) +
    geom_histogram(bins = 50, fill = "steelblue", color = "white") +
    theme_minimal() +
    labs(title = paste("Distribution of", cluster), x = "Polygenic Score", y = "Count")
  
  plot_list[[cluster]] <- p
}

# Save Histograms
pdf_file <- file.path(output_dir, "pPS_histograms.pdf")
pdf(pdf_file, width = 12, height = 10)
do.call(grid.arrange, c(plot_list, ncol = 3))
dev.off()
cat("Saved histograms to", pdf_file, "\n")

# B. Correlation Matrix
cat("Creating correlation matrix...\n")
if(length(cluster_names) > 1) {
  cor_mat <- cor(prs_plot_data %>% dplyr::select(all_of(cluster_names)), use = "complete.obs")
  
  # Melt for ggplot
  melted_cormat <- reshape2::melt(cor_mat)
  
  p_cor <- ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), name="Pearson\nCorrelation") +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1)) +
    labs(title = "Correlation between Partitioned Scores", x="", y="")
  
  cor_file <- file.path(output_dir, "pPS_correlation.pdf")
  ggsave(cor_file, plot = p_cor, width = 8, height = 7)
  cat("Saved correlation plot to", cor_file, "\n")
}

cat("Summary generation complete.\n")