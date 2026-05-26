#!/usr/bin/env Rscript
# cluster_weights_setup_generalized.R
# Usage: Rscript cluster_weights_setup_generalized.R <weights_dir> <cutoff> <chain_file> <genome_build> <output_dir>

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
})

# 1. Parse Arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) stop("Usage: script.R <weights_dir> <input_weights_file> <cutoff> <chain_file> <genome_build> <output_dir>")

weights_dir <- args[1]
input_weights_file <- args[2]  # <-- NEW
cutoff <- as.numeric(args[3])
chain_file <- args[4]
genome_build <- args[5]
output_dir <- args[6]

cat("------------------------------------------------\n")
cat("STEP 1: Setting up Cluster Weights (Single File Mode)\n")
cat("------------------------------------------------\n")

# 3. Read The Single Consolidated Weights File
# 3. Read The Single Consolidated Weights File
weights_file <- file.path(weights_dir, input_weights_file)
if (!file.exists(weights_file)) stop("Weights file not found: ", weights_file)

cat("Reading weights from:", weights_file, "\n")
df <- fread(weights_file)

# Handle BETA_aligned if it exists (some weight files use this column name)
if ("BETA_aligned" %in% colnames(df) && !"BETA" %in% colnames(df)) {
  colnames(df)[colnames(df) == "BETA_aligned"] <- "BETA"
}

# Basic Validation
req_cols <- c("CHR", "POS", "REF", "ALT", "Risk_Allele", "BETA")
missing_cols <- setdiff(req_cols, colnames(df))
if(length(missing_cols) > 0) stop("Missing required columns: ", paste(missing_cols, collapse=", "))

# 2. Dynamically Generate Cluster Names (Replacing the need for an external txt file)
standard_cols <- c("CHR", "POS", "REF", "ALT", "Risk_Allele", "BETA", "BETA_aligned", "Total_GRS", "VAR_ID", "SNP")
cluster_names <- setdiff(colnames(df), standard_cols)

cat("Detected clusters:", paste(cluster_names, collapse=", "), "\n")

# Save it to the output directory so downstream scripts (like Step 5) can read it
write_lines(cluster_names, file.path(output_dir, "cluster_names.txt"))


# 4. Handle Coordinates (LiftOver hg19 -> hg38 if needed)
df$POS <- as.integer(df$POS)

if (grepl("38", genome_build) && !grepl("native", genome_build, ignore.case = TRUE)) {
  cat("Performing LiftOver to hg38...\n")
  if (!file.exists(chain_file)) stop("Chain file not found: ", chain_file)

  chain <- import.chain(chain_file)
  gr <- GRanges(seqnames = paste0("chr", df$CHR),
                ranges = IRanges(start = df$POS, width = 1))

  lift_result <- liftOver(gr, chain)

  # Keep only variants that mapped successfully
  kept_indices <- which(elementNROWS(lift_result) > 0)
  df <- df[kept_indices, ]

  mapped_gr <- unlist(lift_result)
  df$POS_Final <- start(mapped_gr)
  df$CHR_Final <- gsub("chr", "", as.character(seqnames(mapped_gr))) # Remove 'chr' prefix for now

  cat("LiftOver complete. Dropped", length(lift_result) - length(kept_indices), "variants.\n")
} else {
  if (grepl("38", genome_build)) {
    cat("Keeping coordinates as provided (already hg38)...\n")
  } else {
    cat("Keeping coordinates as provided (hg19)...\n")
  }
  df$POS_Final <- df$POS
  df$CHR_Final <- as.character(df$CHR)
}

# 5. formatting for PLINK
# Ensure Risk_Allele matches REF or ALT
df <- df %>%
  mutate(
    valid_risk = Risk_Allele == REF | Risk_Allele == ALT,
    # Construct distinct VAR_ID for tracking
    VAR_ID = paste(CHR_Final, POS_Final, REF, ALT, sep = ":")
  )

if(any(!df$valid_risk)) {
  warning(sum(!df$valid_risk), " variants have Risk_Allele matching neither REF nor ALT. Removing them.")
  df <- df %>% filter(valid_risk)
}

# 6. Generate Output Files
final_df <- df %>%
  dplyr::select(VAR_ID, Risk_Allele, Total_GRS = BETA, all_of(cluster_names)) %>%
  # Apply cutoff to cluster weights
  mutate(across(all_of(cluster_names), ~ ifelse(abs(.) < cutoff, 0, .)))

# Write PLINK format
output_weights_file <- file.path(output_dir, "pPS_weights_plink_format.tsv")
write_tsv(final_df, output_weights_file)

# Write BCFtools input (BED format)
bcftools_input <- df %>%
  dplyr::select(Chr = CHR_Final, start = POS_Final) %>%
  mutate(Chr = paste0("chr", Chr), end = start) %>%
  dplyr::select(Chr, start, end) %>%
  arrange(Chr, start)

write_delim(bcftools_input, file.path(output_dir, "bcftools_input.txt"), delim = "\t", col_names = FALSE)

cat("Done. Processed", nrow(final_df), "variants.\n")