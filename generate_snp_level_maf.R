
#!/usr/bin/env Rscript
# generate_snp_level_maf.R
# Final robust version: Handles CHR type mismatches (char vs int)

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggplot2)
})

# 1. Parse Arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript generate_snp_level_maf.R <maf_file> <weights_file> <output_dir>")
}

maf_file <- args[1]
weights_file <- args[2]
output_dir <- args[3]

cat("Reading MAF data from:", maf_file, "\n")
cat("Reading Weights from:", weights_file, "\n")

# 2. Read MAF Data
# Expecting columns: CHROM, POS, REF, ALT, AF (Allele Frequency)

maf_data <- fread(maf_file, header = FALSE, col.names = c("CHR", "POS", "REF", "ALT", "AF")) %>%
  # [FIX] Robustly handle chr prefixes in VCF output
  mutate(CHR = as.integer(gsub("chr", "", as.character(CHR))))

# FIX: Normalize MAF Data Types (Strip 'chr' and ensure integer)
maf_data <- maf_data %>%
  mutate(CHR = gsub("^chr", "", CHR)) %>%  # Remove 'chr' prefix
  mutate(CHR = as.integer(CHR)) %>%         # Convert to integer
  mutate(POS = as.integer(POS)) %>%         # Ensure POS is integer
  mutate(MAF = ifelse(AF > 0.5, 1 - AF, AF)) %>%
  mutate(VAR_ID = paste(CHR, POS, REF, ALT, sep = ":"))

# 3. Read Weights
weights <- fread(weights_file) 

# Normalize Weights Column Names
if("VAR_ID_updated" %in% colnames(weights)) {
  weights <- weights %>% dplyr::rename(VAR_ID = VAR_ID_updated)
}



# Ensure Weights have CHR/POS columns as integers
if(!"CHR" %in% colnames(weights)) {
  var_parts <- tstrsplit(weights[["VAR_ID"]], ":", fixed = TRUE)
  weights <- weights %>%
    mutate(CHR = as.integer(gsub("^chr", "", var_parts[[1]])),
           POS = as.integer(var_parts[[2]]),
           REF = var_parts[[3]],
           ALT = var_parts[[4]])
} else {
  # If CHR exists, ensure it is integer
  weights <- weights %>%
    mutate(CHR = as.integer(gsub("^chr", "", CHR))) %>%
    mutate(POS = as.integer(POS))
  # Extract REF/ALT from VAR_ID if not already present (needed for allele-specific join)
  if (!"REF" %in% colnames(weights)) {
    var_parts <- tstrsplit(weights[["VAR_ID"]], ":", fixed = TRUE)
    weights <- weights %>%
      mutate(REF = var_parts[[3]], ALT = var_parts[[4]])
  }
}

# 4. Merge Data
cat("Merging datasets...\n")
merged_data <- weights %>%
  dplyr::select(VAR_ID, CHR, POS, REF, ALT, risk_allele) %>%
  inner_join(maf_data %>% dplyr::select(CHR, POS, REF, ALT, AF, MAF), by = c("CHR", "POS", "REF", "ALT")) %>%
  mutate(risk_allele_freq = case_when(risk_allele==ALT ~ AF,
                                      risk_allele==REF ~ 1-AF,
                                      TRUE ~ NA_real_)) %>%
  dplyr::select(VAR_ID, risk_allele, risk_allele_freq, everything())

cat("Merged", nrow(merged_data), "variants with MAF data.\n")

# 5. Save Table
output_table <- file.path(output_dir, "snp_level_maf.txt")
write_tsv(merged_data, output_table)
cat("Saved RAF table to:", output_table, "\n")

# 6. Generate QC Plot
p <- ggplot(merged_data, aes(x = risk_allele_freq)) +
  geom_histogram(binwidth = 0.01, fill = "forestgreen", color = "white") +
  theme_minimal() +
  labs(title = "Risk Allele Frequency Distribution",
       subtitle = paste("Cohort:", basename(output_dir)),
       x = "RAF",
       y = "Count")

output_plot <- file.path(output_dir, "snp_level_maf_distribution.png")
ggsave(output_plot, p, width = 10, height = 8, units='in',dpi = 300)
cat("Saved MAF plot to:", output_plot, "\n")