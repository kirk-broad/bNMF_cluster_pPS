
#!/usr/bin/env Rscript
# update_weights_allele_alignment.R
# Diagnostic Version 2 - Deep Dive on Risk Alleles

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
weights_file <- args[1]
vcf_allele_file <- args[2]
output_file <- ifelse(length(args) >= 3, args[3], "updated_weights.txt")

cat("Reading weights:", weights_file, "\n")
cat("Reading VCF info:", vcf_allele_file, "\n")

# 1. READ WEIGHTS
weights <- fread(weights_file, header = TRUE, sep="\t")

# Dynamically determine cluster names from the columns
# We ignore standard columns to isolate just the cluster/score columns
standard_cols <- c("VAR_ID", "VAR_ID_hg38", "chr", "pos", "REF", "ALT", "risk_allele", "Risk_Allele", "w_a1", "w_a2", "a1", "a2", "CHR", "POS", "BETA", "SNP", "Proxy_Flag", "orig_snp")
cluster_names <- setdiff(colnames(weights), standard_cols)

# Ensure Total_GRS is included if it exists in the file
if("Total_GRS" %in% colnames(weights) && !"Total_GRS" %in% cluster_names) {
  cluster_names <- unique(c("Total_GRS", cluster_names))
}

cat("Detected cluster columns:", paste(cluster_names, collapse=", "), "\n")

# --- INSERT FIX HERE (Must be before diagnostics) ---
# Check for capitalization mismatch (common if Proxy Step was skipped)
if ("Risk_Allele" %in% names(weights) && !"risk_allele" %in% names(weights)) {
  cat("Renaming 'Risk_Allele' to 'risk_allele'...\n")
  weights <- weights %>% dplyr::rename(risk_allele = Risk_Allele)
}
# Also handle VAR_ID mismatch if present
if ("VAR_ID_hg38" %in% names(weights) && !"VAR_ID" %in% names(weights)) {
  weights <- weights %>% dplyr::rename(VAR_ID = VAR_ID_hg38)
}
# ----------------------------------------------------

# DIAGNOSTIC: Check Risk Allele immediately
cat("Raw Weights Rows:", nrow(weights), "\n")
cat("Risk Allele Column Sample (First 5):\n")
print(head(weights$risk_allele, 5))
na_risk <- sum(is.na(weights$risk_allele))
cat("Number of NA risk alleles:", na_risk, "\n")

# Clean
var_parts <- tstrsplit(weights[["VAR_ID"]], ":", fixed = TRUE)
weights <- weights %>%
  mutate(risk_allele = toupper(trimws(risk_allele)),
         chr = var_parts[[1]], pos = as.integer(var_parts[[2]]),
         w_a1 = toupper(trimws(var_parts[[3]])),
         w_a2 = toupper(trimws(var_parts[[4]])))

# 2. READ VCF INFO
vcf_info <- fread(vcf_allele_file, header = FALSE, sep = "\t", col.names = c("VAR_ID_VCF"))
vcf_parts <- tstrsplit(vcf_info[["VAR_ID_VCF"]], ":", fixed = TRUE)
vcf_info <- vcf_info %>%
  mutate(chr = vcf_parts[[1]], pos = as.integer(vcf_parts[[2]]),
         REF = toupper(trimws(vcf_parts[[3]])),
         ALT = toupper(trimws(vcf_parts[[4]])))

cat("VCF entries:", nrow(vcf_info), "\n")

# 3. NORMALIZE FOR JOIN
weights$chr_norm <- gsub("^chr", "", as.character(weights$chr), ignore.case = TRUE)
vcf_info$chr_norm <- gsub("^chr", "", as.character(vcf_info$chr), ignore.case = TRUE)
weights$pos <- as.integer(weights$pos)
vcf_info$pos <- as.integer(vcf_info$pos)

# 4. MERGE
merged <- inner_join(weights, vcf_info, by = c("chr_norm", "pos"))
cat("Variants matched by CHR:POS:", nrow(merged), "\n")

if(nrow(merged) == 0) {
  stop("No variants matched by position!")
}

# 5. FILTER
merged_final <- merged %>%
  mutate(
    # Check if alleles match (either direct or flipped)
    match_direct = (w_a1 == REF & w_a2 == ALT),
    match_flip   = (w_a1 == ALT & w_a2 == REF),
    valid_alleles = coalesce(match_direct | match_flip, FALSE),
    
    # Check if risk allele exists in VCF (Handle NA explicitly)
    risk_matches_ref = (risk_allele == REF),
    risk_matches_alt = (risk_allele == ALT),
    risk_valid = coalesce(risk_matches_ref | risk_matches_alt, FALSE)
  )

# 6. DETAILED DIAGNOSTICS
cat("--- FILTERING DIAGNOSTICS ---\n")
cat("1. Allele Pairs Mismatch (Weights vs VCF):", sum(!merged_final$valid_alleles), "\n")
cat("2. Risk Allele Invalid (Not REF or ALT):  ", sum(!merged_final$risk_valid), "\n")

# Debug the failures
if(sum(!merged_final$risk_valid) > 0) {
  cat("\nExample FAILED Risk Alleles:\n")
  bad_rows <- merged_final %>% filter(!risk_valid) %>% head(5)
  print(bad_rows %>% select(VAR_ID, REF, ALT, risk_allele, w_a1, w_a2))
}

# 7. FINALIZE
final_output <- merged_final %>%
  filter(valid_alleles & risk_valid) %>%
  mutate(
    chr_final = ifelse(!is.na(chr.y), chr.y, chr.x),
    VAR_ID_updated = paste(chr_final, pos, REF, ALT, sep = ":")
  ) %>%
  dplyr::select(VAR_ID_updated, risk_allele, all_of(cluster_names))

cat("Writing", nrow(final_output), "aligned variants to", output_file, "\n")
fwrite(final_output, file = output_file, sep = "\t")