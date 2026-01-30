
#!/usr/bin/env Rscript
# cluster_weights_setup_generalized.R
# Usage: Rscript cluster_weights_setup_generalized.R <weights_dir> <cutoff> <chain_file> <genome_build> <output_dir>

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(readxl)
  library(rtracklayer)
  library(GenomicRanges)
})

# 1. Parse Arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("Usage: script.R <weights_dir> <cutoff> <chain_file> <genome_build> <output_dir>", call.=FALSE)
}

weights_dir <- args[1]
cutoff <- as.numeric(args[2])
chain_file <- args[3]
genome_build <- args[4] # "hg19" or "hg38"
output_dir <- args[5]

cat("------------------------------------------------\n")
cat("STEP 1: Setting up Cluster Weights\n")
cat("Weights Dir:  ", weights_dir, "\n")
cat("Target Build: ", genome_build, "\n")
cat("------------------------------------------------\n")

# 2. Read Cluster Names & Weights
cluster_names_file <- file.path(weights_dir, 'cluster_names.txt')
if (!file.exists(cluster_names_file)) stop("Cluster names file not found: ", cluster_names_file)

cluster_names <- make.unique(read_lines(cluster_names_file))
k <- length(cluster_names)

# Dynamically find the Excel file based on K
final_col <- LETTERS[k+3]
weights_filename <- sprintf("sorted_cluster_weights_K%i_rev.xlsx", k)
weights_path <- file.path(weights_dir, weights_filename)

if (!file.exists(weights_path)) stop("Weights Excel file not found: ", weights_path)

cat("Reading raw weights from:", weights_filename, "\n")
df_weights <- read_excel(weights_path, range = sprintf("A1:%s2000", final_col))
names(df_weights)[4:(k+3)] <- cluster_names

# Clean and Prep Data (Assuming raw input is ALWAYS hg19)
df_weights <- df_weights %>%
  separate(VAR_ID, into=c("CHR","POS_hg19","REF","ALT"), sep = "_", remove = F) %>%
  mutate(ChrPos = paste(CHR, POS_hg19, sep = ":")) %>%
  mutate(group = seq.int(nrow(.))) %>%
  mutate(POS_hg19 = as.integer(POS_hg19)) %>%
  dplyr::select(group, ChrPos, VAR_ID_hg19=VAR_ID, CHR, POS_hg19, REF, ALT, rsID, all_of(cluster_names)) %>%
  drop_na()

# 3. Load GWAS Summary Stats
gwas_stats_path <- file.path(weights_dir, "alignment_GWAS_summStats.csv")
if (!file.exists(gwas_stats_path)) stop("GWAS stats file not found: ", gwas_stats_path)

input_snps <- fread(gwas_stats_path, stringsAsFactors = FALSE, data.table = F)

# Handle GWAS column variations
if (!"Other_Allele" %in% colnames(input_snps)) {
  input_snps <- input_snps %>%
    mutate(Other_Allele = case_when(
      Risk_Allele == REF ~ ALT,
      Risk_Allele == ALT ~ REF,
      TRUE ~ NA_character_
    ))
}
if (!"ChrPos" %in% names(input_snps)) {
  input_snps <- input_snps %>% dplyr::rename(ChrPos = SNP)
}

input_snps <- input_snps %>%
  dplyr::select(ChrPos, Risk_Allele, Other_Allele, REF, ALT, P_VALUE, BETA_aligned) %>%
  mutate(BETA_aligned = ifelse(BETA_aligned == Inf, max(BETA_aligned[BETA_aligned != Inf], na.rm=T), BETA_aligned))

# 4. Handle Coordinates (Conditional LiftOver)
df_weights_rev2 <- df_weights %>% mutate_at(.vars = c('CHR','POS_hg19'), as.integer)

if (grepl("38", genome_build)) {
  # --- CASE A: HG38 (Perform LiftOver) ---
  cat("Performing LiftOver to hg38...\n")
  
  if (!file.exists(chain_file)) stop("Chain file required for hg38 but not found: ", chain_file)
  chain <- import.chain(chain_file)
  
  gr <- GRanges(seqnames = paste0("chr", df_weights_rev2$CHR),
                strand = "*",
                ranges = IRanges(start = df_weights_rev2$POS_hg19, width = 1))
  
  lift_result <- liftOver(gr, chain)
  mapped_gr <- unlist(lift_result)
  mapped_gr$group <- rep(seq_along(lift_result), elementNROWS(lift_result))
  
  df_coord_map <- as.data.frame(mapped_gr) %>%
    dplyr::select(group, seqnames, Pos_Final = start) %>%
    mutate(CHR_Final = gsub("chr", "", seqnames))
    
} else {
  # --- CASE B: HG19 (Keep Original Coordinates) ---
  cat("Keeping coordinates in hg19 (No LiftOver)...\n")
  
  df_coord_map <- df_weights_rev2 %>%
    dplyr::select(group, CHR, POS_hg19) %>%
    mutate(CHR_Final = as.character(CHR), 
           Pos_Final = POS_hg19)
}

# 5. Merge and Align
input_aligned <- df_weights_rev2 %>%
  inner_join(df_coord_map, by = "group") %>%
  inner_join(input_snps, by = "ChrPos") %>% # Join by hg19 position to get alleles
  mutate(
    weights_risk_allele = case_when(
      ALT.x == Risk_Allele ~ ALT.x,
      REF.x == Risk_Allele ~ REF.x,
      ALT.x == Other_Allele ~ REF.x,
      REF.x == Other_Allele ~ ALT.x,
      TRUE ~ NA_character_
    ),
    # Construct Final VAR_ID
    VAR_ID_Final = paste(CHR_Final, Pos_Final, 
                         ifelse(weights_risk_allele == REF.x, REF.x, ALT.x), # Final REF
                         ifelse(weights_risk_allele == REF.x, ALT.x, REF.x), # Final ALT
                         sep = ":"),
    alignment_status = case_when(
      is.na(weights_risk_allele) ~ "failed",
      weights_risk_allele == Risk_Allele ~ "direct_match",
      TRUE ~ "strand_issue"
    )
  ) %>%
  filter(alignment_status == "direct_match") %>%
  mutate(Risk_Allele_final = weights_risk_allele, BETA_final = BETA_aligned)

# 6. Generate Output Files
cat("Successfully aligned:", nrow(input_aligned), "variants.\n")

final_df <- input_aligned %>%
  dplyr::select(VAR_ID = VAR_ID_Final, 
                Risk_Allele = Risk_Allele_final, 
                Total_GRS = BETA_final, 
                all_of(cluster_names)) %>%
  mutate(across(all_of(cluster_names), ~ ifelse(. < cutoff, 0, .)))

# Save PLINK format weights
output_weights_file <- file.path(output_dir, "pPS_weights_plink_format.tsv")
write_tsv(final_df, output_weights_file)

# Create BCFtools input
bcftools_input <- final_df %>%
  separate(VAR_ID, into = c("chr", "pos", "ref", "alt"), sep = ":", remove = FALSE) %>%
  arrange(as.integer(chr), as.integer(pos)) %>%
  mutate(Chr = paste0("chr", chr), start = pos, end = pos) %>%
  dplyr::select(Chr, start, end)

output_bed_file <- file.path(output_dir, "bcftools_input.txt")
write_delim(bcftools_input, output_bed_file, delim = "\t", col_names = FALSE)
write_lines(cluster_names, file.path(output_dir, "cluster_names.txt"))

cat("Done. Weights and bed file saved to:", output_dir, "\n")