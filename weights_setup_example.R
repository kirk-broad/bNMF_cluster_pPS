#!/usr/bin/env Rscript
# prep_bNMF_inputs.R
# Purpose: Converts raw bNMF Excel weights and GWAS stats into the standard
#          "cluster_weights.tsv" format required by the PPS pipeline.

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(readxl)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript prep_bNMF_inputs.R <weights_dir>")

weights_dir <- args[1]
cat("------------------------------------------------\n")
cat("PREP: Converting bNMF outputs to Standard Input\n")
cat("------------------------------------------------\n")

# 1. Read Cluster Names
cluster_names_file <- file.path(weights_dir, 'cluster_names.txt')
if (!file.exists(cluster_names_file)) stop("cluster_names.txt not found")
cluster_names <- read_lines(cluster_names_file)
k <- length(cluster_names)

# 2. Find and Read Excel Weights
weights_filename <- sprintf("sorted_cluster_weights_K%i_rev.xlsx", k)
weights_path <- file.path(weights_dir, weights_filename)
if (!file.exists(weights_path)) stop("Excel weights file not found: ", weights_filename)

cat("Reading Excel weights:", weights_filename, "\n")
final_col <- LETTERS[k+3] 
df_weights <- read_excel(weights_path, range = sprintf("A1:%s2000", final_col))
names(df_weights)[4:(k+3)] <- cluster_names

# Clean Weights Data
df_weights_clean <- df_weights %>%
  separate(VAR_ID, into=c("CHR","POS","REF","ALT"), sep = "_", remove = FALSE) %>%
  mutate(POS = as.integer(POS)) %>%
  mutate(ChrPos = paste(CHR, POS, sep = ":")) %>%
  dplyr::select(ChrPos, CHR, POS, REF, ALT, all_of(cluster_names)) %>%
  drop_na()

# 3. Read GWAS Stats & PERFORM ALIGNMENT
gwas_file <- file.path(weights_dir, "alignment_GWAS_summStats.csv")
if (!file.exists(gwas_file)) stop("alignment_GWAS_summStats.csv not found")

cat("Reading GWAS stats:", gwas_file, "\n")
df_gwas <- fread(gwas_file)

# Handle column naming variations
if(!"ChrPos" %in% names(df_gwas)) df_gwas$ChrPos <- df_gwas$SNP

# --- YOUR ALIGNMENT LOGIC HERE ---
cat("Aligning Betas to Risk Allele...\n")
df_gwas_clean <- df_gwas %>%
  mutate(
    # Ensure Other_Allele is consistently set
    Other_Allele = case_when(
      Risk_Allele == ALT ~ REF,
      Risk_Allele == REF ~ ALT,
      TRUE ~ NA_character_
    ),
    # Flip BETA if Risk_Allele is REF (assuming input BETA is relative to ALT)
    BETA_aligned = case_when(
      Risk_Allele == ALT ~ BETA,
      Risk_Allele == REF ~ -BETA,
      TRUE ~ NA_real_
    )
  ) %>%
  dplyr::select(ChrPos, Risk_Allele, Other_Allele, BETA_aligned)

# Logic Check: Ensure all weights are positive (risk-increasing)
if (all(df_gwas_clean$BETA_aligned > 0, na.rm = TRUE)) {
  cat("SUCCESS: All aligned betas are positive.\n")
} else {
  warning("WARNING: Some aligned betas are negative. Check input GWAS stats.")
  print(head(df_gwas_clean %>% filter(BETA_aligned < 0)))
}

# 4. Merge and Create Standard Output
merged <- df_weights_clean %>%
  inner_join(df_gwas_clean, by = "ChrPos") %>%
  # Rename BETA_aligned to BETA for the standard format
  dplyr::select(CHR, POS, REF, ALT, Risk_Allele, BETA = BETA_aligned, all_of(cluster_names))

output_file <- file.path(weights_dir, "cluster_weights.tsv")
write_tsv(merged, output_file)

cat("Created standardized input:", output_file, "\n")



#---- PLOS ----

summ_stats <- fread("../../../PLOS/PLOS_94_SNPs_MVP.DIAMANTE_SummStats.txt")
snps <- fread("../../../PLOS/PLOS_94_SNPs.csv")

weights <- fread("../../../PLOS/PLOS_weights_combined.csv")
clusters_plos <- unique(weights$Cluster)

merged_plos <- snps %>%
  left_join(weights, by=c('variant'='SNP')) %>%
  pivot_wider(id_cols = c("VAR_ID_hg19","RISK_ALLELE"),names_from = "Cluster",values_from = "weight") %>%
  inner_join(summ_stats %>% select(VAR_ID_hg19=VAR_ID, ODDS_RATIO), by="VAR_ID_hg19") %>%
  separate(VAR_ID_hg19, into=c("CHR","POS","REF","ALT"),sep = "_") %>%
  rename(Risk_Allele=RISK_ALLELE) %>%
  
  mutate(    # Ensure Other_Allele is consistently set
    Other_Allele = case_when(
      Risk_Allele == ALT ~ REF,
      Risk_Allele == REF ~ ALT,
      TRUE ~ NA_character_
    ),
    # Flip BETA if Risk_Allele is REF (assuming input BETA is relative to ALT)
    BETA_aligned = case_when(
      Risk_Allele == ALT ~ log(ODDS_RATIO),
      Risk_Allele == REF ~ -log(ODDS_RATIO),
      TRUE ~ NA_real_
    )
  ) %>%
  mutate(across(clusters_plos, ~replace_na(.x, 0))) %>%

  select(CHR, POS, REF, ALT, Risk_Allele, BETA_aligned, all_of(clusters_plos))
names(merged_plos) <- gsub("\\/", "\\_", names(merged_plos))

if (all(merged_plos$BETA_aligned > 0, na.rm = TRUE)) {
  cat("SUCCESS: All aligned betas are positive.\n")
} else {
  warning("WARNING: Some aligned betas are negative. Check input GWAS stats.")
  print(head(merged_plos %>% filter(BETA_aligned < 0)))
}



output_file_plos <- file.path("../../../PLOS/", "cluster_weights.tsv")
write_tsv(merged_plos, output_file_plos)

#---- Suzuki ----
library(strex)
df_suzuki_clusters <- fread("../../../Collaborations/T2DGGI/DATA/GGI_ST6_snp_clusters.csv")
suzuki_clusters <- unique(df_suzuki_clusters %>% filter(`Cluster assignment` != "") %>% pull(`Cluster assignment`))
suzuki_clusters
suzuki_rsIDs <- fread("../../../Collaborations/T2DGGI/DATA/GGI_1289snps_rsID_map.txt")
suzuki_gwas <- fread("../../../Collaborations/T2DGGI/DATA/Suzuki_alignment_GWAS_summStats.csv")

merged_suzuki <- suzuki_rsIDs %>%
  separate(VAR_ID, into=c("CHR","POS","REF","ALT"),sep = "_") %>%
  mutate(ChrPos = paste(CHR, POS, sep=":")) %>%
  inner_join(suzuki_gwas, by='ChrPos') %>%
  filter((ALT.x==ALT.y & REF.x==REF.y) | (ALT.x==REF.y & REF.x==ALT.y)) %>%
  inner_join(df_suzuki_clusters %>% select(rsID=`Index SNV`, cluster=`Cluster assignment`), by='rsID') %>%
  mutate(BETA_aligned = case_when(
                                Risk_Allele == ALT.y ~ BETA,
                                Risk_Allele == REF.y ~ -BETA,
                                TRUE ~ NA_real_
                              ),
         weight = BETA_aligned,
  ) %>%
  pivot_wider(id_cols = c("CHR","POS","REF.y","ALT.y","Risk_Allele","BETA_aligned"),names_from = "cluster",values_from = "weight") %>%
  mutate(across(suzuki_clusters, ~replace_na(.x, 0))) %>%
  select(CHR, POS, REF=REF.y, ALT=ALT.y, Risk_Allele, BETA_aligned, all_of(suzuki_clusters))
names(merged_suzuki) <- gsub("\\ ", "\\_", names(merged_suzuki))

  

if (all(merged_suzuki$BETA_aligned > 0, na.rm = TRUE)) {
  cat("SUCCESS: All aligned betas are positive.\n")
} else {
  warning("WARNING: Some aligned betas are negative. Check input GWAS stats.")
  print(head(merged_suzuki %>% filter(BETA_aligned < 0)))
}


output_file_suzuki <- file.path("../../../Collaborations/T2DGGI/DATA", "cluster_weights.tsv")
write_tsv(merged_suzuki, output_file_suzuki)



