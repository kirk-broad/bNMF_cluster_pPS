#!/usr/bin/env Rscript
# update_weights_with_proxies.R
# Final Debug Version: Uses coalesce() and traces risk_allele integrity

library(LDlinkR)
library(data.table)
library(strex)
library(tidyverse)
library(dplyr)

# --- Parameters ---
weights_file <- "./pPS_files/pPS_weights_plink_format.tsv"
output_dir <- "./pPS_files"
ldlink_token <- "YOUR_LDLINK_TOKEN_HERE" 
genome_build_input <- "grch38"
vcf_file_template <- "VCF_FILE_TEMPLATE_PLACEHOLDER"
perform_search <- TRUE  # <--- NEW TOGGLE (Default to TRUE)

# --- Normalize Build ---
api_build <- ifelse(grepl("37|19", genome_build_input), "grch37", "grch38")
cat("Using Genome Build for LDlink:", api_build, "\n")

# --- Read the weights file ---
cat("Reading weights file from:", weights_file, "\n")
weights <- fread(weights_file)

# --- CHECKPOINT 1: Load ---
cat("DEBUG [1]: Loaded weights. Rows:", nrow(weights), "\n")
risk_col_idx <- grep("^risk_allele$", names(weights), ignore.case = TRUE)
if(length(risk_col_idx) > 0) {
  cat("DEBUG [1]: Found risk column:", names(weights)[risk_col_idx], "\n")
  # Rename immediately
  names(weights)[risk_col_idx] <- "risk_allele"
  # Force Character
  weights$risk_allele <- as.character(weights$risk_allele)
  cat("DEBUG [1]: Sample:", paste(head(weights$risk_allele, 3), collapse=","), "\n")
} else {
  stop("CRITICAL: No Risk Allele column found in input file!")
}

# --- Column Normalization (VAR_ID) ---
if ("VAR_ID" %in% names(weights)) {
  # do nothing
} else if (api_build == "grch37" && "VAR_ID_hg19" %in% names(weights)) {
  weights <- weights %>% dplyr::rename(VAR_ID = VAR_ID_hg19)
} else if ("VAR_ID_hg38" %in% names(weights)) {
  weights <- weights %>% dplyr::rename(VAR_ID = VAR_ID_hg38)
} else {
  stop("Could not find VAR_ID column!")
}

# --- Identify missing SNPs ---
missing_snps_file <- "./pPS_files/missing_snps.txt"
if(file.exists(missing_snps_file) && file.size(missing_snps_file) > 0) {
  df_missing_snps <- fread(missing_snps_file, col.names = c("CHR","POS1","POS2"), data.table = FALSE) %>%
    mutate(SNP = paste0(CHR,":",POS1))
  cat("Missing SNPs:", nrow(df_missing_snps), "\n")
  missing_snps <- df_missing_snps$SNP
} else {
  missing_snps <- character(0)
}

# --- HELPER: Safe Query ---
query_ldlink_safe <- function(snp, attempt=1) {
  if(attempt > 5) return(NULL)
  if(attempt > 1) {
    cat("  [Retry", attempt, "] Waiting", attempt*5, "s...\n")
    Sys.sleep(attempt * 5)
  } else {
    Sys.sleep(3) # Base delay
  }
  
  res <- tryCatch({
    LDproxy(snp, pop = "EUR", r2d = "r2", token = ldlink_token, genome_build = api_build)
  }, error = function(e) return(NULL))
  
  if(is.null(res)) return(query_ldlink_safe(snp, attempt + 1))
  
  if(is.data.frame(res) && nrow(res) >= 1) {
    first_cell <- as.character(res[1,1])
    if(grepl("error", first_cell, ignore.case=TRUE)) {
      if(grepl("concurrent", first_cell, ignore.case=TRUE)) return(query_ldlink_safe(snp, attempt + 1))
      return(NULL)
    }
  }
  return(res)
}

# --- Query Loop ---
# --- Query Loop ---

valid_proxies <- list()

# [CHANGE] Add check for perform_search
if(perform_search && length(missing_snps) > 0) {
  cat("Searching for proxies enabled...\n")

  for(snp in missing_snps) {
    cat("Processing:", snp, "...")
    res <- query_ldlink_safe(snp)
    
    if(!is.null(res) && is.data.frame(res) && "R2" %in% names(res)) {
      res$R2 <- as.numeric(as.character(res$R2))
      proxies <- res[res$R2 >= 0.5 & res$Coord != snp, ]
      
      if(nrow(proxies) > 0) {
        proxies <- proxies[order(-proxies$R2), ]
        
        # 1. Determine correct VCF file path
        chr_raw <- str_before_first(snp, ":")
        chr_num <- gsub("chr", "", chr_raw, ignore.case = TRUE)
        vcf_file <- sprintf(vcf_file_template, chr_num)
        
        # 2. DETECT VCF FORMAT (hg19/hg38 check)
        # We check the VCF header once per SNP to be safe
        vcf_has_chr <- tryCatch({
          # Peek at first variant
          header_chk <- system(paste("bcftools view -H", vcf_file, "| head -n 1 | cut -f1"), intern = TRUE)
          grepl("chr", header_chk)
        }, error = function(e) FALSE)
        
        # 3. Check proxies against VCF
        found_proxy <- FALSE
        for(i in 1:min(nrow(proxies), 15)) {
          cand <- proxies[i, ]
          cand$orig_snp <- snp
          
          # ADJUST COORDINATE: Match LDlink output (chr1:123) to VCF format
          check_coord <- cand$Coord
          if(vcf_has_chr && !grepl("chr", check_coord)) {
            check_coord <- paste0("chr", check_coord)
          } else if(!vcf_has_chr && grepl("chr", check_coord)) {
            check_coord <- gsub("chr", "", check_coord)
          }
          
          check_cmd <- sprintf("bcftools view -r %s %s | grep -v '^#' | wc -l", check_coord, vcf_file)
          n <- tryCatch({ as.integer(system(check_cmd, intern = TRUE)) }, error = function(e) 0)
          
          if(n > 0) {
            cat(" Replaced with", cand$RS_Number, "\n")
            valid_proxies[[length(valid_proxies)+1]] <- cand
            found_proxy <- TRUE
            break
          }
        }
        if(!found_proxy) cat(" No matching proxy in VCF\n")
        
      } else { cat(" No proxies > R2 0.5\n") }
    } else { cat(" LDlink query failed\n") }
  }
} else {
  cat("Proxy search skipped (Disabled or no missing SNPs). Missing SNPs will be flagged as 'No_proxy'.\n")
}



# --- Process Results ---
if(length(valid_proxies) > 0) {
  df_my_proxies <- rbindlist(valid_proxies) %>%
    mutate(
      proxy_a1 = str_extract(Alleles, "(?<=\\()[^/]+"),
      proxy_a2 = str_extract(Alleles, "(?<=/)[^\\)]+"),
      proxy_VAR_ID = paste(Coord, proxy_a1, proxy_a2, sep = ":")
    )
  
  weights_with_snp <- weights %>% mutate(orig_snp = str_before_nth(VAR_ID, ":", 2))
  
  df_proxy_map <- weights_with_snp %>%
    inner_join(df_my_proxies, by = 'orig_snp') %>%
    mutate(proxy_pos = as.integer(str_after_first(Coord, ":")), Proxy_Flag = "Replaced") 
  
  get_proxy_risk_allele <- function(correlated_str, orig_risk) {
    if(is.na(correlated_str) || is.na(orig_risk)) return(NA_character_)
    pairs <- str_split(correlated_str, ",")[[1]]
    for(pair in pairs) {
      if(!grepl("=", pair)) next
      split_pair <- str_split(pair, "=")[[1]]
      if(length(split_pair) == 2 && str_trim(split_pair[1]) == orig_risk) {
        val <- toupper(trimws(split_pair[2]))
        if(val %in% c("NULL","NA","","."))
          return(NA_character_)
        return(val)
      }   
    }
    return(NA_character_)
  }
  
  df_proxy_map$proxy_risk_allele <- mapply(get_proxy_risk_allele, df_proxy_map$Correlated_Alleles, df_proxy_map$risk_allele)
  
  df_proxy_map <- df_proxy_map %>%
    dplyr::select(VAR_ID, proxy_VAR_ID, proxy_pos, proxy_a1, proxy_a2, proxy_risk_allele, Proxy_Flag)
  
} else {
  df_proxy_map <- data.frame(VAR_ID=character(), proxy_VAR_ID=character(), proxy_pos=integer(), proxy_a1=character(), proxy_a2=character(), proxy_risk_allele=character(), Proxy_Flag=character())
}

# --- FINAL MERGE (THE FIX) ---
if(!"pos" %in% names(weights)) {
  var_parts <- tstrsplit(weights[["VAR_ID"]], ":", fixed = TRUE)
  weights <- weights %>%
    mutate(chr = var_parts[[1]], pos = as.integer(var_parts[[2]]),
           a1 = var_parts[[3]], a2 = var_parts[[4]])
}

# 1. Prepare base status
all_snps_status <- weights %>%
  mutate(orig_snp = str_before_nth(VAR_ID, ":", 2)) %>%
  mutate(Proxy_Flag = ifelse(orig_snp %in% missing_snps, 
                             ifelse(VAR_ID %in% df_proxy_map$VAR_ID, "Replaced", "No_proxy"), 
                             "Original"))

cat("DEBUG [2]: Before Merge. Risk Allele NAs:", sum(is.na(all_snps_status$risk_allele)), "\n")

# 2. Join
new_weights <- all_snps_status %>%
  left_join(df_proxy_map, by = "VAR_ID") 

cat("DEBUG [3]: After Join. Risk Allele NAs:", sum(is.na(new_weights$risk_allele)), "\n")

# 3. Mutate with COALESCE (Safer than ifelse)
new_weights <- new_weights %>%
  mutate(
    VAR_ID      = coalesce(proxy_VAR_ID, VAR_ID),
    pos         = coalesce(proxy_pos, pos),
    a1          = coalesce(proxy_a1, a1),
    a2          = coalesce(proxy_a2, a2),
    
    proxy_risk_allele_clean = toupper(trimws(as.character(proxy_risk_allele))),
    proxy_risk_allele_clean = na_if(proxy_risk_allele_clean, "NULL"),
    proxy_risk_allele_clean = na_if(proxy_risk_allele_clean, ""),
    proxy_risk_allele_clean = na_if(proxy_risk_allele_clean, "NA"),
    proxy_risk_allele_clean = na_if(proxy_risk_allele_clean, "."),
    
    risk_allele_clean = toupper(trimws(as.character(risk_allele))),
    risk_allele_clean = na_if(risk_allele_clean, "NULL"),
    risk_allele_clean = na_if(risk_allele_clean, ""),
    risk_allele_clean = na_if(risk_allele_clean, "NA"),
    risk_allele_clean = na_if(risk_allele_clean, "."),
    
    risk_allele = coalesce(proxy_risk_allele_clean, risk_allele_clean),
    
    Proxy_Flag  = coalesce(Proxy_Flag.y, Proxy_Flag.x)
  ) %>%
  mutate(chr = paste0("chr", chr), SNP = paste(chr, pos, sep = ":")) %>%
  dplyr::select(-proxy_risk_allele_clean, -risk_allele_clean) %>%
  dplyr::select(all_of(c(names(weights), "SNP", "Proxy_Flag")))

cat("DEBUG [4]: Final. Risk Allele NAs:", sum(is.na(new_weights$risk_allele)), "\n")
cat("DEBUG [4]: Sample:", paste(head(new_weights$risk_allele, 3), collapse=","), "\n")
cat("DEBUG [X]: NULL-string risk alleles:", sum(new_weights$risk_allele %in% c("NULL","NA","")), "\n")

write_delim(new_weights, file.path(output_dir, "updated_cluster_weights.txt"), delim = "\t")
for_VCF <- new_weights %>% dplyr::select(Chr = chr, start = pos, end = pos)
write_delim(for_VCF, file.path(output_dir, "updated_bcftools_input.txt"), delim = "\t", col_names = FALSE)

cat("Done!\n")