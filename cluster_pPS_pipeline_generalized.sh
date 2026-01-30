#!/bin/bash
# cluster_pPS_pipeline_generalized.sh

# --- ARGUMENT PARSING ---
# Default to Step 1 if no argument provided.
# Run as "bash cluster_pPS_pipeline_generalized.sh 4" to start at Step 4
START_STEP=${1:-1}

#============================================================================
# SETUP AND CONFIGURATION
#============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -f "${SCRIPT_DIR}/cluster_pPS_config.sh" ]; then
    source "${SCRIPT_DIR}/cluster_pPS_config.sh"
else
    echo "ERROR: Config file not found."
    exit 1
fi

validate_config

# Setup environment
if [ "$USE_BROAD_MODULES" = true ]; then
    source /broad/software/scripts/useuse
    use Bcftools || true
    use Tabix || true
    use R-4.3 || true
elif [ "$USE_SLURM_MODULES" = true ]; then
    if [ -n "$CUSTOM_SOFTWARE_SETUP" ]; then
        eval "$CUSTOM_SOFTWARE_SETUP"
    fi
elif [ -n "$CUSTOM_SOFTWARE_SETUP" ]; then
    eval "$CUSTOM_SOFTWARE_SETUP"
fi

set -e 

print_config
mkdir -p "$OUTPUT_DIR" "$TEMP_DIR"

# Define this GLOBALLY so it is available even if Step 2 is skipped
# Define these GLOBALLY so they are available even if Step 2 is skipped
COMBINED_VCF="${OUTPUT_DIR}/${FINAL_VCF_NAME}"
VCF_LIST="${TEMP_DIR}/vcf_filelist.txt"

LOGFILE="${OUTPUT_DIR}/pipeline_${STUDY_NAME}_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOGFILE") 2>&1

echo "Pipeline started at: $(date)"
echo "Starting at Step: $START_STEP"

#============================================================================
# STEP 1: SETUP CLUSTER WEIGHTS
#============================================================================

if [ "$START_STEP" -le 1 ]; then
    echo "============================================================================"
    echo "STEP 1: Setting up cluster weights (Target Build: ${GENOME_BUILD})"
    echo "============================================================================"

    WEIGHTS_OUTPUT="${OUTPUT_DIR}/pPS_weights_plink_format.tsv"

    if [ -f "$CHAIN_FILE" ]; then
        CHAIN_FILE_PATH="$CHAIN_FILE"
    else
        CHAIN_FILE_PATH="${CLUSTER_WEIGHTS_DIR}/${CHAIN_FILE}"
    fi

    # 1. ENSURE WEIGHTS FILE EXISTS
    if [ -f "$WEIGHTS_OUTPUT" ] && [ -s "$WEIGHTS_OUTPUT" ]; then
        echo "Found existing weights file at: $WEIGHTS_OUTPUT"
        echo "Skipping setup script."
    else
        echo "Running setup script..."
        Rscript "${SCRIPT_DIR}/cluster_weights_setup_generalized.R" "$CLUSTER_WEIGHTS_DIR" "$CUTOFF" "$CHAIN_FILE_PATH" "$GENOME_BUILD" "$OUTPUT_DIR"
        
        if [ ! -f "$WEIGHTS_OUTPUT" ]; then 
            echo "ERROR: Weights file generation failed."
            exit 1
        fi
    fi

    # 2. DETECT VCF FORMAT (Auto-detect hg19 vs hg38 naming)
    # We do this every time to ensure compatibility with the current VCFs
    
    # Find a test VCF (chr22 or chr1)
    TEST_VCF=$(echo "$VCF_PATH_TEMPLATE" | sed 's/{CHR}/22/')
    if [ ! -f "$TEST_VCF" ]; then TEST_VCF=$(echo "$VCF_PATH_TEMPLATE" | sed 's/{CHR}/1/'); fi

    echo "Checking VCF naming convention in: $TEST_VCF"
    
    # Default to 'chr' prefix (hg38 standard)
    R_PASTE_CMD="paste0('chr', chr)"
    
    if [ -f "$TEST_VCF" ]; then
        # Peek at the first variant's CHROM column
        VCF_CHR_FORMAT=$($BCFTOOLS_PATH view -H "$TEST_VCF" | head -n 1 | cut -f1)
        
        if [[ "$VCF_CHR_FORMAT" == chr* ]]; then
            echo "Detected 'chr' prefix in VCF (e.g., $VCF_CHR_FORMAT). Format: chr1, chr22"
            R_PASTE_CMD="paste0('chr', chr)"
        else
            echo "Detected NO 'chr' prefix in VCF (e.g., $VCF_CHR_FORMAT). Format: 1, 22"
            R_PASTE_CMD="chr"
        fi
    else
        echo "WARNING: Could not find test VCF. Defaulting to 'chr' prefix."
    fi

    # 3. GENERATE BCFTOOLS INPUT FILE (Using detected format)
    echo "Generating bcftools input file..."
    Rscript -e "
    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(dplyr))
    suppressPackageStartupMessages(library(tidyr))
    
    weights <- fread('${WEIGHTS_OUTPUT}')
    
    # Handle column name variations
    if('VAR_ID_hg38' %in% names(weights)) { weights <- weights %>% rename(VAR_ID = VAR_ID_hg38) }

    bcftools_input <- weights %>% 
        separate(VAR_ID, into = c('chr', 'pos', 'a1', 'a2'), sep = ':', remove = FALSE) %>% 
        mutate(chr = gsub('chr', '', chr)) %>% 
        arrange(as.integer(chr), as.integer(pos)) %>% 
        # Apply Dynamic Prefix
        mutate(Chr = ${R_PASTE_CMD}, start = pos, end = pos) %>% 
        select(Chr, start, end)
        
    fwrite(bcftools_input, '${OUTPUT_DIR}/bcftools_input.txt', sep = '\t', col.names = FALSE)
    "
fi

#============================================================================
# STEP 2: EXTRACT VARIANTS
#============================================================================

if [ "$START_STEP" -le 2 ]; then
    echo "============================================================================"
    echo "STEP 2: Extracting variants"
    echo "============================================================================"

    VCF_LIST="${TEMP_DIR}/vcf_filelist.txt"
    rm -f "$VCF_LIST"
    touch "$VCF_LIST"

    for chr in {1..22}; do
        echo "Processing chromosome $chr..."
        VCF_FILE=$(echo "$VCF_PATH_TEMPLATE" | sed "s/{CHR}/$chr/")
        OUTPUT_VCF="${TEMP_DIR}/$(echo "$VCF_OUTPUT_TEMPLATE" | sed "s/{CHR}/$chr/")"
        
        if [ ! -f "$VCF_FILE" ]; then echo "WARNING: VCF missing for chr $chr"; continue; fi
        
        $BCFTOOLS_PATH view -R "${OUTPUT_DIR}/bcftools_input.txt" -Oz -o "$OUTPUT_VCF" "$VCF_FILE"
        $BCFTOOLS_PATH index "$OUTPUT_VCF" --tbi
        
        num_variants=$($BCFTOOLS_PATH view -H "$OUTPUT_VCF" | wc -l)
        if [ "$num_variants" -gt 0 ]; then echo "$OUTPUT_VCF" >> "$VCF_LIST"; fi
    done

    echo "Combining VCFs..."
    $BCFTOOLS_PATH concat -f "$VCF_LIST" -Oz -o "$COMBINED_VCF"
    $BCFTOOLS_PATH index "$COMBINED_VCF" --tbi
    rm -f ${TEMP_DIR}/extracted_chr*.vcf.gz*
fi

#============================================================================
# STEP 3: CHECK MISSING
#============================================================================

if [ "$START_STEP" -le 3 ]; then
    echo "============================================================================"
    echo "STEP 3: Checking missing SNPs"
    echo "============================================================================"

    $BCFTOOLS_PATH query -f '%CHROM\t%POS\t%POS\n' "$COMBINED_VCF" > "${OUTPUT_DIR}/present_snps.txt"
    sort "${OUTPUT_DIR}/bcftools_input.txt" > "${OUTPUT_DIR}/expected_snps.txt"
    sort "${OUTPUT_DIR}/present_snps.txt" > "${OUTPUT_DIR}/sorted_present_snps.txt"
    comm -23 "${OUTPUT_DIR}/expected_snps.txt" "${OUTPUT_DIR}/sorted_present_snps.txt" > "${OUTPUT_DIR}/missing_snps.txt"

    MISSING_COUNT=$(wc -l < "${OUTPUT_DIR}/missing_snps.txt")
    echo "Missing SNPs: $MISSING_COUNT"
else
    # If skipping Step 3, we MUST assume missing_snps.txt exists and read it
    # otherwise Step 4 (Proxy Search) will silently fail or skip.
    if [ -f "${OUTPUT_DIR}/missing_snps.txt" ]; then
        MISSING_COUNT=$(wc -l < "${OUTPUT_DIR}/missing_snps.txt")
        echo "Skipping Step 3: Loaded $MISSING_COUNT missing SNPs from existing file."
    else
        echo "Skipping Step 3: No missing SNPs file found. Assuming 0 missing."
        MISSING_COUNT=0
    fi
fi

#============================================================================
# STEP 4: PROXY SEARCH
#============================================================================

if [ "$START_STEP" -le 4 ]; then

    if [ "$ENABLE_PROXY_SEARCH" = true ] && [ "$MISSING_COUNT" -gt 0 ]; then
        echo "============================================================================"
        echo "STEP 4: Running Proxy Search"
        echo "============================================================================"
        
        # PASS GENOME_BUILD TO R SCRIPT HERE
    # Prepare the VCF path pattern for R (escape slashes for sed)
        # This converts /path/to/chr{CHR}.vcf to /path/to/chr%s.vcf
    # Prepare the VCF path pattern (converts {CHR} to %s)
        VCF_PATTERN=$(echo "$VCF_PATH_TEMPLATE" | sed 's/{CHR}/%s/')
        
        # Configure the R script by replacing placeholders
        sed -e "s/ldlink_token <- \"YOUR_LDLINK_TOKEN_HERE\"/ldlink_token <- \"$LDLINK_TOKEN\"/" \
            -e "s|weights_file <- \"./pPS_files/pPS_weights_plink_format.tsv\"|weights_file <- \"${OUTPUT_DIR}/pPS_weights_plink_format.tsv\"|" \
            -e "s|output_dir <- \"./pPS_files\"|output_dir <- \"${OUTPUT_DIR}\"|" \
            -e "s|missing_snps_file <- \"./pPS_files/missing_snps.txt\"|missing_snps_file <- \"${OUTPUT_DIR}/missing_snps.txt\"|" \
            -e "s/pop = \"EUR\"/pop = \"$POPULATION\"/" \
            -e "s/r2d = \"r2\"/r2d = \"r2\"/" \
            -e "s/genome_build_input <- \"grch38\"/genome_build_input <- \"$GENOME_BUILD\"/" \
            -e "s|vcf_file_template <- \"VCF_FILE_TEMPLATE_PLACEHOLDER\"|vcf_file_template <- \"$VCF_PATTERN\"|" \
            "${SCRIPT_DIR}/update_weights_with_proxies.R" > "${TEMP_DIR}/update_weights_with_proxies_configured.R"
        
        Rscript "${TEMP_DIR}/update_weights_with_proxies_configured.R"
        
    echo "Re-extracting SNPs with updated proxy list..."
        
        # Reset file list
        rm -f "$VCF_LIST"
        touch "$VCF_LIST"
        
        # Loop through chromosomes again with the NEW input file
        for chr in {1..22}; do
            VCF_FILE=$(echo "$VCF_PATH_TEMPLATE" | sed "s/{CHR}/$chr/")
            OUTPUT_VCF="${TEMP_DIR}/extracted_proxy_chr${chr}.vcf.gz"
            
            if [ ! -f "$VCF_FILE" ]; then
                continue
            fi
            
            # Use updated_bcftools_input.txt which contains proxies
            $BCFTOOLS_PATH view -R "${OUTPUT_DIR}/updated_bcftools_input.txt" -Oz -o "$OUTPUT_VCF" "$VCF_FILE"
            # ADDED -f here
            $BCFTOOLS_PATH index -f "$OUTPUT_VCF" --tbi
            
            num_variants=$($BCFTOOLS_PATH view -H "$OUTPUT_VCF" | wc -l)
            if [ "$num_variants" -gt 0 ]; then
                echo "$OUTPUT_VCF" >> "$VCF_LIST"
                # Ensure IDs are standardized to CHR:POS:REF:ALT to help PLINK matching
                $BCFTOOLS_PATH annotate --set-id '%CHROM:%POS:%REF:%ALT' "$OUTPUT_VCF" -Oz -o "${OUTPUT_VCF}.tmp"
                mv "${OUTPUT_VCF}.tmp" "$OUTPUT_VCF"
                # ADDED -f here too
                $BCFTOOLS_PATH index -f "$OUTPUT_VCF" --tbi
            else
                rm -f "$OUTPUT_VCF" "${OUTPUT_VCF}.tbi"
            fi
        done
        
        # Re-concatenate
        echo "Creating final combined VCF with proxies..."
        $BCFTOOLS_PATH concat -f "$VCF_LIST" -Oz -o "$COMBINED_VCF"
        $BCFTOOLS_PATH index -f "$COMBINED_VCF" --tbi
        
        WEIGHT_FILE="${OUTPUT_DIR}/updated_cluster_weights.txt"
    else
        WEIGHT_FILE="${OUTPUT_DIR}/pPS_weights_plink_format.tsv"
    fi
fi

# --- HANDLE SKIPPING STEP 4 ---
# If we skipped step 4, we need to determine which weights file to use
if [ -z "$WEIGHT_FILE" ]; then
    if [ -f "${OUTPUT_DIR}/updated_cluster_weights.txt" ]; then
        echo "Found updated weights with proxies. Using that."
        WEIGHT_FILE="${OUTPUT_DIR}/updated_cluster_weights.txt"
    else
        echo "Using original weights (no proxies found or step skipped)."
        WEIGHT_FILE="${OUTPUT_DIR}/pPS_weights_plink_format.tsv"
    fi
fi

#============================================================================
# STEP 5: ALLELE ALIGNMENT
#============================================================================


echo "============================================================================"
echo "STEP 5: Allele Alignment"
echo "============================================================================"

$BCFTOOLS_PATH query -f '%CHROM:%POS:%REF:%ALT\n' "$COMBINED_VCF" > "${OUTPUT_DIR}/vcf_variant_ids.txt"
Rscript "${SCRIPT_DIR}/update_weights_allele_alignment.R" \
    "$WEIGHT_FILE" \
    "${OUTPUT_DIR}/vcf_variant_ids.txt" \
    "${OUTPUT_DIR}/cluster_names.txt" \
    "${OUTPUT_DIR}/updated_cluster_weights_aligned.txt"

WEIGHT_FILE="${OUTPUT_DIR}/updated_cluster_weights_aligned.txt"

#============================================================================
# STEP 6: CALCULATE SCORES
#============================================================================

echo "============================================================================"
echo "STEP 6: Calculating Scores"
echo "============================================================================"

num_cols=$(awk -F'\t' 'NR==1 {print NF-2}' "$WEIGHT_FILE")
score_cols=""
for ((i=3; i<=$((num_cols+2)); i++)); do score_cols+="$i,"; done
score_cols="${score_cols%,}"

$PLINK2_PATH --vcf "$COMBINED_VCF" \
    --score "$WEIGHT_FILE" 1 2 header ignore-dup-ids list-variants cols=scoresums \
    --score-col-nums "$score_cols" \
    --out "${OUTPUT_DIR}/prs_scores"

#============================================================================
# FINISH
#============================================================================
echo "Pipeline completed successfully."