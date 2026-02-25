
#!/bin/bash
# cluster_pPS_config.sh
# Configuration file for cluster-based polygenic score pipeline

#============================================================================
# COLLABORATOR-SPECIFIC SETTINGS - MODIFY THESE FOR YOUR COHORT
#============================================================================

# Study Information
STUDY_NAME="T2D_CLUSTERS"
COHORT_DESCRIPTION="MY_COHORT"

# Paths to Input Data
# Use {CHR} as placeholder for chromosome number
VCF_PATH_TEMPLATE="/PATH/TO/VCF/FILE/FILE_NAME.chr{CHR}.vcf.gz"

# Cluster Weights Directory
# MUST CONTAIN: 
# 1. cluster_names.txt
# 2. sorted_cluster_weights_K*_rev.xlsx (The raw weights)
# 3. alignment_GWAS_summStats.csv
CLUSTER_WEIGHTS_DIR="/PATH/TO/CLUSTER_WEIGHTS"
INPUT_WEIGHTS_FILE="NatMed_cluster_weights.tsv"  # <-- UPDATE FILE NAME

# Analysis Parameters
GENOME_BUILD="hg38"                           # Options: "hg19" or "hg38"
CUTOFF=0.5                                    # Weight cutoff (set to 0 to keep all)
CHAIN_FILE="/PATH/TO/CHAIN_FILE/hg19ToHg38.over.chain"            # Required if GENOME_BUILD="hg38"
ENABLE_PROXY_SEARCH=false    # Change this to false if you want to skip searching
ENABLE_DOSAGE_EXTRACTION=true
ENABLE_SUMMARY_PLOTS=true

# Software Paths
PLINK2_PATH="/PATH/TO/PLINK2/plink2"
BCFTOOLS_PATH="bcftools"
TABIX_PATH="tabix"

# Environment Setup
USE_BROAD_MODULES=true
USE_SLURM_MODULES=false
CUSTOM_SOFTWARE_SETUP=""

# LDproxy Configuration  
LDLINK_TOKEN="YOUR_LDLINK_TOKEN_HERE"
POPULATION="EUR"
R2_THRESHOLD=0.7

# Output Settings
OUTPUT_DIR="./pPS_results_${STUDY_NAME}"
TEMP_DIR="./temp_${STUDY_NAME}"

#============================================================================
# EXPORTS (DO NOT MODIFY)
#============================================================================
export STUDY_NAME COHORT_DESCRIPTION VCF_PATH_TEMPLATE CLUSTER_WEIGHTS_DIR
export PLINK2_PATH BCFTOOLS_PATH TABIX_PATH
export LDLINK_TOKEN POPULATION R2_THRESHOLD GENOME_BUILD CUTOFF CHAIN_FILE
export OUTPUT_DIR TEMP_DIR
export USE_BROAD_MODULES USE_SLURM_MODULES CUSTOM_SOFTWARE_SETUP
export ENABLE_PROXY_SEARCH ENABLE_DOSAGE_EXTRACTION ENABLE_SUMMARY_PLOTS
export MAX_MISSING_RATE=0.1 MIN_SNP_COUNT=1
export VCF_OUTPUT_TEMPLATE="extracted_chr{CHR}.vcf.gz"
export FINAL_VCF_NAME="combined_${STUDY_NAME}.vcf.gz"

validate_config() {
    echo "Validating configuration..."
# Check if CHAIN_FILE exists directly, OR inside the weights dir
    if [ "$GENOME_BUILD" == "hg38" ]; then
        if [ ! -f "$CHAIN_FILE" ] && [ ! -f "${CLUSTER_WEIGHTS_DIR}/${CHAIN_FILE}" ]; then
            echo "ERROR: hg38 requested but chain file not found at: $CHAIN_FILE"
            echo "Checked both absolute path and inside CLUSTER_WEIGHTS_DIR"
            exit 1
        fi
    fi
}

print_config() {
    echo "============================================================================"
    echo "CONFIGURATION SUMMARY"
    echo "============================================================================"
    echo "Study Name:       $STUDY_NAME"
    echo "Cohort Desc:      $COHORT_DESCRIPTION"
    echo "Genome Build:     $GENOME_BUILD"
    echo "Chain File:       $CHAIN_FILE"
    echo "Weights Dir:      $CLUSTER_WEIGHTS_DIR"
    echo "Output Dir:       $OUTPUT_DIR"
    echo "VCF Template:     $VCF_PATH_TEMPLATE"
    echo "Cutoff:           $CUTOFF"
    echo "============================================================================"
}