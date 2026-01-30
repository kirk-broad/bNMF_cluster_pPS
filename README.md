This revised **README.md** reflects your unified pipeline, where Step 1 now handles the conversion from raw bNMF Excel files and performs LiftOver. It is designed to be clear for both your future self and external collaborators at the Broad.

---

# PPS Pipeline: bNMF Cluster-Based Polygenic Scoring

This pipeline automates the generation of partitioned Polygenic Scores (pPS) directly from bNMF clustering results. It is designed for use on Broad Institute infrastructure but is flexible enough for any SLURM-based HPC environment.

## 🌟 Key Features

* **Raw Data Input**: Directly accepts raw bNMF Excel weight files.
* **Automated LiftOver**: Converts hg19 coordinates to hg38 automatically during setup.
* **GWAS Alignment**: Ensures risk alleles in your clusters match the underlying GWAS summary statistics.
* **Proxy Search**: Integrated LDlink API support to replace missing SNPs with high-LD proxies ().
* **Broad-Ready**: Pre-configured to load Broad Institute modules (PLINK2, Bcftools, R-4.3).

---

## 🛠 Prerequisites

### Required Software

* **PLINK2**: For score calculation.
* **Bcftools & Tabix**: For VCF manipulation and indexing.
* **R ( 4.3)**: With the following packages:
* **CRAN**: `tidyverse`, `data.table`, `readxl`, `LDlinkR`, `ggplot2`, `gridExtra`, `reshape2`
* **Bioconductor**: `rtracklayer`, `GenomicRanges`



### Installation (R Dependencies)

```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("rtracklayer", "GenomicRanges"))
install.packages(c("tidyverse", "data.table", "readxl", "LDlinkR", "ggplot2", "gridExtra", "reshape2"))

```

---

## 📁 Required Input Files

Place these in your `CLUSTER_WEIGHTS_DIR`:

1. **bNMF Weights**: `sorted_cluster_weights_K[X]_rev.xlsx` (Raw bNMF output).
2. **GWAS Stats**: `alignment_GWAS_summStats.csv` (For allele alignment).
3. **Cluster Names**: `cluster_names.txt` (One name per line).
4. **Chain File**: `hg19ToHg38.over.chain` (Required for LiftOver to hg38).

---

## 🚀 Getting Started

### 1. Configuration

Edit `cluster_pPS_config.sh` to point to your data:

```bash
# Path to your cohort VCFs (use {CHR} as placeholder)
VCF_PATH_TEMPLATE="/path/to/genotypes/MGBB.chr{CHR}.vcf.gz"

# Your LDlink token for proxy searches
LDLINK_TOKEN="your_token_here"

# Set to true for Broad Institute users
USE_BROAD_MODULES=true

```

### 2. Execution

Run the full pipeline from your project root:

```bash
chmod +x cluster_pPS_pipeline_generalized.sh
./cluster_pPS_pipeline_generalized.sh

```

*Note: You can skip to a specific step by providing the step number as an argument (e.g., `./cluster_pPS_pipeline_generalized.sh 6` to just re-run scoring).*

---

## 📊 Pipeline Workflow

1. **Step 1 (Setup)**: `cluster_weights_setup_generalized.R` parses the Excel file, performs LiftOver, and aligns alleles.
2. **Step 2 (Extraction)**: `bcftools` extracts variants from cohort VCFs based on the new weights.
3. **Step 4 (Proxies)**: (Optional) If SNPs are missing, the pipeline queries LDlink and re-extracts proxies.
4. **Step 5 (Alignment)**: Performs a final check against the VCF's actual REF/ALT alleles.
5. **Step 6 (Scoring)**: Runs `PLINK2 --score` for every cluster simultaneously.
6. **Step 7 (QC)**: Generates MAF distributions and cluster correlation plots.

---

## 📂 Outputs

All results are stored in `pPS_results_[STUDY_NAME]/`:

* **`prs_scores.sscore`**: The final polygenic scores for each individual.
* **`pPS_summary_stats.txt`**: Mean, SD, and quantiles for each cluster.
* **`pPS_histograms.pdf`**: Distribution plots of your cluster scores.
* **`snp_level_maf.txt`**: Detailed table of Risk Allele Frequencies (RAF) in your cohort.

---

## 🆘 Troubleshooting

* **Module Load Errors**: Ensure `source /broad/software/scripts/useuse` is accessible on your node.
* **LiftOver Failures**: Verify your `hg19ToHg38.over.chain` path is correct in the config.
* **Missing Scripts**: Ensure `generate_snp_level_maf.R` and `pPS_summary_plots.R` are in the same folder as the main script.

---

Would you like me to help you **write the `git commit` message** for this first push to your new repository?