# PPS Pipeline: bNMF Cluster-Based Polygenic Scoring

![Version](https://img.shields.io/badge/version-2.1-blue)

This pipeline automates the generation of partitioned Polygenic Scores (pPS) from bNMF clustering results. It is designed to be flexible, accepting either standardized flat files or raw bNMF outputs, and runs on both Broad Institute infrastructure and standard SLURM clusters.

## 🚀 Quick Start

1.  **Prepare Inputs**:
    * **Option A (Standard)**: Place your formatted `cluster_weights.tsv` in your weights directory.
    * **Option B (bNMF Users)**: Use the helper script to convert your Excel files (see [Input Data](#input-data-requirements)).

2.  **Configure**:
    ```bash
    nano cluster_pPS_config.sh
    # Set VCF paths, LDlink token, and Environment variables
    ```

3.  **Run**:
    ```bash
    ./cluster_pPS_pipeline_generalized.sh
    ```

---

## 🛠 Prerequisites

### Required Software
* **PLINK2** & **Bcftools**: Must be in your `$PATH` or loaded via modules.
* **R (≥ 4.3)**: Required packages:
    ```r
    # CRAN
    install.packages(c("tidyverse", "data.table", "readxl", "LDlinkR", "ggplot2", "gridExtra", "reshape2"))
    
    # Bioconductor
    if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(c("rtracklayer", "GenomicRanges"))
    ```

### Required Data
* **VCF Files**: One per chromosome (hg38), bgzip compressed.
* **Chain File**: `hg19ToHg38.over.chain` (Required for automated LiftOver).

---

## 📁 Input Data Requirements

Place the following files in your `CLUSTER_WEIGHTS_DIR`.

### Option 1: Standardized Flat File (Recommended)
If you have your own weights, provide a single tab-separated file named `cluster_weights.tsv`.

**Format:**
```text
CHR  POS    REF  ALT  Risk_Allele  BETA  Cluster_BetaCell  Cluster_Lipodystrophy
1    12345  G    A    A            0.05  0.05              0.00
2    67890  T    C    T            0.12  0.00              0.12
```

---

## 📋 Changelog

### v2.1 (2026-02-25)
- Added `bcftools norm -m -any` to Steps 2 and 4 to correctly handle multi-allelic sites in WGS VCFs — variants queried by position are split into biallelic records before allele matching
- Added `cluster_names.txt` support for flexible cluster naming
- Fixed HPC SIGILL crash in Step 7 (`generate_snp_level_maf.R`) on nodes with restricted CPU instruction sets
- Replaced `tidyr::separate()` with `data.table::tstrsplit()` across all R scripts to avoid SIGILL on HPC nodes

### v2.0 (2026-01-30)
- Initial public release of the generalized bNMF pPS pipeline
- Support for Broad Institute and SLURM cluster environments
- Automated LiftOver (hg19 → hg38), LDlink proxy search, allele alignment, and PLINK2 scoring
