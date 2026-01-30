# Quick Start Guide

## 1. Configure the Pipeline
```bash
# Edit configuration with your settings
nano cluster_pPS_config.sh
```

Key settings to change:
- `STUDY_NAME`: Your study name
- `VCF_PATH_TEMPLATE`: Path to your VCF files with {CHR} placeholder
- `CLUSTER_WEIGHTS_DIR`: Directory with cluster weights from Kirk
- `LDLINK_TOKEN`: Your LDlink API token
- `POPULATION`: Your population (EUR, AFR, AMR, EAS, SAS)

## 2. Validate Setup
```bash
./validate_setup.sh
```

## 3. Run Pipeline
```bash
./cluster_pPS_pipeline_generalized.sh
```

## 4. Check Results
Results will be in `pPS_results_[STUDY_NAME]/`

For detailed instructions, see `README.md`.
