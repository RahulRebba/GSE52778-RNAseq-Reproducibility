#!/usr/bin/env bash
set -euo pipefail

echo "1/5 Validating metadata..."
python scripts/01_validate_metadata.py

echo "2/5 Running primary DESeq2 analysis..."
Rscript scripts/02_airway_deseq2.R

echo "3/5 Running paired-vs-unpaired sensitivity analysis..."
Rscript scripts/03_sensitivity_analysis.R

echo "4/5 Running GO Biological Process enrichment..."
Rscript scripts/04_go_enrichment.R

echo "5/5 Building SQLite reporting database..."
python scripts/05_load_sqlite.py

echo
echo "Pipeline finished successfully."
echo "Generated outputs: results/generated/"
echo "Power BI dashboard: powerbi/GSE52778_PowerBI_Dashboard_Final.pbix"