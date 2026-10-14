$ErrorActionPreference = "Stop"

Write-Host "1/5 Validating metadata..."
python scripts/01_validate_metadata.py

Write-Host "2/5 Running primary DESeq2 analysis..."
Rscript scripts/02_airway_deseq2.R

Write-Host "3/5 Running paired-vs-unpaired sensitivity analysis..."
Rscript scripts/03_sensitivity_analysis.R

Write-Host "4/5 Running GO Biological Process enrichment..."
Rscript scripts/04_go_enrichment.R

Write-Host "5/5 Building SQLite reporting database..."
python scripts/05_load_sqlite.py

Write-Host ""
Write-Host "Pipeline finished successfully."
Write-Host "Generated outputs: results/generated/"
Write-Host "Power BI dashboard: powerbi/GSE52778_PowerBI_Dashboard_Final.pbix"