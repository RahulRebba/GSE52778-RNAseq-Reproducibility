nextflow.enable.dsl=2

/*
Prototype Nextflow workflow for the GSE52778 project.

The validated reproduction route is documented in README.md and
docs/REPRODUCE_ON_WINDOWS.md.

This workflow is retained as orchestration scaffolding and has not
been validated as the final production execution path.
*/

process VALIDATE_METADATA {
    output:
    path "metadata_validation_report.txt"

    script:
    """
    python scripts/01_validate_metadata.py
    cp results/generated/metadata_validation_report.txt .
    """
}

process DESEQ2_ANALYSIS {
    output:
    path "de_done.txt"

    script:
    """
    Rscript scripts/02_airway_deseq2.R
    echo done > de_done.txt
    """
}

process SENSITIVITY_ANALYSIS {
    output:
    path "sensitivity_done.txt"

    script:
    """
    Rscript scripts/03_sensitivity_analysis.R
    echo done > sensitivity_done.txt
    """
}

process GO_ENRICHMENT {
    output:
    path "go_done.txt"

    script:
    """
    Rscript scripts/04_go_enrichment.R
    echo done > go_done.txt
    """
}

process BUILD_SQLITE {
    output:
    path "sql_done.txt"

    script:
    """
    python scripts/05_load_sqlite.py
    echo done > sql_done.txt
    """
}

workflow {
    VALIDATE_METADATA()
    DESEQ2_ANALYSIS()
    SENSITIVITY_ANALYSIS()
    GO_ENRICHMENT()
    BUILD_SQLITE()
}