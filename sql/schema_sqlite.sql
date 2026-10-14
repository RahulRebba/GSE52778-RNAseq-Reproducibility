-- ============================================================
-- GSE52778 Gene Expression Project
-- Final SQLite relational schema
-- ============================================================

PRAGMA foreign_keys = ON;


-- ------------------------------------------------------------
-- 1. STUDY
-- One row per public gene-expression study
-- ------------------------------------------------------------

CREATE TABLE IF NOT EXISTS study (
    study_id TEXT PRIMARY KEY,
    title TEXT NOT NULL,
    organism TEXT,
    platform TEXT,
    source_url TEXT
);


-- ------------------------------------------------------------
-- 2. SAMPLE
-- One row per biological RNA-seq sample
-- ------------------------------------------------------------

CREATE TABLE IF NOT EXISTS sample (
    sample_id TEXT PRIMARY KEY,
    study_id TEXT NOT NULL,
    geo_sample TEXT NOT NULL UNIQUE,
    donor TEXT NOT NULL,
    treatment TEXT NOT NULL,
    dex TEXT NOT NULL,
    include_primary_analysis INTEGER NOT NULL
        CHECK (include_primary_analysis IN (0, 1)),

    FOREIGN KEY (study_id)
        REFERENCES study(study_id)
        ON DELETE CASCADE
);


-- ------------------------------------------------------------
-- 3. SAMPLE QC
-- Sample-level QC and normalization information
-- ------------------------------------------------------------

CREATE TABLE IF NOT EXISTS sample_qc (
    sample_id TEXT PRIMARY KEY,
    metadata_status TEXT NOT NULL,
    raw_library_millions REAL,
    size_factor REAL,
    pc1 REAL,
    pc2 REAL,
    outlier_flag TEXT,
    note TEXT,

    FOREIGN KEY (sample_id)
        REFERENCES sample(sample_id)
        ON DELETE CASCADE
);


-- ------------------------------------------------------------
-- 4. GENE
-- Gene annotation stored once, independent of analysis
-- ------------------------------------------------------------

CREATE TABLE IF NOT EXISTS gene (
    ensembl_gene_id TEXT PRIMARY KEY,
    gene_symbol TEXT,
    gene_name TEXT,
    gene_biotype TEXT
);


-- ------------------------------------------------------------
-- 5. ANALYSIS RUN
-- Describes each statistical analysis/model
-- ------------------------------------------------------------

CREATE TABLE IF NOT EXISTS analysis_run (
    analysis_id INTEGER PRIMARY KEY AUTOINCREMENT,
    study_id TEXT NOT NULL,
    analysis_name TEXT NOT NULL,
    analysis_type TEXT NOT NULL,
    is_primary INTEGER NOT NULL
        CHECK (is_primary IN (0, 1)),
    design_formula TEXT,
    filter_rule TEXT,
    alpha REAL,
    effect_size_threshold REAL,
    result_provenance TEXT,
    created_at TEXT DEFAULT CURRENT_TIMESTAMP,

    FOREIGN KEY (study_id)
        REFERENCES study(study_id)
        ON DELETE CASCADE
);


-- ------------------------------------------------------------
-- 6. ANALYSIS RESULT
-- Differential-expression statistics for each gene
-- within each analysis
-- ------------------------------------------------------------

CREATE TABLE IF NOT EXISTS analysis_result (
    analysis_id INTEGER NOT NULL,
    ensembl_gene_id TEXT NOT NULL,
    baseMean REAL,
    log2FoldChange REAL,
    lfcSE REAL,
    stat REAL,
    pvalue REAL,
    padj REAL,
    result_class TEXT,

    PRIMARY KEY (
        analysis_id,
        ensembl_gene_id
    ),

    FOREIGN KEY (analysis_id)
        REFERENCES analysis_run(analysis_id)
        ON DELETE CASCADE,

    FOREIGN KEY (ensembl_gene_id)
        REFERENCES gene(ensembl_gene_id)
        ON DELETE CASCADE
);


-- ------------------------------------------------------------
-- 7. GO ENRICHMENT
-- Simplified GO Biological Process enrichment results
-- ------------------------------------------------------------

CREATE TABLE IF NOT EXISTS go_enrichment (
    enrichment_id INTEGER PRIMARY KEY AUTOINCREMENT,
    analysis_id INTEGER NOT NULL,
    direction TEXT NOT NULL,
    ontology TEXT NOT NULL DEFAULT 'BP',
    go_id TEXT NOT NULL,
    description TEXT NOT NULL,
    gene_ratio TEXT,
    background_ratio TEXT,
    pvalue REAL,
    adjusted_p_value REAL,
    q_value REAL,
    gene_count INTEGER,
    gene_ids TEXT,

    FOREIGN KEY (analysis_id)
        REFERENCES analysis_run(analysis_id)
        ON DELETE CASCADE,

    UNIQUE (
        analysis_id,
        direction,
        go_id
    )
);


-- ------------------------------------------------------------
-- 8. ANALYSIS SUMMARY
-- Headline project-level metrics for reporting / Power BI
-- ------------------------------------------------------------

CREATE TABLE IF NOT EXISTS analysis_summary (
    summary_id INTEGER PRIMARY KEY AUTOINCREMENT,
    study_id TEXT NOT NULL,
    metric TEXT NOT NULL,
    numeric_value REAL,
    text_value TEXT,
    unit TEXT,

    FOREIGN KEY (study_id)
        REFERENCES study(study_id)
        ON DELETE CASCADE,

    UNIQUE (
        study_id,
        metric
    )
);


-- ============================================================
-- INDEXES
-- Improve filtering, joins and dashboard queries
-- ============================================================

CREATE INDEX IF NOT EXISTS idx_sample_study
ON sample(study_id);


CREATE INDEX IF NOT EXISTS idx_sample_donor
ON sample(donor);


CREATE INDEX IF NOT EXISTS idx_analysis_study
ON analysis_run(study_id);


CREATE INDEX IF NOT EXISTS idx_result_gene
ON analysis_result(ensembl_gene_id);


CREATE INDEX IF NOT EXISTS idx_result_padj
ON analysis_result(analysis_id, padj);


CREATE INDEX IF NOT EXISTS idx_result_class
ON analysis_result(analysis_id, result_class);


CREATE INDEX IF NOT EXISTS idx_go_direction
ON go_enrichment(
    analysis_id,
    direction,
    adjusted_p_value
);