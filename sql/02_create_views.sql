-- ============================================================
-- GSE52778 Reporting Views
-- Reusable analytical layer for SQL and Power BI
-- ============================================================


-- ------------------------------------------------------------
-- 1. SAMPLE QC VIEW
-- Combines sample metadata with QC metrics
-- ------------------------------------------------------------

DROP VIEW IF EXISTS vw_sample_qc;

CREATE VIEW vw_sample_qc AS
SELECT
    s.sample_id,
    s.geo_sample,
    s.donor,
    s.treatment,
    s.dex,
    q.metadata_status,
    q.raw_library_millions,
    q.size_factor,
    q.pc1,
    q.pc2,
    q.outlier_flag
FROM sample AS s
INNER JOIN sample_qc AS q
    ON s.sample_id = q.sample_id;


-- ------------------------------------------------------------
-- 2. PRIMARY DIFFERENTIAL EXPRESSION VIEW
-- Combines primary DESeq2 results with gene annotation
-- ------------------------------------------------------------

DROP VIEW IF EXISTS vw_primary_de_results;

CREATE VIEW vw_primary_de_results AS
SELECT
    g.ensembl_gene_id,
    g.gene_symbol,
    g.gene_name,
    g.gene_biotype,
    ar.baseMean,
    ar.log2FoldChange,
    ar.lfcSE,
    ar.stat,
    ar.pvalue,
    ar.padj,
    ar.result_class
FROM analysis_result AS ar
INNER JOIN gene AS g
    ON ar.ensembl_gene_id = g.ensembl_gene_id
WHERE ar.analysis_id = 1;


-- ------------------------------------------------------------
-- 3. SENSITIVITY COMPARISON VIEW
-- Compares paired and unpaired models gene by gene
-- ------------------------------------------------------------

DROP VIEW IF EXISTS vw_sensitivity_comparison;

CREATE VIEW vw_sensitivity_comparison AS
SELECT
    g.gene_symbol,
    p.ensembl_gene_id,
    p.log2FoldChange AS paired_log2fc,
    u.log2FoldChange AS unpaired_log2fc,
    p.padj AS paired_padj,
    u.padj AS unpaired_padj,
    p.result_class AS paired_result_class,
    u.result_class AS unpaired_result_class,

    CASE
        WHEN p.padj < 0.05
             AND u.padj < 0.05
            THEN 'Both'

        WHEN p.padj < 0.05
             AND (
                 u.padj >= 0.05
                 OR u.padj IS NULL
             )
            THEN 'Paired only'

        WHEN (
                 p.padj >= 0.05
                 OR p.padj IS NULL
             )
             AND u.padj < 0.05
            THEN 'Unpaired only'

        ELSE 'Neither'
    END AS overlap_class

FROM analysis_result AS p

INNER JOIN analysis_result AS u
    ON p.ensembl_gene_id = u.ensembl_gene_id

INNER JOIN gene AS g
    ON p.ensembl_gene_id = g.ensembl_gene_id

WHERE p.analysis_id = 1
  AND u.analysis_id = 2;


-- ------------------------------------------------------------
-- 4. GO ENRICHMENT VIEW
-- Simplified GO Biological Process results
-- from the primary analysis
-- ------------------------------------------------------------

DROP VIEW IF EXISTS vw_go_enrichment;

CREATE VIEW vw_go_enrichment AS
SELECT
    ge.direction,
    ge.ontology,
    ge.go_id,
    ge.description,
    ge.gene_ratio,
    ge.background_ratio,
    ge.adjusted_p_value,
    ge.q_value,
    ge.gene_count,
    ge.gene_ids
FROM go_enrichment AS ge

INNER JOIN analysis_run AS ar
    ON ge.analysis_id = ar.analysis_id

WHERE ar.is_primary = 1;