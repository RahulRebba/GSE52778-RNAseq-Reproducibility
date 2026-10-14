-- ============================================================
-- GSE52778 Example Analytical Queries
-- Demonstrates validation, filtering, aggregation, joins,
-- CASE logic, sensitivity analysis and pathway interpretation
-- ============================================================


-- ------------------------------------------------------------
-- 1. VERIFY EXPERIMENTAL DESIGN
-- Confirm four donors, each with untreated and treated samples
-- ------------------------------------------------------------

SELECT
    donor,
    treatment,
    COUNT(*) AS sample_count
FROM sample
GROUP BY
    donor,
    treatment
ORDER BY
    donor,
    treatment;


-- ------------------------------------------------------------
-- 2. SAMPLE-LEVEL QC
-- View sample metadata together with normalization/QC metrics
-- ------------------------------------------------------------

SELECT
    sample_id,
    donor,
    treatment,
    raw_library_millions,
    size_factor,
    pc1,
    pc2,
    outlier_flag
FROM vw_sample_qc
ORDER BY
    donor,
    treatment;


-- ------------------------------------------------------------
-- 3. PRIMARY DE RESULT-CLASS SUMMARY
-- Count genes in each significance/effect-size category
-- ------------------------------------------------------------

SELECT
    result_class,
    COUNT(*) AS gene_count
FROM vw_primary_de_results
GROUP BY
    result_class
ORDER BY
    gene_count DESC;


-- ------------------------------------------------------------
-- 4. STRONGEST SIGNIFICANT DIFFERENTIAL-EXPRESSION RESULTS
-- Rank significant genes by absolute log2 fold-change
-- ------------------------------------------------------------

SELECT
    gene_symbol,
    ensembl_gene_id,
    log2FoldChange,
    padj,
    result_class
FROM vw_primary_de_results
WHERE
    padj < 0.05
ORDER BY
    ABS(log2FoldChange) DESC
LIMIT 20;


-- ------------------------------------------------------------
-- 5. CLASSIFY SIGNIFICANT GENES BY DIRECTION USING CASE
-- Demonstrates conditional logic
-- ------------------------------------------------------------

SELECT
    CASE
        WHEN log2FoldChange > 0 THEN 'Upregulated'
        WHEN log2FoldChange < 0 THEN 'Downregulated'
        ELSE 'No change'
    END AS direction,
    COUNT(*) AS gene_count
FROM vw_primary_de_results
WHERE
    padj < 0.05
GROUP BY
    CASE
        WHEN log2FoldChange > 0 THEN 'Upregulated'
        WHEN log2FoldChange < 0 THEN 'Downregulated'
        ELSE 'No change'
    END
ORDER BY
    gene_count DESC;


-- ------------------------------------------------------------
-- 6. GENES WITH >= 2-FOLD CHANGE AND PADJ < 0.05
-- log2FoldChange >= 1  = at least 2-fold increase
-- log2FoldChange <= -1 = at least 2-fold decrease
-- ------------------------------------------------------------

SELECT
    gene_symbol,
    ensembl_gene_id,
    log2FoldChange,
    padj,
    result_class
FROM vw_primary_de_results
WHERE
    padj < 0.05
    AND ABS(log2FoldChange) >= 1
ORDER BY
    ABS(log2FoldChange) DESC;


-- ------------------------------------------------------------
-- 7. PAIRED VS UNPAIRED SENSITIVITY OVERLAP
-- Compare significance stability between the two models
-- ------------------------------------------------------------

SELECT
    overlap_class,
    COUNT(*) AS gene_count
FROM vw_sensitivity_comparison
GROUP BY
    overlap_class
ORDER BY
    gene_count DESC;


-- ------------------------------------------------------------
-- 8. INSPECT GENES SIGNIFICANT ONLY IN THE PAIRED MODEL
-- Shows genes benefiting from donor adjustment
-- ------------------------------------------------------------

SELECT
    gene_symbol,
    ensembl_gene_id,
    paired_log2fc,
    unpaired_log2fc,
    paired_padj,
    unpaired_padj
FROM vw_sensitivity_comparison
WHERE
    overlap_class = 'Paired only'
ORDER BY
    paired_padj ASC
LIMIT 20;


-- ------------------------------------------------------------
-- 9. TOP GO BIOLOGICAL-PROCESS TERMS
-- Show strongest enrichment terms from the primary analysis
-- ------------------------------------------------------------

SELECT
    direction,
    go_id,
    description,
    gene_ratio,
    adjusted_p_value,
    gene_count
FROM vw_go_enrichment
ORDER BY
    adjusted_p_value ASC
LIMIT 20;


-- ------------------------------------------------------------
-- 10. GO TERM COUNTS BY DIRECTION
-- Confirm upregulated and downregulated enrichment sets
-- ------------------------------------------------------------

SELECT
    direction,
    COUNT(*) AS term_count
FROM vw_go_enrichment
GROUP BY
    direction
ORDER BY
    term_count DESC;


-- ------------------------------------------------------------
-- 11. TRACE ANALYSIS PROVENANCE
-- Identify which model generated each analysis set
-- ------------------------------------------------------------

SELECT
    st.study_id,
    ar.analysis_id,
    ar.analysis_name,
    ar.analysis_type,
    ar.is_primary,
    ar.design_formula,
    ar.filter_rule,
    ar.alpha,
    ar.effect_size_threshold,
    ar.result_provenance
FROM study AS st
INNER JOIN analysis_run AS ar
    ON st.study_id = ar.study_id
ORDER BY
    ar.analysis_id;


-- ------------------------------------------------------------
-- 12. TRACE SAMPLE RECORDS BACK TO THE STUDY
-- Demonstrates table relationships and data lineage
-- ------------------------------------------------------------

SELECT
    st.study_id,
    s.sample_id,
    s.geo_sample,
    s.donor,
    s.treatment,
    q.metadata_status,
    q.outlier_flag
FROM study AS st
INNER JOIN sample AS s
    ON s.study_id = st.study_id
LEFT JOIN sample_qc AS q
    ON q.sample_id = s.sample_id
ORDER BY
    s.donor,
    s.treatment;