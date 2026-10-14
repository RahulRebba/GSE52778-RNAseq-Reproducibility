#!/usr/bin/env python3

from pathlib import Path
import sqlite3
import pandas as pd


# ============================================================
# PROJECT PATHS
# ============================================================

ROOT = Path(__file__).resolve().parents[1]

DB_PATH = ROOT / "results" / "generated" / "gse52778_project.db"
SCHEMA_PATH = ROOT / "sql" / "schema_sqlite.sql"
VIEWS_PATH = ROOT / "sql" / "02_create_views.sql"

META_PATH = ROOT / "data" / "metadata" / "sample_metadata_cleaned.csv"

NORM_PATH = ROOT / "results" / "generated" / "normalization_qc_manual.csv"
PCA_PATH = ROOT / "results" / "generated" / "pca_coordinates_manual.csv"

PRIMARY_DE_PATH = (
    ROOT
    / "results"
    / "generated"
    / "deseq2_results_annotated_manual.csv"
)

UNPAIRED_DE_PATH = (
    ROOT
    / "results"
    / "generated"
    / "deseq2_results_unpaired_manual.csv"
)

GO_UP_PATH = (
    ROOT
    / "results"
    / "generated"
    / "go_bp_upregulated_simplified.csv"
)

GO_DOWN_PATH = (
    ROOT
    / "results"
    / "generated"
    / "go_bp_downregulated_simplified.csv"
)

SUMMARY_PATH = (
    ROOT
    / "results"
    / "generated"
    / "analysis_summary_manual.csv"
)


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def sql_value(value):
    """
    Convert pandas / NumPy values into SQLite-safe Python values.
    Missing values become NULL.
    """

    if pd.isna(value):
        return None

    if hasattr(value, "item"):
        return value.item()

    return value


def boolean_to_integer(value):
    """
    Convert TRUE/FALSE-style values to SQLite 1/0.
    """

    text = str(value).strip().lower()

    if text in {"1", "true", "yes", "y"}:
        return 1

    if text in {"0", "false", "no", "n"}:
        return 0

    raise ValueError(
        f"Cannot interpret boolean value: {value}"
    )


def split_summary_value(value):
    """
    Store numerical summary values as REAL.
    Store non-numerical values as TEXT.
    """

    if pd.isna(value):
        return None, None

    text = str(value).strip()

    try:
        return float(text), None
    except ValueError:
        return None, text


# ============================================================
# VERIFY REQUIRED INPUT FILES
# ============================================================

required_files = [
    SCHEMA_PATH,
    VIEWS_PATH,
    META_PATH,
    NORM_PATH,
    PCA_PATH,
    PRIMARY_DE_PATH,
    UNPAIRED_DE_PATH,
    GO_UP_PATH,
    GO_DOWN_PATH,
    SUMMARY_PATH,
]

missing_files = [
    path
    for path in required_files
    if not path.exists()
]

if missing_files:
    raise FileNotFoundError(
        "Missing required files:\n"
        + "\n".join(str(path) for path in missing_files)
    )


# ============================================================
# READ SOURCE DATA
# ============================================================

metadata = pd.read_csv(META_PATH)

normalization = pd.read_csv(NORM_PATH)

pca = pd.read_csv(PCA_PATH)

primary_de = pd.read_csv(PRIMARY_DE_PATH)

unpaired_de = pd.read_csv(UNPAIRED_DE_PATH)

go_up = pd.read_csv(GO_UP_PATH)

go_down = pd.read_csv(GO_DOWN_PATH)

analysis_summary = pd.read_csv(SUMMARY_PATH)


# ============================================================
# VALIDATE SOURCE DATA
# ============================================================

if len(metadata) != 8:
    raise ValueError(
        f"Expected 8 samples, found {len(metadata)}"
    )

if metadata["sra_run"].duplicated().any():
    raise ValueError("Duplicate SRA sample IDs found.")

if metadata["geo_sample"].duplicated().any():
    raise ValueError("Duplicate GEO sample IDs found.")

if len(normalization) != 8:
    raise ValueError(
        f"Expected 8 normalization rows, found {len(normalization)}"
    )

if len(pca) != 8:
    raise ValueError(
        f"Expected 8 PCA rows, found {len(pca)}"
    )

if primary_de["gene_id"].duplicated().any():
    raise ValueError(
        "Duplicate gene IDs found in primary DE results."
    )

if unpaired_de["gene_id"].duplicated().any():
    raise ValueError(
        "Duplicate gene IDs found in unpaired DE results."
    )

primary_genes = set(primary_de["gene_id"])

unpaired_genes = set(unpaired_de["gene_id"])

if primary_genes != unpaired_genes:
    raise ValueError(
        "Primary and unpaired analyses contain different gene sets."
    )


# ============================================================
# CREATE / OPEN DATABASE
# ============================================================

DB_PATH.parent.mkdir(
    parents=True,
    exist_ok=True
)

connection = sqlite3.connect(DB_PATH)

connection.execute(
    "PRAGMA foreign_keys = ON;"
)

connection.executescript(
    SCHEMA_PATH.read_text(
        encoding="utf-8"
    )
)

connection.executescript(
    VIEWS_PATH.read_text(
        encoding="utf-8"
    )
)

# ============================================================
# CLEAR EXISTING DATA
# ============================================================

connection.execute(
    "DELETE FROM go_enrichment"
)

connection.execute(
    "DELETE FROM analysis_result"
)

connection.execute(
    "DELETE FROM analysis_run"
)

connection.execute(
    "DELETE FROM sample_qc"
)

connection.execute(
    "DELETE FROM sample"
)

connection.execute(
    "DELETE FROM gene"
)

connection.execute(
    "DELETE FROM analysis_summary"
)

connection.execute(
    "DELETE FROM study"
)

connection.execute(
    """
    DELETE FROM sqlite_sequence
    WHERE name IN (
        'analysis_run',
        'go_enrichment',
        'analysis_summary'
    )
    """
)


# ============================================================
# INSERT STUDY
# ============================================================

connection.execute(
    """
    INSERT INTO study (
        study_id,
        title,
        organism,
        platform,
        source_url
    )
    VALUES (?, ?, ?, ?, ?)
    """,
    (
        "GSE52778",
        "Human Airway Smooth Muscle Transcriptome Changes "
        "in Response to Asthma Medications",
        "Homo sapiens",
        "Illumina HiSeq 2000",
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778",
    ),
)


# ============================================================
# INSERT SAMPLES
# ============================================================

sample_rows = []

for _, row in metadata.iterrows():

    sample_rows.append(
        (
            row["sra_run"],
            row["study"],
            row["geo_sample"],
            row["donor"],
            row["treatment"],
            row["dex"],
            boolean_to_integer(
                row["include_primary_analysis"]
            ),
        )
    )

connection.executemany(
    """
    INSERT INTO sample (
        sample_id,
        study_id,
        geo_sample,
        donor,
        treatment,
        dex,
        include_primary_analysis
    )
    VALUES (?, ?, ?, ?, ?, ?, ?)
    """,
    sample_rows,
)


# ============================================================
# BUILD SAMPLE QC TABLE
# ============================================================

pca_qc = (
    pca[
        [
            "name",
            "PC1",
            "PC2",
        ]
    ]
    .rename(
        columns={
            "name": "sample_id"
        }
    )
)

sample_qc = normalization.merge(
    pca_qc,
    on="sample_id",
    how="inner",
    validate="one_to_one",
)

if len(sample_qc) != 8:
    raise ValueError(
        "Sample QC merge did not produce 8 rows."
    )

metadata_ids = set(
    metadata["sra_run"]
)

qc_ids = set(
    sample_qc["sample_id"]
)

if metadata_ids != qc_ids:
    raise ValueError(
        "Sample IDs do not match between metadata and QC files."
    )

sample_qc_rows = []

for _, row in sample_qc.iterrows():

    sample_qc_rows.append(
        (
            row["sample_id"],
            "PASS",
            sql_value(
                row["raw_library_millions"]
            ),
            sql_value(
                row["size_factor"]
            ),
            sql_value(
                row["PC1"]
            ),
            sql_value(
                row["PC2"]
            ),
            "NOT_FLAGGED",
            "Retained in primary analysis after QC review",
        )
    )

connection.executemany(
    """
    INSERT INTO sample_qc (
        sample_id,
        metadata_status,
        raw_library_millions,
        size_factor,
        pc1,
        pc2,
        outlier_flag,
        note
    )
    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
    """,
    sample_qc_rows,
)


# ============================================================
# INSERT GENE ANNOTATION
# ============================================================

gene_rows = []

for _, row in primary_de.iterrows():

    gene_rows.append(
        (
            row["gene_id"],
            sql_value(
                row["symbol"]
            ),
            sql_value(
                row["gene_name"]
            ),
            sql_value(
                row["gene_biotype"]
            ),
        )
    )

connection.executemany(
    """
    INSERT INTO gene (
        ensembl_gene_id,
        gene_symbol,
        gene_name,
        gene_biotype
    )
    VALUES (?, ?, ?, ?)
    """,
    gene_rows,
)


# ============================================================
# INSERT ANALYSIS RUNS
# ============================================================

cursor = connection.execute(
    """
    INSERT INTO analysis_run (
        study_id,
        analysis_name,
        analysis_type,
        is_primary,
        design_formula,
        filter_rule,
        alpha,
        effect_size_threshold,
        result_provenance
    )
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
    """,
    (
        "GSE52778",
        "Dexamethasone vs untreated - paired donor model",
        "differential_expression",
        1,
        "~ cell + dex",
        "count >= 10 in at least 4 samples",
        0.05,
        1.0,
        "Locally generated DESeq2 primary analysis",
    ),
)

primary_analysis_id = cursor.lastrowid


cursor = connection.execute(
    """
    INSERT INTO analysis_run (
        study_id,
        analysis_name,
        analysis_type,
        is_primary,
        design_formula,
        filter_rule,
        alpha,
        effect_size_threshold,
        result_provenance
    )
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
    """,
    (
        "GSE52778",
        "Dexamethasone vs untreated - unpaired sensitivity model",
        "sensitivity",
        0,
        "~ dex",
        "count >= 10 in at least 4 samples",
        0.05,
        1.0,
        "Locally generated DESeq2 sensitivity analysis",
    ),
)

unpaired_analysis_id = cursor.lastrowid


# ============================================================
# FUNCTION TO INSERT DE RESULTS
# ============================================================

def insert_de_results(
    dataframe,
    analysis_id,
):
    rows = []

    for _, row in dataframe.iterrows():

        rows.append(
            (
                analysis_id,
                row["gene_id"],
                sql_value(
                    row["baseMean"]
                ),
                sql_value(
                    row["log2FoldChange"]
                ),
                sql_value(
                    row["lfcSE"]
                ),
                sql_value(
                    row["stat"]
                ),
                sql_value(
                    row["pvalue"]
                ),
                sql_value(
                    row["padj"]
                ),
                sql_value(
                    row["result_class"]
                ),
            )
        )

    connection.executemany(
        """
        INSERT INTO analysis_result (
            analysis_id,
            ensembl_gene_id,
            baseMean,
            log2FoldChange,
            lfcSE,
            stat,
            pvalue,
            padj,
            result_class
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        rows,
    )


# ============================================================
# INSERT PRIMARY + SENSITIVITY RESULTS
# ============================================================

insert_de_results(
    primary_de,
    primary_analysis_id,
)

insert_de_results(
    unpaired_de,
    unpaired_analysis_id,
)


# ============================================================
# INSERT GO ENRICHMENT RESULTS
# ============================================================

def insert_go_results(
    dataframe,
    direction,
):
    rows = []

    for _, row in dataframe.iterrows():

        rows.append(
            (
                primary_analysis_id,
                direction,
                "BP",
                row["ID"],
                row["Description"],
                row["GeneRatio"],
                row["BgRatio"],
                sql_value(
                    row["pvalue"]
                ),
                sql_value(
                    row["p.adjust"]
                ),
                sql_value(
                    row["qvalue"]
                ),
                int(
                    row["Count"]
                ),
                sql_value(
                    row["geneID"]
                ),
            )
        )

    connection.executemany(
        """
        INSERT INTO go_enrichment (
            analysis_id,
            direction,
            ontology,
            go_id,
            description,
            gene_ratio,
            background_ratio,
            pvalue,
            adjusted_p_value,
            q_value,
            gene_count,
            gene_ids
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        rows,
    )


insert_go_results(
    go_up,
    "upregulated",
)

insert_go_results(
    go_down,
    "downregulated",
)


# ============================================================
# INSERT ANALYSIS SUMMARY
# ============================================================

summary_rows = []

for _, row in analysis_summary.iterrows():

    numeric_value, text_value = split_summary_value(
        row["value"]
    )

    summary_rows.append(
        (
            "GSE52778",
            row["metric"],
            numeric_value,
            text_value,
            None,
        )
    )

connection.executemany(
    """
    INSERT INTO analysis_summary (
        study_id,
        metric,
        numeric_value,
        text_value,
        unit
    )
    VALUES (?, ?, ?, ?, ?)
    """,
    summary_rows,
)


# ============================================================
# COMMIT DATABASE
# ============================================================

connection.commit()


# ============================================================
# VALIDATION QUERIES
# ============================================================

table_counts = {}

for table in [
    "study",
    "sample",
    "sample_qc",
    "gene",
    "analysis_run",
    "analysis_result",
    "go_enrichment",
    "analysis_summary",
]:

    table_counts[table] = connection.execute(
        f"SELECT COUNT(*) FROM {table}"
    ).fetchone()[0]


primary_significant = connection.execute(
    """
    SELECT COUNT(*)
    FROM analysis_result
    WHERE analysis_id = ?
      AND padj < 0.05
    """,
    (
        primary_analysis_id,
    ),
).fetchone()[0]


primary_twofold = connection.execute(
    """
    SELECT COUNT(*)
    FROM analysis_result
    WHERE analysis_id = ?
      AND padj < 0.05
      AND ABS(log2FoldChange) >= 1
    """,
    (
        primary_analysis_id,
    ),
).fetchone()[0]


unpaired_significant = connection.execute(
    """
    SELECT COUNT(*)
    FROM analysis_result
    WHERE analysis_id = ?
      AND padj < 0.05
    """,
    (
        unpaired_analysis_id,
    ),
).fetchone()[0]


foreign_key_errors = connection.execute(
    "PRAGMA foreign_key_check"
).fetchall()


print(
    f"\nDatabase: {DB_PATH}"
)

print(
    "\nTable row counts:"
)

for table, count in table_counts.items():
    print(
        f"  {table}: {count}"
    )

print(
    f"\nPrimary significant genes: {primary_significant}"
)

print(
    f"Primary >=2-fold significant genes: {primary_twofold}"
)

print(
    f"Unpaired significant genes: {unpaired_significant}"
)

print(
    f"Foreign-key errors: {len(foreign_key_errors)}"
)


connection.close()