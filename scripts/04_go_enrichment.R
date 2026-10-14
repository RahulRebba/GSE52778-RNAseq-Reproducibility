# ============================================================
# GSE52778 GO Biological Process Enrichment
# ============================================================
# Purpose:
# Perform GO Biological Process enrichment separately for
# strongly upregulated and strongly downregulated genes.
#
# Input:
# results/generated/deseq2_results_annotated_manual.csv
#
# Definition of genes used for enrichment:
# adjusted p-value < 0.05 AND |log2FoldChange| >= 1
#
# Background universe:
# all genes tested in the primary DESeq2 analysis
# ============================================================


suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
})


# ------------------------------------------------------------
# 1. File paths
# ------------------------------------------------------------

input_file <- "results/generated/deseq2_results_annotated_manual.csv"
output_dir <- "results/generated"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# ------------------------------------------------------------
# 2. Read the primary annotated DESeq2 results
# ------------------------------------------------------------

de <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Genes in annotated DESeq2 table:", nrow(de), "\n")


# ------------------------------------------------------------
# 3. Define enrichment gene sets
# ------------------------------------------------------------
# The project definition of a >=2-fold significant change is:
#
# padj < 0.05
# AND
# |log2FoldChange| >= 1
#
# We use the numerical thresholds directly rather than relying
# on the text stored in result_class.

up_genes <- de$gene_id[
  !is.na(de$padj) &
  de$padj < 0.05 &
  de$log2FoldChange >= 1
]

down_genes <- de$gene_id[
  !is.na(de$padj) &
  de$padj < 0.05 &
  de$log2FoldChange <= -1
]

universe_genes <- de$gene_id

cat(">=2-fold upregulated genes:", length(up_genes), "\n")
cat(">=2-fold downregulated genes:", length(down_genes), "\n")
cat("Background genes:", length(universe_genes), "\n")


# ------------------------------------------------------------
# 4. Remove possible Ensembl version suffixes
# ------------------------------------------------------------
# Example:
# ENSG00000123456.7 -> ENSG00000123456

clean_ensembl <- function(x) {
  sub("\\..*$", "", x)
}

up_genes <- clean_ensembl(up_genes)
down_genes <- clean_ensembl(down_genes)
universe_genes <- clean_ensembl(universe_genes)


# ------------------------------------------------------------
# 5. Convert Ensembl IDs to Entrez IDs
# ------------------------------------------------------------

universe_map <- bitr(
  universe_genes,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

up_map <- bitr(
  up_genes,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

down_map <- bitr(
  down_genes,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

universe_entrez <- unique(universe_map$ENTREZID)
up_entrez <- unique(up_map$ENTREZID)
down_entrez <- unique(down_map$ENTREZID)

cat("Mapped background genes:", length(universe_entrez), "\n")
cat("Mapped upregulated genes:", length(up_entrez), "\n")
cat("Mapped downregulated genes:", length(down_entrez), "\n")


# ------------------------------------------------------------
# 6. GO Biological Process enrichment - upregulated genes
# ------------------------------------------------------------

go_up <- enrichGO(
  gene = up_entrez,
  universe = universe_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)


# ------------------------------------------------------------
# 7. GO Biological Process enrichment - downregulated genes
# ------------------------------------------------------------

go_down <- enrichGO(
  gene = down_entrez,
  universe = universe_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)


# ------------------------------------------------------------
# 8. Save full enrichment results
# ------------------------------------------------------------

go_up_df <- as.data.frame(go_up)
go_down_df <- as.data.frame(go_down)

write.csv(
  go_up_df,
  file.path(output_dir, "go_bp_upregulated_full.csv"),
  row.names = FALSE
)

write.csv(
  go_down_df,
  file.path(output_dir, "go_bp_downregulated_full.csv"),
  row.names = FALSE
)

cat("Full GO BP upregulated terms:", nrow(go_up_df), "\n")
cat("Full GO BP downregulated terms:", nrow(go_down_df), "\n")


# ------------------------------------------------------------
# 9. Reduce redundant GO terms
# ------------------------------------------------------------
# GO is hierarchical, so many significant terms describe
# closely related biological concepts.
#
# simplify() removes highly similar terms while retaining the
# term with the smallest adjusted p-value.

go_up_simplified <- simplify(
  go_up,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang"
)

go_down_simplified <- simplify(
  go_down,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang"
)


# ------------------------------------------------------------
# 10. Save simplified enrichment results
# ------------------------------------------------------------

go_up_simplified_df <- as.data.frame(go_up_simplified)
go_down_simplified_df <- as.data.frame(go_down_simplified)

write.csv(
  go_up_simplified_df,
  file.path(output_dir, "go_bp_upregulated_simplified.csv"),
  row.names = FALSE
)

write.csv(
  go_down_simplified_df,
  file.path(output_dir, "go_bp_downregulated_simplified.csv"),
  row.names = FALSE
)

cat(
  "Simplified GO BP upregulated terms:",
  nrow(go_up_simplified_df),
  "\n"
)

cat(
  "Simplified GO BP downregulated terms:",
  nrow(go_down_simplified_df),
  "\n"
)


# ------------------------------------------------------------
# 11. Create GO enrichment dot plots
# ------------------------------------------------------------

p_up <- dotplot(
  go_up_simplified,
  showCategory = 10
) +
  ggtitle("GO Biological Process - Upregulated genes") +
  theme_minimal(base_size = 11)

ggsave(
  filename = file.path(
    output_dir,
    "09_go_bp_upregulated_dotplot.png"
  ),
  plot = p_up,
  width = 10,
  height = 7,
  dpi = 300
)


p_down <- dotplot(
  go_down_simplified,
  showCategory = 10
) +
  ggtitle("GO Biological Process - Downregulated genes") +
  theme_minimal(base_size = 11)

ggsave(
  filename = file.path(
    output_dir,
    "10_go_bp_downregulated_dotplot.png"
  ),
  plot = p_down,
  width = 10,
  height = 7,
  dpi = 300
)


cat("\nGO enrichment analysis completed successfully.\n")