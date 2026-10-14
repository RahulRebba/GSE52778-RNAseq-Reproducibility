# GSE52778 reproducible RNA-seq analysis
# Real data: Bioconductor "airway" package
# Biological question: effect of dexamethasone while controlling for donor/cell line.
# Primary design: ~ cell + dex
#
# IMPORTANT:
# DESeq2 requires raw/estimated counts, not FPKM/TPM.
# This script loads the count data packaged with Bioconductor airway.

options(stringsAsFactors = FALSE)

root <- normalizePath(file.path(getwd()), winslash="/", mustWork=FALSE)
if (basename(root) == "scripts") root <- dirname(root)
outdir <- file.path(root, "results", "generated")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
bioc_pkgs <- c("airway","DESeq2")
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, ask=FALSE, update=FALSE)
}
cran_pkgs <- c("ggplot2","pheatmap")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(airway)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

data("airway", package="airway")
se <- airway
counts <- assay(se)
meta <- as.data.frame(colData(se))

# ------------------------
# 1. DATA INTEGRITY CHECKS
# ------------------------
stopifnot(ncol(counts) == nrow(meta))
stopifnot(all(colnames(counts) == rownames(meta)))
stopifnot(all(counts >= 0))
stopifnot(all(abs(counts - round(counts)) < 1e-8))

meta$sample_id <- rownames(meta)
meta$dex <- relevel(factor(meta$dex), ref="untrt")
meta$cell <- factor(meta$cell)

write.csv(meta, file.path(outdir, "sample_metadata_from_airway.csv"), row.names=FALSE)
write.csv(data.frame(
  metric=c("genes","samples","missing_counts","duplicate_sample_ids"),
  value=c(nrow(counts),ncol(counts),sum(is.na(counts)),sum(duplicated(colnames(counts))))
), file.path(outdir, "data_integrity_summary.csv"), row.names=FALSE)

# Library sizes
lib <- data.frame(
  sample_id=colnames(counts),
  donor=meta$cell,
  treatment=meta$dex,
  library_size=colSums(counts)
)
write.csv(lib, file.path(outdir, "library_sizes.csv"), row.names=FALSE)

png(file.path(outdir, "01_library_sizes.png"), width=1600, height=900, res=150)
ggplot(lib, aes(x=sample_id, y=library_size/1e6, fill=treatment)) +
  geom_col() +
  labs(title="GSE52778 library size by sample",
       x="SRA run", y="Library size (millions of counts)") +
  theme_minimal(base_size=13) +
  theme(axis.text.x=element_text(angle=45,hjust=1))
dev.off()

# ------------------------
# 2. LOW-EXPRESSION FILTER
# ------------------------
keep <- rowSums(counts >= 10) >= 4
counts_f <- counts[keep,,drop=FALSE]

filter_summary <- data.frame(
  genes_before=nrow(counts),
  genes_after=sum(keep),
  genes_removed=sum(!keep),
  rule="raw count >= 10 in at least 4 samples"
)
write.csv(filter_summary, file.path(outdir, "filter_summary.csv"), row.names=FALSE)

# ------------------------
# 3. PRIMARY PAIRED DESEQ2 MODEL
# ------------------------
dds <- DESeqDataSetFromMatrix(
  countData=round(counts_f),
  colData=meta,
  design=~ cell + dex
)
dds <- DESeq(dds)

norm_counts <- counts(dds, normalized=TRUE)
write.csv(as.data.frame(norm_counts), file.path(outdir, "normalized_counts.csv"))

# ------------------------
# 4. EXPLORATORY QC: VST, PCA, CORRELATION
# ------------------------
vsd <- vst(dds, blind=FALSE)

pca_data <- plotPCA(vsd, intgroup=c("dex","cell"), returnData=TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
write.csv(pca_data, file.path(outdir, "pca_coordinates.csv"), row.names=FALSE)

p <- ggplot(pca_data, aes(PC1, PC2, color=dex, shape=cell, label=name)) +
  geom_point(size=4) +
  labs(title="GSE52778 PCA: dexamethasone vs untreated",
       x=paste0("PC1: ",percentVar[1],"% variance"),
       y=paste0("PC2: ",percentVar[2],"% variance")) +
  theme_minimal(base_size=14)
ggsave(file.path(outdir, "02_pca.png"), p, width=10, height=6, dpi=180)

cor_mat <- cor(assay(vsd), method="pearson")
write.csv(cor_mat, file.path(outdir, "sample_correlation_matrix.csv"))

png(file.path(outdir, "03_sample_correlation_heatmap.png"),
    width=1500,height=1300,res=150)
pheatmap(cor_mat, main="GSE52778 sample correlation")
dev.off()

# ------------------------
# 5. DIFFERENTIAL EXPRESSION
# ------------------------
res <- results(dds, contrast=c("dex","trt","untrt"), alpha=0.05)
res_df <- as.data.frame(res)
res_df$ensembl_gene_id <- rownames(res_df)

# Try to carry symbols from airway row metadata without inventing mappings.
anno <- as.data.frame(rowData(se))
sym_candidates <- c("symbol","gene_name","Symbol","gene_symbol")
sym_col <- sym_candidates[sym_candidates %in% colnames(anno)]
if (length(sym_col) > 0) {
  map <- anno[,c(sym_col[1]),drop=FALSE]
  map$ensembl_gene_id <- rownames(anno)
  colnames(map)[1] <- "gene_symbol"
  res_df <- merge(res_df,map,by="ensembl_gene_id",all.x=TRUE,sort=FALSE)
} else {
  res_df$gene_symbol <- NA_character_
}

res_df$direction <- "Not significant"
res_df$direction[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange >= 1] <- "Up"
res_df$direction[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange <= -1] <- "Down"

res_df <- res_df[order(res_df$padj),]
write.csv(res_df, file.path(outdir, "differential_expression_full.csv"), row.names=FALSE)
write.csv(subset(res_df,direction=="Up"),
          file.path(outdir, "significant_upregulated_genes.csv"), row.names=FALSE)
write.csv(subset(res_df,direction=="Down"),
          file.path(outdir, "significant_downregulated_genes.csv"), row.names=FALSE)

summary_out <- data.frame(
  metric=c("tested_rows","FDR_lt_0.05","FDR_and_absLFC_ge_1","up","down"),
  value=c(nrow(res_df),
          sum(!is.na(res_df$padj) & res_df$padj<0.05),
          sum(!is.na(res_df$padj) & res_df$padj<0.05 & abs(res_df$log2FoldChange)>=1),
          sum(res_df$direction=="Up"),
          sum(res_df$direction=="Down"))
)
write.csv(summary_out,file.path(outdir,"de_summary.csv"),row.names=FALSE)

# MA plot
png(file.path(outdir, "04_MA_plot.png"), width=1600,height=1000,res=160)
plotMA(res, ylim=c(-6,6), main="GSE52778 DESeq2 MA plot")
dev.off()

# Volcano plot
vol <- res_df
vol$minuslog10padj <- -log10(pmax(vol$padj, .Machine$double.xmin))
pvol <- ggplot(vol, aes(log2FoldChange, minuslog10padj, color=direction)) +
  geom_point(alpha=0.55, size=1.5, na.rm=TRUE) +
  geom_vline(xintercept=c(-1,1), linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  labs(title="GSE52778 volcano plot",
       x="log2 fold change (Dex / untreated)",
       y="-log10 adjusted p-value") +
  theme_minimal(base_size=14)
ggsave(file.path(outdir, "05_volcano_plot.png"), pvol, width=10, height=7, dpi=180)

# Top 30 heatmap
top30 <- head(res_df$ensembl_gene_id[!is.na(res_df$padj)], 30)
top30 <- intersect(top30, rownames(assay(vsd)))
if (length(top30) > 1) {
  mat <- assay(vsd)[top30,,drop=FALSE]
  mat <- mat - rowMeans(mat)
  ann <- data.frame(treatment=meta$dex, donor=meta$cell)
  rownames(ann) <- rownames(meta)
  png(file.path(outdir, "06_top30_heatmap.png"), width=1600,height=1500,res=150)
  pheatmap(mat, annotation_col=ann, show_rownames=FALSE,
           main="Top 30 differential genes")
  dev.off()
}

# ------------------------
# 6. SESSION / REPRODUCIBILITY
# ------------------------
capture.output(sessionInfo(), file=file.path(outdir,"sessionInfo.txt"))

cat("\nDONE\n")
cat("Generated outputs in:", outdir, "\n")
cat("Primary design: ~ cell + dex\n")
cat("Important: compare your de_summary.csv with results/reference/reference_summary.csv.\n")
