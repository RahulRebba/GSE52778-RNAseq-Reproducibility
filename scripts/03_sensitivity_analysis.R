# GSE52778 sensitivity analysis
# Purpose: test whether the treatment conclusions are sensitive to model specification.
# No sample is removed: the verified reference analysis did not flag an outlier.
# Compare paired model (~cell + dex) vs unpaired (~dex).
options(stringsAsFactors=FALSE)
root <- normalizePath(getwd(), winslash="/", mustWork=FALSE)
if (basename(root)=="scripts") root <- dirname(root)
outdir <- file.path(root,"results","generated")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

suppressPackageStartupMessages({
  library(airway); library(DESeq2)
})
data("airway")
cts <- assay(airway)
meta <- as.data.frame(colData(airway))
meta$dex <- relevel(factor(meta$dex), ref="untrt")
meta$cell <- factor(meta$cell)
keep <- rowSums(cts >= 10) >= 4
cts <- round(cts[keep,,drop=FALSE])

fit_one <- function(design_formula) {
  d <- DESeqDataSetFromMatrix(cts, meta, design_formula)
  d <- DESeq(d)
  r <- as.data.frame(results(d, contrast=c("dex","trt","untrt"), alpha=0.05))
  r$gene <- rownames(r)
  r
}

paired <- fit_one(~cell + dex)
unpaired <- fit_one(~dex)

sig_p <- paired$gene[!is.na(paired$padj) & paired$padj < .05 & abs(paired$log2FoldChange)>=1]
sig_u <- unpaired$gene[!is.na(unpaired$padj) & unpaired$padj < .05 & abs(unpaired$log2FoldChange)>=1]

summary <- data.frame(
  metric=c("paired_significant","unpaired_significant","overlap","paired_only","unpaired_only"),
  value=c(length(sig_p),length(sig_u),length(intersect(sig_p,sig_u)),
          length(setdiff(sig_p,sig_u)),length(setdiff(sig_u,sig_p)))
)
write.csv(summary,file.path(outdir,"sensitivity_paired_vs_unpaired_summary.csv"),row.names=FALSE)

merged <- merge(
  paired[,c("gene","log2FoldChange","padj")],
  unpaired[,c("gene","log2FoldChange","padj")],
  by="gene",suffixes=c("_paired","_unpaired"),all=TRUE
)
write.csv(merged,file.path(outdir,"sensitivity_paired_vs_unpaired_full.csv"),row.names=FALSE)
cat("Sensitivity analysis complete.\n")
