#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
})

input_file <- "data/geo/GSE268951/supplementary/GSE268951_all.fpkm_mRNA_anno.txt.gz"
outdir <- "results/mouse_gse268951_count_de"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

control_samples <- c("A-1", "A-2", "A-3", "B-1", "B-2", "B-3")
oe_samples <- c("C-1", "C-2", "C-3", "D-1", "D-2", "D-3")
samples <- c(control_samples, oe_samples)
count_cols <- paste0(samples, "_expected_count")

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

df <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE)
counts <- as.matrix(df[, count_cols])
storage.mode(counts) <- "numeric"
counts[is.na(counts)] <- 0
counts <- round(counts)
rownames(counts) <- make.unique(df$gene_id)
colnames(counts) <- samples

group <- factor(c(rep("CON", length(control_samples)), rep("OE", length(oe_samples))), levels = c("CON", "OE"))
dge <- DGEList(counts = counts, group = group)
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge)
etest <- exactTest(dge, pair = c("CON", "OE"))
edge <- topTags(etest, n = Inf, sort.by = "none")$table
edge$gene_id <- rownames(edge)
edge <- edge[, c("gene_id", "logFC", "logCPM", "PValue", "FDR")]

annot <- df[, c("gene_id", "GeneID", "Chr", "Start", "End", "Strand", "GeneSymbol", "Description", "KO", "GO_BP", "GO_CC", "GO_MF")]
colnames(annot)[colnames(annot) == "GeneID"] <- "geneID"
out <- merge(edge, annot, by = "gene_id", all.x = TRUE, sort = FALSE)
out <- out[order(out$FDR, out$PValue), ]
write_tsv(out, file.path(outdir, "gse268951_transcript_edgeR_exacttest_public_counts.tsv"))

message("Wrote ", file.path(outdir, "gse268951_transcript_edgeR_exacttest_public_counts.tsv"))
