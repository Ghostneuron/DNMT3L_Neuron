#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
})

input_file <- "data/geo/GSE268951/supplementary/GSE268951_all.fpkm_mRNA_anno.txt.gz"
outdir <- "results/mouse_gse268951_count_de"
note_file <- "notes/gse268951_count_de_first_pass.md"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(note_file), recursive = TRUE, showWarnings = FALSE)

control_samples <- c("A-1", "A-2", "A-3", "B-1", "B-2", "B-3")
oe_samples <- c("C-1", "C-2", "C-3", "D-1", "D-2", "D-3")
samples <- c(control_samples, oe_samples)
count_cols <- paste0(samples, "_expected_count")
fpkm_cols <- paste0(samples, "_FPKM")
alpha <- 0.05
lfc_relaxed <- 0.58
lfc_large <- 1.0

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

safe_padj <- function(x) {
  x[is.na(x)] <- 1
  x
}

collapse_text <- function(x) {
  values <- unique(x[!is.na(x) & x != "" & x != "-"])
  if (length(values) == 0) {
    return("")
  }
  paste(values, collapse = ";")
}

df <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE)
df <- df[df$GeneID != "" & df$GeneID != "-" & df$GeneSymbol != "" & df$GeneSymbol != "-", , drop = FALSE]

for (col in c(count_cols, fpkm_cols)) {
  df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
  df[[col]][is.na(df[[col]])] <- 0
}

numeric_agg <- aggregate(df[, c(count_cols, fpkm_cols)], by = list(gene_id = df$GeneID), FUN = sum)
first_idx <- !duplicated(df$GeneID)
annot <- df[first_idx, c("GeneID", "GeneSymbol", "Description", "Chr", "Start", "End", "Strand", "KO"), drop = FALSE]
colnames(annot)[colnames(annot) == "GeneID"] <- "gene_id"

go_bp <- aggregate(df$GO_BP, by = list(gene_id = df$GeneID), FUN = collapse_text)
go_cc <- aggregate(df$GO_CC, by = list(gene_id = df$GeneID), FUN = collapse_text)
go_mf <- aggregate(df$GO_MF, by = list(gene_id = df$GeneID), FUN = collapse_text)
colnames(go_bp)[2] <- "GO_BP"
colnames(go_cc)[2] <- "GO_CC"
colnames(go_mf)[2] <- "GO_MF"

gene_df <- merge(annot, numeric_agg, by = "gene_id", all.x = TRUE, sort = FALSE)
gene_df <- merge(gene_df, go_bp, by = "gene_id", all.x = TRUE, sort = FALSE)
gene_df <- merge(gene_df, go_cc, by = "gene_id", all.x = TRUE, sort = FALSE)
gene_df <- merge(gene_df, go_mf, by = "gene_id", all.x = TRUE, sort = FALSE)

counts <- as.matrix(gene_df[, count_cols])
storage.mode(counts) <- "numeric"
counts <- round(counts)
rownames(counts) <- make.unique(gene_df$gene_id)
colnames(counts) <- samples

sample_info <- data.frame(
  sample = samples,
  sex = factor(c(rep("male", 3), rep("female", 3), rep("male", 3), rep("female", 3))),
  condition = factor(c(rep("control", length(control_samples)), rep("Dnmt3l_OE", length(oe_samples))),
                     levels = c("control", "Dnmt3l_OE")),
  row.names = samples
)

keep <- rowSums(counts >= 10) >= 3
dds <- DESeqDataSetFromMatrix(countData = counts[keep, , drop = FALSE], colData = sample_info, design = ~ sex + condition)
dds <- DESeq(dds, quiet = TRUE)
deseq_res <- results(dds, contrast = c("condition", "Dnmt3l_OE", "control"), alpha = alpha)
deseq_df <- as.data.frame(deseq_res)
deseq_df$gene_rowname <- rownames(deseq_df)

gene_df$gene_rowname <- make.unique(gene_df$gene_id)
merged <- merge(gene_df, deseq_df, by = "gene_rowname", all.x = TRUE, sort = FALSE)
merged$tested_by_deseq2 <- merged$gene_rowname %in% rownames(dds)
merged$deseq2_baseMean <- merged$baseMean
merged$deseq2_log2fc <- merged$log2FoldChange
merged$deseq2_lfcSE <- merged$lfcSE
merged$deseq2_stat <- merged$stat
merged$deseq2_pvalue <- merged$pvalue
merged$deseq2_padj <- safe_padj(merged$padj)
merged$control_mean_count <- rowMeans(as.matrix(merged[, paste0(control_samples, "_expected_count")]))
merged$dnmt3l_oe_mean_count <- rowMeans(as.matrix(merged[, paste0(oe_samples, "_expected_count")]))
merged$control_mean_fpkm <- rowMeans(as.matrix(merged[, paste0(control_samples, "_FPKM")]))
merged$dnmt3l_oe_mean_fpkm <- rowMeans(as.matrix(merged[, paste0(oe_samples, "_FPKM")]))
merged$fpkm_log2fc <- log2(merged$dnmt3l_oe_mean_fpkm + 0.1) - log2(merged$control_mean_fpkm + 0.1)
merged$deseq2_padj05 <- merged$deseq2_padj <= alpha
merged$deseq2_padj05_abslog2fc058 <- merged$deseq2_padj <= alpha & abs(merged$deseq2_log2fc) >= lfc_relaxed
merged$deseq2_padj05_abslog2fc1 <- merged$deseq2_padj <= alpha & abs(merged$deseq2_log2fc) >= lfc_large
merged$deseq2_nominal_p05_abslog2fc1 <- merged$deseq2_pvalue <= alpha & abs(merged$deseq2_log2fc) >= lfc_large

dge <- DGEList(counts = counts[keep, , drop = FALSE], group = sample_info$condition, genes = gene_df[keep, c("gene_id", "GeneSymbol")])
dge <- calcNormFactors(dge)
design <- model.matrix(~ sex + condition, data = sample_info)
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, coef = "conditionDnmt3l_OE")
edger <- topTags(qlf, n = Inf, sort.by = "none")$table
edger$gene_rowname <- rownames(edger)
edger_small <- edger[, c("gene_rowname", "logFC", "logCPM", "F", "PValue", "FDR")]
colnames(edger_small) <- c("gene_rowname", "edger_log2fc", "edger_logCPM", "edger_F", "edger_pvalue", "edger_fdr")
merged <- merge(merged, edger_small, by = "gene_rowname", all.x = TRUE, sort = FALSE)
merged$edger_fdr <- safe_padj(merged$edger_fdr)
merged$edger_fdr05_abslog2fc1 <- merged$edger_fdr <= alpha & abs(merged$edger_log2fc) >= lfc_large

output_cols <- c(
  "gene_id", "GeneSymbol", "Description", "Chr", "Start", "End", "Strand",
  count_cols, fpkm_cols,
  "control_mean_count", "dnmt3l_oe_mean_count", "control_mean_fpkm", "dnmt3l_oe_mean_fpkm",
  "fpkm_log2fc", "tested_by_deseq2", "deseq2_baseMean", "deseq2_log2fc", "deseq2_lfcSE",
  "deseq2_stat", "deseq2_pvalue", "deseq2_padj", "deseq2_padj05",
  "deseq2_padj05_abslog2fc058", "deseq2_padj05_abslog2fc1", "deseq2_nominal_p05_abslog2fc1",
  "edger_log2fc", "edger_logCPM", "edger_F", "edger_pvalue", "edger_fdr", "edger_fdr05_abslog2fc1",
  "GO_BP", "GO_CC", "GO_MF", "KO"
)

gene_results <- merged[, output_cols]
gene_results <- gene_results[order(gene_results$deseq2_padj, -abs(gene_results$deseq2_log2fc)), ]
write_tsv(gene_results, file.path(outdir, "gse268951_deseq2_edger_count_de_results.tsv"))
write_tsv(gene_results[gene_results$deseq2_padj05_abslog2fc058, ], file.path(outdir, "gse268951_deseq2_fdr05_abslog2fc058_genes.tsv"))
write_tsv(gene_results[gene_results$deseq2_padj05_abslog2fc1, ], file.path(outdir, "gse268951_deseq2_fdr05_abslog2fc1_genes.tsv"))
write_tsv(gene_results[gene_results$deseq2_nominal_p05_abslog2fc1, ], file.path(outdir, "gse268951_deseq2_nominal_p05_abslog2fc1_genes.tsv"))
write_tsv(
  gene_results[gene_results$deseq2_padj05_abslog2fc1 & gene_results$edger_fdr05_abslog2fc1, ],
  file.path(outdir, "gse268951_consensus_deseq2_edger_fdr05_abslog2fc1_genes.tsv")
)

summary_rows <- data.frame(
  metric = c(
    "input_transcript_rows",
    "collapsed_gene_rows",
    "samples",
    "design",
    "tested_by_deseq2",
    "deseq2_fdr05",
    "deseq2_fdr05_abslog2fc058",
    "deseq2_fdr05_abslog2fc1",
    "deseq2_nominal_p05_abslog2fc1",
    "edger_fdr05_abslog2fc1",
    "consensus_deseq2_edger_fdr05_abslog2fc1"
  ),
  value = c(
    nrow(df),
    nrow(gene_df),
    length(samples),
    "~ sex + condition",
    sum(merged$tested_by_deseq2, na.rm = TRUE),
    sum(merged$deseq2_padj <= alpha, na.rm = TRUE),
    sum(merged$deseq2_padj05_abslog2fc058, na.rm = TRUE),
    sum(merged$deseq2_padj05_abslog2fc1, na.rm = TRUE),
    sum(merged$deseq2_nominal_p05_abslog2fc1, na.rm = TRUE),
    sum(merged$edger_fdr05_abslog2fc1, na.rm = TRUE),
    sum(merged$deseq2_padj05_abslog2fc1 & merged$edger_fdr05_abslog2fc1, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)
write_tsv(summary_rows, file.path(outdir, "gse268951_count_de_summary.tsv"))

target_symbols <- c("Dnmt3l", "Dnmt3a", "Dnmt3b", "Stat1", "Stat3", "Hes3", "Ascl1", "Neurog2", "Neurod2", "Cdkn1a")
target_rows <- gene_results[gene_results$GeneSymbol %in% target_symbols, ]
target_rows <- target_rows[match(target_symbols[target_symbols %in% target_rows$GeneSymbol], target_rows$GeneSymbol), ]
write_tsv(target_rows, file.path(outdir, "gse268951_target_marker_count_de.tsv"))

note_lines <- c(
  "# GSE268951 Count-Based DE First Pass",
  "",
  "Date: 2026-06-11",
  "",
  "This public-data reanalysis uses the GEO supplementary table `GSE268951_all.fpkm_mRNA_anno.txt.gz`, which contains both per-sample `expected_count` and FPKM columns despite the filename. Transcript rows were collapsed to gene level by summing `expected_count` and FPKM values across rows with the same `GeneID`.",
  "",
  "These public expected-count DESeq2/edgeR results are the manuscript-facing `GSE268951` RNA-seq effects. FPKM columns are retained only for abundance/context.",
  "",
  "The sensitivity DESeq2 model used design `~ sex + condition`, comparing Dnmt3l-OE samples (`C-1` to `C-3`, `D-1` to `D-3`) against controls (`A-1` to `A-3`, `B-1` to `B-3`). edgeR quasi-likelihood results are included as an additional sensitivity check.",
  "",
  "## Summary",
  "",
  paste0("- Genes tested by DESeq2: ", sum(merged$tested_by_deseq2, na.rm = TRUE), "."),
  paste0("- DESeq2 FDR <= 0.05 genes: ", sum(merged$deseq2_padj <= alpha, na.rm = TRUE), "."),
  paste0("- DESeq2 FDR <= 0.05 and |log2FC| >= 0.58 genes: ", sum(merged$deseq2_padj05_abslog2fc058, na.rm = TRUE), "."),
  paste0("- DESeq2 FDR <= 0.05 and |log2FC| >= 1 genes: ", sum(merged$deseq2_padj05_abslog2fc1, na.rm = TRUE), "."),
  paste0("- edgeR FDR <= 0.05 and |log2FC| >= 1 genes: ", sum(merged$edger_fdr05_abslog2fc1, na.rm = TRUE), "."),
  paste0("- DESeq2/edgeR consensus FDR <= 0.05 and |log2FC| >= 1 genes: ", sum(merged$deseq2_padj05_abslog2fc1 & merged$edger_fdr05_abslog2fc1, na.rm = TRUE), "."),
  "",
  "## Outputs",
  "",
  "- `results/mouse_gse268951_count_de/gse268951_count_de_summary.tsv`",
  "- `results/mouse_gse268951_count_de/gse268951_deseq2_edger_count_de_results.tsv`",
  "- `results/mouse_gse268951_count_de/gse268951_deseq2_fdr05_abslog2fc058_genes.tsv`",
  "- `results/mouse_gse268951_count_de/gse268951_deseq2_fdr05_abslog2fc1_genes.tsv`",
  "- `results/mouse_gse268951_count_de/gse268951_target_marker_count_de.tsv`"
)
writeLines(note_lines, note_file)

message("Wrote ", file.path(outdir, "gse268951_deseq2_edger_count_de_results.tsv"))
message("Wrote ", file.path(outdir, "gse268951_count_de_summary.tsv"))
message("Wrote ", note_file)
