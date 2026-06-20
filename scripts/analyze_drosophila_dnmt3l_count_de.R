#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
})

input_file <- "data/geo/GSE263069/supplementary/GSE263069_all.fpkm_anno.txt.gz"
outdir <- "results/drosophila_dnmt3l_count_de"
figdir <- "figures/drosophila_dnmt3l_count_de"
note_file <- "notes/drosophila_dnmt3l_count_de_first_pass.md"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(note_file), recursive = TRUE, showWarnings = FALSE)

control_samples <- c("C1", "C2", "C3")
oe_samples <- c("D1", "D2", "D3")
count_cols <- c(control_samples, oe_samples)
fpkm_cols <- paste0(count_cols, "_FPKM")
alpha <- 0.05
lfc_relaxed <- 0.58
lfc_large <- 1.0

modules <- list(
  neural_development = "neuro|nervous system|axon|dendrit|mushroom body|brain|synap|glia|photoreceptor",
  synapse_neurotransmission = "synap|neurotrans|cholinergic|serotonergic|glutamate|GABA|neuromuscular|axon guidance",
  ion_transport_calcium = "calcium|voltage[- ]gated|ion transport|ion channel|potassium|sodium|calmodulin|transmembrane transport",
  chromatin_transcription = "chromatin|histone|DNA binding|transcription|RNA polymerase|enhancer|nucleus",
  cell_cycle_growth = "cell cycle|mitotic|proliferation|growth|cell fate|apoptotic|differentiation",
  stress_immune = "stress|immune|defense|heat|oxidative|antimicrobial|response to",
  metabolism_detox = "metabolic|metabolism|drug|oxidation|monooxygenase|cytochrome|glutathione|hydrolase"
)

selected_markers <- data.frame(
  marker_group = c(
    rep("JAK/STAT", 7),
    rep("neurogenic bHLH/Notch", 14),
    "axon/synapse"
  ),
  requested_symbol = c(
    "Stat92E", "hop", "upd1", "upd2", "upd3", "Socs36E", "Socs44A",
    "ac", "sc", "ase", "ato", "da", "sens", "pros", "elav", "dpn", "N", "Dl", "E(spl)m8-HLH", "H", "neur",
    "Abl"
  ),
  stringsAsFactors = FALSE
)

read_input <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

safe_padj <- function(x) {
  x[is.na(x)] <- 1
  x
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

bh_adjust <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  out[ok] <- p.adjust(x[ok], method = "BH")
  out
}

clean_text <- function(x) {
  x <- gsub("%3B", ";", x)
  x[is.na(x)] <- ""
  x
}

hyper_enrichment <- function(results, query_col) {
  bg <- results[results$tested_by_deseq2, , drop = FALSE]
  query <- bg[[query_col]]
  rows <- list()
  searchable <- paste(
    clean_text(bg$Description),
    clean_text(bg$GO_BP),
    clean_text(bg$GO_CC),
    clean_text(bg$GO_MF),
    sep = " "
  )
  for (module in names(modules)) {
    hits <- grepl(modules[[module]], searchable, ignore.case = TRUE, perl = TRUE)
    tab <- matrix(
      c(
        sum(query & hits, na.rm = TRUE),
        sum(query & !hits, na.rm = TRUE),
        sum(!query & hits, na.rm = TRUE),
        sum(!query & !hits, na.rm = TRUE)
      ),
      nrow = 2,
      byrow = TRUE
    )
    ft <- fisher.test(tab, alternative = "greater")
    module_query <- bg[query & hits, , drop = FALSE]
    module_query <- module_query[order(module_query$deseq2_pvalue, -abs(module_query$deseq2_log2fc)), , drop = FALSE]
    rows[[length(rows) + 1]] <- data.frame(
      query_set = query_col,
      module = module,
      background_genes = nrow(bg),
      module_background_genes = sum(hits, na.rm = TRUE),
      query_genes = sum(query, na.rm = TRUE),
      module_query_genes = sum(query & hits, na.rm = TRUE),
      module_query_up = sum(query & hits & bg$deseq2_log2fc > 0, na.rm = TRUE),
      module_query_down = sum(query & hits & bg$deseq2_log2fc < 0, na.rm = TRUE),
      fold_enrichment = (sum(query & hits, na.rm = TRUE) / max(sum(query, na.rm = TRUE), 1)) /
        (sum(hits, na.rm = TRUE) / max(nrow(bg), 1)),
      odds_ratio = unname(ft$estimate),
      pvalue = ft$p.value,
      example_genes = paste(head(module_query$GeneSymbol, 25), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

plot_png_pdf <- function(prefix, width = 7, height = 5, code) {
  png(paste0(prefix, ".png"), width = width, height = height, units = "in", res = 300)
  force(code)
  dev.off()
  pdf(paste0(prefix, ".pdf"), width = width, height = height)
  force(code)
  dev.off()
}

df <- read_input(input_file)
counts <- as.matrix(df[, count_cols])
storage.mode(counts) <- "numeric"
counts <- round(counts)
rownames(counts) <- make.unique(df$gene_id)

sample_info <- data.frame(
  sample = count_cols,
  condition = factor(c(rep("control", length(control_samples)), rep("DNMT3L_OE", length(oe_samples))),
                     levels = c("control", "DNMT3L_OE")),
  row.names = count_cols
)

keep <- rowSums(counts >= 10) >= 3
dds <- DESeqDataSetFromMatrix(countData = counts[keep, , drop = FALSE], colData = sample_info, design = ~ condition)
dds <- DESeq(dds, quiet = TRUE)
deseq_res <- results(dds, contrast = c("condition", "DNMT3L_OE", "control"), alpha = alpha)
deseq_df <- as.data.frame(deseq_res)
deseq_df$gene_rowname <- rownames(deseq_df)

annot <- df
annot$gene_rowname <- make.unique(df$gene_id)
merged <- merge(annot, deseq_df, by = "gene_rowname", all.x = TRUE, sort = FALSE)
merged$tested_by_deseq2 <- merged$gene_rowname %in% rownames(dds)
merged$deseq2_baseMean <- merged$baseMean
merged$deseq2_log2fc <- merged$log2FoldChange
merged$deseq2_lfcSE <- merged$lfcSE
merged$deseq2_stat <- merged$stat
merged$deseq2_pvalue <- merged$pvalue
merged$deseq2_padj <- safe_padj(merged$padj)
merged$control_mean_count <- rowMeans(as.matrix(merged[, control_samples]))
merged$dnmt3l_oe_mean_count <- rowMeans(as.matrix(merged[, oe_samples]))
merged$control_mean_fpkm <- rowMeans(as.matrix(merged[, paste0(control_samples, "_FPKM")]))
merged$dnmt3l_oe_mean_fpkm <- rowMeans(as.matrix(merged[, paste0(oe_samples, "_FPKM")]))
merged$fpkm_log2fc <- log2(merged$dnmt3l_oe_mean_fpkm + 0.1) - log2(merged$control_mean_fpkm + 0.1)
merged$deseq2_padj05 <- merged$deseq2_padj <= alpha
merged$deseq2_padj05_abslog2fc058 <- merged$deseq2_padj <= alpha & abs(merged$deseq2_log2fc) >= lfc_relaxed
merged$deseq2_padj05_abslog2fc1 <- merged$deseq2_padj <= alpha & abs(merged$deseq2_log2fc) >= lfc_large
merged$deseq2_nominal_p05_abslog2fc1 <- merged$deseq2_pvalue <= alpha & abs(merged$deseq2_log2fc) >= lfc_large

dge <- DGEList(counts = counts[keep, , drop = FALSE], group = sample_info$condition, genes = df[keep, c("gene_id", "GeneSymbol")])
dge <- calcNormFactors(dge)
design <- model.matrix(~ condition, data = sample_info)
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, coef = "conditionDNMT3L_OE")
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
write_tsv(gene_results, file.path(outdir, "gse263069_deseq2_edger_count_de_results.tsv"))
write_tsv(gene_results[gene_results$deseq2_padj05_abslog2fc058, ], file.path(outdir, "gse263069_deseq2_fdr05_abslog2fc058_genes.tsv"))
write_tsv(gene_results[gene_results$deseq2_padj05_abslog2fc1, ], file.path(outdir, "gse263069_deseq2_fdr05_abslog2fc1_genes.tsv"))
write_tsv(gene_results[gene_results$deseq2_nominal_p05_abslog2fc1, ], file.path(outdir, "gse263069_deseq2_nominal_p05_abslog2fc1_genes.tsv"))
write_tsv(
  gene_results[gene_results$deseq2_padj05_abslog2fc1 & gene_results$edger_fdr05_abslog2fc1, ],
  file.path(outdir, "gse263069_consensus_deseq2_edger_fdr05_abslog2fc1_genes.tsv")
)

module_tables <- list(
  hyper_enrichment(merged, "deseq2_padj05_abslog2fc058"),
  hyper_enrichment(merged, "deseq2_padj05_abslog2fc1"),
  hyper_enrichment(merged, "deseq2_nominal_p05_abslog2fc1")
)
module_enrichment <- do.call(rbind, module_tables)
module_enrichment$fdr <- ave(module_enrichment$pvalue, module_enrichment$query_set, FUN = function(x) p.adjust(x, method = "BH"))
module_enrichment <- module_enrichment[order(module_enrichment$query_set, module_enrichment$fdr, module_enrichment$pvalue), ]
write_tsv(module_enrichment, file.path(outdir, "gse263069_deseq2_module_enrichment.tsv"))

marker_rows <- merge(selected_markers, gene_results, by.x = "requested_symbol", by.y = "GeneSymbol", all.x = TRUE, sort = FALSE)
marker_rows <- marker_rows[order(marker_rows$marker_group, marker_rows$requested_symbol), ]
write_tsv(marker_rows, file.path(outdir, "gse263069_deseq2_selected_marker_effects.tsv"))

summary_rows <- data.frame(
  metric = c(
    "input_genes",
    "deseq2_tested_genes_after_count_filter",
    "deseq2_fdr05_genes",
    "deseq2_fdr05_abslog2fc058_genes",
    "deseq2_fdr05_abslog2fc1_genes",
    "deseq2_nominal_p05_abslog2fc1_genes",
    "edger_fdr05_abslog2fc1_genes",
    "deseq2_and_edger_fdr05_abslog2fc1_overlap"
  ),
  value = c(
    nrow(df),
    nrow(dds),
    sum(merged$deseq2_padj05, na.rm = TRUE),
    sum(merged$deseq2_padj05_abslog2fc058, na.rm = TRUE),
    sum(merged$deseq2_padj05_abslog2fc1, na.rm = TRUE),
    sum(merged$deseq2_nominal_p05_abslog2fc1, na.rm = TRUE),
    sum(merged$edger_fdr05_abslog2fc1, na.rm = TRUE),
    sum(merged$deseq2_padj05_abslog2fc1 & merged$edger_fdr05_abslog2fc1, na.rm = TRUE)
  )
)
write_tsv(summary_rows, file.path(outdir, "gse263069_count_de_summary.tsv"))

plot_png_pdf(file.path(figdir, "gse263069_deseq2_ma_plot"), width = 7, height = 5, {
  x <- log10(merged$deseq2_baseMean + 1)
  y <- merged$deseq2_log2fc
  ok <- is.finite(x) & is.finite(y) & merged$tested_by_deseq2
  sig <- merged$deseq2_padj05_abslog2fc1
  plot(x[ok], y[ok], pch = 16, cex = 0.35, col = ifelse(sig[ok], "#B33A3A", "#8A8A8A"),
       xlab = "log10 mean normalized count + 1", ylab = "DESeq2 log2FC (DNMT3L-OE vs control)",
       main = "GSE263069 count-based MA plot")
  abline(h = 0, col = "grey35", lwd = 0.8)
  abline(h = c(-lfc_large, lfc_large), lty = 2, col = "grey55")
  legend("topright", legend = c("FDR <= 0.05 & |log2FC| >= 1", "other"),
         col = c("#B33A3A", "#8A8A8A"), pch = 16, bty = "n", cex = 0.8)
})

plot_png_pdf(file.path(figdir, "gse263069_deseq2_volcano"), width = 7, height = 5.5, {
  x <- merged$deseq2_log2fc
  y <- -log10(merged$deseq2_padj)
  y[!is.finite(y)] <- NA
  sig <- merged$deseq2_padj05_abslog2fc1
  plot(x, y, pch = 16, cex = 0.35, col = ifelse(sig, "#B33A3A", "#8A8A8A"),
       xlab = "DESeq2 log2FC (DNMT3L-OE vs control)", ylab = "-log10(FDR)",
       main = "GSE263069 count-based DESeq2")
  abline(v = c(-lfc_large, lfc_large), lty = 2, col = "grey45")
  abline(h = -log10(alpha), lty = 2, col = "grey45")
  legend("topright", legend = c("FDR <= 0.05 & |log2FC| >= 1", "other"),
         col = c("#B33A3A", "#8A8A8A"), pch = 16, bty = "n", cex = 0.8)
})

plot_png_pdf(file.path(figdir, "gse263069_deseq2_vs_fpkm_log2fc"), width = 6, height = 6, {
  tested <- merged$tested_by_deseq2 & is.finite(merged$deseq2_log2fc) & is.finite(merged$fpkm_log2fc)
  plot(merged$fpkm_log2fc[tested], merged$deseq2_log2fc[tested],
       pch = 16, cex = 0.35, col = "#4C78A8",
       xlab = "FPKM first-pass log2FC", ylab = "DESeq2 count log2FC",
       main = "Effect-size concordance")
  abline(0, 1, lty = 2, col = "grey45")
  abline(h = 0, v = 0, col = "grey85")
})

plot_png_pdf(file.path(figdir, "gse263069_deseq2_module_enrichment"), width = 8, height = 5, {
  sub <- module_enrichment[module_enrichment$query_set == "deseq2_padj05_abslog2fc1", ]
  if (nrow(sub) == 0 || sum(sub$module_query_genes, na.rm = TRUE) == 0) {
    sub <- module_enrichment[module_enrichment$query_set == "deseq2_padj05_abslog2fc058", ]
  }
  sub <- sub[is.finite(sub$fold_enrichment), , drop = FALSE]
  if (nrow(sub) == 0) {
    plot.new()
    text(0.5, 0.5, "No finite module enrichment values", cex = 1.1)
  } else {
    sub <- sub[order(sub$fold_enrichment), ]
    vals <- pmin(sub$fold_enrichment, 8)
    cols <- ifelse(sub$fdr <= 0.05, "#B33A3A", "#4C78A8")
    label_map <- c(
      metabolism_detox = "metabolism",
      ion_transport_calcium = "ion transport",
      stress_immune = "stress/immune",
      synapse_neurotransmission = "synapse",
      neural_development = "neural dev.",
      cell_cycle_growth = "cell cycle",
      chromatin_transcription = "chromatin"
    )
    labels <- unname(label_map[sub$module])
    labels[is.na(labels)] <- gsub("_", " ", sub$module[is.na(labels)])
    par(mar = c(4.5, 10, 3.5, 1.5))
    barplot(vals, names.arg = labels, horiz = TRUE, las = 1, col = cols,
            xlab = "Fold enrichment", main = "DESeq2 module enrichment: FDR <= 0.05, |log2FC| >= 1",
            cex.names = 0.9)
    legend("bottomright", legend = c("FDR <= 0.05", "not FDR"), fill = c("#B33A3A", "#4C78A8"), bty = "n", cex = 0.8)
  }
})

plot_png_pdf(file.path(figdir, "gse263069_deseq2_selected_marker_effects"), width = 8, height = 6, {
  sub <- marker_rows[!is.na(marker_rows$deseq2_log2fc), , drop = FALSE]
  if (nrow(sub) == 0) {
    plot.new()
    text(0.5, 0.5, "Selected markers were not detected in count-DE table", cex = 1.0)
  } else {
    sub <- sub[order(sub$marker_group, sub$deseq2_log2fc), , drop = FALSE]
    colors <- ifelse(
      sub$marker_group == "JAK/STAT",
      "#59A14F",
      ifelse(sub$marker_group == "axon/synapse", "#8E6C8A", "#4C78A8")
    )
    labels <- paste0(sub$requested_symbol, " (", sub$marker_group, ")")
    par(mar = c(4.5, 14, 3.5, 2))
    barplot(
      sub$deseq2_log2fc,
      names.arg = labels,
      horiz = TRUE,
      las = 1,
      col = colors,
      border = NA,
      xlab = "DESeq2 log2FC (DNMT3L-OE vs control)",
      main = "Selected fly JAK/STAT and neurogenic marker effects",
      cex.names = 0.72
    )
    abline(v = 0, col = "grey20", lwd = 0.8)
    abline(v = c(-1, 1), col = "grey55", lty = 2, lwd = 0.8)
  }
})

top_deseq <- head(gene_results[gene_results$deseq2_padj05_abslog2fc1, c("GeneSymbol", "deseq2_log2fc", "deseq2_padj")], 15)
top_consensus <- head(gene_results[gene_results$deseq2_padj05_abslog2fc1 & gene_results$edger_fdr05_abslog2fc1, c("GeneSymbol", "deseq2_log2fc", "deseq2_padj", "edger_fdr")], 15)
top_text <- if (nrow(top_deseq) == 0) {
  "No genes passed DESeq2 FDR <= 0.05 and |log2FC| >= 1."
} else {
  paste(apply(top_deseq, 1, function(x) paste0("`", x[["GeneSymbol"]], "` log2FC=", signif(as.numeric(x[["deseq2_log2fc"]]), 3), ", FDR=", signif(as.numeric(x[["deseq2_padj"]]), 3))), collapse = "; ")
}
consensus_text <- if (nrow(top_consensus) == 0) {
  "No genes passed both DESeq2 and edgeR FDR <= 0.05 with |log2FC| >= 1."
} else {
  paste(apply(top_consensus, 1, function(x) paste0("`", x[["GeneSymbol"]], "` DESeq2 log2FC=", signif(as.numeric(x[["deseq2_log2fc"]]), 3), ", DESeq2 FDR=", signif(as.numeric(x[["deseq2_padj"]]), 3), ", edgeR FDR=", signif(as.numeric(x[["edger_fdr"]]), 3))), collapse = "; ")
}

note_lines <- c(
  "# GSE263069 Count-Based DNMT3L-OE Differential Expression",
  "",
  "Date: 2026-06-11",
  "",
  "## Method",
  "",
  "The `C1-C3` and `D1-D3` expected-count columns from `GSE263069_all.fpkm_anno.txt.gz` were rounded to integer counts and analyzed with DESeq2 using design `~ condition` (`DNMT3L_OE` vs `control`). Genes were retained for DESeq2 testing if they had at least 10 counts in at least 3 samples. edgeR quasi-likelihood testing was run as a sensitivity analysis using the same filtered count matrix.",
  "",
  "This is stronger than the earlier Welch test on `log2(FPKM + 0.1)`, but it is still based on public processed expected counts rather than raw FASTQ re-quantification.",
  "",
  "## Summary",
  "",
  paste0("- Input genes: ", nrow(df), "."),
  paste0("- DESeq2-tested genes after count filter: ", nrow(dds), "."),
  paste0("- DESeq2 FDR <= 0.05 genes: ", sum(merged$deseq2_padj05, na.rm = TRUE), "."),
  paste0("- DESeq2 FDR <= 0.05 and |log2FC| >= 0.58 genes: ", sum(merged$deseq2_padj05_abslog2fc058, na.rm = TRUE), "."),
  paste0("- DESeq2 FDR <= 0.05 and |log2FC| >= 1 genes: ", sum(merged$deseq2_padj05_abslog2fc1, na.rm = TRUE), "."),
  paste0("- DESeq2 nominal P <= 0.05 and |log2FC| >= 1 genes: ", sum(merged$deseq2_nominal_p05_abslog2fc1, na.rm = TRUE), "."),
  paste0("- edgeR QL FDR <= 0.05 and |log2FC| >= 1 genes: ", sum(merged$edger_fdr05_abslog2fc1, na.rm = TRUE), "."),
  paste0("- DESeq2/edgeR consensus FDR <= 0.05 and |log2FC| >= 1 genes: ", sum(merged$deseq2_padj05_abslog2fc1 & merged$edger_fdr05_abslog2fc1, na.rm = TRUE), "."),
  "",
  "## Top DESeq2 FDR/large-effect genes",
  "",
  top_text,
  "",
  "## DESeq2/edgeR consensus large-effect genes",
  "",
  consensus_text,
  "",
  "## Interpretation",
  "",
  "The count-based analysis should replace the FPKM-Welch screen for gene-level claims. If a gene is significant only in the FPKM screen but not in the count-based DESeq2/edgeR tables, it should remain exploratory.",
  "",
  "For the manuscript, the Drosophila data remain a methylation-system-independent comparator. They can support the plausibility that human DNMT3L perturbs neural/developmental/stress transcription outside mammalian DNMT3A/DNMT3B CpG methylation, but the claim strength should follow the count-based FDR results and the edgeR sensitivity results.",
  ""
)
writeLines(note_lines, note_file)

message("Wrote ", outdir)
message("Wrote ", figdir)
message("Wrote ", note_file)
