#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default) {
  key <- paste0("--", name, "=")
  hit <- args[startsWith(args, key)]
  if (length(hit) == 0) {
    return(default)
  }
  sub(key, "", hit[[1]], fixed = TRUE)
}

gse301886_counts <- arg_value("gse301886-counts", "data/geo/GSE301886/supplementary/GSE301886_T21_raw_counts.txt.gz")
gse301886_meta <- arg_value("gse301886-meta", "data/geo/GSE301886/metadata/sample_metadata.tsv")
gse305139_counts <- arg_value("gse305139-counts", "data/geo/GSE305139/supplementary/GSE305139_raw_counts_tissue_bulk.tsv.gz")
gse305139_meta <- arg_value("gse305139-meta", "data/geo/GSE305139/metadata/sample_metadata.tsv")
gencode_gtf <- arg_value("gencode-gtf", "data/reference/gencode/gencode.v39.annotation.gtf.gz")
outdir <- arg_value("outdir", "results/human_ds_bulk_count_de")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

parse_characteristics <- function(x) {
  parts <- trimws(strsplit(x, ";", fixed = TRUE)[[1]])
  parts <- parts[nzchar(parts) & grepl("=", parts, fixed = TRUE)]
  out <- setNames(rep("", length(parts)), rep("", length(parts)))
  for (part in parts) {
    kv <- strsplit(part, "=", fixed = TRUE)[[1]]
    key <- tolower(trimws(kv[[1]]))
    value <- trimws(paste(kv[-1], collapse = "="))
    out[[key]] <- value
  }
  out
}

clean_ensembl <- function(x) sub("\\..*$", "", x)

normalize_group <- function(x) {
  y <- tolower(trimws(x))
  ifelse(y %in% c("t21", "trisomy 21", "ds", "down syndrome"), "DS",
    ifelse(y %in% c("control", "con", "disomic", "euploid"), "Control", trimws(x))
  )
}

read_gencode_symbols <- function(path) {
  if (!file.exists(path)) {
    return(data.frame(gene_id = character(), gene_symbol = character()))
  }
  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)
  gene_id <- character()
  gene_symbol <- character()
  repeat {
    lines <- readLines(con, n = 50000, warn = FALSE)
    if (length(lines) == 0) {
      break
    }
    lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
    if (length(lines) == 0) {
      next
    }
    parts <- strsplit(lines, "\t", fixed = TRUE)
    gene_lines <- vapply(parts, function(x) length(x) >= 9 && x[[3]] == "gene", logical(1))
    if (!any(gene_lines)) {
      next
    }
    attrs <- vapply(parts[gene_lines], function(x) x[[9]], character(1))
    ids <- sub('.*gene_id "([^"]+)".*', "\\1", attrs)
    syms <- sub('.*gene_name "([^"]+)".*', "\\1", attrs)
    ids <- clean_ensembl(ids)
    keep <- nzchar(ids) & nzchar(syms) & ids != attrs & syms != attrs
    gene_id <- c(gene_id, ids[keep])
    gene_symbol <- c(gene_symbol, syms[keep])
  }
  unique(data.frame(gene_id = gene_id, gene_symbol = gene_symbol, stringsAsFactors = FALSE))
}

metadata_gse301886 <- function(path) {
  raw <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  rows <- lapply(seq_len(nrow(raw)), function(i) {
    chars <- parse_characteristics(raw$characteristics[[i]])
    title <- raw$title[[i]]
    sample <- trimws(sub("^.*,", "", title))
    data.frame(
      sample = sample,
      geo_sample_accession = raw$accession[[i]],
      dataset = "GSE301886",
      context = chars[["tissue"]],
      group = normalize_group(chars[["karyotype"]]),
      age_pcw = suppressWarnings(as.numeric(chars[["gestational age"]])),
      sex = chars[["sex"]],
      final_analysis = TRUE,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

metadata_gse305139 <- function(path) {
  raw <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  rows <- lapply(seq_len(nrow(raw)), function(i) {
    chars <- parse_characteristics(raw$characteristics[[i]])
    data.frame(
      sample = chars[["library igf_id"]],
      geo_sample_accession = raw$accession[[i]],
      dataset = "GSE305139",
      context = "fetal_brain_final",
      group = normalize_group(ifelse(nzchar(chars[["group"]]), chars[["group"]], chars[["genotype"]])),
      age_pcw = suppressWarnings(as.numeric(chars[["dev stage_pcw"]])),
      sex = chars[["sex"]],
      final_analysis = toupper(chars[["in final_analysis"]]) != "FALSE",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

read_gse301886_counts <- function(path, gencode) {
  tab <- read.delim(gzfile(path), check.names = FALSE, stringsAsFactors = FALSE)
  gene_id <- clean_ensembl(tab[[1]])
  idx <- match(gene_id, gencode$gene_id)
  gene_symbol <- ifelse(is.na(idx), "", gencode$gene_symbol[idx])
  counts <- as.matrix(tab[, -1, drop = FALSE])
  storage.mode(counts) <- "numeric"
  rownames(counts) <- make.unique(gene_id)
  list(counts = counts, gene_info = data.frame(gene_id = gene_id, gene_symbol = gene_symbol, stringsAsFactors = FALSE))
}

read_gse305139_counts <- function(path) {
  tab <- read.delim(gzfile(path), check.names = FALSE, stringsAsFactors = FALSE)
  gene_id <- clean_ensembl(tab[[1]])
  gene_symbol <- tab[[2]]
  counts <- as.matrix(tab[, -(1:2), drop = FALSE])
  storage.mode(counts) <- "numeric"
  rownames(counts) <- make.unique(gene_id)
  list(counts = counts, gene_info = data.frame(gene_id = gene_id, gene_symbol = gene_symbol, stringsAsFactors = FALSE))
}

fit_edgeR_context <- function(dataset, context, counts_all, gene_info, meta_all) {
  meta <- meta_all[meta_all$context == context & meta_all$group %in% c("Control", "DS") & meta_all$final_analysis, , drop = FALSE]
  meta <- meta[match(intersect(colnames(counts_all), meta$sample), meta$sample), , drop = FALSE]
  counts <- counts_all[, meta$sample, drop = FALSE]
  group_counts <- table(meta$group)
  if (!all(c("Control", "DS") %in% names(group_counts)) || any(group_counts[c("Control", "DS")] < 2)) {
    stop(sprintf("Not enough samples for %s %s", dataset, context))
  }

  meta$group <- factor(meta$group, levels = c("Control", "DS"))
  meta$sex <- factor(ifelse(nzchar(meta$sex), meta$sex, "unknown"))
  y0 <- DGEList(counts = counts)
  cpm_values <- cpm(y0, log = FALSE)
  log2_cpm_plus_1 <- log2(cpm_values + 1)
  ds_idx <- meta$group == "DS"
  control_idx <- meta$group == "Control"

  design_full <- model.matrix(~ age_pcw + sex + group, data = meta)
  design_group <- model.matrix(~ group, data = meta)
  design <- design_full
  model_formula <- "~ age_pcw + sex + group"
  if (qr(design_full)$rank < ncol(design_full)) {
    design <- design_group
    model_formula <- "~ group"
  }
  rownames(design) <- meta$sample

  y <- DGEList(counts = counts)
  keep <- filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  design_filtered <- design[colnames(y), , drop = FALSE]
  y <- estimateDisp(y, design_filtered, robust = TRUE)
  fit <- glmQLFit(y, design_filtered, robust = TRUE)
  coef_name <- grep("^groupDS$", colnames(design_filtered), value = TRUE)
  if (length(coef_name) != 1) {
    stop(sprintf("Could not find groupDS coefficient for %s %s", dataset, context))
  }
  qlf <- glmQLFTest(fit, coef = coef_name)
  de <- qlf$table
  de$PValue <- qlf$table$PValue
  de$FDR <- p.adjust(de$PValue, method = "BH")
  de$gene_id <- clean_ensembl(rownames(de))

  result <- data.frame(
    dataset = dataset,
    context = context,
    gene_id = gene_info$gene_id,
    gene_symbol = gene_info$gene_symbol,
    n_ds = sum(ds_idx),
    n_control = sum(control_idx),
    model_formula = model_formula,
    final_analysis_only = TRUE,
    mean_ds_log2_cpm_plus_1 = rowMeans(log2_cpm_plus_1[, ds_idx, drop = FALSE]),
    mean_control_log2_cpm_plus_1 = rowMeans(log2_cpm_plus_1[, control_idx, drop = FALSE]),
    delta_ds_minus_control_log2_cpm_plus_1 =
      rowMeans(log2_cpm_plus_1[, ds_idx, drop = FALSE]) -
      rowMeans(log2_cpm_plus_1[, control_idx, drop = FALSE]),
    mean_ds_cpm = rowMeans(cpm_values[, ds_idx, drop = FALSE]),
    mean_control_cpm = rowMeans(cpm_values[, control_idx, drop = FALSE]),
    ds_detected_fraction = rowMeans(counts[, ds_idx, drop = FALSE] > 0),
    control_detected_fraction = rowMeans(counts[, control_idx, drop = FALSE] > 0),
    tested_by_edgeR = FALSE,
    logFC = NA_real_,
    logCPM = NA_real_,
    F = NA_real_,
    PValue = NA_real_,
    FDR = NA_real_,
    stringsAsFactors = FALSE
  )

  m <- match(result$gene_id, de$gene_id)
  hit <- !is.na(m)
  result$tested_by_edgeR[hit] <- TRUE
  result$logFC[hit] <- de$logFC[m[hit]]
  result$logCPM[hit] <- de$logCPM[m[hit]]
  result$F[hit] <- de$F[m[hit]]
  result$PValue[hit] <- de$PValue[m[hit]]
  result$FDR[hit] <- de$FDR[m[hit]]
  result$effect_log2fc_for_claim <- ifelse(result$tested_by_edgeR, result$logFC, result$delta_ds_minus_control_log2_cpm_plus_1)
  result$effect_source <- ifelse(result$tested_by_edgeR, "edgeR_QLF", "mean_log2CPM_low_expression_not_tested")
  result$direction <- ifelse(
    is.na(result$effect_log2fc_for_claim) | abs(result$effect_log2fc_for_claim) < 0.25,
    "flat",
    ifelse(result$effect_log2fc_for_claim > 0, "up", "down")
  )
  result <- result[order(result$FDR, -abs(result$logFC), na.last = TRUE), ]
  rownames(result) <- NULL

  list(result = result, sample_manifest = meta)
}

write_table <- function(df, path) {
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

gencode <- read_gencode_symbols(gencode_gtf)
gse301886 <- read_gse301886_counts(gse301886_counts, gencode)
gse305139 <- read_gse305139_counts(gse305139_counts)
meta301886 <- metadata_gse301886(gse301886_meta)
meta305139 <- metadata_gse305139(gse305139_meta)

jobs <- list(
  list(dataset = "GSE301886", context = "DLPFC", data = gse301886, meta = meta301886),
  list(dataset = "GSE301886", context = "Hippocampus", data = gse301886, meta = meta301886),
  list(dataset = "GSE305139", context = "fetal_brain_final", data = gse305139, meta = meta305139)
)

all_results <- list()
all_samples <- list()
summary_rows <- list()
for (job in jobs) {
  fit <- fit_edgeR_context(job$dataset, job$context, job$data$counts, job$data$gene_info, job$meta)
  result <- fit$result
  file_context <- gsub("[^A-Za-z0-9]+", "_", job$context)
  out_file <- file.path(outdir, sprintf("%s_%s_edgeR_de.tsv", job$dataset, file_context))
  write_table(result, out_file)
  all_results[[length(all_results) + 1]] <- result
  samples <- fit$sample_manifest
  samples$context <- job$context
  all_samples[[length(all_samples) + 1]] <- samples
  tested <- result[result$tested_by_edgeR, , drop = FALSE]
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    dataset = job$dataset,
    context = job$context,
    n_ds = unique(result$n_ds)[[1]],
    n_control = unique(result$n_control)[[1]],
    model_formula = unique(result$model_formula)[[1]],
    total_gene_rows = nrow(result),
    edgeR_tested_gene_rows = nrow(tested),
    fdr05_genes = sum(tested$FDR <= 0.05, na.rm = TRUE),
    fdr05_abslogfc058_genes = sum(tested$FDR <= 0.05 & abs(tested$logFC) >= 0.58, na.rm = TRUE),
    fdr05_abslogfc1_genes = sum(tested$FDR <= 0.05 & abs(tested$logFC) >= 1, na.rm = TRUE),
    output_file = out_file,
    stringsAsFactors = FALSE
  )
}

combined <- do.call(rbind, all_results)
samples <- do.call(rbind, all_samples)
summary <- do.call(rbind, summary_rows)
write_table(combined, file.path(outdir, "human_ds_bulk_count_de_all_contexts.tsv"))
write_table(samples, file.path(outdir, "human_ds_bulk_count_de_sample_manifest.tsv"))
write_table(summary, file.path(outdir, "human_ds_bulk_count_de_summary.tsv"))

message("Wrote count-based human DS bulk DE results to ", outdir)
