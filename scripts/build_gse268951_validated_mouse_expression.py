#!/usr/bin/env python3
"""Build audit-only GSE268951 transcript and gene-level expression tables.

This legacy helper preserves the local CON-VS-OE edgeR-style table comparison.
Its outputs are not manuscript-facing because the local table is not a public
GEO supplementary file. Use analyze_gse268951_count_de.R for manuscript-facing
public expected-count DESeq2/edgeR results.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_USER_TABLE = "data/geo/GSE268951/supplementary/GSE268951_CON-VS-OE_DE_anno.xls"
DEFAULT_PRESERVED_COMPARISON = "results/mouse_gse268951_count_de/gse268951_correct_table_vs_public_transcript_edger_comparison.tsv"
DEFAULT_RAW = "data/geo/GSE268951/supplementary/GSE268951_all.fpkm_mRNA_anno.txt.gz"
DEFAULT_OUTDIR = "results/mouse_gse268951_count_de"

SAMPLES = ["A-1", "A-2", "A-3", "B-1", "B-2", "B-3", "C-1", "C-2", "C-3", "D-1", "D-2", "D-3"]
CONTROL = ["A-1", "A-2", "A-3", "B-1", "B-2", "B-3"]
OE = ["C-1", "C-2", "C-3", "D-1", "D-2", "D-3"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--user-table", default=DEFAULT_USER_TABLE)
    parser.add_argument("--preserved-comparison", default=DEFAULT_PRESERVED_COMPARISON)
    parser.add_argument("--raw", default=DEFAULT_RAW)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    return parser.parse_args()


def safe_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def load_primary_stats(user_table: Path, preserved_comparison: Path) -> tuple[pd.DataFrame, str]:
    if user_table.exists():
        table = pd.read_csv(user_table, sep="\t", low_memory=False, na_values=["NA", ""])
        stats = table.rename(columns={"geneID": "gene_id_gene", "gene_id": "transcript_id"}).copy()
        stats = stats[["transcript_id", "logFC", "logCPM", "PValue", "FDR", "gene_id_gene"]]
        return stats, "user_CON-VS-OE_edgeR_transcript_table"

    comp = pd.read_csv(preserved_comparison, sep="\t", low_memory=False, na_values=["NA", ""])
    stats = pd.DataFrame(
        {
            "transcript_id": comp["gene_id"],
            "logFC": comp["logFC_user"],
            "logCPM": comp["logCPM_user"],
            "PValue": comp["PValue_user"],
            "FDR": comp["FDR_user"],
            "gene_id_gene": comp["geneID_user"],
        }
    )
    return stats, "preserved_user_CON-VS-OE_edgeR_transcript_stats"


def fmt(value: object, digits: int = 6) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return ""
    if not math.isfinite(number):
        return ""
    return f"{number:.{digits}g}"


def direction(value: float, threshold: float = 0.25) -> str:
    if not math.isfinite(value):
        return ""
    if value >= threshold:
        return "up"
    if value <= -threshold:
        return "down"
    return "flat"


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    raw = pd.read_csv(args.raw, sep="\t", low_memory=False)
    raw = raw.rename(columns={"gene_id": "transcript_id", "GeneID": "gene_id"})

    primary, source = load_primary_stats(Path(args.user_table), Path(args.preserved_comparison))
    primary = primary.rename(columns={"gene_id_gene": "gene_id_from_primary"})
    for col in ["logFC", "logCPM", "PValue", "FDR"]:
        primary[col] = safe_numeric(primary[col])

    count_cols = [f"{sample}_expected_count" for sample in SAMPLES]
    fpkm_cols = [f"{sample}_FPKM" for sample in SAMPLES]
    for col in count_cols + fpkm_cols:
        raw[col] = safe_numeric(raw[col]).fillna(0.0)

    merged = raw.merge(primary, on="transcript_id", how="left")
    merged["gene_id"] = merged["gene_id"].fillna(merged["gene_id_from_primary"])
    merged = merged[
        merged["gene_id"].notna()
        & merged["gene_id"].ne("")
        & merged["gene_id"].ne("-")
        & merged["GeneSymbol"].notna()
        & merged["GeneSymbol"].ne("")
        & merged["GeneSymbol"].ne("-")
    ].copy()
    merged["control_mean_count"] = merged[[f"{sample}_expected_count" for sample in CONTROL]].mean(axis=1)
    merged["dnmt3l_oe_mean_count"] = merged[[f"{sample}_expected_count" for sample in OE]].mean(axis=1)
    merged["control_mean_fpkm"] = merged[[f"{sample}_FPKM" for sample in CONTROL]].mean(axis=1)
    merged["dnmt3l_oe_mean_fpkm"] = merged[[f"{sample}_FPKM" for sample in OE]].mean(axis=1)
    merged["fpkm_log2fc"] = np.log2((merged["dnmt3l_oe_mean_fpkm"] + 0.1) / (merged["control_mean_fpkm"] + 0.1))
    merged["gse268951_effect_source"] = source
    merged["gse268951_de_p05_fc2"] = (merged["PValue"].le(0.05) & merged["logFC"].abs().ge(1.0))
    merged["gse268951_de_fdr05_fc2"] = (merged["FDR"].le(0.05) & merged["logFC"].abs().ge(1.0))
    merged["gse268951_direction"] = merged["logFC"].map(direction)

    transcript_cols = [
        "transcript_id",
        "gene_id",
        "GeneSymbol",
        "Description",
        "Chr",
        "Start",
        "End",
        "Strand",
        "logFC",
        "logCPM",
        "PValue",
        "FDR",
        "control_mean_count",
        "dnmt3l_oe_mean_count",
        "control_mean_fpkm",
        "dnmt3l_oe_mean_fpkm",
        "fpkm_log2fc",
        "gse268951_de_p05_fc2",
        "gse268951_de_fdr05_fc2",
        "gse268951_direction",
        "gse268951_effect_source",
        "KO",
        "GO_BP",
        "GO_CC",
        "GO_MF",
    ]
    transcript = merged[transcript_cols].copy()
    transcript.to_csv(outdir / "gse268951_validated_transcript_edgeR_de.tsv", sep="\t", index=False)

    gene_sort = merged.copy()
    gene_sort["FDR_sort"] = gene_sort["FDR"].fillna(1.0)
    gene_sort["PValue_sort"] = gene_sort["PValue"].fillna(1.0)
    gene_sort["abs_logFC_sort"] = gene_sort["logFC"].abs().fillna(-1.0)
    gene_sort["logCPM_sort"] = gene_sort["logCPM"].fillna(-999.0)
    gene_sort = gene_sort.sort_values(
        ["gene_id", "FDR_sort", "PValue_sort", "abs_logFC_sort", "logCPM_sort"],
        ascending=[True, True, True, False, False],
    )
    gene = gene_sort.drop_duplicates("gene_id", keep="first").copy()
    tx_counts = merged.groupby("gene_id")["transcript_id"].nunique().rename("transcript_count_for_gene")
    gene = gene.merge(tx_counts, on="gene_id", how="left")

    gene_out = pd.DataFrame(
        {
            "gene_id": gene["gene_id"],
            "representative_transcript_id": gene["transcript_id"],
            "transcript_count_for_gene": gene["transcript_count_for_gene"],
            "gene_symbol": gene["GeneSymbol"],
            "chrom": gene["Chr"],
            "start": gene["Start"],
            "end": gene["End"],
            "strand": gene["Strand"],
            "description": gene["Description"],
            "gse268951_control_mean_count": gene["control_mean_count"].map(fmt),
            "gse268951_dnmt3l_oe_mean_count": gene["dnmt3l_oe_mean_count"].map(fmt),
            "gse268951_control_mean_fpkm": gene["control_mean_fpkm"].map(fmt),
            "gse268951_dnmt3l_oe_mean_fpkm": gene["dnmt3l_oe_mean_fpkm"].map(fmt),
            "gse268951_log2fc_oe_vs_control": gene["logFC"].map(fmt),
            "gse268951_primary_pvalue": gene["PValue"].map(fmt),
            "gse268951_primary_fdr": gene["FDR"].map(fmt),
            "gse268951_edgeR_logFC": gene["logFC"].map(fmt),
            "gse268951_edgeR_logCPM": gene["logCPM"].map(fmt),
            "gse268951_edgeR_pvalue": gene["PValue"].map(fmt),
            "gse268951_edgeR_fdr": gene["FDR"].map(fmt),
            "gse268951_fpkm_log2fc": gene["fpkm_log2fc"].map(fmt),
            "gse268951_effect_source": gene["gse268951_effect_source"],
            "gse268951_direction": gene["gse268951_direction"],
            "gse268951_de_p05_fc2": np.where(gene["gse268951_de_p05_fc2"], "yes", "no"),
            "gse268951_de_fdr05_fc2": np.where(gene["gse268951_de_fdr05_fc2"], "yes", "no"),
            "go_bp_terms": gene["GO_BP"].fillna("").astype(str),
            "go_cc_terms": gene["GO_CC"].fillna("").astype(str),
            "go_mf_terms": gene["GO_MF"].fillna("").astype(str),
            "KO": gene["KO"].fillna("").astype(str),
        }
    )
    gene_out.to_csv(outdir / "gse268951_validated_gene_edgeR_best_transcript.tsv", sep="\t", index=False)

    targets = ["Dnmt3l", "Dnmt3a", "Dnmt3b", "Stat1", "Stat3", "Hes3", "Ascl1", "Neurog2", "Neurod2", "Cdkn1a"]
    target_out = gene_out[gene_out["gene_symbol"].isin(targets)].copy()
    target_out["target_order"] = target_out["gene_symbol"].map({gene: i for i, gene in enumerate(targets)})
    target_out = target_out.sort_values("target_order").drop(columns=["target_order"])
    target_out.to_csv(outdir / "gse268951_validated_target_marker_gene_effects.tsv", sep="\t", index=False)

    summary = pd.DataFrame(
        [
            ("primary_source", source),
            ("transcript_rows", len(transcript)),
            ("gene_rows", len(gene_out)),
            ("transcript_fdr05_abslogfc1", int(merged["gse268951_de_fdr05_fc2"].sum())),
            ("gene_best_transcript_fdr05_abslogfc1", int((gene_out["gse268951_de_fdr05_fc2"] == "yes").sum())),
            ("gene_best_transcript_p05_abslogfc1", int((gene_out["gse268951_de_p05_fc2"] == "yes").sum())),
        ],
        columns=["metric", "value"],
    )
    summary.to_csv(outdir / "gse268951_validated_expression_summary.tsv", sep="\t", index=False)

    print(f"Wrote {outdir / 'gse268951_validated_transcript_edgeR_de.tsv'}")
    print(f"Wrote {outdir / 'gse268951_validated_gene_edgeR_best_transcript.tsv'}")
    print(f"Wrote {outdir / 'gse268951_validated_target_marker_gene_effects.tsv'}")
    print(f"Wrote {outdir / 'gse268951_validated_expression_summary.tsv'}")


if __name__ == "__main__":
    main()
