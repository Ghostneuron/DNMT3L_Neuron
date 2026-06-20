#!/usr/bin/env python3
"""Legacy comparison of GSE268951 public count-DE results with an older external table.

This script is audit-only. Current manuscript analyses should use public
processed expected counts via analyze_gse268951_count_de.R, not local or
external CON-VS-OE tables.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_USER_DE = "data/geo/GSE268951/supplementary/GSE268951_legacy_Group_Control-VS-Dnmt3l-overexpression_DE_anno.xls"
DEFAULT_OUR_DE = "results/mouse_gse268951_count_de/gse268951_deseq2_edger_count_de_results.tsv"
DEFAULT_OUTDIR = "results/mouse_gse268951_count_de"
DEFAULT_NOTE = "notes/gse268951_legacy_wrong_table_comparison.md"

TARGETS = ["Dnmt3l", "Dnmt3a", "Dnmt3b", "Stat1", "Stat3", "Hes3", "Ascl1", "Neurog2", "Neurod2", "Cdkn1a"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--user-de", default=DEFAULT_USER_DE)
    parser.add_argument("--our-de", default=DEFAULT_OUR_DE)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--note", default=DEFAULT_NOTE)
    return parser.parse_args()


def safe_float(value: object) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return math.nan
    return number if math.isfinite(number) else math.nan


def fmt(value: object, digits: int = 6) -> str:
    number = safe_float(value)
    if not math.isfinite(number):
        return ""
    return f"{number:.{digits}g}"


def sig_flag(log2fc: pd.Series, padj: pd.Series, lfc: float = 1.0) -> pd.Series:
    return padj.le(0.05) & log2fc.abs().ge(lfc)


def direction(values: pd.Series, threshold: float = 0.0) -> pd.Series:
    return np.where(values > threshold, "up", np.where(values < -threshold, "down", "flat"))


def pick_best_user_rows(user: pd.DataFrame) -> pd.DataFrame:
    user = user.copy()
    user["padj_for_sort"] = user["padj"].fillna(1.0)
    user["baseMean_for_sort"] = user["baseMean"].fillna(-1.0)
    user["abs_lfc_for_sort"] = user["log2FoldChange"].abs().fillna(-1.0)
    user = user.sort_values(
        ["gene_id", "padj_for_sort", "baseMean_for_sort", "abs_lfc_for_sort"],
        ascending=[True, True, False, False],
    )
    return user.drop_duplicates("gene_id", keep="first")


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    user = pd.read_csv(args.user_de, sep="\t", low_memory=False, na_values=["NA", ""])
    ours = pd.read_csv(args.our_de, sep="\t", low_memory=False, na_values=["NA", ""])

    numeric_user = ["baseMean", "foldchange", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
    for col in numeric_user:
        user[col] = pd.to_numeric(user[col], errors="coerce")
    numeric_ours = [
        "deseq2_baseMean",
        "deseq2_log2fc",
        "deseq2_lfcSE",
        "deseq2_stat",
        "deseq2_pvalue",
        "deseq2_padj",
        "control_mean_count",
        "dnmt3l_oe_mean_count",
        "control_mean_fpkm",
        "dnmt3l_oe_mean_fpkm",
        "fpkm_log2fc",
        "edger_log2fc",
        "edger_fdr",
    ]
    for col in numeric_ours:
        ours[col] = pd.to_numeric(ours[col], errors="coerce")

    user_best = pick_best_user_rows(user)
    merged = user_best.merge(
        ours,
        on="gene_id",
        how="inner",
        suffixes=("_user", "_ours"),
    )

    merged["user_sig_fdr05_abslog2fc1"] = sig_flag(merged["log2FoldChange"], merged["padj"])
    merged["ours_sig_fdr05_abslog2fc1"] = sig_flag(merged["deseq2_log2fc"], merged["deseq2_padj"])
    merged["user_direction"] = direction(merged["log2FoldChange"])
    merged["ours_direction"] = direction(merged["deseq2_log2fc"])
    finite = merged["log2FoldChange"].notna() & merged["deseq2_log2fc"].notna()
    merged["delta_ours_minus_user_log2fc"] = merged["deseq2_log2fc"] - merged["log2FoldChange"]
    merged["abs_delta_log2fc"] = merged["delta_ours_minus_user_log2fc"].abs()

    comparison_cols = [
        "gene_id",
        "GeneSymbol_user",
        "GeneSymbol_ours",
        "baseMean",
        "deseq2_baseMean",
        "log2FoldChange",
        "deseq2_log2fc",
        "delta_ours_minus_user_log2fc",
        "abs_delta_log2fc",
        "pvalue",
        "deseq2_pvalue",
        "padj",
        "deseq2_padj",
        "user_sig_fdr05_abslog2fc1",
        "ours_sig_fdr05_abslog2fc1",
        "user_direction",
        "ours_direction",
        "control_mean_count",
        "dnmt3l_oe_mean_count",
        "control_mean_fpkm",
        "dnmt3l_oe_mean_fpkm",
        "fpkm_log2fc",
    ]
    comparison = merged[comparison_cols].sort_values("abs_delta_log2fc", ascending=False)
    comparison.to_csv(outdir / "gse268951_user_vs_rederived_count_de_comparison.tsv", sep="\t", index=False)

    target = comparison[
        comparison["GeneSymbol_user"].isin(TARGETS) | comparison["GeneSymbol_ours"].isin(TARGETS)
    ].copy()
    target["target_order"] = target["GeneSymbol_user"].map({gene: i for i, gene in enumerate(TARGETS)})
    target = target.sort_values(["target_order", "GeneSymbol_user"], na_position="last").drop(columns=["target_order"])
    target.to_csv(outdir / "gse268951_user_vs_rederived_target_marker_comparison.tsv", sep="\t", index=False)

    user_sig_ids = set(user_best.loc[sig_flag(user_best["log2FoldChange"], user_best["padj"]), "gene_id"])
    our_sig_ids = set(ours.loc[sig_flag(ours["deseq2_log2fc"], ours["deseq2_padj"]), "gene_id"])
    overlap_sig = user_sig_ids & our_sig_ids
    union_sig = user_sig_ids | our_sig_ids

    finite_merged = merged.loc[finite].copy()
    pearson = finite_merged["log2FoldChange"].corr(finite_merged["deseq2_log2fc"], method="pearson")
    spearman = finite_merged["log2FoldChange"].rank().corr(finite_merged["deseq2_log2fc"].rank(), method="pearson")
    sign_agree = (finite_merged["user_direction"] == finite_merged["ours_direction"]).mean()

    summary_rows = [
        ("user_rows", len(user)),
        ("user_unique_gene_ids", user["gene_id"].nunique()),
        ("user_duplicate_gene_id_rows", len(user) - user["gene_id"].nunique()),
        ("our_rows", len(ours)),
        ("our_unique_gene_ids", ours["gene_id"].nunique()),
        ("intersect_gene_ids_after_user_dedup", len(merged)),
        ("finite_log2fc_overlap_genes", int(finite.sum())),
        ("pearson_log2fc_overlap", pearson),
        ("spearman_log2fc_overlap", spearman),
        ("direction_agreement_fraction", sign_agree),
        ("user_fdr05_abslog2fc1_genes", len(user_sig_ids)),
        ("ours_fdr05_abslog2fc1_genes", len(our_sig_ids)),
        ("shared_fdr05_abslog2fc1_genes", len(overlap_sig)),
        ("union_fdr05_abslog2fc1_genes", len(union_sig)),
        ("jaccard_fdr05_abslog2fc1", len(overlap_sig) / len(union_sig) if union_sig else math.nan),
        ("user_only_fdr05_abslog2fc1_genes", len(user_sig_ids - our_sig_ids)),
        ("ours_only_fdr05_abslog2fc1_genes", len(our_sig_ids - user_sig_ids)),
        ("median_abs_log2fc_delta_overlap", finite_merged["abs_delta_log2fc"].median()),
        ("genes_abs_log2fc_delta_gt_1", int((finite_merged["abs_delta_log2fc"] > 1).sum())),
    ]
    summary = pd.DataFrame(summary_rows, columns=["metric", "value"])
    summary.to_csv(outdir / "gse268951_user_vs_rederived_count_de_summary.tsv", sep="\t", index=False)

    note = [
        "# GSE268951 User DE Table Comparison",
        "",
        "Date: 2026-06-11",
        "",
        f"User table: `{args.user_de}`",
        f"Redervied table: `{args.our_de}`",
        "",
        "## Summary",
        "",
    ]
    for metric, value in summary_rows:
        note.append(f"- {metric}: {fmt(value) if isinstance(value, float) else value}")
    note.extend(
        [
            "",
            "## Interpretation",
            "",
            "The user-provided table is a DESeq2-style gene annotation result and should be treated as the preferred count-based result if it is from the original analysis/report. The public GEO `expected_count` columns do not reproduce it closely in this first reanalysis.",
            "",
            "The strongest practical consequence is that downstream manuscript analyses should use the user-provided DE table as the primary mouse `GSE268951` expression effect source, while retaining the public processed-count reanalysis only as a sensitivity/audit table.",
            "",
            "## Outputs",
            "",
            "- `results/mouse_gse268951_count_de/gse268951_user_vs_rederived_count_de_summary.tsv`",
            "- `results/mouse_gse268951_count_de/gse268951_user_vs_rederived_count_de_comparison.tsv`",
            "- `results/mouse_gse268951_count_de/gse268951_user_vs_rederived_target_marker_comparison.tsv`",
        ]
    )
    Path(args.note).write_text("\n".join(note) + "\n")

    print(f"Wrote {outdir / 'gse268951_user_vs_rederived_count_de_summary.tsv'}")
    print(f"Wrote {outdir / 'gse268951_user_vs_rederived_count_de_comparison.tsv'}")
    print(f"Wrote {outdir / 'gse268951_user_vs_rederived_target_marker_comparison.tsv'}")
    print(f"Wrote {args.note}")


if __name__ == "__main__":
    main()
