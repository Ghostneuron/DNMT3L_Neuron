#!/usr/bin/env python3
"""Integrate GSE305146 peak-level DA betas with candidate gene windows."""

from __future__ import annotations

import argparse
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analyze_gse305146_marker_atac_motif_deviation import CLASS_LABELS, CLASS_ORDER
from analyze_gse305146_marker_atac_motif_instance_peaksets import FAMILY_ORDER, MOTIF_FAMILY_LABELS
from infer_gse305146_marker_cell_classes_and_targets import CORE_GENES


MECHANISM_LABELS = {
    "cross_species_dnmt3a_linked_methylation_plus_developmental_cpg_tf_context": "DNMT3A methyl. + dev/CpG",
    "cross_species_dnmt3a_linked_methylation_plus_developmental_tf_context": "DNMT3A methyl. + dev TF",
    "cross_species_dnmt3a_linked_methylation_without_selected_tf_instance": "DNMT3A methyl. only",
    "mouse_dnmt3a_linked_methylation_with_tf_context_human_unresolved": "Mouse methyl. + TF; human unresolved",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gene-windows",
        default="results/mechanism_question/gse305146_marker_atac_window_gene_windows.tsv",
    )
    parser.add_argument(
        "--candidate-axis",
        default="results/mechanism_question/final_candidate_mechanism_axis_summary.tsv",
    )
    parser.add_argument(
        "--tested-regions",
        default="results/mechanism_question/gse305146_marker_atac_peak_da_tested_regions.tsv",
    )
    parser.add_argument(
        "--peak-effects",
        default="results/mechanism_question/gse305146_marker_atac_peak_da_effects.tsv.gz",
    )
    parser.add_argument(
        "--motif-hits",
        default="results/mechanism_question/gse305146_marker_atac_motif_deviation_region_hits.tsv",
    )
    parser.add_argument(
        "--ranked-motif",
        default="results/mechanism_question/gse305146_peak_da_ranked_motif_shift.tsv",
    )
    parser.add_argument("--outdir", default="results/mechanism_question")
    parser.add_argument("--figdir", default="figures/mechanism_question")
    parser.add_argument("--note", default="notes/gse305146_candidate_peak_da_integration.md")
    parser.add_argument("--class-beta-threshold", type=float, default=0.10)
    parser.add_argument("--local-motif-beta-threshold", type=float, default=0.02)
    return parser.parse_args()


def fmt_float(value: float) -> str:
    return f"{value:.3g}" if math.isfinite(value) else "NA"


def load_gene_metadata(candidate_axis_path: Path) -> tuple[pd.DataFrame, list[str]]:
    candidate = pd.read_csv(candidate_axis_path, sep="\t")
    rows = []
    for gene in CORE_GENES:
        rows.append(
            {
                "gene_symbol": gene,
                "gene_role": "core_target",
                "tier": "",
                "candidate_group": "",
                "final_mechanism_axis_call": "",
                "final_mechanism_axis_label": "",
                "short_interpretation": "",
                "gene_order": len(rows),
            }
        )
    for _, row in candidate.iterrows():
        rows.append(
            {
                "gene_symbol": row["gene_symbol"],
                "gene_role": "final_candidate",
                "tier": row.get("tier", ""),
                "candidate_group": row.get("candidate_group", ""),
                "final_mechanism_axis_call": row.get("final_mechanism_axis_call", ""),
                "final_mechanism_axis_label": MECHANISM_LABELS.get(
                    row.get("final_mechanism_axis_call", ""),
                    row.get("final_mechanism_axis_call", ""),
                ),
                "short_interpretation": row.get("short_interpretation", ""),
                "gene_order": len(rows),
            }
        )
    meta = pd.DataFrame(rows)
    order = meta["gene_symbol"].tolist()
    return meta, order


def collapse_windows(windows_path: Path, gene_meta: pd.DataFrame) -> pd.DataFrame:
    windows = pd.read_csv(windows_path, sep="\t")
    grouped = (
        windows.groupby(["gene_symbol", "chrom"], dropna=False, observed=True)
        .agg(
            gene_start=("gene_start", "min"),
            gene_end=("gene_end", "max"),
            window_start=("window_start", "min"),
            window_end=("window_end", "max"),
            mean_window_peak_features=("n_window_peak_features", "mean"),
            max_window_peak_features=("n_window_peak_features", "max"),
            n_libraries=("library_prefix", "nunique"),
        )
        .reset_index()
    )
    grouped = grouped.merge(gene_meta, on="gene_symbol", how="left")
    grouped["gene_order"] = grouped["gene_order"].fillna(999).astype(int)
    grouped = grouped.sort_values(["gene_order", "gene_symbol"])
    return grouped


def overlap_peaks_to_genes(regions: pd.DataFrame, windows: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for _, gene in windows.iterrows():
        sub = regions[
            (regions["chrom"] == gene["chrom"])
            & (regions["end"] >= int(gene["window_start"]))
            & (regions["start"] <= int(gene["window_end"]))
        ].copy()
        for _, peak in sub.iterrows():
            rows.append(
                {
                    "gene_symbol": gene["gene_symbol"],
                    "gene_role": gene.get("gene_role", ""),
                    "tier": gene.get("tier", ""),
                    "candidate_group": gene.get("candidate_group", ""),
                    "final_mechanism_axis_call": gene.get("final_mechanism_axis_call", ""),
                    "final_mechanism_axis_label": gene.get("final_mechanism_axis_label", ""),
                    "chrom": gene["chrom"],
                    "gene_window_start": int(gene["window_start"]),
                    "gene_window_end": int(gene["window_end"]),
                    "peak_id": peak["peak_id"],
                    "peak_chrom": peak["chrom"],
                    "peak_start": int(peak["start"]),
                    "peak_end": int(peak["end"]),
                    "background_rank": int(peak["background_rank"]),
                    "peak_delta_frequency_ds_minus_control": float(peak["delta_frequency_ds_minus_control"]),
                    "peak_max_frequency": float(peak["max_frequency"]),
                }
            )
    return pd.DataFrame(rows)


def read_peak_motif_hits(motif_hits_path: Path, local_peak_ids: set[str]) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    usecols = ["peak_id", "motif_family", "motif_family_label"]
    for chunk in pd.read_csv(motif_hits_path, sep="\t", usecols=usecols, chunksize=200_000):
        chunk = chunk[chunk["peak_id"].isin(local_peak_ids)].drop_duplicates()
        if not chunk.empty:
            rows.append(chunk)
    if not rows:
        return pd.DataFrame(columns=usecols)
    return pd.concat(rows, ignore_index=True).drop_duplicates()


def load_local_effects(peak_effects_path: Path, peak_gene_map: pd.DataFrame) -> pd.DataFrame:
    local_peak_ids = set(peak_gene_map["peak_id"])
    usecols = [
        "inferred_cell_class",
        "peak_id",
        "n_ds_donors",
        "n_control_donors",
        "n_ds_cells",
        "n_control_cells",
        "mean_ds_log2_peak_cpm",
        "mean_control_log2_peak_cpm",
        "delta_mean_ds_minus_control_log2_peak_cpm",
        "pcw_adjusted_group_beta",
        "pcw_adjusted_group_p",
        "fdr_bh",
        "significant_da",
        "direction",
    ]
    frames: list[pd.DataFrame] = []
    for chunk in pd.read_csv(peak_effects_path, sep="\t", usecols=usecols, chunksize=200_000):
        chunk = chunk[chunk["peak_id"].isin(local_peak_ids)]
        if not chunk.empty:
            frames.append(chunk)
    if not frames:
        return pd.DataFrame(columns=usecols)
    effects = pd.concat(frames, ignore_index=True)
    effects = effects.merge(peak_gene_map, on="peak_id", how="left")
    effects["significant_da"] = effects["significant_da"].astype(str).str.lower().isin({"true", "1"})
    return effects


def summarize_gene_class(local: pd.DataFrame, threshold: float) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (gene, cell_class), sub in local.groupby(["gene_symbol", "inferred_cell_class"], dropna=False, observed=True):
        beta = pd.to_numeric(sub["pcw_adjusted_group_beta"], errors="coerce")
        delta = pd.to_numeric(sub["delta_mean_ds_minus_control_log2_peak_cpm"], errors="coerce")
        sig = sub[sub["significant_da"]]
        top_pos = sub.sort_values("pcw_adjusted_group_beta", ascending=False).iloc[0]
        top_neg = sub.sort_values("pcw_adjusted_group_beta", ascending=True).iloc[0]
        rows.append(
            {
                "gene_symbol": gene,
                "inferred_cell_class": cell_class,
                "gene_role": sub["gene_role"].iloc[0],
                "tier": sub["tier"].iloc[0],
                "candidate_group": sub["candidate_group"].iloc[0],
                "final_mechanism_axis_call": sub["final_mechanism_axis_call"].iloc[0],
                "final_mechanism_axis_label": sub["final_mechanism_axis_label"].iloc[0],
                "chrom": sub["chrom"].iloc[0],
                "gene_window_start": int(sub["gene_window_start"].iloc[0]),
                "gene_window_end": int(sub["gene_window_end"].iloc[0]),
                "n_local_recurrent_peaks": int(sub["peak_id"].nunique()),
                "mean_pcw_adjusted_beta": float(beta.mean()),
                "median_pcw_adjusted_beta": float(beta.median()),
                "mean_abs_pcw_adjusted_beta": float(beta.abs().mean()),
                "mean_delta_log2_peak_cpm": float(delta.mean()),
                "n_beta_ge_threshold": int((beta >= threshold).sum()),
                "n_beta_le_negative_threshold": int((beta <= -threshold).sum()),
                "n_nominal_p_lt_0_05": int((pd.to_numeric(sub["pcw_adjusted_group_p"], errors="coerce") < 0.05).sum()),
                "n_strict_da": int(sig["peak_id"].nunique()),
                "n_strict_da_ds_open": int((sig["direction"] == "DS_open").sum()),
                "n_strict_da_ds_closed": int((sig["direction"] == "DS_closed").sum()),
                "top_positive_peak_id": top_pos["peak_id"],
                "top_positive_beta": float(top_pos["pcw_adjusted_group_beta"]),
                "top_negative_peak_id": top_neg["peak_id"],
                "top_negative_beta": float(top_neg["pcw_adjusted_group_beta"]),
                "n_ds_donors": int(sub["n_ds_donors"].max()),
                "n_control_donors": int(sub["n_control_donors"].max()),
            }
        )
    out = pd.DataFrame(rows)
    return out


def summarize_motif_class(
    local: pd.DataFrame,
    motif_hits: pd.DataFrame,
    ranked: pd.DataFrame,
    beta_threshold: float,
) -> pd.DataFrame:
    if motif_hits.empty:
        return pd.DataFrame()
    annotated = local.merge(motif_hits, on="peak_id", how="inner")
    rows: list[dict[str, object]] = []
    local_counts = (
        local.groupby(["gene_symbol", "inferred_cell_class"], dropna=False, observed=True)
        .agg(n_local_recurrent_peaks=("peak_id", "nunique"), mean_beta_all_local=("pcw_adjusted_group_beta", "mean"))
        .reset_index()
    )
    for (gene, cell_class, family), sub in annotated.groupby(
        ["gene_symbol", "inferred_cell_class", "motif_family"], dropna=False, observed=True
    ):
        beta = pd.to_numeric(sub["pcw_adjusted_group_beta"], errors="coerce")
        row = {
            "gene_symbol": gene,
            "inferred_cell_class": cell_class,
            "motif_family": family,
            "motif_family_label": MOTIF_FAMILY_LABELS.get(str(family), str(family)),
            "gene_role": sub["gene_role"].iloc[0],
            "tier": sub["tier"].iloc[0],
            "candidate_group": sub["candidate_group"].iloc[0],
            "final_mechanism_axis_call": sub["final_mechanism_axis_call"].iloc[0],
            "final_mechanism_axis_label": sub["final_mechanism_axis_label"].iloc[0],
            "n_local_motif_peaks": int(sub["peak_id"].nunique()),
            "mean_beta_local_motif_peaks": float(beta.mean()),
            "median_beta_local_motif_peaks": float(beta.median()),
            "n_local_motif_beta_ge_threshold": int((beta >= beta_threshold).sum()),
            "n_local_motif_beta_le_negative_threshold": int((beta <= -beta_threshold).sum()),
        }
        rows.append(row)
    out = pd.DataFrame(rows)
    out = out.merge(local_counts, on=["gene_symbol", "inferred_cell_class"], how="left")
    out["delta_local_motif_beta_minus_all_local"] = out["mean_beta_local_motif_peaks"] - out["mean_beta_all_local"]
    ranked_cols = [
        "inferred_cell_class",
        "motif_family",
        "delta_mean_beta_motif_minus_nonmotif",
        "rank_biserial",
        "fdr_bh",
        "shift_direction",
        "significant_ranked_shift",
    ]
    out = out.merge(ranked[ranked_cols], on=["inferred_cell_class", "motif_family"], how="left")
    out["local_global_sign_concordant"] = np.sign(out["mean_beta_local_motif_peaks"]) == np.sign(
        out["delta_mean_beta_motif_minus_nonmotif"]
    )
    out["local_ranked_context_support"] = (
        out["significant_ranked_shift"].astype(str).str.lower().isin({"true", "1"})
        & out["local_global_sign_concordant"]
        & (out["mean_beta_local_motif_peaks"].abs() >= beta_threshold)
    )
    return out


def summarize_gene(local_summary: pd.DataFrame, motif_summary: pd.DataFrame, gene_order: list[str], threshold: float) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    order_map = {gene: idx for idx, gene in enumerate(gene_order)}
    for gene, sub in local_summary.groupby("gene_symbol", dropna=False, observed=True):
        top_pos = sub.sort_values("mean_pcw_adjusted_beta", ascending=False).iloc[0]
        top_neg = sub.sort_values("mean_pcw_adjusted_beta", ascending=True).iloc[0]
        motif_sub = motif_summary[motif_summary["gene_symbol"] == gene] if not motif_summary.empty else pd.DataFrame()
        support = motif_sub[motif_sub.get("local_ranked_context_support", False)] if not motif_sub.empty else pd.DataFrame()
        pos_support = support[support["mean_beta_local_motif_peaks"] > 0] if not support.empty else pd.DataFrame()
        neg_support = support[support["mean_beta_local_motif_peaks"] < 0] if not support.empty else pd.DataFrame()
        rows.append(
            {
                "gene_symbol": gene,
                "gene_role": sub["gene_role"].iloc[0],
                "tier": sub["tier"].iloc[0],
                "candidate_group": sub["candidate_group"].iloc[0],
                "final_mechanism_axis_call": sub["final_mechanism_axis_call"].iloc[0],
                "final_mechanism_axis_label": sub["final_mechanism_axis_label"].iloc[0],
                "chrom": sub["chrom"].iloc[0],
                "n_local_recurrent_peaks": int(sub["n_local_recurrent_peaks"].max()),
                "n_classes_tested": int(sub["inferred_cell_class"].nunique()),
                "mean_beta_across_classes": float(sub["mean_pcw_adjusted_beta"].mean()),
                "n_classes_mean_beta_ge_threshold": int((sub["mean_pcw_adjusted_beta"] >= threshold).sum()),
                "n_classes_mean_beta_le_negative_threshold": int((sub["mean_pcw_adjusted_beta"] <= -threshold).sum()),
                "max_positive_class": top_pos["inferred_cell_class"],
                "max_positive_mean_beta": float(top_pos["mean_pcw_adjusted_beta"]),
                "max_negative_class": top_neg["inferred_cell_class"],
                "max_negative_mean_beta": float(top_neg["mean_pcw_adjusted_beta"]),
                "n_strict_da_total": int(sub["n_strict_da"].sum()),
                "strict_da_classes": ";".join(sorted(sub.loc[sub["n_strict_da"] > 0, "inferred_cell_class"].unique())),
                "n_local_ranked_motif_context_rows": int(len(support)),
                "positive_ranked_motif_contexts": ";".join(
                    sorted(
                        {
                            f"{row.inferred_cell_class}:{row.motif_family_label}"
                            for row in pos_support.itertuples(index=False)
                        }
                    )
                ),
                "negative_ranked_motif_contexts": ";".join(
                    sorted(
                        {
                            f"{row.inferred_cell_class}:{row.motif_family_label}"
                            for row in neg_support.itertuples(index=False)
                        }
                    )
                ),
                "gene_order": order_map.get(gene, 999),
            }
        )
    out = pd.DataFrame(rows).sort_values(["gene_order", "gene_symbol"]).drop(columns=["gene_order"])
    return out


def add_missing_gene_rows(gene_summary: pd.DataFrame, windows: pd.DataFrame, gene_order: list[str]) -> pd.DataFrame:
    present = set(gene_summary["gene_symbol"]) if not gene_summary.empty else set()
    missing_rows: list[dict[str, object]] = []
    order_map = {gene: idx for idx, gene in enumerate(gene_order)}
    for _, row in windows.iterrows():
        gene = row["gene_symbol"]
        if gene in present:
            continue
        missing_rows.append(
            {
                "gene_symbol": gene,
                "gene_role": row.get("gene_role", ""),
                "tier": row.get("tier", ""),
                "candidate_group": row.get("candidate_group", ""),
                "final_mechanism_axis_call": row.get("final_mechanism_axis_call", ""),
                "final_mechanism_axis_label": row.get("final_mechanism_axis_label", ""),
                "chrom": row.get("chrom", ""),
                "n_local_recurrent_peaks": 0,
                "n_classes_tested": 0,
                "mean_beta_across_classes": math.nan,
                "n_classes_mean_beta_ge_threshold": 0,
                "n_classes_mean_beta_le_negative_threshold": 0,
                "max_positive_class": "",
                "max_positive_mean_beta": math.nan,
                "max_negative_class": "",
                "max_negative_mean_beta": math.nan,
                "n_strict_da_total": 0,
                "strict_da_classes": "",
                "n_local_ranked_motif_context_rows": 0,
                "positive_ranked_motif_contexts": "",
                "negative_ranked_motif_contexts": "",
                "gene_order": order_map.get(gene, 999),
            }
        )
    if missing_rows:
        gene_summary = pd.concat([gene_summary, pd.DataFrame(missing_rows)], ignore_index=True)
    gene_summary["gene_order"] = gene_summary["gene_symbol"].map(order_map).fillna(999).astype(int)
    return gene_summary.sort_values(["gene_order", "gene_symbol"]).drop(columns=["gene_order"])


def mechanism_axis_summary(local_summary: pd.DataFrame) -> pd.DataFrame:
    candidate = local_summary[local_summary["gene_role"] == "final_candidate"].copy()
    rows: list[dict[str, object]] = []
    for (axis, cell_class), sub in candidate.groupby(
        ["final_mechanism_axis_label", "inferred_cell_class"], dropna=False, observed=True
    ):
        rows.append(
            {
                "final_mechanism_axis_label": axis,
                "inferred_cell_class": cell_class,
                "n_genes": int(sub["gene_symbol"].nunique()),
                "mean_local_peak_beta": float(sub["mean_pcw_adjusted_beta"].mean()),
                "median_local_peak_beta": float(sub["mean_pcw_adjusted_beta"].median()),
                "n_genes_positive_mean_beta": int((sub["mean_pcw_adjusted_beta"] > 0).sum()),
                "n_genes_negative_mean_beta": int((sub["mean_pcw_adjusted_beta"] < 0).sum()),
            }
        )
    return pd.DataFrame(rows)


def plot_gene_class_heatmap(local_summary: pd.DataFrame, gene_order: list[str], figdir: Path) -> None:
    matrix = np.full((len(gene_order), len(CLASS_ORDER)), np.nan)
    lookup = {
        (row.gene_symbol, row.inferred_cell_class): float(row.mean_pcw_adjusted_beta)
        for row in local_summary.itertuples(index=False)
    }
    for i, gene in enumerate(gene_order):
        for j, cell_class in enumerate(CLASS_ORDER):
            matrix[i, j] = lookup.get((gene, cell_class), np.nan)
    finite = matrix[np.isfinite(matrix)]
    vmax = max(0.1, float(np.nanpercentile(np.abs(finite), 95))) if finite.size else 0.2
    vmax = min(vmax, 0.6)
    fig, ax = plt.subplots(figsize=(11, max(7.0, 0.28 * len(gene_order))))
    im = ax.imshow(matrix, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(np.arange(len(CLASS_ORDER)))
    ax.set_xticklabels([CLASS_LABELS.get(name, name) for name in CLASS_ORDER], rotation=45, ha="right")
    ax.set_yticks(np.arange(len(gene_order)))
    ax.set_yticklabels(gene_order)
    ax.set_title("Candidate-window recurrent peak DA betas")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Mean PCW-adjusted DS beta across local peaks")
    ax.set_xticks(np.arange(-0.5, len(CLASS_ORDER), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(gene_order), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.7)
    ax.tick_params(which="minor", bottom=False, left=False)
    fig.tight_layout()
    fig.savefig(figdir / "gse305146_candidate_peak_da_local_beta_heatmap.png", dpi=300)
    fig.savefig(figdir / "gse305146_candidate_peak_da_local_beta_heatmap.pdf")
    plt.close(fig)


def plot_mechanism_heatmap(axis_summary: pd.DataFrame, figdir: Path) -> None:
    if axis_summary.empty:
        return
    axes = [label for label in MECHANISM_LABELS.values() if label in set(axis_summary["final_mechanism_axis_label"])]
    matrix = np.full((len(axes), len(CLASS_ORDER)), np.nan)
    lookup = {
        (row.final_mechanism_axis_label, row.inferred_cell_class): float(row.mean_local_peak_beta)
        for row in axis_summary.itertuples(index=False)
    }
    for i, axis in enumerate(axes):
        for j, cell_class in enumerate(CLASS_ORDER):
            matrix[i, j] = lookup.get((axis, cell_class), np.nan)
    finite = matrix[np.isfinite(matrix)]
    vmax = max(0.05, float(np.nanpercentile(np.abs(finite), 95))) if finite.size else 0.1
    vmax = min(vmax, 0.4)
    fig, ax = plt.subplots(figsize=(10.5, 3.4))
    im = ax.imshow(matrix, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(np.arange(len(CLASS_ORDER)))
    ax.set_xticklabels([CLASS_LABELS.get(name, name) for name in CLASS_ORDER], rotation=45, ha="right")
    ax.set_yticks(np.arange(len(axes)))
    ax.set_yticklabels(axes)
    ax.set_title("Mechanism-axis local recurrent peak betas")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Mean local peak beta")
    ax.set_xticks(np.arange(-0.5, len(CLASS_ORDER), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(axes), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.7)
    ax.tick_params(which="minor", bottom=False, left=False)
    fig.tight_layout()
    fig.savefig(figdir / "gse305146_candidate_peak_da_mechanism_axis_heatmap.png", dpi=300)
    fig.savefig(figdir / "gse305146_candidate_peak_da_mechanism_axis_heatmap.pdf")
    plt.close(fig)


def write_note(path: Path, gene_summary: pd.DataFrame, local_summary: pd.DataFrame, motif_summary: pd.DataFrame) -> None:
    top_pos = gene_summary.sort_values("max_positive_mean_beta", ascending=False).head(10)
    top_neg = gene_summary.sort_values("max_negative_mean_beta", ascending=True).head(10)
    strict = gene_summary[gene_summary["n_strict_da_total"] > 0]
    motif_support = (
        motif_summary[motif_summary["local_ranked_context_support"]]
        if not motif_summary.empty and "local_ranked_context_support" in motif_summary
        else pd.DataFrame()
    )
    lines = [
        "# GSE305146 Candidate-Local Peak DA Integration",
        "",
        "Question: do final candidate/core gene windows contain recurrent peaks whose PCW-adjusted DS accessibility betas and motif-family context match the broader GSE305146 regulatory signals?",
        "",
        "## Inputs",
        "",
        "- Gene windows: core genes plus final 24 candidates, using the same gene body +/-100 kb windows from the marker-inferred ATAC target-window pass.",
        "- Peak effects: donor-level PCW-adjusted recurrent-peak DA betas from the strict peak DA model.",
        "- Motifs: HOMER motif-family hits from the recurrent-background motif-deviation scan and ranked motif beta shifts.",
        "",
        "## Local Strict DA",
        "",
    ]
    if strict.empty:
        lines.append("- No candidate/core gene windows contained strict FDR-significant local DA peaks.")
    else:
        for _, row in strict.iterrows():
            lines.append(
                f"- {row['gene_symbol']}: {int(row['n_strict_da_total'])} strict DA local peak-class rows; "
                f"classes {row['strict_da_classes'] or 'NA'}."
            )
    lines.extend(["", "## Strongest Positive Local Peak-Beta Contexts", ""])
    for _, row in top_pos.iterrows():
        lines.append(
            f"- {row['gene_symbol']} in {row['max_positive_class']}: max class mean beta "
            f"{fmt_float(float(row['max_positive_mean_beta']))}; local peaks {int(row['n_local_recurrent_peaks'])}."
        )
    lines.extend(["", "## Strongest Negative Local Peak-Beta Contexts", ""])
    for _, row in top_neg.iterrows():
        lines.append(
            f"- {row['gene_symbol']} in {row['max_negative_class']}: min class mean beta "
            f"{fmt_float(float(row['max_negative_mean_beta']))}; local peaks {int(row['n_local_recurrent_peaks'])}."
        )
    lines.extend(["", "## Ranked Motif Context At Candidate Windows", ""])
    if motif_support.empty:
        lines.append("- No local motif rows met ranked-shift concordance criteria.")
    else:
        counts = (
            motif_support.groupby(["motif_family_label"], observed=True)
            .agg(n_gene_class_rows=("gene_symbol", "size"), n_genes=("gene_symbol", "nunique"))
            .reset_index()
            .sort_values(["n_gene_class_rows", "n_genes"], ascending=False)
            .head(12)
        )
        for _, row in counts.iterrows():
            lines.append(
                f"- {row['motif_family_label']}: {int(row['n_gene_class_rows'])} gene-class rows across "
                f"{int(row['n_genes'])} genes."
            )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This analysis localizes the global GSE305146 peak-DA and motif-shift layers to candidate/core gene windows. It does not prove that a local peak regulates the nearby gene, nor that a TF directly binds there. It is a prioritization layer: local recurrent peaks with positive or negative DS betas and concordant ranked motif context strengthen candidate-specific regulatory plausibility.",
            "",
            "The result should be read together with methylation support, expression shifts, target-window accessibility, motif-deviation, and strict DA. In this pass, strict FDR peak DA was not found inside the candidate/core windows; autosomal candidate interpretation should therefore rely more on distributed local betas and motif context than on strict individual-peak FDR.",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    figdir = Path(args.figdir)
    outdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    gene_meta, gene_order = load_gene_metadata(Path(args.candidate_axis))
    windows = collapse_windows(Path(args.gene_windows), gene_meta)
    windows.to_csv(
        outdir / "gse305146_candidate_peak_da_gene_windows.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    regions = pd.read_csv(args.tested_regions, sep="\t")
    peak_gene_map = overlap_peaks_to_genes(regions, windows)
    peak_gene_map.to_csv(
        outdir / "gse305146_candidate_peak_da_peak_gene_map.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    motif_hits = read_peak_motif_hits(Path(args.motif_hits), set(peak_gene_map["peak_id"]))
    motif_hits.to_csv(
        outdir / "gse305146_candidate_peak_da_local_motif_hits.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    local = load_local_effects(Path(args.peak_effects), peak_gene_map)
    local_summary = summarize_gene_class(local, args.class_beta_threshold)
    local_summary["gene_symbol"] = pd.Categorical(local_summary["gene_symbol"], categories=gene_order, ordered=True)
    local_summary["inferred_cell_class"] = pd.Categorical(local_summary["inferred_cell_class"], categories=CLASS_ORDER, ordered=True)
    local_summary = local_summary.sort_values(["gene_symbol", "inferred_cell_class"])
    local_summary["gene_symbol"] = local_summary["gene_symbol"].astype(str)
    local_summary["inferred_cell_class"] = local_summary["inferred_cell_class"].astype(str)
    local_summary.to_csv(
        outdir / "gse305146_candidate_peak_da_gene_class_summary.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    ranked = pd.read_csv(args.ranked_motif, sep="\t")
    motif_summary = summarize_motif_class(local, motif_hits, ranked, args.local_motif_beta_threshold)
    if not motif_summary.empty:
        motif_summary["gene_symbol"] = pd.Categorical(motif_summary["gene_symbol"], categories=gene_order, ordered=True)
        motif_summary["inferred_cell_class"] = pd.Categorical(
            motif_summary["inferred_cell_class"], categories=CLASS_ORDER, ordered=True
        )
        motif_summary["motif_family"] = pd.Categorical(motif_summary["motif_family"], categories=FAMILY_ORDER, ordered=True)
        motif_summary = motif_summary.sort_values(["gene_symbol", "inferred_cell_class", "motif_family"])
        motif_summary["gene_symbol"] = motif_summary["gene_symbol"].astype(str)
        motif_summary["inferred_cell_class"] = motif_summary["inferred_cell_class"].astype(str)
        motif_summary["motif_family"] = motif_summary["motif_family"].astype(str)
    motif_summary.to_csv(
        outdir / "gse305146_candidate_peak_da_gene_motif_class_summary.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    gene_summary = summarize_gene(local_summary, motif_summary, gene_order, args.class_beta_threshold)
    gene_summary = add_missing_gene_rows(gene_summary, windows, gene_order)
    gene_summary.to_csv(
        outdir / "gse305146_candidate_peak_da_gene_summary.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    axis_summary = mechanism_axis_summary(local_summary)
    axis_summary.to_csv(
        outdir / "gse305146_candidate_peak_da_mechanism_axis_summary.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    plot_gene_class_heatmap(local_summary, gene_order, figdir)
    plot_mechanism_heatmap(axis_summary, figdir)
    write_note(Path(args.note), gene_summary, local_summary, motif_summary)
    print("Wrote candidate-local GSE305146 peak DA integration")
    print(outdir / "gse305146_candidate_peak_da_gene_summary.tsv")
    print(outdir / "gse305146_candidate_peak_da_gene_class_summary.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
