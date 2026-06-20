#!/usr/bin/env python3
"""Summarize GSE162170 fetal cortex chromVAR activity for mechanism TF families."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


FAMILY_PATTERNS: list[tuple[str, str]] = [
    ("STAT", r"STAT"),
    ("LHX_HOMEBOX", r"\bLHX|\bLhx"),
    ("DLX_HOMEBOX", r"\bDLX|\bDlx"),
    ("GSX_HOMEBOX", r"\bGSX|\bGsx"),
    ("EMX_HOMEBOX", r"\bEMX|\bEmx"),
    ("OLIG_BHLH", r"\bOLIG|\bOlig"),
    ("NEUROGENIC_BHLH", r"NEUROD|NeuroD|ATOH|Atoh|TCF4|NeuroG|NEUROG"),
    ("RFX_XBOX", r"\bRFX|\bRfx"),
    ("MEF2", r"\bMEF2|\bMef2"),
    ("SOX", r"\bSOX|\bSox"),
    ("NFIA_NFIX", r"\bNFIA|\bNFIA\b|\bNFIX|\bNfix"),
    ("FOXP", r"\bFOXP|\bFoxp"),
    ("TEAD", r"\bTEAD|\bTead"),
    ("ZEB", r"\bZEB|\bZeb"),
    ("EOMES", r"\bEOMES|\bEomes"),
]

FAMILY_ORDER = [
    "LHX_HOMEBOX",
    "DLX_HOMEBOX",
    "GSX_HOMEBOX",
    "EMX_HOMEBOX",
    "OLIG_BHLH",
    "NEUROGENIC_BHLH",
    "NFIA_NFIX",
    "RFX_XBOX",
    "MEF2",
    "SOX",
    "STAT",
    "FOXP",
    "TEAD",
    "ZEB",
    "EOMES",
]

AGE_ORDER = ["pcw16", "pcw20", "pcw21", "pcw24"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--chromvar",
        default="data/geo/GSE162170/supplementary/GSE162170_atac_chromVAR_TF_activities.tsv.gz",
    )
    parser.add_argument(
        "--metadata",
        default="data/geo/GSE162170/supplementary/GSE162170_atac_cell_metadata.txt.gz",
    )
    parser.add_argument("--outdir", default="results/mechanism_question")
    parser.add_argument("--figdir", default="figures/mechanism_question")
    parser.add_argument("--min-group-cells", type=int, default=50)
    return parser.parse_args()


def motif_tf_name(motif_id: str) -> str:
    return motif_id.split("_", 1)[1] if "_" in motif_id else motif_id


def motif_families(motif_id: str) -> list[str]:
    name = motif_tf_name(motif_id)
    families = [family for family, pattern in FAMILY_PATTERNS if re.search(pattern, name)]
    return families


def fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        return f"{value:.6g}"
    return str(value)


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fieldnames})


def load_metadata(path: Path, cell_ids: list[str]) -> pd.DataFrame:
    metadata = pd.read_csv(path, sep="\t", compression="gzip", dtype=str)
    metadata = metadata.set_index("Cell.ID").reindex(cell_ids)
    metadata["pcw"] = metadata["Age"].str.extract(r"pcw(\d+)").astype(float)
    metadata["age_order"] = metadata["Age"].map({age: i for i, age in enumerate(AGE_ORDER)})
    metadata["sample_id"] = metadata["Sample.ID"]
    metadata["tissue_id"] = metadata["Tissue.ID"]
    metadata["cluster_id"] = metadata["Iterative.LSI.Clusters"]
    metadata["cell_id"] = metadata.index
    missing = int(metadata["Age"].isna().sum())
    if missing:
        raise ValueError(f"{missing} chromVAR cells were not found in the ATAC metadata")
    return metadata.reset_index(drop=True)


def load_selected_chromvar(path: Path) -> tuple[list[str], list[dict[str, object]]]:
    rows: list[dict[str, object]] = []
    with gzip.open(path, "rt") as handle:
        cell_ids = next(handle).strip().split()
        for line in handle:
            fields = line.strip().split()
            if not fields:
                continue
            motif_id = fields[0]
            families = motif_families(motif_id)
            if not families:
                continue
            values = np.asarray(fields[1:], dtype=np.float32)
            if values.size != len(cell_ids):
                raise ValueError(f"{motif_id} has {values.size} values but expected {len(cell_ids)}")
            rows.append(
                {
                    "motif_id": motif_id,
                    "motif_tf": motif_tf_name(motif_id),
                    "families": families,
                    "values": values,
                }
            )
    return cell_ids, rows


def summarize_group(
    values: np.ndarray,
    metadata: pd.DataFrame,
    group_cols: list[str],
    min_cells: int = 1,
) -> pd.DataFrame:
    selected_cols = list(dict.fromkeys(group_cols + ["pcw", "sample_id", "tissue_id"]))
    frame = metadata[selected_cols].copy()
    frame["activity"] = values
    grouped = frame.groupby(group_cols, dropna=False)
    summary = grouped.agg(
        n_cells=("activity", "size"),
        mean_activity=("activity", "mean"),
        median_activity=("activity", "median"),
        sd_activity=("activity", "std"),
        n_samples=("sample_id", "nunique"),
        n_tissues=("tissue_id", "nunique"),
        pcw=("pcw", "median"),
    ).reset_index()
    summary = summary[summary["n_cells"] >= min_cells].copy()
    return summary


def add_motif_columns(summary: pd.DataFrame, motif_row: dict[str, object], family: str) -> pd.DataFrame:
    summary.insert(0, "motif_family", family)
    summary.insert(0, "motif_tf", motif_row["motif_tf"])
    summary.insert(0, "motif_id", motif_row["motif_id"])
    return summary


def spearman_or_empty(x: pd.Series, y: pd.Series) -> tuple[float, float]:
    if len(x) < 3 or x.nunique() < 2 or y.nunique() < 2:
        return (float("nan"), float("nan"))
    rho, p_value = stats.spearmanr(x, y)
    return (float(rho), float(p_value))


def make_family_tables(
    motif_age: pd.DataFrame,
    motif_cluster: pd.DataFrame,
    motif_age_cluster: pd.DataFrame,
    motif_sample: pd.DataFrame,
    min_group_cells: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    family_age = (
        motif_age.groupby(["motif_family", "Age", "pcw"], dropna=False)
        .agg(
            n_motifs=("motif_id", "nunique"),
            mean_activity=("mean_activity", "mean"),
            median_of_motif_medians=("median_activity", "median"),
            mean_n_cells=("n_cells", "mean"),
        )
        .reset_index()
    )
    family_age["age_order"] = family_age["Age"].map({age: i for i, age in enumerate(AGE_ORDER)})
    family_age = family_age.sort_values(["motif_family", "age_order"]).drop(columns=["age_order"])

    family_cluster = (
        motif_cluster.groupby(["motif_family", "cluster_id"], dropna=False)
        .agg(
            n_motifs=("motif_id", "nunique"),
            mean_activity=("mean_activity", "mean"),
            median_of_motif_medians=("median_activity", "median"),
            mean_n_cells=("n_cells", "mean"),
            n_samples=("n_samples", "max"),
            n_tissues=("n_tissues", "max"),
        )
        .reset_index()
    )
    family_cluster["cluster_num"] = family_cluster["cluster_id"].str.extract(r"c(\d+)").astype(float)
    family_cluster = family_cluster.sort_values(["motif_family", "cluster_num"]).drop(columns=["cluster_num"])

    family_age_cluster = (
        motif_age_cluster.groupby(["motif_family", "Age", "pcw", "cluster_id"], dropna=False)
        .agg(
            n_motifs=("motif_id", "nunique"),
            mean_activity=("mean_activity", "mean"),
            median_of_motif_medians=("median_activity", "median"),
            mean_n_cells=("n_cells", "mean"),
            n_samples=("n_samples", "max"),
            n_tissues=("n_tissues", "max"),
        )
        .reset_index()
    )
    family_age_cluster = family_age_cluster[family_age_cluster["mean_n_cells"] >= min_group_cells].copy()
    family_age_cluster["age_order"] = family_age_cluster["Age"].map({age: i for i, age in enumerate(AGE_ORDER)})
    family_age_cluster["cluster_num"] = family_age_cluster["cluster_id"].str.extract(r"c(\d+)").astype(float)
    family_age_cluster = family_age_cluster.sort_values(["motif_family", "age_order", "cluster_num"]).drop(
        columns=["age_order", "cluster_num"]
    )

    trends = []
    for family, family_df in motif_sample.groupby("motif_family"):
        sample_mean = (
            family_df.groupby(["sample_id", "Age", "pcw"], dropna=False)
            .agg(mean_activity=("mean_activity", "mean"), n_motifs=("motif_id", "nunique"))
            .reset_index()
        )
        rho, p_value = spearman_or_empty(sample_mean["pcw"], sample_mean["mean_activity"])
        trends.append(
            {
                "motif_family": family,
                "n_samples": len(sample_mean),
                "n_motifs": sample_mean["n_motifs"].max(),
                "spearman_rho_pcw_vs_activity": rho,
                "spearman_p_value_descriptive": p_value,
                "youngest_age_mean": sample_mean.loc[sample_mean["pcw"].idxmin(), "mean_activity"],
                "oldest_age_mean": sample_mean.loc[sample_mean["pcw"].idxmax(), "mean_activity"],
            }
        )
    family_trend = pd.DataFrame(trends)
    return family_age, family_cluster, family_age_cluster, family_trend


def plot_family_age(family_age: pd.DataFrame, outpath: Path) -> None:
    plot_families = [
        "LHX_HOMEBOX",
        "DLX_HOMEBOX",
        "GSX_HOMEBOX",
        "OLIG_BHLH",
        "NEUROGENIC_BHLH",
        "RFX_XBOX",
        "MEF2",
        "STAT",
        "SOX",
    ]
    colors = {
        "LHX_HOMEBOX": "#1f77b4",
        "DLX_HOMEBOX": "#2ca02c",
        "GSX_HOMEBOX": "#9467bd",
        "OLIG_BHLH": "#17becf",
        "NEUROGENIC_BHLH": "#d62728",
        "RFX_XBOX": "#7f7f7f",
        "MEF2": "#8c564b",
        "STAT": "#e377c2",
        "SOX": "#bcbd22",
    }
    fig, ax = plt.subplots(figsize=(8, 4.8))
    for family in plot_families:
        sub = family_age[family_age["motif_family"] == family].copy()
        if sub.empty:
            continue
        sub["age_order"] = sub["Age"].map({age: i for i, age in enumerate(AGE_ORDER)})
        sub = sub.sort_values("age_order")
        ax.plot(
            sub["pcw"],
            sub["mean_activity"],
            marker="o",
            linewidth=1.8,
            markersize=5,
            color=colors.get(family, "black"),
            label=family,
        )
    ax.axhline(0, color="#333333", linewidth=0.8, alpha=0.6)
    ax.set_xlabel("Post-conception week")
    ax.set_ylabel("Mean chromVAR activity")
    ax.set_title("GSE162170 fetal cortex TF motif activity")
    ax.legend(frameon=False, fontsize=8, ncol=3)
    fig.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def plot_family_cluster_heatmap(family_cluster: pd.DataFrame, outpath: Path) -> None:
    plot_families = [family for family in FAMILY_ORDER if family in set(family_cluster["motif_family"])]
    matrix = family_cluster.pivot_table(index="motif_family", columns="cluster_id", values="mean_activity", aggfunc="mean")
    matrix = matrix.reindex(plot_families)
    cluster_order = sorted(matrix.columns, key=lambda value: int(re.sub(r"\D", "", value) or 0))
    matrix = matrix[cluster_order]

    fig_width = max(8, 0.45 * len(cluster_order) + 3)
    fig_height = max(5, 0.35 * len(plot_families) + 2)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    image = ax.imshow(matrix.values, aspect="auto", cmap="RdBu_r", vmin=-1.5, vmax=1.5)
    ax.set_xticks(np.arange(len(cluster_order)))
    ax.set_xticklabels(cluster_order, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(plot_families)))
    ax.set_yticklabels(plot_families, fontsize=8)
    ax.set_xlabel("ATAC cluster ID")
    ax.set_title("GSE162170 motif family activity by ATAC cluster")
    cbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Mean chromVAR activity")
    fig.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    figdir = Path(args.figdir)

    cell_ids, chromvar_rows = load_selected_chromvar(Path(args.chromvar))
    metadata = load_metadata(Path(args.metadata), cell_ids)

    inventory_rows = []
    motif_age_frames = []
    motif_cluster_frames = []
    motif_age_cluster_frames = []
    motif_sample_frames = []

    for motif_row in chromvar_rows:
        values = motif_row["values"]
        assert isinstance(values, np.ndarray)
        for family in motif_row["families"]:
            inventory_rows.append(
                {
                    "motif_id": motif_row["motif_id"],
                    "motif_tf": motif_row["motif_tf"],
                    "motif_family": family,
                    "n_cells": len(cell_ids),
                    "overall_mean_activity": float(np.mean(values)),
                    "overall_median_activity": float(np.median(values)),
                    "overall_sd_activity": float(np.std(values, ddof=1)),
                }
            )
            motif_age_frames.append(add_motif_columns(summarize_group(values, metadata, ["Age"]), motif_row, family))
            motif_cluster_frames.append(
                add_motif_columns(summarize_group(values, metadata, ["cluster_id"], args.min_group_cells), motif_row, family)
            )
            motif_age_cluster_frames.append(
                add_motif_columns(
                    summarize_group(values, metadata, ["Age", "cluster_id"], args.min_group_cells),
                    motif_row,
                    family,
                )
            )
            motif_sample_frames.append(add_motif_columns(summarize_group(values, metadata, ["sample_id", "Age"]), motif_row, family))

    motif_age = pd.concat(motif_age_frames, ignore_index=True)
    motif_cluster = pd.concat(motif_cluster_frames, ignore_index=True)
    motif_age_cluster = pd.concat(motif_age_cluster_frames, ignore_index=True)
    motif_sample = pd.concat(motif_sample_frames, ignore_index=True)
    family_age, family_cluster, family_age_cluster, family_trend = make_family_tables(
        motif_age, motif_cluster, motif_age_cluster, motif_sample, args.min_group_cells
    )

    top_contexts = (
        family_age_cluster.sort_values(["motif_family", "mean_activity"], ascending=[True, False])
        .groupby("motif_family", as_index=False)
        .head(5)
        .copy()
    )
    top_contexts["rank_within_family"] = top_contexts.groupby("motif_family")["mean_activity"].rank(
        ascending=False, method="first"
    )
    top_contexts = top_contexts.sort_values(["motif_family", "rank_within_family"])

    write_tsv(
        outdir / "gse162170_chromvar_target_motif_inventory.tsv",
        inventory_rows,
        [
            "motif_id",
            "motif_tf",
            "motif_family",
            "n_cells",
            "overall_mean_activity",
            "overall_median_activity",
            "overall_sd_activity",
        ],
    )
    motif_age.to_csv(outdir / "gse162170_chromvar_target_motif_age_summary.tsv", sep="\t", index=False)
    motif_cluster.to_csv(outdir / "gse162170_chromvar_target_motif_cluster_summary.tsv", sep="\t", index=False)
    motif_age_cluster.to_csv(outdir / "gse162170_chromvar_target_motif_age_cluster_summary.tsv", sep="\t", index=False)
    motif_sample.to_csv(outdir / "gse162170_chromvar_target_motif_sample_summary.tsv", sep="\t", index=False)
    family_age.to_csv(outdir / "gse162170_chromvar_family_age_summary.tsv", sep="\t", index=False)
    family_cluster.to_csv(outdir / "gse162170_chromvar_family_cluster_summary.tsv", sep="\t", index=False)
    family_age_cluster.to_csv(outdir / "gse162170_chromvar_family_age_cluster_summary.tsv", sep="\t", index=False)
    family_trend.to_csv(outdir / "gse162170_chromvar_family_age_trends_descriptive.tsv", sep="\t", index=False)
    top_contexts.to_csv(outdir / "gse162170_chromvar_family_top_contexts.tsv", sep="\t", index=False)

    plot_family_age(family_age, figdir / "gse162170_chromvar_motif_family_age.png")
    plot_family_cluster_heatmap(family_cluster, figdir / "gse162170_chromvar_motif_family_cluster_heatmap.png")

    n_families = len({row["motif_family"] for row in inventory_rows})
    print(f"Selected {len(chromvar_rows)} motif rows across {n_families} families and {len(cell_ids)} ATAC cells")
    print(f"Wrote GSE162170 chromVAR mechanism tables to {outdir}")
    print(f"Wrote figures to {figdir}")


if __name__ == "__main__":
    main()
