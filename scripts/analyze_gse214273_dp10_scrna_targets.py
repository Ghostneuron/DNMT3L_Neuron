#!/usr/bin/env python3
"""Analyze Dp(10)2Yey hippocampus scRNA-seq target-gene expression from GSE214273.

The GEO supplementary file contains raw 10x-style Matrix Market files but no
cell-type annotation table. This script therefore uses two conservative views:
sample-level pseudobulk expression and per-sample target-gene detection
fractions across all recovered cells.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import itertools
import math
import tarfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


SAMPLES = [
    {
        "accession": "GSM6602935",
        "sample": "WT_1",
        "prefix": "GSM6602935_CLE413A1",
        "group": "WT",
        "replicate": "1",
    },
    {
        "accession": "GSM6602936",
        "sample": "Dp2Yey_1",
        "prefix": "GSM6602936_CLE413A2",
        "group": "Dp(10)2Yey",
        "replicate": "1",
    },
    {
        "accession": "GSM6602937",
        "sample": "WT_2",
        "prefix": "GSM6602937_CLE413A3",
        "group": "WT",
        "replicate": "2",
    },
    {
        "accession": "GSM6602938",
        "sample": "Dp2Yey_2",
        "prefix": "GSM6602938_CLE413A4",
        "group": "Dp(10)2Yey",
        "replicate": "2",
    },
]

TARGET_ORDER = ["Dnmt3l", "Dnmt3a", "Dnmt3b", "Stat1", "Stat3"]
TARGET_ALIASES = {
    "Dnmt3l": {"Dnmt3l", "DNMT3L", "ENSMUSG00000000730"},
    "Dnmt3a": {"Dnmt3a", "DNMT3A", "ENSMUSG00000020661"},
    "Dnmt3b": {"Dnmt3b", "DNMT3B", "ENSMUSG00000027478"},
    "Stat1": {"Stat1", "STAT1", "ENSMUSG00000026104"},
    "Stat3": {"Stat3", "STAT3", "ENSMUSG00000004040"},
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-tar", default="data/geo/GSE214273/supplementary/GSE214273_RAW.tar")
    parser.add_argument("--raw-dir", default="data/geo/GSE214273/supplementary/raw_mtx")
    parser.add_argument("--outdir", default="results/gse214273_dp10_scrna")
    parser.add_argument("--figdir", default="figures/gse214273_dp10_scrna")
    parser.add_argument("--force-extract", action="store_true")
    return parser.parse_args()


def fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, (float, np.floating)):
        number = float(value)
        return f"{number:.6g}" if math.isfinite(number) else ""
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    return str(value)


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})


def write_tsv_gz(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})


def ensure_extracted(raw_tar: Path, raw_dir: Path, force: bool) -> None:
    expected = [raw_dir / f"{sample['prefix']}_matrix.mtx.gz" for sample in SAMPLES]
    if not force and all(path.exists() for path in expected):
        return
    raw_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(raw_tar) as tar:
        tar.extractall(raw_dir)


def read_features(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with gzip.open(path, "rt") as handle:
        for index, line in enumerate(handle):
            parts = line.rstrip("\n").split("\t")
            gene_id = parts[0] if parts else ""
            symbol = parts[1] if len(parts) > 1 else gene_id
            feature_type = parts[2] if len(parts) > 2 else ""
            rows.append(
                {
                    "feature_index_1based": index + 1,
                    "gene_id": gene_id,
                    "gene_symbol": symbol,
                    "feature_type": feature_type,
                }
            )
    return rows


def target_indices(features: list[dict[str, object]]) -> dict[int, str]:
    mapping: dict[int, str] = {}
    for row in features:
        gene_id = str(row["gene_id"])
        symbol = str(row["gene_symbol"])
        for target, aliases in TARGET_ALIASES.items():
            if gene_id in aliases or symbol in aliases:
                mapping[int(row["feature_index_1based"]) - 1] = target
    missing = [gene for gene in TARGET_ORDER if gene not in set(mapping.values())]
    if missing:
        raise RuntimeError(f"Missing target genes in feature file: {', '.join(missing)}")
    return mapping


def read_barcode_count(path: Path) -> int:
    count = 0
    with gzip.open(path, "rt") as handle:
        for count, _line in enumerate(handle, start=1):
            pass
    return count


def read_matrix_counts(
    matrix_path: Path,
    n_features: int,
    target_by_index: dict[int, str],
) -> dict[str, object]:
    gene_counts = np.zeros(n_features, dtype=np.float64)
    detected_cells = np.zeros(n_features, dtype=np.int32)
    target_cell_counts: dict[str, np.ndarray] = {}
    n_cells = None
    nnz_expected = None
    with gzip.open(matrix_path, "rt") as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            n_rows, n_cols, nnz_expected = [int(value) for value in line.split()]
            if n_rows != n_features:
                raise RuntimeError(f"{matrix_path} has {n_rows} rows but feature file has {n_features}")
            n_cells = n_cols
            target_cell_counts = {gene: np.zeros(n_cells, dtype=np.float32) for gene in TARGET_ORDER}
            break
        if n_cells is None:
            raise RuntimeError(f"No Matrix Market dimension line found in {matrix_path}")

        nnz_seen = 0
        for line in handle:
            row_s, col_s, value_s = line.split()
            feature_index = int(row_s) - 1
            cell_index = int(col_s) - 1
            count = float(value_s)
            gene_counts[feature_index] += count
            if count > 0:
                detected_cells[feature_index] += 1
            gene = target_by_index.get(feature_index)
            if gene is not None:
                target_cell_counts[gene][cell_index] += count
            nnz_seen += 1

    return {
        "n_cells": int(n_cells),
        "nnz_expected": int(nnz_expected or 0),
        "nnz_seen": nnz_seen,
        "gene_counts": gene_counts,
        "detected_cells": detected_cells,
        "target_cell_counts": target_cell_counts,
        "total_umi": float(gene_counts.sum()),
    }


def mean(values: list[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    return sum(finite) / len(finite) if finite else math.nan


def exact_permutation_pvalue(wt_values: list[float], dp_values: list[float]) -> float:
    observed = mean(dp_values) - mean(wt_values)
    pooled = wt_values + dp_values
    dp_size = len(dp_values)
    if not pooled or dp_size == 0:
        return math.nan
    total = 0
    extreme = 0
    indices = range(len(pooled))
    for dp_indices in itertools.combinations(indices, dp_size):
        dp_set = set(dp_indices)
        perm_dp = [pooled[idx] for idx in indices if idx in dp_set]
        perm_wt = [pooled[idx] for idx in indices if idx not in dp_set]
        delta = mean(perm_dp) - mean(perm_wt)
        total += 1
        if abs(delta) >= abs(observed) - 1e-12:
            extreme += 1
    return extreme / total if total else math.nan


def direction(delta: float, threshold: float = 0.1) -> str:
    if not math.isfinite(delta) or abs(delta) < threshold:
        return "flat"
    return "up" if delta > 0 else "down"


def analyze(raw_dir: Path) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    features = read_features(raw_dir / f"{SAMPLES[0]['prefix']}_features.tsv.gz")
    target_by_index = target_indices(features)
    n_features = len(features)

    sample_gene_rows: list[dict[str, object]] = []
    target_cell_rows: list[dict[str, object]] = []
    sample_metrics: dict[str, dict[str, object]] = {}
    sample_gene_log2: dict[str, np.ndarray] = {}
    sample_gene_cpm: dict[str, np.ndarray] = {}
    sample_gene_detected_fraction: dict[str, np.ndarray] = {}

    for sample in SAMPLES:
        matrix_path = raw_dir / f"{sample['prefix']}_matrix.mtx.gz"
        barcode_path = raw_dir / f"{sample['prefix']}_barcodes.tsv.gz"
        barcode_count = read_barcode_count(barcode_path)
        counts = read_matrix_counts(matrix_path, n_features=n_features, target_by_index=target_by_index)
        n_cells = int(counts["n_cells"])
        if barcode_count != n_cells:
            raise RuntimeError(f"{sample['sample']} barcode count {barcode_count} != matrix cells {n_cells}")

        gene_counts = counts["gene_counts"]
        detected_cells = counts["detected_cells"]
        total_umi = float(counts["total_umi"])
        cpm = gene_counts / total_umi * 1_000_000 if total_umi else np.full(n_features, np.nan)
        log2_cpm = np.log2(cpm + 1.0)
        detected_fraction = detected_cells / n_cells if n_cells else np.full(n_features, np.nan)

        sample_metrics[sample["sample"]] = {**sample, "n_cells": n_cells, "total_umi": total_umi}
        sample_gene_cpm[sample["sample"]] = cpm
        sample_gene_log2[sample["sample"]] = log2_cpm
        sample_gene_detected_fraction[sample["sample"]] = detected_fraction

        for idx, feature in enumerate(features):
            sample_gene_rows.append(
                {
                    "dataset": "GSE214273",
                    "accession": sample["accession"],
                    "sample": sample["sample"],
                    "group": sample["group"],
                    "replicate": sample["replicate"],
                    "n_cells": n_cells,
                    "total_umi": total_umi,
                    "gene_id": feature["gene_id"],
                    "gene_symbol": feature["gene_symbol"],
                    "count": gene_counts[idx],
                    "detected_cells": detected_cells[idx],
                    "detected_fraction": detected_fraction[idx],
                    "cpm": cpm[idx],
                    "log2_cpm_plus_1": log2_cpm[idx],
                    "target_gene": "yes" if str(feature["gene_symbol"]) in TARGET_ORDER else "no",
                }
            )

        for gene, cell_counts in counts["target_cell_counts"].items():
            for cell_index, cell_count in enumerate(cell_counts, start=1):
                target_cell_rows.append(
                    {
                        "dataset": "GSE214273",
                        "accession": sample["accession"],
                        "sample": sample["sample"],
                        "group": sample["group"],
                        "replicate": sample["replicate"],
                        "cell_index": cell_index,
                        "gene_symbol": gene,
                        "count": float(cell_count),
                        "detected": int(cell_count > 0),
                    }
                )

    effect_rows: list[dict[str, object]] = []
    wt_samples = [sample["sample"] for sample in SAMPLES if sample["group"] == "WT"]
    dp_samples = [sample["sample"] for sample in SAMPLES if sample["group"] == "Dp(10)2Yey"]
    for idx, feature in enumerate(features):
        wt_log2 = [float(sample_gene_log2[sample][idx]) for sample in wt_samples]
        dp_log2 = [float(sample_gene_log2[sample][idx]) for sample in dp_samples]
        wt_cpm = [float(sample_gene_cpm[sample][idx]) for sample in wt_samples]
        dp_cpm = [float(sample_gene_cpm[sample][idx]) for sample in dp_samples]
        wt_detect = [float(sample_gene_detected_fraction[sample][idx]) for sample in wt_samples]
        dp_detect = [float(sample_gene_detected_fraction[sample][idx]) for sample in dp_samples]
        delta_log2 = mean(dp_log2) - mean(wt_log2)
        delta_detect = mean(dp_detect) - mean(wt_detect)
        symbol = str(feature["gene_symbol"])
        effect_rows.append(
            {
                "dataset": "GSE214273",
                "comparison": "Dp(10)2Yey_vs_WT_pseudobulk",
                "gene_id": feature["gene_id"],
                "gene_symbol": symbol,
                "target_gene": "yes" if symbol in TARGET_ORDER else "no",
                "wt_mean_cpm": mean(wt_cpm),
                "dp10_mean_cpm": mean(dp_cpm),
                "wt_mean_log2_cpm_plus_1": mean(wt_log2),
                "dp10_mean_log2_cpm_plus_1": mean(dp_log2),
                "delta_log2_cpm_plus_1": delta_log2,
                "wt_mean_detected_fraction": mean(wt_detect),
                "dp10_mean_detected_fraction": mean(dp_detect),
                "delta_detected_fraction": delta_detect,
                "permutation_pvalue_log2_cpm": exact_permutation_pvalue(wt_log2, dp_log2),
                "direction": direction(delta_log2),
            }
        )

    target_effect_rows = [row for row in effect_rows if row["target_gene"] == "yes"]
    target_effect_rows.sort(key=lambda row: TARGET_ORDER.index(str(row["gene_symbol"])))
    return sample_gene_rows, effect_rows, target_effect_rows, target_cell_rows


def make_target_plots(
    sample_gene_rows: list[dict[str, object]],
    target_effect_rows: list[dict[str, object]],
    effect_rows: list[dict[str, object]],
    figdir: Path,
) -> None:
    figdir.mkdir(parents=True, exist_ok=True)
    target_sample_rows = [
        row for row in sample_gene_rows if row["gene_symbol"] in TARGET_ORDER
    ]
    colors = {"WT": "#3f7f93", "Dp(10)2Yey": "#b04a5a"}
    x_positions = np.arange(len(TARGET_ORDER))
    offsets = {"WT_1": -0.24, "WT_2": -0.08, "Dp2Yey_1": 0.08, "Dp2Yey_2": 0.24}

    fig, ax = plt.subplots(figsize=(7.2, 4.1))
    for row in target_sample_rows:
        gene = str(row["gene_symbol"])
        x = TARGET_ORDER.index(gene) + offsets[str(row["sample"])]
        ax.scatter(
            x,
            float(row["log2_cpm_plus_1"]),
            color=colors[str(row["group"])],
            edgecolor="black",
            linewidth=0.25,
            s=52,
            label=str(row["group"]),
            alpha=0.9,
        )
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), frameon=False, loc="upper left")
    ax.set_xticks(x_positions)
    ax.set_xticklabels(TARGET_ORDER)
    ax.set_ylabel("log2(CPM + 1)")
    ax.set_title("GSE214273 hippocampus scRNA-seq pseudobulk")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(figdir / f"target_gene_log2cpm_by_sample.{suffix}", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    genes = [str(row["gene_symbol"]) for row in target_effect_rows]
    deltas = [float(row["delta_log2_cpm_plus_1"]) for row in target_effect_rows]
    bar_colors = ["#b04a5a" if value > 0 else "#3f7f93" for value in deltas]
    ax.barh(genes, deltas, color=bar_colors)
    ax.axvline(0, color="#222222", linewidth=0.8)
    ax.set_xlabel("Dp(10)2Yey - WT delta log2(CPM + 1)")
    ax.set_title("Target-gene expression effect")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(figdir / f"target_gene_delta_log2cpm.{suffix}", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.1))
    for row in target_sample_rows:
        gene = str(row["gene_symbol"])
        x = TARGET_ORDER.index(gene) + offsets[str(row["sample"])]
        ax.scatter(
            x,
            float(row["detected_fraction"]) * 100.0,
            color=colors[str(row["group"])],
            edgecolor="black",
            linewidth=0.25,
            s=52,
            alpha=0.9,
        )
    ax.set_xticks(x_positions)
    ax.set_xticklabels(TARGET_ORDER)
    ax.set_ylabel("Detected cells (%)")
    ax.set_title("Target-gene detection across hippocampal cells")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(figdir / f"target_gene_detection_fraction_by_sample.{suffix}", dpi=220)
    plt.close(fig)

    xs = []
    ys = []
    labels = []
    for row in effect_rows:
        wt = float(row["wt_mean_log2_cpm_plus_1"])
        dp = float(row["dp10_mean_log2_cpm_plus_1"])
        xs.append((wt + dp) / 2.0)
        ys.append(float(row["delta_log2_cpm_plus_1"]))
        labels.append(str(row["gene_symbol"]))
    fig, ax = plt.subplots(figsize=(6.6, 4.5))
    ax.scatter(xs, ys, s=4, color="#6b6b6b", alpha=0.35, linewidth=0)
    for row in target_effect_rows:
        wt = float(row["wt_mean_log2_cpm_plus_1"])
        dp = float(row["dp10_mean_log2_cpm_plus_1"])
        x = (wt + dp) / 2.0
        y = float(row["delta_log2_cpm_plus_1"])
        ax.scatter([x], [y], s=48, color="#d49a2a", edgecolor="black", linewidth=0.3)
        ax.text(x + 0.04, y, str(row["gene_symbol"]), fontsize=8, va="center")
    ax.axhline(0, color="#222222", linewidth=0.8)
    ax.set_xlabel("Mean log2(CPM + 1)")
    ax.set_ylabel("Dp(10)2Yey - WT delta log2(CPM + 1)")
    ax.set_title("Genome-wide pseudobulk expression shift")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(figdir / f"genomewide_delta_log2cpm_target_labels.{suffix}", dpi=220)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    raw_tar = Path(args.raw_tar)
    raw_dir = Path(args.raw_dir)
    outdir = Path(args.outdir)
    figdir = Path(args.figdir)

    ensure_extracted(raw_tar, raw_dir, force=args.force_extract)
    sample_gene_rows, effect_rows, target_effect_rows, target_cell_rows = analyze(raw_dir)

    sample_gene_fields = [
        "dataset",
        "accession",
        "sample",
        "group",
        "replicate",
        "n_cells",
        "total_umi",
        "gene_id",
        "gene_symbol",
        "count",
        "detected_cells",
        "detected_fraction",
        "cpm",
        "log2_cpm_plus_1",
        "target_gene",
    ]
    effect_fields = [
        "dataset",
        "comparison",
        "gene_id",
        "gene_symbol",
        "target_gene",
        "wt_mean_cpm",
        "dp10_mean_cpm",
        "wt_mean_log2_cpm_plus_1",
        "dp10_mean_log2_cpm_plus_1",
        "delta_log2_cpm_plus_1",
        "wt_mean_detected_fraction",
        "dp10_mean_detected_fraction",
        "delta_detected_fraction",
        "permutation_pvalue_log2_cpm",
        "direction",
    ]
    target_cell_fields = [
        "dataset",
        "accession",
        "sample",
        "group",
        "replicate",
        "cell_index",
        "gene_symbol",
        "count",
        "detected",
    ]
    write_tsv_gz(outdir / "gse214273_sample_gene_pseudobulk.tsv.gz", sample_gene_rows, sample_gene_fields)
    write_tsv(outdir / "gse214273_gene_effect_summary.tsv", effect_rows, effect_fields)
    write_tsv(outdir / "gse214273_target_gene_effect_summary.tsv", target_effect_rows, effect_fields)
    write_tsv_gz(outdir / "gse214273_target_cell_values.tsv.gz", target_cell_rows, target_cell_fields)
    make_target_plots(sample_gene_rows, target_effect_rows, effect_rows, figdir)

    print(
        "Wrote "
        f"{len(effect_rows)} gene effects, {len(target_effect_rows)} target effects, "
        f"{len(target_cell_rows)} target cell rows, and figures to {figdir}."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
