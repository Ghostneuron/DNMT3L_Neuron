#!/usr/bin/env python3
"""Make cross-species candidate-gene summary figures."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DEFAULT_RESULTS = "results/cross_species_candidate_expression"
DEFAULT_FIGDIR = "figures/cross_species_candidate_expression"
DEFAULT_GENES = [
    "FZD7",
    "TGIF2",
    "PLPP3",
    "PRDM8",
    "CRABP2",
    "FOXP1",
    "BCL11B",
    "RSPO1",
    "CYP26A1",
    "RXRA",
    "SCUBE3",
    "TDGF1",
]
DEFAULT_CELL_TYPES = [
    "IPC",
    "RG",
    "NEU_CUX2",
    "NEU_TLE4",
    "NEU_RELN",
    "NEU_RORB",
    "OPC",
    "MIC",
    "AST",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", default=DEFAULT_RESULTS)
    parser.add_argument("--figdir", default=DEFAULT_FIGDIR)
    parser.add_argument("--genes", default=",".join(DEFAULT_GENES))
    parser.add_argument("--cell-types", default=",".join(DEFAULT_CELL_TYPES))
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows.append({key.strip(): (value.strip() if isinstance(value, str) else value) for key, value in row.items()})
    return rows


def safe_float(value: object) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return math.nan
    return number if math.isfinite(number) else math.nan


def get_bulk_delta(rows: list[dict[str, str]], gene: str, dataset: str, context: str) -> float:
    for row in rows:
        if row["human_symbol"] == gene and row["dataset"] == dataset and row["context"] == context:
            return safe_float(row["delta_ds_minus_control_log2_cpm_plus_1"])
    return math.nan


def get_cell_delta(rows: list[dict[str, str]], gene: str, label: str, label_type: str = "author_cell_type") -> float:
    for row in rows:
        if row["human_symbol"] == gene and row["label_type"] == label_type and row["label"] == label:
            return safe_float(row["delta_ds_minus_control_log2_cpm_plus_1"])
    return math.nan


def clipped_matrix(values: np.ndarray, limit: float = 1.5) -> np.ndarray:
    return np.clip(values, -limit, limit)


def heatmap(
    matrix: np.ndarray,
    row_labels: list[str],
    col_labels: list[str],
    title: str,
    output_base: Path,
    vlim: float = 1.5,
) -> None:
    masked = np.ma.masked_invalid(clipped_matrix(matrix, vlim))
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#e5e5e5")
    height = max(4.0, 0.38 * len(row_labels) + 1.6)
    width = max(7.0, 0.48 * len(col_labels) + 2.2)
    fig, ax = plt.subplots(figsize=(width, height))
    im = ax.imshow(masked, aspect="auto", cmap=cmap, vmin=-vlim, vmax=vlim)
    ax.set_xticks(np.arange(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=9)
    ax.set_title(title, fontsize=11, pad=10)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("DS/control or model/control delta", fontsize=8)
    cbar.ax.tick_params(labelsize=8)
    fig.tight_layout()
    fig.savefig(output_base.with_suffix(".png"), dpi=300)
    fig.savefig(output_base.with_suffix(".pdf"))
    plt.close(fig)


def regional_barplot(
    bulk_rows: list[dict[str, str]],
    genes: list[str],
    output_base: Path,
) -> None:
    contexts = [
        ("GSE301886", "DLPFC", "DLPFC"),
        ("GSE301886", "Hippocampus", "Hippocampus"),
        ("GSE305139", "fetal_brain_final", "Fetal cortex"),
    ]
    values = np.array([[get_bulk_delta(bulk_rows, gene, dataset, context) for dataset, context, _ in contexts] for gene in genes])
    x = np.arange(len(genes))
    width = 0.24
    colors = ["#4778b3", "#d47f35", "#5f9b68"]
    fig, ax = plt.subplots(figsize=(8.8, 4.2))
    for idx, (_dataset, _context, label) in enumerate(contexts):
        ax.bar(x + (idx - 1) * width, values[:, idx], width=width, label=label, color=colors[idx])
    ax.axhline(0, color="#333333", linewidth=0.8)
    finite_values = values[np.isfinite(values)]
    if finite_values.size:
        lower = min(-0.1, float(finite_values.min()) * 1.18)
        upper = max(0.1, float(finite_values.max()) * 1.18)
        ax.set_ylim(lower, upper)
    ax.set_xticks(x)
    ax.set_xticklabels(genes, rotation=35, ha="right", fontsize=9)
    ax.set_ylabel("DS - control delta log2(CPM + 1)")
    ax.set_title("Human DS Fetal Bulk Candidate-Gene Effects", fontsize=11, pad=10)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(output_base.with_suffix(".png"), dpi=300)
    fig.savefig(output_base.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    args = parse_args()
    results_dir = Path(args.results_dir)
    figdir = Path(args.figdir)
    figdir.mkdir(parents=True, exist_ok=True)
    genes = [gene.strip().upper() for gene in args.genes.split(",") if gene.strip()]
    cell_types = [label.strip() for label in args.cell_types.split(",") if label.strip()]

    priority_rows = read_tsv(results_dir / "cross_species_candidate_priority.tsv")
    bulk_rows = read_tsv(results_dir / "human_bulk_candidate_effects.tsv")
    cell_rows = read_tsv(results_dir / "gse305153_candidate_celltype_effects.tsv")
    priority_by_gene = {row["human_symbol"]: row for row in priority_rows}

    heatmap_columns = [
        "Mouse OE",
        "Dp10",
        "DLPFC",
        "Hipp",
        "FetalCtx",
        "H5AD all",
        *cell_types,
    ]
    matrix = []
    for gene in genes:
        priority = priority_by_gene.get(gene, {})
        matrix.append(
            [
                safe_float(priority.get("mouse_gse268951_log2fc")),
                safe_float(priority.get("mouse_dp10_delta")),
                get_bulk_delta(bulk_rows, gene, "GSE301886", "DLPFC"),
                get_bulk_delta(bulk_rows, gene, "GSE301886", "Hippocampus"),
                get_bulk_delta(bulk_rows, gene, "GSE305139", "fetal_brain_final"),
                safe_float(priority.get("gse305153_all_cells_delta")),
                *[get_cell_delta(cell_rows, gene, label) for label in cell_types],
            ]
        )
    heatmap(
        np.array(matrix, dtype=float),
        genes,
        heatmap_columns,
        "Cross-Species Candidate-Gene Effect Summary",
        figdir / "cross_species_candidate_effect_heatmap",
    )

    top_genes = ["FZD7", "TGIF2", "PLPP3", "PRDM8", "CRABP2", "FOXP1"]
    regional_barplot(bulk_rows, top_genes, figdir / "human_fetal_bulk_top_candidate_region_effects")

    cell_matrix = np.array([[get_cell_delta(cell_rows, gene, label) for label in cell_types] for gene in top_genes])
    heatmap(
        cell_matrix,
        top_genes,
        cell_types,
        "GSE305153 Fetal Cortex Cell-Type Effects",
        figdir / "gse305153_top_candidate_celltype_heatmap",
    )


if __name__ == "__main__":
    main()
