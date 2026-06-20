#!/usr/bin/env python3
"""Make figures for human methylation and mouse-candidate convergence."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DEFAULT_RESULTS = "results/human_methylation_convergence"
DEFAULT_FIGDIR = "figures/human_methylation_convergence"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", default=DEFAULT_RESULTS)
    parser.add_argument("--figdir", default=DEFAULT_FIGDIR)
    parser.add_argument("--top-n", type=int, default=24)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [{key.strip(): value.strip() for key, value in row.items()} for row in reader]


def safe_float(value: object) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return math.nan
    return number if math.isfinite(number) else math.nan


def heatmap(matrix: np.ndarray, row_labels: list[str], col_labels: list[str], output_base: Path) -> None:
    display = np.ma.masked_invalid(matrix)
    cmap = plt.get_cmap("YlGnBu").copy()
    cmap.set_bad("#e5e5e5")
    height = max(5, 0.32 * len(row_labels) + 1.6)
    fig, ax = plt.subplots(figsize=(8.5, height))
    im = ax.imshow(display, aspect="auto", cmap=cmap)
    ax.set_xticks(np.arange(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=8)
    ax.set_title("Human Methylation and Mouse-Candidate Convergence", fontsize=11, pad=10)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Column-scaled support", fontsize=8)
    cbar.ax.tick_params(labelsize=8)
    fig.tight_layout()
    fig.savefig(output_base.with_suffix(".png"), dpi=300)
    fig.savefig(output_base.with_suffix(".pdf"))
    plt.close(fig)


def barplot(rows: list[dict[str, str]], output_base: Path) -> None:
    genes = [row["gene_symbol"] for row in rows]
    strict = np.array([safe_float(row["n_human_triple_phase3_fetal_mean_positive_probes"]) for row in rows])
    expression = np.array([safe_float(row["human_expression_support_score"]) for row in rows])
    mouse = np.array([safe_float(row["mouse_evidence_score"]) for row in rows])
    x = np.arange(len(genes))
    width = 0.26
    fig, ax = plt.subplots(figsize=(10, 4.8))
    ax.bar(x - width, strict, width=width, color="#5d9a63", label="Human methylation probes")
    ax.bar(x, expression, width=width, color="#c7793a", label="Human expr concordance with mouse OE")
    ax.bar(x + width, mouse, width=width, color="#4d79ad", label="Mouse score")
    ax.set_xticks(x)
    ax.set_xticklabels(genes, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Support count or score")
    ax.set_title("Top Human-Mouse Convergent Candidates", fontsize=11, pad=10)
    ax.legend(frameon=False, fontsize=8, ncols=3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(output_base.with_suffix(".png"), dpi=300)
    fig.savefig(output_base.with_suffix(".pdf"))
    plt.close(fig)


def scale_columns(matrix: np.ndarray) -> np.ndarray:
    out = matrix.copy().astype(float)
    for col in range(out.shape[1]):
        values = out[:, col]
        finite = np.isfinite(values)
        if not np.any(finite):
            continue
        max_value = np.max(values[finite])
        if max_value > 0:
            out[finite, col] = values[finite] / max_value
    return out


def main() -> None:
    args = parse_args()
    results_dir = Path(args.results_dir)
    figdir = Path(args.figdir)
    figdir.mkdir(parents=True, exist_ok=True)

    rows = read_tsv(results_dir / "human_mouse_convergent_candidate_genes.tsv")
    strict_rows = [
        row for row in rows if safe_float(row["n_human_triple_phase3_fetal_mean_positive_probes"]) > 0
    ][: args.top_n]
    if not strict_rows:
        return

    columns = [
        ("Strict fetal methyl probes", "n_human_triple_phase3_fetal_mean_positive_probes"),
        ("Fetal-majority probes", "n_human_triple_phase3_fetal_ds_majority_probes"),
        ("Any DS methyl probes", "n_human_triple_phase3_any_ds_probes"),
        ("Mouse score", "mouse_evidence_score"),
        ("Human expr concordance\nwith mouse OE", "human_expression_support_score"),
        ("Combined score", "combined_methylation_expression_score"),
    ]
    matrix = np.array([[safe_float(row[field]) for _label, field in columns] for row in strict_rows])
    heatmap(
        scale_columns(matrix),
        [row["gene_symbol"] for row in strict_rows],
        [label for label, _field in columns],
        figdir / "human_mouse_convergent_support_heatmap",
    )
    barplot(strict_rows[:16], figdir / "top_human_mouse_convergent_candidates")


if __name__ == "__main__":
    main()
