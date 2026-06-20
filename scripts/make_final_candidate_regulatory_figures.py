#!/usr/bin/env python3
"""Make figures for final candidate evidence and regulatory context."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch


DEFAULT_MATRIX = "results/final_candidate_evidence/final_candidate_evidence_matrix_with_regulatory_context.tsv"
DEFAULT_FIGDIR = "figures/final_candidate_evidence"


TIER_COLORS = {
    "1_primary": "#2f6f8f",
    "2_secondary_balanced": "#6a9f58",
    "3_methylation_heavy": "#b68a3a",
    "4_expression_pathway": "#9a5c83",
}


TIER_LABELS = {
    "1_primary": "Primary methylation-convergent",
    "2_secondary_balanced": "Secondary methylation + expression",
    "3_methylation_heavy": "Methylation-heavy developmental",
    "4_expression_pathway": "Expression/pathway unresolved",
}


HEATMAP_COLUMNS = [
    ("Mouse methylation score", "mouse_evidence_score"),
    ("Human methylation score", "human_methylation_support_score"),
    ("Human expr concordance\nwith mouse OE", "human_expression_support_score"),
    ("Mouse-human expr\npriority score", "cross_species_expression_priority_score"),
    ("Human strict regulatory probes", "human_strict_regulatory_probe_count"),
    ("Human regulatory probes", "human_any_regulatory_probe_count"),
    ("Mouse regulatory DMRs", "mouse_any_regulatory_dmr_count"),
    ("Mouse ATAC DMRs", "mouse_e14_5_forebrain_atac_dmr_count"),
    ("Mouse H3K27ac DMRs", "mouse_embryonic_forebrain_h3k27ac_dmr_count"),
    ("Mouse H3K27me3 DMRs", "mouse_embryonic_forebrain_h3k27me3_dmr_count"),
]


BAR_COLUMNS = [
    ("Human strict regulatory probes", "human_strict_regulatory_probe_count", "#4f88b5"),
    ("Human other regulatory probes", "human_other_regulatory_probe_count", "#9ecae1"),
    ("Mouse ATAC DMRs", "mouse_e14_5_forebrain_atac_dmr_count", "#5f9b68"),
    ("Mouse H3K27ac DMRs", "mouse_embryonic_forebrain_h3k27ac_dmr_count", "#e3a347"),
    ("Mouse H3K27me3 DMRs", "mouse_embryonic_forebrain_h3k27me3_dmr_count", "#8c6bb1"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", default=DEFAULT_MATRIX)
    parser.add_argument("--figdir", default=DEFAULT_FIGDIR)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [
            {key.strip(): (value.strip() if isinstance(value, str) else "") for key, value in row.items()}
            for row in reader
        ]


def safe_float(value: object) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return math.nan
    return number if math.isfinite(number) else math.nan


def sort_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    tier_order = {
        "1_primary": 0,
        "2_secondary_balanced": 1,
        "3_methylation_heavy": 2,
        "4_expression_pathway": 3,
    }
    return sorted(
        rows,
        key=lambda row: (
            tier_order.get(row["tier"], 99),
            -safe_float(row.get("regulatory_context_score")),
            row["gene_symbol"],
        ),
    )


def scale_columns(matrix: np.ndarray) -> np.ndarray:
    out = matrix.copy().astype(float)
    for col_idx in range(out.shape[1]):
        values = out[:, col_idx]
        finite = np.isfinite(values)
        if not np.any(finite):
            continue
        max_value = np.max(values[finite])
        if max_value > 0:
            out[finite, col_idx] = values[finite] / max_value
    return out


def add_tier_strip(ax: plt.Axes, rows: list[dict[str, str]]) -> None:
    for idx, row in enumerate(rows):
        ax.add_patch(
            plt.Rectangle(
                (-0.78, idx - 0.5),
                0.14,
                1.0,
                color=TIER_COLORS.get(row["tier"], "#888888"),
                clip_on=True,
                linewidth=0,
            )
        )


def make_evidence_heatmap(rows: list[dict[str, str]], figdir: Path) -> None:
    matrix = np.array([[safe_float(row.get(field)) for _label, field in HEATMAP_COLUMNS] for row in rows])
    display = np.ma.masked_invalid(scale_columns(matrix))

    cmap = plt.get_cmap("YlGnBu").copy()
    cmap.set_bad("#eeeeee")

    fig, ax = plt.subplots(figsize=(11.8, 8.9))
    im = ax.imshow(display, aspect="auto", cmap=cmap, vmin=0, vmax=1)
    add_tier_strip(ax, rows)
    ax.set_xlim(-0.86, len(HEATMAP_COLUMNS) - 0.5)

    ax.set_yticks(np.arange(len(rows)))
    ax.set_yticklabels([row["gene_symbol"] for row in rows], fontsize=8)
    ax.set_xticks(np.arange(len(HEATMAP_COLUMNS)))
    ax.set_xticklabels([label for label, _field in HEATMAP_COLUMNS], rotation=40, ha="right", fontsize=8)
    ax.set_title("Final Candidate Evidence and Regulatory Context", fontsize=12, pad=12)
    ax.tick_params(length=0)
    ax.tick_params(axis="y", pad=10)
    for spine in ax.spines.values():
        spine.set_visible(False)

    legend_handles = [Patch(color=color, label=TIER_LABELS[tier]) for tier, color in TIER_COLORS.items()]
    fig.legend(
        handles=legend_handles,
        frameon=False,
        fontsize=7.5,
        loc="lower center",
        bbox_to_anchor=(0.48, 0.015),
        ncols=2,
    )
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Column-scaled support", fontsize=8)
    cbar.ax.tick_params(labelsize=8)
    fig.tight_layout(rect=(0.03, 0.08, 1, 1))
    fig.savefig(figdir / "final_candidate_evidence_regulatory_heatmap.png", dpi=300)
    fig.savefig(figdir / "final_candidate_evidence_regulatory_heatmap.pdf")
    plt.close(fig)


def make_regulatory_barplot(rows: list[dict[str, str]], figdir: Path) -> None:
    plot_rows = sorted(rows, key=lambda row: (-safe_float(row.get("regulatory_context_score")), row["gene_symbol"]))
    genes = [row["gene_symbol"] for row in plot_rows]
    for row in plot_rows:
        any_reg = safe_float(row.get("human_any_regulatory_probe_count"))
        strict = safe_float(row.get("human_strict_regulatory_probe_count"))
        row["human_other_regulatory_probe_count"] = str(max(0.0, any_reg - strict))

    x = np.arange(len(plot_rows))
    fig, ax = plt.subplots(figsize=(11.5, 5.2))
    bottom = np.zeros(len(plot_rows))
    for label, field, color in BAR_COLUMNS:
        values = np.array([safe_float(row.get(field)) for row in plot_rows])
        values = np.nan_to_num(values, nan=0.0)
        ax.bar(x, values, bottom=bottom, label=label, color=color, width=0.78)
        bottom += values

    ax.set_xticks(x)
    ax.set_xticklabels(genes, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Overlapping methylation intervals")
    ax.set_title("Candidate Regulatory Support from Fetal Brain ATAC and Histone Marks", fontsize=12, pad=10)
    ax.legend(frameon=False, fontsize=8, ncols=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#d9d9d9", linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(figdir / "final_candidate_regulatory_support_barplot.png", dpi=300)
    fig.savefig(figdir / "final_candidate_regulatory_support_barplot.pdf")
    plt.close(fig)


def make_class_count_plot(rows: list[dict[str, str]], figdir: Path) -> None:
    counts: dict[str, int] = {}
    for row in rows:
        label = row["regulatory_context_class"]
        counts[label] = counts.get(label, 0) + 1

    labels = [
        "human_mouse_regulatory_supported",
        "human_mouse_regulatory_possible",
        "human_regulatory_supported",
        "mouse_regulatory_supported",
        "no_regulatory_overlap_detected",
    ]
    labels = [label for label in labels if label in counts]
    values = [counts[label] for label in labels]
    display_labels = [label.replace("_", "\n") for label in labels]
    colors = ["#2f6f8f", "#7aa6c2", "#9bbf67", "#c99b4e", "#bdbdbd"][: len(labels)]

    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    ax.bar(np.arange(len(labels)), values, color=colors, width=0.7)
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(display_labels, fontsize=8)
    ax.set_ylabel("Candidate genes")
    ax.set_title("Regulatory-Context Classes Across Final Candidates", fontsize=11, pad=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for idx, value in enumerate(values):
        ax.text(idx, value + 0.25, str(value), ha="center", va="bottom", fontsize=9)
    ax.set_ylim(0, max(values) + 2)
    fig.tight_layout()
    fig.savefig(figdir / "final_candidate_regulatory_class_counts.png", dpi=300)
    fig.savefig(figdir / "final_candidate_regulatory_class_counts.pdf")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    rows = sort_rows(read_tsv(Path(args.matrix)))
    figdir = Path(args.figdir)
    figdir.mkdir(parents=True, exist_ok=True)

    make_evidence_heatmap(rows, figdir)
    make_regulatory_barplot(rows, figdir)
    make_class_count_plot(rows, figdir)
    print(f"Wrote final candidate regulatory figures to {figdir}")


if __name__ == "__main__":
    main()
