#!/usr/bin/env python3
"""Plot paper-facing summaries from GSE305146 marker-cell-class output tables."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


CORE_GENES = ["DNMT3L", "DNMT3A", "DNMT3B", "STAT1", "STAT3"]
CLASS_ORDER = [
    "RG",
    "Cycling_prog",
    "IPC",
    "NEU_TLE4_deep",
    "NEU_CUX2",
    "NEU_RORB",
    "NEU_RELN",
    "GABA_SST",
    "AST",
    "OPC_OL",
    "MIC",
    "VASC",
]
CLASS_LABELS = {
    "RG": "RG",
    "Cycling_prog": "Cycling",
    "IPC": "IPC",
    "NEU_TLE4_deep": "Deep EN",
    "NEU_CUX2": "CUX2 EN",
    "NEU_RORB": "RORB EN",
    "NEU_RELN": "RELN/CALB2",
    "GABA_SST": "GABA/SST",
    "AST": "AST",
    "OPC_OL": "OPC/OL",
    "MIC": "MIC",
    "VASC": "VASC",
}
MECHANISM_LABELS = {
    "cross_species_dnmt3a_linked_methylation_plus_developmental_cpg_tf_context": "DNMT3A methyl. + dev/CpG",
    "cross_species_dnmt3a_linked_methylation_plus_developmental_tf_context": "DNMT3A methyl. + dev TF",
    "cross_species_dnmt3a_linked_methylation_without_selected_tf_instance": "DNMT3A methyl. only",
    "mouse_dnmt3a_linked_methylation_with_tf_context_human_unresolved": "Mouse methyl. + TF; human unresolved",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--effects", default="results/mechanism_question/gse305146_marker_cell_class_effects.tsv")
    parser.add_argument("--mechanism", default="results/mechanism_question/gse305146_marker_mechanism_axis_cell_class_summary.tsv")
    parser.add_argument("--axis", default="results/mechanism_question/final_candidate_mechanism_axis_summary.tsv")
    parser.add_argument("--figdir", default="figures/mechanism_question")
    return parser.parse_args()


def make_heatmap(matrix: pd.DataFrame, path_prefix: Path, title: str, vmax: float) -> None:
    matrix = matrix.dropna(how="all").dropna(axis=1, how="all")
    if matrix.empty:
        return
    path_prefix.parent.mkdir(parents=True, exist_ok=True)
    values = matrix.to_numpy(dtype=float)
    fig_w = max(7.5, 0.58 * matrix.shape[1] + 3.0)
    fig_h = max(4.2, 0.34 * matrix.shape[0] + 2.2)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(values, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_title(title, fontsize=13, pad=12)
    ax.set_xticks(range(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(matrix.shape[0]))
    ax.set_yticklabels(matrix.index, fontsize=9)
    ax.set_xticks(np.arange(-0.5, matrix.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, matrix.shape[0], 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.032, pad=0.02)
    cbar.set_label("DS - control delta log2(CPM + 1)")
    fig.tight_layout()
    fig.savefig(path_prefix.with_suffix(".png"), dpi=300)
    fig.savefig(path_prefix.with_suffix(".pdf"))
    plt.close(fig)


def main() -> int:
    args = parse_args()
    figdir = Path(args.figdir)
    effects = pd.read_csv(args.effects, sep="\t")
    mechanism = pd.read_csv(args.mechanism, sep="\t")
    axis = pd.read_csv(args.axis, sep="\t")
    class_sizes = axis.groupby("final_mechanism_axis_call")["gene_symbol"].nunique().to_dict()

    core = effects[effects["gene_symbol"].isin(CORE_GENES)].copy()
    core_matrix = core.pivot_table(
        index="gene_symbol",
        columns="inferred_cell_class",
        values="delta_ds_minus_control_log2_cpm",
        aggfunc="first",
    ).reindex(index=CORE_GENES, columns=CLASS_ORDER)
    core_matrix.columns = [CLASS_LABELS.get(col, col) for col in core_matrix.columns]
    make_heatmap(
        core_matrix,
        figdir / "gse305146_marker_core_gene_cell_class_heatmap_curated",
        "GSE305146 core genes by inferred cell class",
        vmax=1.0,
    )

    candidates = effects[effects["gene_role"] == "final_candidate"].copy()
    gene_order = axis["gene_symbol"].tolist()
    candidate_matrix = candidates.pivot_table(
        index="gene_symbol",
        columns="inferred_cell_class",
        values="delta_ds_minus_control_log2_cpm",
        aggfunc="first",
    ).reindex(index=gene_order, columns=CLASS_ORDER)
    candidate_matrix.columns = [CLASS_LABELS.get(col, col) for col in candidate_matrix.columns]
    make_heatmap(
        candidate_matrix,
        figdir / "gse305146_marker_final_candidate_cell_class_heatmap_curated",
        "GSE305146 final candidates by inferred cell class",
        vmax=1.2,
    )

    mechanism_matrix = mechanism.pivot_table(
        index="final_mechanism_axis_call",
        columns="inferred_cell_class",
        values="mean_delta_log2_cpm",
        aggfunc="first",
    ).reindex(columns=CLASS_ORDER)
    mechanism_matrix.index = [
        f"{MECHANISM_LABELS.get(idx, idx)} (n={class_sizes.get(idx, 0)})" for idx in mechanism_matrix.index
    ]
    mechanism_matrix.columns = [CLASS_LABELS.get(col, col) for col in mechanism_matrix.columns]
    make_heatmap(
        mechanism_matrix,
        figdir / "gse305146_marker_mechanism_axis_cell_class_heatmap_curated",
        "GSE305146 mechanism classes by inferred cell class",
        vmax=0.8,
    )

    print("Wrote curated GSE305146 marker heatmaps")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
