#!/usr/bin/env python3
"""Build Extended Data STAT-paper marker reconciliation figure."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


FIGDIR = Path("figures/stat_marker_reconciliation")
INPUT = Path("results/mechanism_question/stat_paper_marker_dmg_audit.tsv")


def fdr_to_score(series: pd.Series) -> pd.Series:
    vals = pd.to_numeric(series, errors="coerce").clip(lower=1e-300)
    return -np.log10(vals)


def main() -> None:
    FIGDIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(INPUT, sep="\t")
    genes = df["human_gene"].tolist()

    human_direct = (
        pd.to_numeric(df["human_dnmt3l_oe_direct_candidate_dmc_count"], errors="coerce").fillna(0)
        + pd.to_numeric(df["human_dnmt3l_oe_direct_candidate_dmr_count"], errors="coerce").fillna(0)
    )
    mouse_nearby = pd.to_numeric(df["mouse_nearby_dnmt3l_dmr_count_100kb"], errors="coerce").fillna(0)
    mouse_promoter_body = (
        pd.to_numeric(df["mouse_gene_body_dnmt3l_dmr_count"], errors="coerce").fillna(0)
        + pd.to_numeric(df["mouse_promoter_2kb_dnmt3l_dmr_count"], errors="coerce").fillna(0)
    )
    mouse_logfc = pd.to_numeric(df["mouse_dnmt3l_oe_log2fc"], errors="coerce")
    mouse_fdr = pd.to_numeric(df["mouse_dnmt3l_oe_fdr"], errors="coerce")
    dnmt3a_fdr = pd.to_numeric(df["mouse_dnmt3a_cko_fdr"], errors="coerce")
    mouse_fdr_score = fdr_to_score(df["mouse_dnmt3l_oe_fdr"]).fillna(0)
    dnmt3a_fdr_score = fdr_to_score(df["mouse_dnmt3a_cko_fdr"]).fillna(0)

    plt.rcParams.update({
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.facecolor": "white",
    })
    fig = plt.figure(figsize=(13.5, 8), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05])

    ax = fig.add_subplot(gs[0, 0])
    x = np.arange(len(genes))
    width = 0.34
    ax.bar(x - width / 2, human_direct, width, color="#4C78A8", label="Human DNMT3L-OE direct DMC/DMR")
    ax.bar(x + width / 2, mouse_promoter_body, width, color="#E45756", label="Mouse promoter/gene-body DMR")
    ax.set_xticks(x, genes, rotation=35, ha="right")
    ax.set_ylabel("Count")
    ax.set_title("A  Direct gene-linked methylation evidence", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=8)

    ax = fig.add_subplot(gs[0, 1])
    ax.bar(x, mouse_nearby, color="#72B7B2")
    ax.set_xticks(x, genes, rotation=35, ha="right")
    ax.set_ylabel("Nearby DMRs within 100 kb")
    ax.set_title("B  Mouse P28 nearby Dnmt3l-OE DMRs", loc="left", fontweight="bold")
    for i, value in enumerate(mouse_nearby):
        ax.text(i, value + max(mouse_nearby) * 0.03, f"{int(value)}", ha="center", fontsize=8)

    ax = fig.add_subplot(gs[1, 0])
    colors = ["#4C78A8" if v < 0 else "#E45756" for v in mouse_logfc.fillna(0)]
    ax.bar(x, mouse_logfc.fillna(0), color=colors)
    ax.axhline(0, color="0.25", lw=0.8)
    ax.set_xticks(x, genes, rotation=35, ha="right")
    ax.set_ylabel("Mouse Dnmt3l-OE log2FC")
    ax.set_title("C  Mouse P28 expression effects are not FDR-significant", loc="left", fontweight="bold")
    vals = mouse_logfc.fillna(0)
    ax.set_ylim(vals.min() - 0.25, vals.max() + 0.2)
    for i, fdr in enumerate(mouse_fdr):
        value = vals.iloc[i]
        y = value + (0.07 if value >= 0 else -0.08)
        va = "bottom" if value >= 0 else "top"
        ax.text(i, y, f"FDR {fdr:.2g}", ha="center", va=va, fontsize=7)

    ax = fig.add_subplot(gs[1, 1])
    matrix = pd.DataFrame(
        {
            "human direct DMC/DMR": human_direct.gt(0).astype(int).to_numpy(),
            "mouse nearby DMR": mouse_nearby.gt(0).astype(int).to_numpy(),
            "mouse P28 DE FDR<0.05": mouse_fdr.lt(0.05).fillna(False).astype(int).to_numpy(),
            "Dnmt3a-cKO DE FDR<0.05": dnmt3a_fdr.lt(0.05).fillna(False).astype(int).to_numpy(),
        },
        index=genes,
    )
    ax.imshow(matrix.T, aspect="auto", cmap=plt.matplotlib.colors.ListedColormap(["#F2F2F2", "#2F6F8F"]))
    ax.set_xticks(np.arange(len(genes)), genes, rotation=35, ha="right")
    ax.set_yticks(np.arange(matrix.shape[1]), matrix.columns)
    ax.set_title("D  Evidence pattern summary", loc="left", fontweight="bold")
    for row in range(matrix.shape[1]):
        for col in range(matrix.shape[0]):
            ax.text(col, row, "yes" if matrix.iloc[col, row] else "no", ha="center", va="center", fontsize=7)

    fig.suptitle(
        "Extended Data Figure 8. STAT-paper marker reconciliation",
        fontweight="bold",
    )
    out = FIGDIR / "extended_data_stat_marker_reconciliation"
    fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    print(out.with_suffix(".png"))
    print(out.with_suffix(".pdf"))


if __name__ == "__main__":
    main()
