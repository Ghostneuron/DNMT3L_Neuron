#!/usr/bin/env python3
"""Plot TcMAC21 P1 forebrain HSA21 expression from Kazuki et al. source data."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SOURCE = Path("data/papers/kazuki2020_tcmac21/elife-56223-fig2-data1-v3.xlsx")
OUTDIR = Path("results/kazuki2020_tcmac21_targets")
FIGDIR = Path("figures/kazuki2020_tcmac21_targets")
SUPP_TABLE = Path("manuscript/Supplementary Tables/ST2e_kazuki2020_tcmac21_p1_forebrain_hsa21_expression.tsv")


def clean_symbol(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).replace("\xa0", " ").strip()


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series.replace("N", np.nan), errors="coerce")


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    FIGDIR.mkdir(parents=True, exist_ok=True)

    xl = pd.ExcelFile(SOURCE)
    pcg_sheet = next(sheet for sheet in xl.sheet_names if "FPKM of PCG" in sheet)
    df = pd.read_excel(SOURCE, sheet_name=pcg_sheet)
    df = df.rename(
        columns={
            "Gene name(verified list)": "gene_symbol",
            "Gene start (bp)": "gene_start_bp",
            "Gene end (bp)": "gene_end_bp",
            "Karyotype band": "karyotype_band",
            "Mouse ortholog": "mouse_ortholog",
            "Similarity (between mouse and human)": "human_mouse_similarity",
            "Location on mouse chromosome": "mouse_chromosome",
            "Eu(HSA21 PCG, FPKM)": "eu_hsa21_fpkm",
            "MAC21(HSA21 PCG, FPKM)": "mac21_hsa21_fpkm",
            "Eu(HSA21 mouse ortholog, FPKM)": "eu_mouse_ortholog_fpkm",
            "MAC21(HSA21 mouse ortholog, FPKM)": "mac21_mouse_ortholog_fpkm",
            "Note": "note",
        }
    )
    for col in [
        "gene_symbol",
        "mouse_ortholog",
        "karyotype_band",
        "human_mouse_similarity",
        "mouse_chromosome",
        "note",
    ]:
        df[col] = df[col].map(clean_symbol)
    for col in [
        "gene_start_bp",
        "gene_end_bp",
        "eu_hsa21_fpkm",
        "mac21_hsa21_fpkm",
        "eu_mouse_ortholog_fpkm",
        "mac21_mouse_ortholog_fpkm",
    ]:
        df[col] = numeric(df[col])

    df.insert(0, "source_dataset", "Kazuki2020_eLife_Figure2_source_data1")
    df.insert(1, "source_context", "TcMAC21_P1_forebrain_RNAseq")
    df.insert(2, "source_sheet", pcg_sheet)
    df["hsa21_detected_gt0"] = df["mac21_hsa21_fpkm"].fillna(0).gt(0)
    df["hsa21_detected_ge1"] = df["mac21_hsa21_fpkm"].fillna(0).ge(1)
    df["mouse_ortholog_detected_ge1"] = df["mac21_mouse_ortholog_fpkm"].fillna(0).ge(1)
    df.to_csv(OUTDIR / "kazuki2020_hsa21_pcg_expression_clean.tsv", sep="\t", index=False)
    SUPP_TABLE.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(SUPP_TABLE, sep="\t", index=False)

    selected = ["POTED", "APP", "TRAPPC10", "DNMT3L", "PRMT2"]
    selected_df = df[df["gene_symbol"].isin(selected)].copy()
    selected_df["plot_order"] = selected_df["gene_symbol"].map({gene: i for i, gene in enumerate(selected)})
    selected_df = selected_df.sort_values("plot_order")
    selected_df.drop(columns=["plot_order"]).to_csv(
        OUTDIR / "kazuki2020_selected_hsa21_targets.tsv", sep="\t", index=False
    )

    summary = pd.DataFrame(
        [
            {
                "metric": "hsa21_pcg_total",
                "value": int(len(df)),
            },
            {
                "metric": "hsa21_pcg_mac21_fpkm_gt0",
                "value": int(df["hsa21_detected_gt0"].sum()),
            },
            {
                "metric": "hsa21_pcg_mac21_fpkm_ge1",
                "value": int(df["hsa21_detected_ge1"].sum()),
            },
            {
                "metric": "dnmt3l_hsa21_mac21_fpkm",
                "value": float(df.loc[df["gene_symbol"] == "DNMT3L", "mac21_hsa21_fpkm"].iloc[0]),
            },
            {
                "metric": "dnmt3l_mouse_ortholog_mac21_fpkm",
                "value": float(
                    df.loc[df["gene_symbol"] == "DNMT3L", "mac21_mouse_ortholog_fpkm"].iloc[0]
                ),
            },
        ]
    )
    summary.to_csv(OUTDIR / "kazuki2020_tcmac21_expression_summary.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(11.4, 8.6),
        gridspec_kw={"height_ratios": [1.05, 1.25], "hspace": 0.38},
    )

    ax = axes[0]
    plot_df = df.sort_values("gene_start_bp").copy()
    x = plot_df["gene_start_bp"] / 1_000_000
    y = np.log10(plot_df["mac21_hsa21_fpkm"].fillna(0) + 0.1)
    colors = np.where(plot_df["hsa21_detected_ge1"], "#4C78A8", "#B8C0CC")
    ax.scatter(x, y, s=22, c=colors, alpha=0.9, linewidths=0)
    highlight_colors = {
        "POTED": "#6E5A9B",
        "APP": "#54A24B",
        "TRAPPC10": "#F58518",
        "DNMT3L": "#D85C4A",
        "PRMT2": "#4C78A8",
    }
    for gene, color in highlight_colors.items():
        row = plot_df.loc[plot_df["gene_symbol"] == gene]
        if row.empty:
            continue
        xx = float(row["gene_start_bp"].iloc[0]) / 1_000_000
        yy = float(np.log10(row["mac21_hsa21_fpkm"].fillna(0).iloc[0] + 0.1))
        ax.scatter([xx], [yy], s=62, c=color, edgecolor="white", linewidth=0.8, zorder=4)
        dy = 0.16 if gene != "DNMT3L" else -0.24
        ax.text(xx, yy + dy, gene, ha="center", va="center", fontsize=8.4, color=color, fontweight="bold")
    ax.axhline(np.log10(1.1), color="#6B7280", linestyle="--", linewidth=1.0)
    ax.text(
        14.0,
        np.log10(1.1) + 0.06,
        "FPKM = 1",
        fontsize=8.2,
        color="#4B5563",
        va="bottom",
    )
    ax.set_title("TcMAC21 P1 forebrain: human HSA21q protein-coding expression", fontsize=13, fontweight="bold")
    ax.set_xlabel("Human chromosome 21 coordinate (Mb)")
    ax.set_ylabel("log10(FPKM + 0.1)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#E5E7EB", linewidth=0.8)

    ax = axes[1]
    measures = [
        ("Eu HSA21 background", "eu_hsa21_fpkm", "#9CA3AF"),
        ("TcMAC21 HSA21 transgene", "mac21_hsa21_fpkm", "#D85C4A"),
        ("Eu mouse ortholog", "eu_mouse_ortholog_fpkm", "#4C78A8"),
        ("TcMAC21 mouse ortholog", "mac21_mouse_ortholog_fpkm", "#F58518"),
    ]
    genes = selected_df["gene_symbol"].tolist()
    xpos = np.arange(len(genes))
    width = 0.18
    for idx, (label, col, color) in enumerate(measures):
        vals = selected_df[col].to_numpy(dtype=float)
        offset = (idx - 1.5) * width
        ax.bar(xpos + offset, vals, width=width, label=label, color=color)
        for xval, val in zip(xpos + offset, vals):
            if np.isnan(val):
                continue
            if val == 0:
                ax.text(xval, 0.012, "0", ha="center", va="bottom", fontsize=7.2, color="#4B5563")
            elif val < 0.25:
                ax.text(xval, val + 0.015, f"{val:.3g}", ha="center", va="bottom", fontsize=7.0, color="#4B5563")
    ax.set_yscale("symlog", linthresh=0.1)
    ax.set_ylim(0, 400)
    ax.set_xticks(xpos)
    ax.set_xticklabels(genes, rotation=0, fontsize=9.5)
    ax.set_ylabel("FPKM")
    ax.set_title("Selected HSA21 genes: human transgene signal and mouse ortholog context", fontsize=12.2, fontweight="bold")
    ax.legend(
        frameon=False,
        ncol=1,
        fontsize=8.4,
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        borderaxespad=0.0,
        handlelength=1.8,
    )
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#E5E7EB", linewidth=0.8)
    axes[0].text(-0.065, 1.06, "a", transform=axes[0].transAxes, fontsize=16, fontweight="bold")
    axes[1].text(-0.065, 1.06, "b", transform=axes[1].transAxes, fontsize=16, fontweight="bold")
    fig.text(
        0.99,
        0.012,
        "Source: Kazuki et al. 2020 eLife Fig. 2 source data 1. P1 forebrain RNA-seq FPKM.",
        ha="right",
        va="bottom",
        fontsize=8.2,
        color="#4B5563",
    )
    fig.subplots_adjust(bottom=0.09, right=0.82, hspace=0.42)

    fig.savefig(FIGDIR / "kazuki2020_tcmac21_p1_forebrain_hsa21_expression.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGDIR / "kazuki2020_tcmac21_p1_forebrain_hsa21_expression.pdf", bbox_inches="tight")
    plt.close(fig)

    print(FIGDIR / "kazuki2020_tcmac21_p1_forebrain_hsa21_expression.png")
    print(OUTDIR / "kazuki2020_selected_hsa21_targets.tsv")
    print(OUTDIR / "kazuki2020_tcmac21_expression_summary.tsv")
    print(SUPP_TABLE)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
