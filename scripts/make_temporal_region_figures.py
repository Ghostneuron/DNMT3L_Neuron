#!/usr/bin/env python3
"""Make temporal and regional summary figures for target genes.

Inputs are the project result tables already generated from BrainSpan, GEO
single-cell, Allen ISH, and Allen ABC single-cell analyses.
"""

from __future__ import annotations

import csv
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/dnmt3l_neuron_mplconfig")
os.environ.setdefault("XDG_CACHE_HOME", "/private/tmp/dnmt3l_neuron_xdg_cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


GENE_ORDER_HUMAN = ["DNMT3L", "DNMT3A", "DNMT3B", "STAT1", "STAT3"]
GENE_ORDER_MOUSE = ["Dnmt3l", "Dnmt3a", "Dnmt3b", "Stat1", "Stat3"]
GENE_COLORS_HUMAN = {
    "DNMT3L": "#7A3E9D",
    "DNMT3A": "#1B6CA8",
    "DNMT3B": "#2E8B57",
    "STAT1": "#C44E52",
    "STAT3": "#D9902F",
}
GENE_COLORS_MOUSE = {
    "Dnmt3l": "#7A3E9D",
    "Dnmt3a": "#1B6CA8",
    "Dnmt3b": "#2E8B57",
    "Stat1": "#C44E52",
    "Stat3": "#D9902F",
}

REGION_LABELS = {
    "A1C": "Auditory ctx",
    "AMY": "Amygdala",
    "CB": "Cerebellum",
    "CBC": "Cerebellar ctx",
    "CGE": "CGE",
    "DFC": "DLPFC",
    "DTH": "Dorsal thal.",
    "HIP": "Hippocampus",
    "IPC": "Inf. parietal ctx",
    "ITC": "Inf. temporal ctx",
    "LGE": "LGE",
    "M1C": "Motor ctx",
    "M1C-S1C": "Motor-sensory ctx",
    "MD": "Mediodorsal thal.",
    "MFC": "Medial PFC",
    "MGE": "MGE",
    "OFC": "Orbital PFC",
    "Ocx": "Occipital neocortex",
    "PCx": "Parietal neocortex",
    "S1C": "Somatosensory ctx",
    "STC": "Sup. temporal ctx",
    "STR": "Striatum",
    "TCx": "Temporal neocortex",
    "URL": "Rostral rhombic lip",
    "V1C": "Visual ctx",
    "VFC": "Ventrolateral PFC",
}

PRENATAL_REGION_ORDER = [
    "MGE",
    "LGE",
    "CGE",
    "DFC",
    "MFC",
    "OFC",
    "VFC",
    "M1C",
    "S1C",
    "A1C",
    "V1C",
    "STC",
    "ITC",
    "HIP",
    "AMY",
    "STR",
    "DTH",
    "MD",
    "CBC",
    "CB",
]

ADULT_REGION_ORDER = [
    "DFC",
    "MFC",
    "OFC",
    "VFC",
    "M1C",
    "S1C",
    "A1C",
    "V1C",
    "STC",
    "ITC",
    "IPC",
    "HIP",
    "AMY",
    "STR",
    "MD",
    "CBC",
]

MOUSE_AGE_ORDER = ["E11.5", "E13.5", "E15.5", "E18.5", "P4", "P14", "P28", "P56"]


def setup_style() -> None:
    sns.set_theme(style="whitegrid", context="paper")
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "figure.dpi": 150,
            "savefig.dpi": 300,
        }
    )


def savefig(fig: plt.Figure, outdir: Path, stem: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"{stem}.png", bbox_inches="tight")
    fig.savefig(outdir / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def pivot_gene_region(df: pd.DataFrame, region_order: list[str]) -> pd.DataFrame:
    subset = df[df["structure_acronym"].isin(region_order)].copy()
    subset["gene_symbol"] = pd.Categorical(
        subset["gene_symbol"],
        categories=GENE_ORDER_HUMAN,
        ordered=True,
    )
    subset["structure_acronym"] = pd.Categorical(
        subset["structure_acronym"],
        categories=region_order,
        ordered=True,
    )
    matrix = subset.pivot_table(
        index="gene_symbol",
        columns="structure_acronym",
        values="mean_rpkm",
        observed=False,
    )
    matrix = matrix.reindex(index=GENE_ORDER_HUMAN, columns=region_order)
    matrix.columns = [REGION_LABELS.get(col, col) for col in matrix.columns]
    return matrix


def plot_human_brainspan_temporal(outdir: Path) -> None:
    df = pd.read_csv("results/allen_brainspan/brainspan_target_gene_age_summary.tsv", sep="\t")
    df = df[df["gene_symbol"].isin(GENE_ORDER_HUMAN)].copy()
    age_table = (
        df[["age_sort_days", "age_name", "developmental_period"]]
        .drop_duplicates()
        .sort_values("age_sort_days")
        .reset_index(drop=True)
    )
    age_to_x = {row.age_name: idx for idx, row in age_table.iterrows()}
    df["x"] = df["age_name"].map(age_to_x)

    fig, axes = plt.subplots(len(GENE_ORDER_HUMAN), 1, figsize=(11, 10), sharex=True)
    for ax, gene in zip(axes, GENE_ORDER_HUMAN):
        gdf = df[df["gene_symbol"] == gene].sort_values("x")
        ax.plot(
            gdf["x"],
            gdf["mean_rpkm"],
            color=GENE_COLORS_HUMAN[gene],
            marker="o",
            linewidth=1.5,
            markersize=3.5,
        )
        ax.fill_between(
            gdf["x"].to_numpy(dtype=float),
            gdf["min_rpkm"].to_numpy(dtype=float),
            gdf["max_rpkm"].to_numpy(dtype=float),
            color=GENE_COLORS_HUMAN[gene],
            alpha=0.12,
            linewidth=0,
        )
        ax.set_ylabel(gene)
        ax.set_ylim(bottom=0)
        ax.grid(axis="y", color="#E6E6E6")
        ax.grid(axis="x", visible=False)

    prenatal_x = [
        idx for idx, row in age_table.iterrows() if row.developmental_period == "prenatal"
    ]
    if prenatal_x:
        last_prenatal = max(prenatal_x)
        for ax in axes:
            ax.axvline(last_prenatal + 0.5, color="#555555", linewidth=0.8, linestyle="--")
        axes[0].text(
            last_prenatal / 2,
            axes[0].get_ylim()[1] * 1.05,
            "prenatal",
            ha="center",
            va="bottom",
            fontsize=8,
        )
        axes[0].text(
            last_prenatal + (len(age_table) - last_prenatal) / 2,
            axes[0].get_ylim()[1] * 1.05,
            "postnatal",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    axes[-1].set_xticks(range(len(age_table)))
    axes[-1].set_xticklabels(age_table["age_name"], rotation=65, ha="right")
    axes[-1].set_xlabel("BrainSpan age")
    fig.suptitle("Human BrainSpan temporal target-gene expression", y=0.995)
    fig.text(0.02, 0.5, "Mean RPKM; shading shows min-max across regions/specimens", rotation=90, va="center")
    savefig(fig, outdir, "human_brainspan_temporal_rpkm")


def plot_human_brainspan_region_heatmaps(outdir: Path) -> None:
    prenatal = pd.read_csv(
        "results/allen_brainspan/brainspan_prenatal_target_gene_region_summary.tsv",
        sep="\t",
    )
    adult = pd.read_csv(
        "results/allen_brainspan/brainspan_adult_target_gene_region_summary.tsv",
        sep="\t",
    )
    prenatal_matrix = pivot_gene_region(prenatal, PRENATAL_REGION_ORDER)
    adult_matrix = pivot_gene_region(adult, ADULT_REGION_ORDER)

    fig, axes = plt.subplots(2, 1, figsize=(13, 6.6), constrained_layout=True)
    for ax, matrix, title in [
        (axes[0], prenatal_matrix, "Prenatal human BrainSpan mean RPKM"),
        (axes[1], adult_matrix, "Adult human BrainSpan mean RPKM"),
    ]:
        log_matrix = np.log10(matrix.astype(float) + 0.01)
        sns.heatmap(
            log_matrix,
            ax=ax,
            cmap="viridis",
            cbar_kws={"label": "log10(mean RPKM + 0.01)"},
            linewidths=0.4,
            linecolor="white",
        )
        ax.set_title(title)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.tick_params(axis="x", rotation=45)
        ax.tick_params(axis="y", rotation=0)
    savefig(fig, outdir, "human_brainspan_region_heatmaps")


def plot_human_gse104276_week(outdir: Path) -> None:
    df = pd.read_csv("results/single_cell/GSE104276/GSE104276_target_gene_by_week.tsv", sep="\t")
    df = df[df["group"].str.startswith("GW", na=False)].copy()
    df["week_num"] = df["group"].str.replace("GW", "", regex=False).astype(int)
    df["gene_symbol"] = pd.Categorical(df["gene_symbol"], categories=GENE_ORDER_HUMAN, ordered=True)
    df = df.sort_values(["gene_symbol", "week_num"])

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for gene in GENE_ORDER_HUMAN:
        gdf = df[df["gene_symbol"] == gene]
        axes[0].plot(
            gdf["week_num"],
            gdf["detected_fraction"],
            marker="o",
            linewidth=1.5,
            markersize=4,
            label=gene,
            color=GENE_COLORS_HUMAN[gene],
        )
        axes[1].plot(
            gdf["week_num"],
            gdf["mean_tpm"],
            marker="o",
            linewidth=1.5,
            markersize=4,
            label=gene,
            color=GENE_COLORS_HUMAN[gene],
        )
    axes[0].set_title("Detection fraction")
    axes[0].set_ylabel("Fraction of cells with TPM > 0")
    axes[0].set_ylim(0, 1)
    axes[1].set_title("Mean expression")
    axes[1].set_ylabel("Mean TPM")
    axes[1].set_yscale("symlog", linthresh=1)
    for ax in axes:
        ax.set_xlabel("Gestational week")
        ax.set_xticks(sorted(df["week_num"].unique()))
        ax.grid(axis="y", color="#E6E6E6")
    axes[1].legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False)
    fig.suptitle("Human fetal PFC single-cell target-gene expression (GSE104276)")
    savefig(fig, outdir, "human_gse104276_fetal_pfc_week_trends")


def plot_mouse_ish_temporal(outdir: Path) -> None:
    df = pd.read_csv("results/allen_ish/allen_ish_delegate_dataset_energy_summary.tsv", sep="\t")
    df = df[df["product_name"] == "Developing Mouse Brain ISH"].copy()
    df = df[df["gene_symbol"].isin(GENE_ORDER_MOUSE)].copy()
    df["age_name"] = pd.Categorical(df["age_name"], categories=MOUSE_AGE_ORDER, ordered=True)
    df["gene_symbol"] = pd.Categorical(df["gene_symbol"], categories=GENE_ORDER_MOUSE, ordered=True)
    df = df.sort_values(["gene_symbol", "age_name"])

    fig, axes = plt.subplots(len(GENE_ORDER_MOUSE), 1, figsize=(8.5, 8.8), sharex=True)
    for ax, gene in zip(axes, GENE_ORDER_MOUSE):
        gdf = df[df["gene_symbol"] == gene]
        x = np.arange(len(MOUSE_AGE_ORDER))
        y = [
            float(gdf[gdf["age_name"] == age]["median_expression_energy"].iloc[0])
            if not gdf[gdf["age_name"] == age].empty
            else np.nan
            for age in MOUSE_AGE_ORDER
        ]
        ax.plot(
            x,
            y,
            color=GENE_COLORS_MOUSE[gene],
            marker="o",
            linewidth=1.5,
            markersize=4,
        )
        ax.set_ylabel(gene)
        ax.set_ylim(bottom=0)
        ax.grid(axis="y", color="#E6E6E6")
        ax.grid(axis="x", visible=False)
    axes[-1].set_xticks(np.arange(len(MOUSE_AGE_ORDER)))
    axes[-1].set_xticklabels(MOUSE_AGE_ORDER)
    axes[-1].set_xlabel("Allen developing mouse brain ISH age")
    fig.suptitle("Mouse temporal expression by Allen developing brain ISH", y=0.995)
    fig.text(0.02, 0.5, "Median expression energy across structures", rotation=90, va="center")
    savefig(fig, outdir, "mouse_allen_ish_temporal_energy")


def plot_mouse_adult_top_regions(outdir: Path) -> None:
    top = pd.read_csv("results/allen_ish/allen_ish_delegate_top_structures.tsv", sep="\t")
    top = top[
        (top["product_name"] == "Mouse Brain adult ISH")
        & (top["age_name"] == "P56")
        & (top["gene_symbol"].isin(GENE_ORDER_MOUSE))
    ].copy()
    top["gene_symbol"] = pd.Categorical(top["gene_symbol"], categories=GENE_ORDER_MOUSE, ordered=True)

    fig, axes = plt.subplots(len(GENE_ORDER_MOUSE), 1, figsize=(9, 11), constrained_layout=True)
    for ax, gene in zip(axes, GENE_ORDER_MOUSE):
        gdf = top[top["gene_symbol"] == gene].sort_values("expression_energy", ascending=True)
        labels = [
            f"{row.structure_acronym}: {row.structure_name[:42]}"
            for row in gdf.itertuples()
        ]
        ax.barh(
            labels,
            gdf["expression_energy"],
            color=GENE_COLORS_MOUSE[gene],
            alpha=0.85,
        )
        ax.set_title(gene)
        ax.set_xlabel("Expression energy")
        ax.grid(axis="x", color="#E6E6E6")
        ax.grid(axis="y", visible=False)
    fig.suptitle("Adult mouse P56 top ISH structures by target gene", y=1.01)
    savefig(fig, outdir, "mouse_adult_ish_top_regions")


def plot_mouse_abc_midbrain_class(outdir: Path) -> None:
    df = pd.read_csv(
        "results/allen_abc/h5ad_targets/WMB-10Xv2-MB-log2/WMB-10Xv2-MB-log2_target_gene_summary.tsv",
        sep="\t",
    )
    df = df[(df["grouping"] == "class") & (df["gene_symbol"].isin(GENE_ORDER_MOUSE))].copy()
    keep_classes = (
        df[["group", "n_cells"]]
        .drop_duplicates()
        .query("n_cells >= 20")
        .sort_values("group")["group"]
        .tolist()
    )
    df = df[df["group"].isin(keep_classes)]
    matrix = df.pivot_table(
        index="group",
        columns="gene_symbol",
        values="detected_fraction",
        observed=False,
    ).reindex(columns=GENE_ORDER_MOUSE)
    n_cells = (
        df[["group", "n_cells"]]
        .drop_duplicates()
        .set_index("group")
        .reindex(matrix.index)["n_cells"]
    )
    matrix.index = [f"{idx} (n={int(n_cells.loc[idx])})" for idx in matrix.index]

    fig, ax = plt.subplots(figsize=(8.5, 7))
    sns.heatmap(
        matrix,
        ax=ax,
        cmap="mako",
        vmin=0,
        vmax=max(0.7, float(np.nanmax(matrix.to_numpy()))),
        linewidths=0.4,
        linecolor="white",
        annot=True,
        fmt=".3f",
        cbar_kws={"label": "Detection fraction"},
    )
    ax.set_title("Allen ABC adult mouse midbrain single-cell detection by class")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=0)
    ax.tick_params(axis="y", rotation=0)
    savefig(fig, outdir, "mouse_abc_midbrain_class_detection")


def write_key_findings(outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    period = pd.read_csv("results/allen_brainspan/brainspan_target_gene_period_summary.tsv", sep="\t")
    human_sc = pd.read_csv("results/single_cell/GSE104276/GSE104276_target_gene_overall_summary.tsv", sep="\t")
    mouse_abc = pd.read_csv(
        "results/allen_abc/h5ad_targets/WMB-10Xv2-MB-log2/WMB-10Xv2-MB-log2_target_gene_summary.tsv",
        sep="\t",
    )
    mouse_abc = mouse_abc[mouse_abc["grouping"] == "overall"]
    mouse_ish = pd.read_csv("results/allen_ish/allen_ish_delegate_dataset_energy_summary.tsv", sep="\t")
    mouse_ish = mouse_ish[mouse_ish["product_name"] == "Developing Mouse Brain ISH"].copy()

    rows: list[dict[str, str]] = []
    for gene in GENE_ORDER_HUMAN:
        p = period[(period["gene_symbol"] == gene) & (period["developmental_period"] == "prenatal")]
        post = period[(period["gene_symbol"] == gene) & (period["developmental_period"] == "postnatal")]
        sc = human_sc[human_sc["gene_symbol"] == gene]
        rows.append(
            {
                "source": "Human BrainSpan + GSE104276",
                "gene": gene,
                "summary": (
                    f"BrainSpan prenatal mean RPKM {p['mean_rpkm'].iloc[0]:.4g}, "
                    f"postnatal mean RPKM {post['mean_rpkm'].iloc[0]:.4g}; "
                    f"GSE104276 fetal PFC detected fraction {sc['detected_fraction'].iloc[0]:.4g}."
                ),
            }
        )
    for gene in GENE_ORDER_MOUSE:
        abc = mouse_abc[mouse_abc["gene_symbol"] == gene]
        ish_gene = mouse_ish[mouse_ish["gene_symbol"] == gene].sort_values("age_name")
        max_age_row = ish_gene.loc[ish_gene["median_expression_energy"].idxmax()]
        rows.append(
            {
                "source": "Mouse Allen ISH + ABC midbrain",
                "gene": gene,
                "summary": (
                    f"ABC adult midbrain detected fraction {abc['detected_fraction'].iloc[0]:.4g}; "
                    f"ISH peak median expression energy at {max_age_row['age_name']} "
                    f"({max_age_row['median_expression_energy']:.4g})."
                ),
            }
        )
    with (outdir / "target_gene_key_findings.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["source", "gene", "summary"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    setup_style()
    figure_dir = Path("figures/target_gene_summary")
    summary_dir = Path("results/summary")

    plot_human_brainspan_temporal(figure_dir)
    plot_human_brainspan_region_heatmaps(figure_dir)
    plot_human_gse104276_week(figure_dir)
    plot_mouse_ish_temporal(figure_dir)
    plot_mouse_adult_top_regions(figure_dir)
    plot_mouse_abc_midbrain_class(figure_dir)
    write_key_findings(summary_dir)

    print(f"Wrote figures to {figure_dir}")
    print(f"Wrote key findings to {summary_dir / 'target_gene_key_findings.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
