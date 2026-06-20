#!/usr/bin/env python3
"""Build an Extended Data methylation/HOMER QC summary figure."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


FIGDIR = Path("figures/methylation_qc")


def pct_to_float(value: object) -> float:
    text = str(value).strip().replace("%", "")
    if not text or text.lower() == "nan":
        return np.nan
    return float(text)


def short_motif(name: str) -> str:
    return str(name).split("/")[0].replace("(Homeobox)", "").replace("(HMG)", "").strip()


def plot_delta_summary(ax: plt.Axes) -> None:
    df = pd.read_csv("results/first_pass_methylation/methylation_delta_summary.tsv", sep="\t")
    keep = [
        ("GSE157078", "human_all_control_vs_DNMT3L_overexpression", "Human DNMT3L-OE"),
        ("GSE271949", "mouse_all_control_vs_Dnmt3l_overexpression", "Mouse Dnmt3l-OE all"),
        ("GSE271949", "mouse_male_control_vs_Dnmt3l_overexpression", "Mouse male"),
        ("GSE271949", "mouse_female_control_vs_Dnmt3l_overexpression", "Mouse female"),
    ]
    rows = []
    for dataset, comparison, label in keep:
        hit = df[(df["dataset"] == dataset) & (df["comparison"] == comparison)]
        if not hit.empty:
            row = hit.iloc[0].copy()
            row["label"] = label
            rows.append(row)
    plot_df = pd.DataFrame(rows)
    y = np.arange(len(plot_df))
    ax.barh(y, plot_df["mean_delta_beta"], color="#4C78A8")
    ax.axvline(0, color="0.25", lw=0.8)
    ax.set_yticks(y, plot_df["label"])
    ax.invert_yaxis()
    ax.set_xlabel("Mean delta beta")
    ax.set_title("A  DNMT3L-OE methylation skew", loc="left", fontweight="bold")
    for i, value in enumerate(plot_df["mean_delta_beta"]):
        ax.text(value + 0.002, i, f"{value:.3f}", va="center", fontsize=8)


def plot_hyper_hypo_counts(ax: plt.Axes) -> None:
    df = pd.read_csv("results/first_pass_methylation/methylation_delta_summary.tsv", sep="\t")
    df = df[df["comparison"].isin([
        "human_all_control_vs_DNMT3L_overexpression",
        "mouse_all_control_vs_Dnmt3l_overexpression",
    ])].copy()
    df["label"] = df["dataset"].map({"GSE157078": "Human DNMT3L-OE", "GSE271949": "Mouse Dnmt3l-OE"})
    x = np.arange(len(df))
    width = 0.36
    ax.bar(x - width / 2, df["hyper_delta_ge_0.05"], width, label="hyper >= 0.05", color="#E45756")
    ax.bar(x + width / 2, df["hypo_delta_le_neg_0.05"], width, label="hypo <= -0.05", color="#4C78A8")
    ax.set_xticks(x, df["label"], rotation=20, ha="right")
    ax.set_yscale("log")
    ax.set_ylabel("CpGs/probes (log scale)")
    ax.set_title("B  Hyper/hypo counts", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=8)


def plot_concordance(ax: plt.Axes) -> None:
    df = pd.read_csv(
        "results/dp10_methylation_dnmt3l_concordance/dp10_dp16_dnmt3l_cpg_concordance_summary.tsv",
        sep="\t",
    )
    x = np.arange(len(df))
    ax.bar(x, df["concordant_fraction"], color=["#72B7B2", "#B279A2"])
    ax.set_xticks(x, df["model"])
    ax.set_ylim(0, 1)
    ax.set_ylabel("Concordant fraction")
    ax.set_title("C  Mouse DS model concordance", loc="left", fontweight="bold")
    for i, value in enumerate(df["concordant_fraction"]):
        ax.text(i, value + 0.03, f"{value:.2f}", ha="center", fontsize=8)


def plot_human_convergence(ax: plt.Axes) -> None:
    df = pd.read_csv("results/human_methylation_convergence/human_methylation_convergence_summary.tsv", sep="\t")
    values = dict(zip(df["metric"], df["value"]))
    metrics = [
        ("DS tested probes", "ds_tested_probes"),
        ("Any DS concordant", "ds_any_context_concordant_probes"),
        ("Fetal DS concordant", "ds_fetal_context_concordant_probes"),
        ("Triple any DS", "human_triple_phase3_any_ds_probes"),
        ("Triple fetal DS", "human_triple_phase3_fetal_ds_probes"),
    ]
    labels = [m[0] for m in metrics]
    counts = [float(values[m[1]]) for m in metrics]
    y = np.arange(len(labels))
    ax.barh(y, counts, color="#59A14F")
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("Probe count")
    ax.set_title("D  Human convergence filters", loc="left", fontweight="bold")
    for i, value in enumerate(counts):
        ax.text(value + max(counts) * 0.015, i, f"{int(value):,}", va="center", fontsize=8)


def plot_homer_top(ax: plt.Axes) -> None:
    oe = pd.read_csv("results/homer_motif_enrichment/homer_known_motif_summary.tsv", sep="\t")
    adult = pd.read_csv(
        "results/adult_human_ds_homer_motif_enrichment_sample25k_s8/homer_known_motif_summary.tsv",
        sep="\t",
    )
    oe = oe[oe["set_id"].eq("GSE157078_human_DNMT3L_OE_DMR_all_hg19")].head(5).copy()
    adult = adult[adult["set_id"].eq("GSE63347_GSE63347_adult_all_brain_DMR_hyper_hg19")].head(5).copy()
    oe["context"] = "Human DNMT3L-OE DMR all"
    adult["context"] = "Adult DS brain hyper DMR"
    df = pd.concat([oe, adult], ignore_index=True)
    df["motif_short"] = df["motif_name"].map(short_motif)
    df["score"] = -np.log10(df["p_value"].astype(float).clip(lower=1e-300))
    df["label"] = df["context"] + " | " + df["motif_short"]
    df = df.sort_values("score")
    colors = df["context"].map({
        "Human DNMT3L-OE DMR all": "#F58518",
        "Adult DS brain hyper DMR": "#9467BD",
    })
    ax.barh(np.arange(len(df)), df["score"], color=colors)
    ax.set_yticks(np.arange(len(df)), df["label"])
    ax.set_xlabel("-log10(HOMER p)")
    ax.set_title("E  Top known HOMER motifs", loc="left", fontweight="bold")


def main() -> None:
    FIGDIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.facecolor": "white",
    })
    fig = plt.figure(figsize=(13, 9), constrained_layout=True)
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.2])
    plot_delta_summary(fig.add_subplot(gs[0, 0]))
    plot_hyper_hypo_counts(fig.add_subplot(gs[0, 1]))
    plot_concordance(fig.add_subplot(gs[1, 0]))
    plot_human_convergence(fig.add_subplot(gs[1, 1]))
    plot_homer_top(fig.add_subplot(gs[2, :]))
    fig.suptitle("Extended Data Figure 4. DNMT3L-linked methylation and HOMER QC", fontweight="bold")
    out = FIGDIR / "extended_data_methylation_homer_qc"
    fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    print(out.with_suffix(".png"))
    print(out.with_suffix(".pdf"))


if __name__ == "__main__":
    main()
