"""Regenerate the Brief Report figures from the frozen summary tables."""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


NAVY = "#154763"
TEAL = "#247F75"
GRAY = "#67757D"
LIGHT = "#C4C9D1"
LEVELS = ["GSE74486", "GSE73747", "GSE73747_frontal", "GSE73747_temporal", "GSE73747_deposited_beta"]
LABELS = ["GSE74486 cerebrum", "GSE73747 cortex", "  Frontal cortex", "  Temporal cortex", "  Deposited beta values"]


def reading(row):
    return f"{row.estimate:.2f} ({row.ci_low:.2f}–{row.ci_high:.2f})"


def base(labels, xlabel, xlimits, ticks, height):
    fig, ax = plt.subplots(figsize=(7.2, height))
    fig.subplots_adjust(left=.30, right=.72, top=.83, bottom=.24)
    ax.set_xscale("log")
    ax.set_xlim(xlimits)
    ax.set_xticks(ticks, labels=[str(t) for t in ticks])
    ax.minorticks_off()
    ax.set_ylim(len(labels) - .5, -.5)
    ax.set_yticks(range(len(labels)), labels=labels)
    ax.tick_params(axis="y", length=0, pad=10)
    ax.tick_params(axis="x", length=3, color=GRAY)
    ax.set_xlabel(xlabel, labelpad=7)
    ax.grid(axis="x", color="#E5E9EC", linewidth=.7)
    ax.axvline(1, linewidth=1, linestyle="--", color=GRAY)
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color(LIGHT)
    return fig, ax


def right_text(ax, heading, y, txt, color):
    if y == 0:
        ax.text(1.06, 1.13, heading, transform=ax.transAxes, weight="bold", fontsize=10)
    ax.text(1.06, y, txt, transform=ax.get_yaxis_transform(), va="center", color=color, fontsize=10)


def save(fig, output, name):
    for kind in ["png", "pdf"]:
        fig.savefig(output / f"{name}.{kind}", dpi=360, facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("tables", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
    estimates = pd.read_csv(args.tables / "estimates.tsv", sep="\t").set_index(["context", "statistic"])
    null = pd.read_csv(args.tables / "null_and_coverage.tsv", sep="\t").set_index("context")

    fig, ax = base(LABELS, "DNMT3L direction odds ratio", (.52, 3.55), [.6, 1, 1.5, 2, 3], 3.12)
    for i, cohort in enumerate(LEVELS):
        row = estimates.loc[cohort, "DNMT3L_OR"]
        color = NAVY if i < 2 else GRAY
        ax.hlines(i, row.ci_low, row.ci_high, color=color, lw=1.8)
        ax.scatter(row.estimate, i, color=color, s=31)
        right_text(ax, "OR (95% CI)", i, reading(row), color)
    save(fig, args.output, "Figure_1_Donor_ORs")

    fig, ax = base(LABELS[:2], "Direction odds ratio", (.72, 1.98), [.8, 1, 1.25, 1.5, 1.8], 2.58)
    for i, cohort in enumerate(LEVELS[:2]):
        lo, med, hi = null.loc[cohort, ["matched_CpG_null_OR_2.5%", "matched_CpG_null_OR_50%", "matched_CpG_null_OR_97.5%"]]
        obs = estimates.loc[(cohort, "DNMT3L_matchable_OR"), "estimate"]
        ax.hlines(i, lo, hi, color=GRAY, lw=2)
        ax.scatter(med, i, facecolor="white", edgecolor=GRAY, s=31, zorder=3)
        ax.scatter(obs, i, marker="D", color=NAVY, s=31, zorder=3)
        right_text(ax, "Observed OR", i, f"{obs:.2f}", NAVY)
    fig.text(.30, .02, "Gray line and circle   Matched CpG 95% range and median\nBlue diamond             Observed DNMT3L OR", fontsize=8, color=GRAY)
    fig.subplots_adjust(bottom=.36)
    save(fig, args.output, "Figure_2_Matched_CpG_Benchmark")

    indices = [LEVELS[0], LEVELS[1], LEVELS[4]]
    fig, ax = base([LABELS[0], LABELS[1], LABELS[4]], "Ratio of ORs (DNMT3L / DNMT3A)", (.66, 4.5), [.8, 1, 1.5, 2, 3, 4], 2.66)
    for i, cohort in enumerate(indices):
        row = estimates.loc[cohort, "DNMT3L_vs_DNMT3A_ratio_of_ORs"]
        color = TEAL if i < 2 else GRAY
        ax.hlines(i, row.ci_low, row.ci_high, color=color, lw=1.8)
        ax.scatter(row.estimate, i, color=color, s=31)
        right_text(ax, "Ratio (95% CI)", i, reading(row), color)
    save(fig, args.output, "Figure_3_DNMT3A_Benchmark")


if __name__ == "__main__":
    main()
