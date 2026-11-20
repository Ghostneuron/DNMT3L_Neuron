"""Draw Brief Report Figure 2 from the donor permutation results."""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path)
    parser.add_argument("figures", type=Path)
    args = parser.parse_args()
    args.figures.mkdir(parents=True, exist_ok=True)
    summary = pd.read_csv(args.results / "followup_permutation.tsv", sep="\t").query('contrast == "primary"')
    draws = pd.read_csv(args.results / "donor_null_draws.tsv.gz", sep="\t").query('contrast == "primary"')
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "pdf.fonttype": 42})
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.58), sharex=True, sharey=True)
    bins = np.arange(-1.9, 1.91, 0.076)
    ticks = np.array([0.25, 0.5, 1, 2, 4])
    for ax, row in zip(axes, summary.rename(columns={"or": "odds_ratio"}).itertuples(), strict=True):
        logged = np.log(draws.loc[draws.cohort.eq(row.cohort), "or"])
        assert len(logged) == row.draws and logged.min() >= bins[0] and logged.max() <= bins[-1]
        ax.hist(logged, bins=bins, color="#AEC5D0", edgecolor="white", linewidth=.25)
        ax.axvline(0, color="#A4AAB0", linewidth=.7)
        for cutoff in (row.critical_low, row.critical_high):
            ax.axvline(np.log(cutoff), color="#5C6770", linestyle="--", linewidth=1.0)
        ax.axvline(np.log(row.odds_ratio), color="#07547A", linewidth=1.8)
        ax.set_title(row.cohort, loc="left", fontweight="bold", pad=7)
        ax.text(1, 1.05, f"OR {row.odds_ratio:.2f}   |   P = {row.p_two_sided:.3f}",
                transform=ax.transAxes, ha="right", color="#07547A")
        ax.set_xticks(np.log(ticks), [f"{x:g}" for x in ticks])
        ax.set_xlim(bins[0], bins[-1])
        ax.set_xlabel("Odds ratio (log scale)", labelpad=4)
        ax.spines[["right", "top"]].set_visible(False)
        ax.spines[["bottom", "left"]].set_color("#A3ADB2")
        ax.tick_params(length=2.5, color="#A3ADB2")
    axes[0].set_ylabel("Donor permutations")
    handles = [Line2D([], [], color="#07547A", linewidth=1.8, label="Observed OR"),
               Line2D([], [], color="#5C6770", linestyle="--", label="Two-sided P < 0.05 thresholds")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(.5, .005))
    fig.subplots_adjust(left=.09, right=.982, bottom=.31, top=.83, wspace=.16)
    for ext, options in (("png", {"dpi": 300}), ("pdf", {"metadata": {"CreationDate": None}})):
        fig.savefig(args.figures / f"Figure_2_Donor_Permutation.{ext}", **options)


if __name__ == "__main__":
    main()
