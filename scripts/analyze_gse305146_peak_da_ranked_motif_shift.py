#!/usr/bin/env python3
"""Ranked motif-family shifts from GSE305146 peak-level DA betas."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from analyze_gse305146_marker_atac_motif_deviation import CLASS_LABELS, CLASS_ORDER
from analyze_gse305146_marker_atac_motif_instance_peaksets import FAMILY_ORDER, MOTIF_FAMILY_LABELS
from analyze_gse305146_marker_atac_peak_da import bh_fdr, read_motif_family_masks


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--effects",
        default="results/mechanism_question/gse305146_marker_atac_peak_da_effects.tsv.gz",
    )
    parser.add_argument(
        "--regions",
        default="results/mechanism_question/gse305146_marker_atac_peak_da_tested_regions.tsv",
    )
    parser.add_argument(
        "--motif-hits",
        default="results/mechanism_question/gse305146_marker_atac_motif_deviation_region_hits.tsv",
    )
    parser.add_argument("--outdir", default="results/mechanism_question")
    parser.add_argument("--figdir", default="figures/mechanism_question")
    parser.add_argument("--note", default="notes/gse305146_peak_da_ranked_motif_shift_first_pass.md")
    parser.add_argument("--fdr-threshold", type=float, default=0.05)
    parser.add_argument("--min-abs-delta-beta", type=float, default=0.02)
    return parser.parse_args()


def fmt_float(value: float) -> str:
    return f"{value:.3g}" if math.isfinite(value) else "NA"


def ranked_motif_shift(args: argparse.Namespace) -> pd.DataFrame:
    regions = pd.read_csv(args.regions, sep="\t")
    motif_masks = read_motif_family_masks(Path(args.motif_hits), regions)
    region_order = regions[["peak_id", "region_index"]].copy()
    effects = pd.read_csv(args.effects, sep="\t")
    effects = effects.merge(region_order, on="peak_id", how="left")
    rows: list[dict[str, object]] = []
    for cell_class in CLASS_ORDER:
        sub = effects[effects["inferred_cell_class"] == cell_class].copy()
        if sub.empty:
            continue
        sub = sub.sort_values("region_index")
        beta = pd.to_numeric(sub["pcw_adjusted_group_beta"], errors="coerce").to_numpy(dtype=float)
        finite = np.isfinite(beta)
        for family in FAMILY_ORDER:
            motif_mask_full = motif_masks.get(family)
            if motif_mask_full is None:
                continue
            motif_mask = motif_mask_full[sub["region_index"].to_numpy(dtype=int)] & finite
            non_mask = (~motif_mask_full[sub["region_index"].to_numpy(dtype=int)]) & finite
            motif_beta = beta[motif_mask]
            non_beta = beta[non_mask]
            if motif_beta.size < 10 or non_beta.size < 10:
                pvalue = 1.0
                u_stat = math.nan
                rank_biserial = math.nan
            else:
                test = mannwhitneyu(motif_beta, non_beta, alternative="two-sided")
                u_stat = float(test.statistic)
                pvalue = float(test.pvalue)
                rank_biserial = (2.0 * u_stat / (motif_beta.size * non_beta.size)) - 1.0
            rows.append(
                {
                    "inferred_cell_class": cell_class,
                    "motif_family": family,
                    "motif_family_label": MOTIF_FAMILY_LABELS.get(family, family),
                    "n_motif_peaks": int(motif_beta.size),
                    "n_nonmotif_peaks": int(non_beta.size),
                    "mean_beta_motif": float(np.nanmean(motif_beta)) if motif_beta.size else math.nan,
                    "mean_beta_nonmotif": float(np.nanmean(non_beta)) if non_beta.size else math.nan,
                    "delta_mean_beta_motif_minus_nonmotif": (
                        float(np.nanmean(motif_beta) - np.nanmean(non_beta)) if motif_beta.size and non_beta.size else math.nan
                    ),
                    "median_beta_motif": float(np.nanmedian(motif_beta)) if motif_beta.size else math.nan,
                    "median_beta_nonmotif": float(np.nanmedian(non_beta)) if non_beta.size else math.nan,
                    "delta_median_beta_motif_minus_nonmotif": (
                        float(np.nanmedian(motif_beta) - np.nanmedian(non_beta)) if motif_beta.size and non_beta.size else math.nan
                    ),
                    "mannwhitney_u": u_stat,
                    "rank_biserial": rank_biserial,
                    "mannwhitney_p": pvalue,
                }
            )
    out = pd.DataFrame(rows)
    out["fdr_bh"] = bh_fdr(pd.to_numeric(out["mannwhitney_p"], errors="coerce").to_numpy())
    out["shift_direction"] = np.where(out["delta_mean_beta_motif_minus_nonmotif"] > 0, "DS_beta_positive", "DS_beta_negative")
    out["significant_ranked_shift"] = (
        (out["fdr_bh"] <= args.fdr_threshold)
        & (out["delta_mean_beta_motif_minus_nonmotif"].abs() >= args.min_abs_delta_beta)
    )
    return out


def plot_delta_heatmap(results: pd.DataFrame, figdir: Path) -> None:
    figdir.mkdir(parents=True, exist_ok=True)
    matrix = np.full((len(FAMILY_ORDER), len(CLASS_ORDER)), np.nan)
    for family_idx, family in enumerate(FAMILY_ORDER):
        for class_idx, cell_class in enumerate(CLASS_ORDER):
            sub = results[(results["motif_family"] == family) & (results["inferred_cell_class"] == cell_class)]
            if not sub.empty:
                matrix[family_idx, class_idx] = float(sub["delta_mean_beta_motif_minus_nonmotif"].iloc[0])
    finite = matrix[np.isfinite(matrix)]
    vmax = max(0.03, float(np.nanpercentile(np.abs(finite), 95))) if finite.size else 0.05
    vmax = min(vmax, 0.15)
    fig, ax = plt.subplots(figsize=(11, 5.6))
    im = ax.imshow(matrix, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(np.arange(len(CLASS_ORDER)))
    ax.set_xticklabels([CLASS_LABELS.get(name, name) for name in CLASS_ORDER], rotation=45, ha="right")
    ax.set_yticks(np.arange(len(FAMILY_ORDER)))
    ax.set_yticklabels([MOTIF_FAMILY_LABELS.get(family, family) for family in FAMILY_ORDER])
    ax.set_title("Peak-level ranked motif beta shifts")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Mean PCW-adjusted beta: motif peaks minus nonmotif peaks")
    ax.set_xticks(np.arange(-0.5, len(CLASS_ORDER), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(FAMILY_ORDER), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)
    fig.tight_layout()
    fig.savefig(figdir / "gse305146_peak_da_ranked_motif_shift_delta_heatmap.png", dpi=300)
    fig.savefig(figdir / "gse305146_peak_da_ranked_motif_shift_delta_heatmap.pdf")
    plt.close(fig)


def plot_signed_fdr_heatmap(results: pd.DataFrame, figdir: Path) -> None:
    matrix = np.zeros((len(FAMILY_ORDER), len(CLASS_ORDER)), dtype=float)
    for family_idx, family in enumerate(FAMILY_ORDER):
        for class_idx, cell_class in enumerate(CLASS_ORDER):
            sub = results[(results["motif_family"] == family) & (results["inferred_cell_class"] == cell_class)]
            if sub.empty:
                continue
            row = sub.iloc[0]
            sign = 1.0 if float(row["delta_mean_beta_motif_minus_nonmotif"]) > 0 else -1.0
            value = min(-math.log10(max(float(row["fdr_bh"]), 1e-300)), 20.0)
            matrix[family_idx, class_idx] = sign * value
    fig, ax = plt.subplots(figsize=(11, 5.6))
    im = ax.imshow(matrix, aspect="auto", cmap="RdBu_r", vmin=-20, vmax=20)
    ax.set_xticks(np.arange(len(CLASS_ORDER)))
    ax.set_xticklabels([CLASS_LABELS.get(name, name) for name in CLASS_ORDER], rotation=45, ha="right")
    ax.set_yticks(np.arange(len(FAMILY_ORDER)))
    ax.set_yticklabels([MOTIF_FAMILY_LABELS.get(family, family) for family in FAMILY_ORDER])
    ax.set_title("Ranked motif shift significance")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("signed -log10(FDR): red positive beta, blue negative beta")
    ax.set_xticks(np.arange(-0.5, len(CLASS_ORDER), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(FAMILY_ORDER), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)
    fig.tight_layout()
    fig.savefig(figdir / "gse305146_peak_da_ranked_motif_shift_fdr_heatmap.png", dpi=300)
    fig.savefig(figdir / "gse305146_peak_da_ranked_motif_shift_fdr_heatmap.pdf")
    plt.close(fig)


def write_note(path: Path, results: pd.DataFrame, args: argparse.Namespace) -> None:
    sig = results[results["significant_ranked_shift"]].copy()
    top_pos = sig.sort_values("delta_mean_beta_motif_minus_nonmotif", ascending=False).head(12)
    top_neg = sig.sort_values("delta_mean_beta_motif_minus_nonmotif", ascending=True).head(12)
    family_summary = (
        sig.groupby(["motif_family", "motif_family_label", "shift_direction"], dropna=False, observed=True)
        .size()
        .reset_index(name="n_classes")
        .sort_values("n_classes", ascending=False)
    )
    lines = [
        "# GSE305146 Peak-Level Ranked Motif Shift",
        "",
        "Question: even when few individual peaks pass strict DA FDR, do motif-family-bearing recurrent peaks show systematic PCW-adjusted DS beta shifts within inferred cell classes?",
        "",
        "## Methods Snapshot",
        "",
        "- Input: PCW-adjusted group betas from the donor-level recurrent-peak DA model.",
        "- Test: Mann-Whitney comparison of beta distributions for motif-bearing peaks versus nonmotif peaks within each inferred cell class.",
        f"- Significant ranked shift: FDR <= {args.fdr_threshold:g} and abs(mean beta difference) >= {args.min_abs_delta_beta:g}.",
        "",
        "## Recurrent Significant Shift Counts",
        "",
    ]
    if family_summary.empty:
        lines.append("- No motif-family ranked shifts passed the reporting threshold.")
    else:
        for _, row in family_summary.head(12).iterrows():
            lines.append(
                f"- {row['motif_family_label']} {row['shift_direction']}: {int(row['n_classes'])} classes."
            )
    lines.extend(["", "## Strongest Positive Motif Beta Shifts", ""])
    if top_pos.empty:
        lines.append("- None.")
    else:
        for _, row in top_pos.iterrows():
            lines.append(
                f"- {row['motif_family_label']} in {row['inferred_cell_class']}: "
                f"delta mean beta {fmt_float(float(row['delta_mean_beta_motif_minus_nonmotif']))}, "
                f"FDR {fmt_float(float(row['fdr_bh']))}."
            )
    lines.extend(["", "## Strongest Negative Motif Beta Shifts", ""])
    if top_neg.empty:
        lines.append("- None.")
    else:
        for _, row in top_neg.iterrows():
            lines.append(
                f"- {row['motif_family_label']} in {row['inferred_cell_class']}: "
                f"delta mean beta {fmt_float(float(row['delta_mean_beta_motif_minus_nonmotif']))}, "
                f"FDR {fmt_float(float(row['fdr_bh']))}."
            )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This ranked analysis complements the strict DA peak table. It is more sensitive to distributed motif-linked accessibility shifts, but less specific than individual peak FDR calls. Therefore it is best used to support broad regulatory-context interpretation, not to claim direct TF binding or methylation independence.",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    figdir = Path(args.figdir)
    outdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)
    results = ranked_motif_shift(args)
    outpath = outdir / "gse305146_peak_da_ranked_motif_shift.tsv"
    results.to_csv(outpath, sep="\t", index=False, lineterminator="\n")
    plot_delta_heatmap(results, figdir)
    plot_signed_fdr_heatmap(results, figdir)
    write_note(Path(args.note), results, args)
    print("Wrote ranked motif shift results")
    print(outpath)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
