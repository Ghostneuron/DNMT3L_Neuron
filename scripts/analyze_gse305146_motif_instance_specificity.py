#!/usr/bin/env python3
"""Motif-hit versus motif-miss accessibility specificity in GSE305146."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


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
DS_FAMILY_ORDER = [
    "LHX_HOMEBOX",
    "DLX_HOMEBOX",
    "GSX_HOMEBOX",
    "ISL_HOMEBOX",
    "EMX_HOMEBOX",
    "NKX_HOMEBOX",
    "OLIG_BHLH",
    "NEURO_BHLH",
    "STAT_CONTROL",
]
CONTROL_FAMILY_ORDER = [
    "RFX_XBOX",
    "MEF2",
    "STAT_CONTROL",
]
SOURCE_ORDER = ["DS_open_recurrent", "Control_open_recurrent"]
FAMILY_ORDER_BY_SOURCE = {
    "DS_open_recurrent": DS_FAMILY_ORDER,
    "Control_open_recurrent": CONTROL_FAMILY_ORDER,
}
SOURCE_LABELS = {
    "DS_open_recurrent": "DS open recurrent",
    "Control_open_recurrent": "Control open recurrent",
}
SOURCE_N_REGIONS = 5000


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--motif-donor",
        default="results/mechanism_question/gse305146_marker_atac_motif_instance_peakset_donor_pseudobulk.tsv",
    )
    parser.add_argument(
        "--source-donor",
        default="results/mechanism_question/gse305146_marker_atac_directional_peakset_donor_pseudobulk.tsv",
    )
    parser.add_argument(
        "--region-summary",
        default="results/mechanism_question/gse305146_marker_atac_motif_instance_region_set_summary.tsv",
    )
    parser.add_argument("--outdir", default="results/mechanism_question")
    parser.add_argument("--figdir", default="figures/mechanism_question")
    parser.add_argument("--note", default="notes/gse305146_marker_atac_motif_instance_specificity.md")
    parser.add_argument("--source-n-regions", type=int, default=SOURCE_N_REGIONS)
    return parser.parse_args()


def fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}" if math.isfinite(value) else ""
    if isinstance(value, np.floating):
        number = float(value)
        return f"{number:.6g}" if math.isfinite(number) else ""
    if isinstance(value, np.integer):
        return str(int(value))
    return str(value)


def write_tsv(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(frame.columns), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in frame.to_dict(orient="records"):
            writer.writerow({field: fmt(row.get(field, "")) for field in frame.columns})


def safe_log2_per_region_cpm(counts: pd.Series, total_counts: pd.Series, n_regions: pd.Series | int) -> pd.Series:
    n = n_regions if isinstance(n_regions, pd.Series) else float(n_regions)
    cpm = counts.astype(float) / n / total_counts.astype(float) * 1_000_000
    cpm = cpm.replace([np.inf, -np.inf], np.nan).fillna(0.0)
    return np.log2(cpm + 1.0)


def order_region_sets(frame: pd.DataFrame) -> list[str]:
    ordered: list[str] = []
    for source in SOURCE_ORDER:
        for family in FAMILY_ORDER_BY_SOURCE.get(source, []):
            for region_set_id in frame.loc[
                (frame["source_peak_set"] == source) & (frame["motif_family"] == family),
                "region_set_id",
            ].drop_duplicates():
                if region_set_id not in ordered:
                    ordered.append(region_set_id)
    for region_set_id in frame["region_set_id"].drop_duplicates():
        if region_set_id not in ordered:
            ordered.append(region_set_id)
    return ordered


def summarize_effects(donor: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    group_cols = [
        "inferred_cell_class",
        "region_set_id",
        "region_set_label",
        "source_peak_set",
        "motif_family",
        "motif_family_label",
        "context_enriched",
        "n_regions_in_set",
        "n_regions_miss",
    ]
    for keys, group in donor.groupby(group_cols, dropna=False):
        record = dict(zip(group_cols, keys))
        ds = group[group["group"] == "DS"]
        control = group[group["group"] == "Control"]
        if ds.empty or control.empty:
            continue
        record.update(
            {
                "n_ds_donors": int(ds["biological_sample_id"].nunique()),
                "n_control_donors": int(control["biological_sample_id"].nunique()),
                "n_ds_cells": int(ds["n_cells"].sum()),
                "n_control_cells": int(control["n_cells"].sum()),
            }
        )
        for metric in [
            "hit_log2_per_region_peak_cpm",
            "miss_log2_per_region_peak_cpm",
            "source_log2_per_region_peak_cpm",
            "hit_minus_miss_log2_per_region_peak_cpm",
        ]:
            ds_mean = float(ds[metric].mean())
            control_mean = float(control[metric].mean())
            record[f"mean_ds_{metric}"] = ds_mean
            record[f"mean_control_{metric}"] = control_mean
            record[f"delta_ds_minus_control_{metric}"] = ds_mean - control_mean
        record["specificity_delta_hit_minus_miss"] = (
            record["delta_ds_minus_control_hit_log2_per_region_peak_cpm"]
            - record["delta_ds_minus_control_miss_log2_per_region_peak_cpm"]
        )
        record["specificity_delta_hit_minus_source"] = (
            record["delta_ds_minus_control_hit_log2_per_region_peak_cpm"]
            - record["delta_ds_minus_control_source_log2_per_region_peak_cpm"]
        )
        rows.append(record)
    return pd.DataFrame(rows)


def summarize_by_region_set(effects: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    group_cols = [
        "region_set_id",
        "region_set_label",
        "source_peak_set",
        "motif_family",
        "motif_family_label",
        "context_enriched",
        "n_regions_in_set",
        "n_regions_miss",
    ]
    for keys, group in effects[effects["inferred_cell_class"].isin(CLASS_ORDER)].groupby(group_cols, dropna=False):
        record = dict(zip(group_cols, keys))
        values = group["specificity_delta_hit_minus_miss"].astype(float)
        hit_deltas = group["delta_ds_minus_control_hit_log2_per_region_peak_cpm"].astype(float)
        miss_deltas = group["delta_ds_minus_control_miss_log2_per_region_peak_cpm"].astype(float)
        record.update(
            {
                "n_inferred_classes": int(values.notna().sum()),
                "mean_specificity_delta_hit_minus_miss": float(values.mean()),
                "median_specificity_delta_hit_minus_miss": float(values.median()),
                "n_classes_specificity_positive": int((values > 0).sum()),
                "n_classes_specificity_negative": int((values < 0).sum()),
                "mean_hit_delta": float(hit_deltas.mean()),
                "mean_miss_delta": float(miss_deltas.mean()),
                "max_specificity_class": str(group.loc[values.idxmax(), "inferred_cell_class"]) if len(values) else "",
                "max_specificity_delta": float(values.max()) if len(values) else math.nan,
                "min_specificity_class": str(group.loc[values.idxmin(), "inferred_cell_class"]) if len(values) else "",
                "min_specificity_delta": float(values.min()) if len(values) else math.nan,
            }
        )
        rows.append(record)
    summary = pd.DataFrame(rows)
    if not summary.empty:
        order = {region_set_id: index for index, region_set_id in enumerate(order_region_sets(summary))}
        summary["_order"] = summary["region_set_id"].map(order).fillna(999).astype(int)
        summary = summary.sort_values(["_order"]).drop(columns=["_order"])
    return summary


def plot_heatmap(effects: pd.DataFrame, figdir: Path) -> None:
    plot_df = effects[effects["inferred_cell_class"].isin(CLASS_ORDER)].copy()
    row_order = order_region_sets(plot_df)
    matrix = (
        plot_df.pivot_table(
            index="region_set_id",
            columns="inferred_cell_class",
            values="specificity_delta_hit_minus_miss",
            aggfunc="mean",
        )
        .reindex(index=row_order, columns=CLASS_ORDER)
    )
    labels = (
        plot_df.drop_duplicates("region_set_id")
        .set_index("region_set_id")
        .reindex(row_order)["region_set_label"]
        .fillna("")
        .tolist()
    )
    values = matrix.to_numpy(dtype=float)
    finite = values[np.isfinite(values)]
    vmax = max(0.25, min(1.0, float(np.nanpercentile(np.abs(finite), 95)) if finite.size else 1.0))

    fig_width = max(11.0, 0.7 * len(CLASS_ORDER) + 4.8)
    fig_height = max(6.0, 0.43 * max(1, len(row_order)) + 2.0)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    image = ax.imshow(matrix, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(CLASS_ORDER)))
    ax.set_xticklabels([CLASS_LABELS.get(item, item) for item in CLASS_ORDER], rotation=45, ha="right")
    ax.set_yticks(range(len(row_order)))
    ax.set_yticklabels(labels)
    ax.set_title("GSE305146 motif-hit specificity within recurrent peak sets")
    cbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Hit delta minus motif-miss delta\nlog2 per-region peak CPM")
    ax.set_xlabel("")
    ax.set_ylabel("")
    fig.tight_layout()
    figdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(figdir / "gse305146_marker_atac_motif_instance_specificity_heatmap.png", dpi=300)
    fig.savefig(figdir / "gse305146_marker_atac_motif_instance_specificity_heatmap.pdf")
    plt.close(fig)


def write_note(note_path: Path, effects: pd.DataFrame, summary: pd.DataFrame, donor: pd.DataFrame) -> None:
    note_path.parent.mkdir(parents=True, exist_ok=True)
    top_pos = effects.sort_values("specificity_delta_hit_minus_miss", ascending=False).head(10)
    top_neg = effects.sort_values("specificity_delta_hit_minus_miss", ascending=True).head(10)
    stat_summary = summary[summary["motif_family"] == "STAT_CONTROL"]
    lines = [
        "# GSE305146 Motif-Instance Specificity Within Peak Sets",
        "",
        "Question: do motif-hit peaks show stronger DS-control accessibility shifts than motif-miss peaks inside the same recurrent source peak set?",
        "",
        "## Inputs",
        "",
        f"- Donor/class motif-hit rows: {len(donor)}.",
        f"- Region sets tested: {summary['region_set_id'].nunique() if not summary.empty else 0}.",
        "- Metric: DS-control delta of log2 per-region peak CPM in motif-hit peaks minus the matched motif-miss complement from the same 5,000-peak source set.",
        "- Sign rule: positive specificity means the motif-hit subset is more DS-positive than the motif-miss complement. For control-open peak sets, this can mean the motif-hit subset is less DS-negative, not truly DS-open.",
        "- Interpretation scope: within-source motif specificity proxy, not chromVAR, not TF binding, and not independent differential accessibility.",
        "",
        "## Top Positive Specificity Effects",
        "",
    ]
    for _, row in top_pos.iterrows():
        lines.append(
            f"- {row['region_set_label']} in {row['inferred_cell_class']}: "
            f"specificity {row['specificity_delta_hit_minus_miss']:.3g}; "
            f"hit delta {row['delta_ds_minus_control_hit_log2_per_region_peak_cpm']:.3g}; "
            f"miss delta {row['delta_ds_minus_control_miss_log2_per_region_peak_cpm']:.3g}."
        )
    lines.extend(["", "## Top Negative Specificity Effects", ""])
    for _, row in top_neg.iterrows():
        lines.append(
            f"- {row['region_set_label']} in {row['inferred_cell_class']}: "
            f"specificity {row['specificity_delta_hit_minus_miss']:.3g}; "
            f"hit delta {row['delta_ds_minus_control_hit_log2_per_region_peak_cpm']:.3g}; "
            f"miss delta {row['delta_ds_minus_control_miss_log2_per_region_peak_cpm']:.3g}."
        )
    lines.extend(["", "## Region-Set Mean Specificity", ""])
    for _, row in summary.iterrows():
        lines.append(
            f"- {row['region_set_label']}: mean specificity "
            f"{row['mean_specificity_delta_hit_minus_miss']:.3g}; "
            f"{int(row['n_classes_specificity_positive'])}/{int(row['n_inferred_classes'])} classes positive; "
            f"context enriched={row['context_enriched']}."
        )
    if not stat_summary.empty:
        lines.extend(["", "## STAT Control Interpretation", ""])
        for _, row in stat_summary.iterrows():
            lines.append(
                f"- {row['region_set_label']}: mean specificity "
                f"{row['mean_specificity_delta_hit_minus_miss']:.3g}; "
                f"context enriched={row['context_enriched']}."
            )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This analysis controls for the direction of the preselected source peak set by comparing motif-hit peaks with the motif-miss complement inside the same DS-open or control-open recurrent peak set. Positive values mean the motif-hit subset is more DS-positive than the complement; negative values mean the complement carries the stronger DS-positive shift.",
            "",
            "The result should be treated as a specificity proxy. It can help prioritize motif families for true chromVAR or differential-accessibility analysis, but it does not prove direct TF activity. The extra motif-hit specificity is modest overall, which argues for a broad altered peak-set program rather than a single motif-family driver. STAT rows remain control rows because they were not enriched in the HOMER foregrounds, and their mean specificity is close to zero.",
        ]
    )
    note_path.write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    figdir = Path(args.figdir)

    motif = pd.read_csv(args.motif_donor, sep="\t")
    source = pd.read_csv(args.source_donor, sep="\t")
    region_summary = pd.read_csv(args.region_summary, sep="\t")

    source = source.rename(
        columns={
            "peak_set": "source_peak_set",
            "peak_set_counts": "source_peak_counts",
            "detected_fraction": "source_detected_fraction",
            "peak_cpm": "source_peak_cpm",
            "log2_peak_cpm": "source_log2_peak_cpm",
        }
    )
    source_cols = [
        "biological_sample_id",
        "sample_name",
        "inferred_cell_class",
        "source_peak_set",
        "source_peak_counts",
        "source_detected_fraction",
        "source_peak_cpm",
        "source_log2_peak_cpm",
    ]
    donor = motif.merge(source[source_cols], on=["biological_sample_id", "sample_name", "inferred_cell_class", "source_peak_set"], how="inner")
    donor = donor.merge(
        region_summary[["region_set_id", "n_regions_in_set"]],
        on="region_set_id",
        how="left",
        suffixes=("", "_summary"),
    )
    if "n_regions_in_set_summary" in donor.columns:
        donor["n_regions_in_set"] = donor["n_regions_in_set"].fillna(donor["n_regions_in_set_summary"])
        donor = donor.drop(columns=["n_regions_in_set_summary"])

    numeric_cols = [
        "n_regions_in_set",
        "region_set_peak_counts",
        "source_peak_counts",
        "total_peak_counts",
        "n_cells",
    ]
    for col in numeric_cols:
        donor[col] = pd.to_numeric(donor[col], errors="coerce")
    donor["n_regions_miss"] = args.source_n_regions - donor["n_regions_in_set"]
    donor["miss_peak_counts"] = donor["source_peak_counts"] - donor["region_set_peak_counts"]
    donor.loc[donor["miss_peak_counts"] < 0, "miss_peak_counts"] = 0.0
    donor.loc[donor["n_regions_miss"] <= 0, "n_regions_miss"] = np.nan

    donor["hit_log2_per_region_peak_cpm"] = safe_log2_per_region_cpm(
        donor["region_set_peak_counts"], donor["total_peak_counts"], donor["n_regions_in_set"]
    )
    donor["miss_log2_per_region_peak_cpm"] = safe_log2_per_region_cpm(
        donor["miss_peak_counts"], donor["total_peak_counts"], donor["n_regions_miss"]
    )
    donor["source_log2_per_region_peak_cpm"] = safe_log2_per_region_cpm(
        donor["source_peak_counts"], donor["total_peak_counts"], args.source_n_regions
    )
    donor["hit_minus_miss_log2_per_region_peak_cpm"] = (
        donor["hit_log2_per_region_peak_cpm"] - donor["miss_log2_per_region_peak_cpm"]
    )

    donor_fields = [
        "biological_sample_id",
        "sample_name",
        "group",
        "genotype",
        "dev_stage_pcw",
        "sex",
        "inferred_cell_class",
        "region_set_id",
        "region_set_label",
        "source_peak_set",
        "motif_family",
        "motif_family_label",
        "context_enriched",
        "n_regions_in_set",
        "n_regions_miss",
        "n_cells",
        "total_peak_counts",
        "region_set_peak_counts",
        "miss_peak_counts",
        "source_peak_counts",
        "hit_log2_per_region_peak_cpm",
        "miss_log2_per_region_peak_cpm",
        "source_log2_per_region_peak_cpm",
        "hit_minus_miss_log2_per_region_peak_cpm",
    ]
    donor = donor[donor_fields]
    effects = summarize_effects(donor)
    summary = summarize_by_region_set(effects)

    write_tsv(outdir / "gse305146_marker_atac_motif_instance_specificity_donor.tsv", donor)
    write_tsv(outdir / "gse305146_marker_atac_motif_instance_specificity_effects.tsv", effects)
    write_tsv(outdir / "gse305146_marker_atac_motif_instance_specificity_summary.tsv", summary)
    plot_heatmap(effects, figdir)
    write_note(Path(args.note), effects, summary, donor)

    print(f"Wrote motif-instance specificity summaries for {len(effects)} class-by-region-set effects")


if __name__ == "__main__":
    main()
