#!/usr/bin/env python3
"""Matched-background motif accessibility validation for GSE305146."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import random
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analyze_gse305146_marker_atac_motif_deviation import CLASS_LABELS, CLASS_ORDER, FAMILY_ORDER
from analyze_gse305146_marker_atac_motif_instance_peaksets import MOTIF_FAMILY_LABELS
from infer_gse305146_marker_cell_classes_and_targets import (
    MARKER_SETS,
    infer_classes,
    parse_features as parse_marker_features,
    selected_libraries,
    stream_matrix,
)


REGION_SIZE = 200
SET_FOREGROUND = "motif_foreground"
SET_MATCHED = "matched_background"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", default="data/geo/GSE305146/metadata/sample_metadata.tsv")
    parser.add_argument("--supplementary", default="data/geo/GSE305146/supplementary")
    parser.add_argument(
        "--background-table",
        default="results/gse305146_scATAC_homer_motif_enrichment/tables/GSE305146_recurrent_peak_background.tsv",
    )
    parser.add_argument(
        "--motif-hits",
        default="results/mechanism_question/gse305146_marker_atac_motif_deviation_region_hits.tsv",
    )
    parser.add_argument("--outdir", default="results/mechanism_question")
    parser.add_argument("--figdir", default="figures/mechanism_question")
    parser.add_argument("--note", default="notes/gse305146_marker_atac_matched_motif_deviation_first_pass.md")
    parser.add_argument("--region-size", type=int, default=REGION_SIZE)
    parser.add_argument("--bin-size", type=int, default=100_000)
    parser.add_argument("--min-total-rna", type=float, default=500.0)
    parser.add_argument("--min-detected-genes", type=int, default=200)
    parser.add_argument("--min-top-score", type=float, default=0.05)
    parser.add_argument("--min-confidence-delta", type=float, default=0.05)
    parser.add_argument("--min-cells-per-donor-class", type=int, default=20)
    parser.add_argument("--seed", type=int, default=305146)
    parser.add_argument("--limit-libraries", type=int, help="Debug: process first N final-analysis libraries.")
    parser.add_argument("--limit-families", type=int, help="Debug: process first N motif families.")
    parser.add_argument("--plot-only", action="store_true", help="Regenerate figures and note from existing output tables.")
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


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})


def safe_float(value: object) -> float:
    try:
        return float(str(value))
    except ValueError:
        return math.nan


def safe_int(value: object) -> int:
    try:
        return int(float(str(value)))
    except ValueError:
        return 0


def read_background(path: Path, region_size: int) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    half = region_size // 2
    frame["center"] = pd.to_numeric(frame["center"], errors="coerce").astype(int)
    frame["chrom"] = frame["chrom"].astype(str)
    frame = frame[frame["chrom"].str.startswith("chr")].copy()
    frame["match_start"] = (frame["center"] - half).clip(lower=0).astype(int)
    frame["match_end"] = frame["match_start"] + region_size
    frame["peak_id"] = (
        "background:"
        + frame["chrom"].astype(str)
        + ":"
        + frame["match_start"].astype(str)
        + "-"
        + frame["match_end"].astype(str)
    )
    for column in ["width", "n_source_peaks", "ds_frequency", "control_frequency", "delta_frequency_ds_minus_control", "max_frequency"]:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame = frame.drop_duplicates("peak_id").reset_index(drop=True)
    return frame


def add_match_bins(background: pd.DataFrame) -> pd.DataFrame:
    frame = background.copy()
    frame["width_bin"] = pd.qcut(frame["width"], q=6, labels=False, duplicates="drop").astype(int)
    frame["source_peak_bin"] = pd.qcut(frame["n_source_peaks"], q=6, labels=False, duplicates="drop").astype(int)
    frame["max_frequency_bin"] = pd.cut(
        frame["max_frequency"],
        bins=[-0.001, 0.1, 0.2, 0.35, 0.5, 0.75, 1.001],
        labels=False,
        include_lowest=True,
    ).astype(int)
    frame["delta_frequency_bin"] = pd.qcut(frame["delta_frequency_ds_minus_control"], q=6, labels=False, duplicates="drop").astype(int)
    return frame


def read_motif_foregrounds(path: Path, families: list[str]) -> dict[str, set[str]]:
    hits = pd.read_csv(path, sep="\t", usecols=["peak_id", "motif_family"])
    foregrounds: dict[str, set[str]] = {}
    for family in families:
        foregrounds[family] = set(hits.loc[hits["motif_family"] == family, "peak_id"].dropna().astype(str))
    return foregrounds


def choose_from_bins(
    foreground: pd.DataFrame,
    candidates: pd.DataFrame,
    rng: random.Random,
) -> tuple[list[dict[str, object]], set[str]]:
    tiers = [
        ["max_frequency_bin", "width_bin", "delta_frequency_bin", "source_peak_bin"],
        ["max_frequency_bin", "width_bin", "delta_frequency_bin"],
        ["max_frequency_bin", "width_bin"],
        ["max_frequency_bin"],
    ]
    pools_by_tier: list[dict[tuple[object, ...], list[str]]] = []
    for tier in tiers:
        pools: dict[tuple[object, ...], list[str]] = defaultdict(list)
        for row in candidates.to_dict(orient="records"):
            pools[tuple(row[column] for column in tier)].append(str(row["peak_id"]))
        for values in pools.values():
            rng.shuffle(values)
        pools_by_tier.append(pools)

    all_candidate_ids = candidates["peak_id"].astype(str).tolist()
    rng.shuffle(all_candidate_ids)
    used: set[str] = set()
    matches: list[dict[str, object]] = []
    candidate_set = set(all_candidate_ids)

    for row in foreground.sample(frac=1.0, random_state=rng.randint(1, 1_000_000_000)).to_dict(orient="records"):
        matched_id = ""
        tier_used = ""
        replacement = "no"
        for tier, pools in zip(tiers, pools_by_tier):
            key = tuple(row[column] for column in tier)
            pool = pools.get(key, [])
            while pool and pool[-1] in used:
                pool.pop()
            if pool:
                matched_id = pool.pop()
                tier_used = "+".join(tier)
                break
        if not matched_id:
            remaining = [peak_id for peak_id in all_candidate_ids if peak_id not in used]
            if remaining:
                matched_id = remaining.pop()
                tier_used = "global_unused"
            elif all_candidate_ids:
                matched_id = rng.choice(all_candidate_ids)
                tier_used = "global_replacement"
                replacement = "yes"
            else:
                continue
        if matched_id not in candidate_set:
            continue
        used.add(matched_id)
        matches.append(
            {
                "motif_family": row["motif_family"],
                "foreground_peak_id": row["peak_id"],
                "matched_peak_id": matched_id,
                "match_tier": tier_used,
                "with_replacement": replacement,
            }
        )
    return matches, used


def build_matched_sets(
    background: pd.DataFrame,
    foregrounds: dict[str, set[str]],
    families: list[str],
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    match_rows: list[dict[str, object]] = []
    qc_rows: list[dict[str, object]] = []
    region_lookup = background.set_index("peak_id").to_dict(orient="index")

    for family_index, family in enumerate(families):
        rng = random.Random(seed + family_index)
        fg_ids = foregrounds.get(family, set()) & set(region_lookup)
        foreground = background[background["peak_id"].isin(fg_ids)].copy()
        foreground["motif_family"] = family
        candidates = background[~background["peak_id"].isin(fg_ids)].copy()
        matches, matched_ids = choose_from_bins(foreground, candidates, rng)
        match_rows.extend(matches)

        for peak_id in sorted(fg_ids):
            region = dict(region_lookup[peak_id])
            rows.append(
                {
                    **region,
                    "peak_id": peak_id,
                    "region_set_id": f"MOTIF_FAMILY|{family}|FOREGROUND",
                    "region_set_type": SET_FOREGROUND,
                    "motif_family": family,
                    "motif_family_label": MOTIF_FAMILY_LABELS.get(family, family),
                }
            )
        for peak_id in sorted(matched_ids):
            region = dict(region_lookup[peak_id])
            rows.append(
                {
                    **region,
                    "peak_id": peak_id,
                    "region_set_id": f"MOTIF_FAMILY|{family}|MATCHED",
                    "region_set_type": SET_MATCHED,
                    "motif_family": family,
                    "motif_family_label": MOTIF_FAMILY_LABELS.get(family, family),
                }
            )

        fg_metrics = foreground[["width", "n_source_peaks", "max_frequency", "delta_frequency_ds_minus_control"]].mean()
        matched = background[background["peak_id"].isin(matched_ids)]
        matched_metrics = matched[["width", "n_source_peaks", "max_frequency", "delta_frequency_ds_minus_control"]].mean()
        qc = {
            "motif_family": family,
            "motif_family_label": MOTIF_FAMILY_LABELS.get(family, family),
            "n_foreground_peaks": len(fg_ids),
            "n_matched_unique_peaks": len(matched_ids),
            "n_match_rows": len(matches),
            "n_replacement_rows": sum(1 for row in matches if row["with_replacement"] == "yes"),
        }
        for metric in ["width", "n_source_peaks", "max_frequency", "delta_frequency_ds_minus_control"]:
            fg_value = float(fg_metrics.get(metric, math.nan))
            bg_value = float(matched_metrics.get(metric, math.nan))
            qc[f"foreground_mean_{metric}"] = fg_value
            qc[f"matched_mean_{metric}"] = bg_value
            qc[f"mean_diff_{metric}"] = fg_value - bg_value
        qc_rows.append(qc)

    region_sets = pd.DataFrame(rows).drop_duplicates(["region_set_id", "peak_id"])
    return region_sets, pd.DataFrame(match_rows), pd.DataFrame(qc_rows)


def build_region_bins(
    region_sets: pd.DataFrame,
    region_set_index: dict[str, int],
    bin_size: int,
) -> dict[str, dict[int, list[tuple[int, int, int]]]]:
    bins: dict[str, dict[int, list[tuple[int, int, int]]]] = defaultdict(lambda: defaultdict(list))
    for row in region_sets.to_dict(orient="records"):
        chrom = str(row["chrom"])
        start = int(row["match_start"])
        end = int(row["match_end"])
        set_index = region_set_index[str(row["region_set_id"])]
        for bin_id in range(start // bin_size, end // bin_size + 1):
            bins[chrom][bin_id].append((start, end, set_index))
    return bins


def map_peak_rows_to_sets(
    feature_path: Path,
    bins: dict[str, dict[int, list[tuple[int, int, int]]]],
    bin_size: int,
) -> tuple[dict[int, list[int]], int, int, dict[str, int]]:
    row_to_sets: dict[int, list[int]] = {}
    last_gene_expression_row = 0
    n_peak_rows = 0
    mapped_by_set: dict[str, int] = defaultdict(int)
    with gzip.open(feature_path, "rt", errors="replace") as handle:
        for row_index, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            feature_type = fields[2]
            if feature_type == "Gene Expression":
                last_gene_expression_row = row_index
                continue
            if feature_type != "Peaks":
                continue
            n_peak_rows += 1
            chrom = fields[3]
            start = safe_int(fields[4])
            end = safe_int(fields[5])
            if end <= start or chrom not in bins:
                continue
            sets: set[int] = set()
            for bin_id in range(start // bin_size, end // bin_size + 1):
                for region_start, region_end, set_index in bins[chrom].get(bin_id, []):
                    if end >= region_start and start <= region_end:
                        sets.add(set_index)
            if sets:
                ordered = sorted(sets)
                row_to_sets[row_index] = ordered
                for set_index in ordered:
                    mapped_by_set[str(set_index)] += 1
    return row_to_sets, last_gene_expression_row, n_peak_rows, dict(mapped_by_set)


def stream_atac_sets(
    matrix_path: Path,
    last_gene_expression_row: int,
    row_to_sets: dict[int, list[int]],
    n_sets: int,
) -> tuple[np.ndarray, np.ndarray]:
    total_peak_counts: np.ndarray | None = None
    set_peak_counts: np.ndarray | None = None
    dims_seen = False
    with gzip.open(matrix_path, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            fields = line.strip().split()
            if not fields:
                continue
            if not dims_seen:
                n_cells = int(fields[1])
                total_peak_counts = np.zeros(n_cells, dtype=np.float64)
                set_peak_counts = np.zeros((n_cells, n_sets), dtype=np.float32)
                dims_seen = True
                continue
            assert total_peak_counts is not None
            assert set_peak_counts is not None
            row_index = int(fields[0])
            if row_index <= last_gene_expression_row:
                continue
            cell_index = int(fields[1]) - 1
            value = float(fields[2])
            total_peak_counts[cell_index] += value
            for set_index in row_to_sets.get(row_index, []):
                set_peak_counts[cell_index, set_index] += value
    if total_peak_counts is None or set_peak_counts is None:
        raise RuntimeError(f"No Matrix Market dimensions found in {matrix_path}")
    return total_peak_counts, set_peak_counts


def region_set_summary(region_sets: pd.DataFrame) -> pd.DataFrame:
    summary = (
        region_sets.groupby(["region_set_id", "region_set_type", "motif_family", "motif_family_label"], dropna=False, observed=True)
        .agg(n_regions=("peak_id", "nunique"))
        .reset_index()
    )
    order = {family: index for index, family in enumerate(FAMILY_ORDER)}
    summary["_order"] = summary["motif_family"].map(order).fillna(999).astype(int)
    summary["_type_order"] = summary["region_set_type"].map({SET_FOREGROUND: 0, SET_MATCHED: 1}).fillna(9).astype(int)
    return summary.sort_values(["_order", "_type_order"]).drop(columns=["_order", "_type_order"])


def summarize_library(
    library: dict[str, str],
    classes: np.ndarray,
    total_peak_counts: np.ndarray,
    set_peak_counts: np.ndarray,
    set_summary: pd.DataFrame,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    class_values = classes.astype(str)
    set_rows = set_summary.to_dict(orient="records")
    for cell_class in sorted(set(class_values)):
        mask = class_values == cell_class
        n_cells = int(mask.sum())
        if n_cells == 0:
            continue
        class_total_peak = float(total_peak_counts[mask].sum())
        for set_index, region_set in enumerate(set_rows):
            set_count = float(set_peak_counts[mask, set_index].sum())
            n_regions = int(region_set["n_regions"])
            per_region_peak_cpm = set_count / n_regions / class_total_peak * 1_000_000 if class_total_peak > 0 and n_regions else math.nan
            rows.append(
                {
                    **library,
                    "inferred_cell_class": cell_class,
                    "region_set_id": region_set["region_set_id"],
                    "region_set_type": region_set["region_set_type"],
                    "motif_family": region_set["motif_family"],
                    "motif_family_label": region_set["motif_family_label"],
                    "n_regions": n_regions,
                    "n_cells": n_cells,
                    "total_peak_counts": class_total_peak,
                    "set_peak_counts": set_count,
                    "per_region_peak_cpm": per_region_peak_cpm,
                    "log2_per_region_peak_cpm": math.log2(per_region_peak_cpm + 1.0) if math.isfinite(per_region_peak_cpm) else math.nan,
                }
            )
    return rows


def combine_donor_pseudobulk(library_rows: list[dict[str, object]]) -> pd.DataFrame:
    frame = pd.DataFrame(library_rows)
    group_cols = [
        "biological_sample_id",
        "sample_name",
        "group",
        "genotype",
        "dev_stage_pcw",
        "sex",
        "inferred_cell_class",
        "region_set_id",
        "region_set_type",
        "motif_family",
        "motif_family_label",
        "n_regions",
    ]
    numeric = ["n_cells", "total_peak_counts", "set_peak_counts"]
    donor = frame.groupby(group_cols, dropna=False, observed=True)[numeric].sum().reset_index()
    donor["per_region_peak_cpm"] = np.where(
        donor["total_peak_counts"] > 0,
        donor["set_peak_counts"] / donor["n_regions"] / donor["total_peak_counts"] * 1_000_000,
        np.nan,
    )
    donor["log2_per_region_peak_cpm"] = np.log2(donor["per_region_peak_cpm"] + 1.0)
    return donor


def matched_ratio_table(donor: pd.DataFrame) -> pd.DataFrame:
    fg = donor[donor["region_set_type"] == SET_FOREGROUND].copy()
    bg = donor[donor["region_set_type"] == SET_MATCHED].copy()
    keys = ["biological_sample_id", "sample_name", "group", "genotype", "dev_stage_pcw", "sex", "inferred_cell_class", "motif_family"]
    merged = fg.merge(
        bg,
        on=keys,
        suffixes=("_foreground", "_matched"),
        how="inner",
    )
    merged["motif_family_label"] = merged["motif_family_label_foreground"]
    merged["matched_log2_ratio"] = merged["log2_per_region_peak_cpm_foreground"] - merged["log2_per_region_peak_cpm_matched"]
    merged["matched_raw_ratio"] = (
        (merged["per_region_peak_cpm_foreground"] + 1e-9) / (merged["per_region_peak_cpm_matched"] + 1e-9)
    )
    return merged


def effect_rows(ratio: pd.DataFrame, min_cells: int) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    filtered = ratio[
        (ratio["n_cells_foreground"] >= min_cells)
        & (ratio["n_cells_matched"] >= min_cells)
    ].copy()
    for (cell_class, family), sub in filtered.groupby(["inferred_cell_class", "motif_family"], dropna=False, observed=True):
        ds = sub[sub["group"] == "DS"]
        control = sub[sub["group"] == "Control"]
        if ds.empty or control.empty:
            continue
        row = {
            "inferred_cell_class": cell_class,
            "motif_family": family,
            "motif_family_label": sub["motif_family_label"].iloc[0],
            "n_foreground_regions": int(sub["n_regions_foreground"].iloc[0]),
            "n_matched_regions": int(sub["n_regions_matched"].iloc[0]),
            "n_ds_donors": int(ds["biological_sample_id"].nunique()),
            "n_control_donors": int(control["biological_sample_id"].nunique()),
            "n_ds_cells": int(ds["n_cells_foreground"].sum()),
            "n_control_cells": int(control["n_cells_foreground"].sum()),
        }
        for metric in [
            "matched_log2_ratio",
            "log2_per_region_peak_cpm_foreground",
            "log2_per_region_peak_cpm_matched",
            "per_region_peak_cpm_foreground",
            "per_region_peak_cpm_matched",
        ]:
            row[f"mean_ds_{metric}"] = float(ds[metric].mean())
            row[f"mean_control_{metric}"] = float(control[metric].mean())
            row[f"delta_ds_minus_control_{metric}"] = float(ds[metric].mean() - control[metric].mean())
        rows.append(row)
    return pd.DataFrame(rows)


def family_summary(effects: pd.DataFrame, qc: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    subset = effects[effects["inferred_cell_class"].isin(CLASS_ORDER)].copy()
    for family, sub in subset.groupby("motif_family", dropna=False, observed=True):
        values = pd.to_numeric(sub["delta_ds_minus_control_matched_log2_ratio"], errors="coerce")
        row = {
            "motif_family": family,
            "motif_family_label": sub["motif_family_label"].iloc[0],
            "n_foreground_regions": int(sub["n_foreground_regions"].iloc[0]),
            "n_matched_regions": int(sub["n_matched_regions"].iloc[0]),
            "n_inferred_classes": int(values.notna().sum()),
            "mean_delta_matched_log2_ratio": float(values.mean()),
            "median_delta_matched_log2_ratio": float(values.median()),
            "n_classes_positive": int((values > 0).sum()),
            "n_classes_negative": int((values < 0).sum()),
            "max_positive_class": str(sub.loc[values.idxmax(), "inferred_cell_class"]) if len(values) else "",
            "max_delta_matched_log2_ratio": float(values.max()) if len(values) else math.nan,
            "max_negative_class": str(sub.loc[values.idxmin(), "inferred_cell_class"]) if len(values) else "",
            "min_delta_matched_log2_ratio": float(values.min()) if len(values) else math.nan,
        }
        qc_sub = qc[qc["motif_family"] == family]
        if not qc_sub.empty:
            for column in qc_sub.columns:
                if column.startswith("mean_diff_") or column in {"n_replacement_rows"}:
                    row[column] = qc_sub[column].iloc[0]
        rows.append(row)
    summary = pd.DataFrame(rows)
    if not summary.empty:
        order = {family: index for index, family in enumerate(FAMILY_ORDER)}
        summary["_order"] = summary["motif_family"].map(order).fillna(999).astype(int)
        summary = summary.sort_values(["_order", "motif_family"]).drop(columns=["_order"])
    return summary


def plot_heatmap(effects: pd.DataFrame, figdir: Path) -> None:
    plot_df = effects[effects["inferred_cell_class"].isin(CLASS_ORDER)].copy()
    row_order = [family for family in FAMILY_ORDER if family in set(plot_df["motif_family"])]
    matrix = (
        plot_df.pivot_table(
            index="motif_family",
            columns="inferred_cell_class",
            values="delta_ds_minus_control_matched_log2_ratio",
            aggfunc="mean",
        )
        .reindex(index=row_order, columns=CLASS_ORDER)
    )
    labels = [MOTIF_FAMILY_LABELS.get(family, family) for family in row_order]
    values = matrix.to_numpy(dtype=float)
    finite = values[np.isfinite(values)]
    vmax = max(0.025, min(0.15, float(np.nanpercentile(np.abs(finite), 95)) if finite.size else 0.1))
    fig_width = max(11.0, 0.7 * len(CLASS_ORDER) + 4.5)
    fig_height = max(4.8, 0.42 * max(1, len(row_order)) + 2.2)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    image = ax.imshow(matrix, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(CLASS_ORDER)))
    ax.set_xticklabels([CLASS_LABELS.get(item, item) for item in CLASS_ORDER], rotation=45, ha="right")
    ax.set_yticks(range(len(row_order)))
    ax.set_yticklabels(labels)
    ax.set_title("GSE305146 matched-background motif accessibility")
    cbar = fig.colorbar(image, ax=ax, fraction=0.028, pad=0.02)
    cbar.set_label("DS - control delta motif/matched log2 ratio")
    fig.tight_layout()
    figdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(figdir / "gse305146_marker_atac_matched_motif_deviation_heatmap.png", dpi=300)
    fig.savefig(figdir / "gse305146_marker_atac_matched_motif_deviation_heatmap.pdf")
    plt.close(fig)


def write_note(path: Path, n_libraries: int, summary: pd.DataFrame, effects: pd.DataFrame) -> None:
    top_pos = effects.sort_values("delta_ds_minus_control_matched_log2_ratio", ascending=False).head(12)
    top_neg = effects.sort_values("delta_ds_minus_control_matched_log2_ratio", ascending=True).head(12)
    lines = [
        "# GSE305146 Matched-Background Motif Accessibility Validation",
        "",
        "Question: do motif-deviation signals survive recurrence/accessibility-like matched background correction?",
        "",
        "## Inputs",
        "",
        f"- Final-analysis libraries processed: {n_libraries}.",
        "- Foregrounds: motif-family peak sets from the recurrent-background HOMER motif scan.",
        "- Matched backgrounds: recurrent peaks without that motif family, matched greedily by max recurrence frequency, width, DS-control recurrence delta, and source-peak count bins.",
        "- Primary readout: DS-control delta of donor/class foreground-minus-matched log2 per-region peak CPM.",
        "- Interpretation scope: matched-background validation proxy; GC sequence matching is still not included.",
        "",
        "## Motif-Family Summary",
        "",
    ]
    for _, row in summary.iterrows():
        lines.append(
            f"- {row['motif_family_label']}: mean delta matched log2 ratio "
            f"{float(row['mean_delta_matched_log2_ratio']):.3g}; "
            f"{int(row['n_classes_positive'])}/{int(row['n_inferred_classes'])} classes positive; "
            f"foreground {int(row['n_foreground_regions'])}, matched {int(row['n_matched_regions'])}."
        )
    lines.extend(["", "## Strongest Positive Matched Effects", ""])
    for _, row in top_pos.iterrows():
        lines.append(
            f"- {row['motif_family_label']} in {row['inferred_cell_class']}: "
            f"delta matched log2 ratio {float(row['delta_ds_minus_control_matched_log2_ratio']):.3g}; "
            f"DS donors {int(row['n_ds_donors'])}, control donors {int(row['n_control_donors'])}."
        )
    lines.extend(["", "## Strongest Negative Matched Effects", ""])
    for _, row in top_neg.iterrows():
        lines.append(
            f"- {row['motif_family_label']} in {row['inferred_cell_class']}: "
            f"delta matched log2 ratio {float(row['delta_ds_minus_control_matched_log2_ratio']):.3g}; "
            f"DS donors {int(row['n_ds_donors'])}, control donors {int(row['n_control_donors'])}."
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This validation asks whether motif-bearing peak accessibility changes more than matched recurrent peaks with similar recurrence and peak-width structure. It is a stronger control than the full-background motif-deviation proxy, but still not full chromVAR because GC/content matching is not included.",
            "",
            "Signals that survive here are better candidates for motif-linked regulatory shifts. Signals that shrink toward zero should be interpreted as broad peak-set or accessibility-context effects rather than motif-specific activity. In the current full run, most developmental homeobox/bHLH families shrink toward small, class-specific effects, while MEF2 remains the clearest negative family and RFX/X-box is modestly negative. STAT control becomes small and balanced across classes.",
            "",
            "Thus, the strongest chromatin-motif conclusion is not uniform STAT or homeobox activation. It is broader fetal DS regulatory-context remodeling with a reproducible MEF2/RFX-linked accessibility-loss component and smaller, class-specific developmental homeobox/bHLH effects.",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def run_plot_only(args: argparse.Namespace) -> int:
    outdir = Path(args.outdir)
    effects = pd.read_csv(outdir / "gse305146_marker_atac_matched_motif_deviation_effects.tsv", sep="\t")
    summary = pd.read_csv(outdir / "gse305146_marker_atac_matched_motif_deviation_summary.tsv", sep="\t")
    library_pb = pd.read_csv(outdir / "gse305146_marker_atac_matched_motif_deviation_library_pseudobulk.tsv", sep="\t")
    n_libraries = int(library_pb["library_prefix"].nunique())
    plot_heatmap(effects, Path(args.figdir))
    write_note(Path(args.note), n_libraries, summary, effects)
    print("Regenerated matched motif-deviation figures")
    return 0


def main() -> int:
    args = parse_args()
    if args.plot_only:
        return run_plot_only(args)

    outdir = Path(args.outdir)
    figdir = Path(args.figdir)
    supplementary = Path(args.supplementary)
    families = [family for family in FAMILY_ORDER if family != "BACKGROUND"]
    if args.limit_families:
        families = families[: args.limit_families]

    background = add_match_bins(read_background(Path(args.background_table), args.region_size))
    foregrounds = read_motif_foregrounds(Path(args.motif_hits), families)
    region_sets, match_map, match_qc = build_matched_sets(background, foregrounds, families, args.seed)
    set_summary = region_set_summary(region_sets)
    set_ids = set_summary["region_set_id"].astype(str).tolist()
    region_set_index = {set_id: index for index, set_id in enumerate(set_ids)}
    bins = build_region_bins(region_sets, region_set_index, args.bin_size)

    region_fields = [
        "peak_id",
        "chrom",
        "match_start",
        "match_end",
        "region_set_id",
        "region_set_type",
        "motif_family",
        "motif_family_label",
        "width",
        "n_source_peaks",
        "ds_frequency",
        "control_frequency",
        "delta_frequency_ds_minus_control",
        "max_frequency",
    ]
    write_tsv(
        outdir / "gse305146_marker_atac_matched_motif_deviation_region_sets.tsv",
        region_sets.to_dict(orient="records"),
        region_fields,
    )
    match_map.to_csv(
        outdir / "gse305146_marker_atac_matched_motif_deviation_match_map.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    match_qc.to_csv(
        outdir / "gse305146_marker_atac_matched_motif_deviation_match_qc.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    set_summary.to_csv(
        outdir / "gse305146_marker_atac_matched_motif_deviation_region_set_summary.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )

    libraries = selected_libraries(Path(args.metadata))
    if args.limit_libraries:
        libraries = libraries[: args.limit_libraries]

    library_rows: list[dict[str, object]] = []
    mapping_rows: list[dict[str, object]] = []
    for index, library in enumerate(libraries, start=1):
        prefix = library["library_prefix"]
        feature_path = supplementary / f"GSE305146_{prefix}_features.tsv.gz"
        matrix_path = supplementary / f"GSE305146_{prefix}_matrix.mtx.gz"
        if not feature_path.exists() or not matrix_path.exists():
            raise FileNotFoundError(f"Missing feature or matrix for {prefix}")
        print(f"[{index}/{len(libraries)}] {prefix}", flush=True)
        row_to_marker, row_to_target, last_gene_row, _present_markers, marker_denominators = parse_marker_features(
            feature_path, []
        )
        row_to_sets, peak_last_gene_row, n_peak_rows, mapped_by_set = map_peak_rows_to_sets(
            feature_path, bins, args.bin_size
        )
        last_gene_row = max(last_gene_row, peak_last_gene_row)
        for set_id, set_index in region_set_index.items():
            mapping_rows.append(
                {
                    **library,
                    "region_set_id": set_id,
                    "n_feature_peak_rows": n_peak_rows,
                    "n_mapped_feature_peak_rows": mapped_by_set.get(str(set_index), 0),
                }
            )
        stats = stream_matrix(
            matrix_path,
            row_to_marker,
            row_to_target,
            last_gene_row,
            n_marker_sets=len(MARKER_SETS),
            n_target_genes=0,
        )
        classes, _confidence, _top_score, _delta = infer_classes(
            stats["total_rna"],
            stats["detected_genes"],
            stats["marker_counts"],
            marker_denominators,
            args,
        )
        total_peak_counts, set_peak_counts = stream_atac_sets(matrix_path, last_gene_row, row_to_sets, len(set_ids))
        library_rows.extend(summarize_library(library, classes, total_peak_counts, set_peak_counts, set_summary))

    library_fields = [
        "library_prefix",
        "geo_sample_accession",
        "biological_sample_id",
        "sample_name",
        "group",
        "genotype",
        "dev_stage_pcw",
        "sex",
        "platform_refs",
        "inferred_cell_class",
        "region_set_id",
        "region_set_type",
        "motif_family",
        "motif_family_label",
        "n_regions",
        "n_cells",
        "total_peak_counts",
        "set_peak_counts",
        "per_region_peak_cpm",
        "log2_per_region_peak_cpm",
    ]
    write_tsv(outdir / "gse305146_marker_atac_matched_motif_deviation_library_pseudobulk.tsv", library_rows, library_fields)
    pd.DataFrame(mapping_rows).to_csv(
        outdir / "gse305146_marker_atac_matched_motif_deviation_mapping_qc.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    donor = combine_donor_pseudobulk(library_rows)
    donor.to_csv(
        outdir / "gse305146_marker_atac_matched_motif_deviation_donor.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    ratio = matched_ratio_table(donor)
    ratio.to_csv(
        outdir / "gse305146_marker_atac_matched_motif_deviation_donor_ratio.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    effects = effect_rows(ratio, args.min_cells_per_donor_class)
    effects.to_csv(
        outdir / "gse305146_marker_atac_matched_motif_deviation_effects.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    summary = family_summary(effects, match_qc)
    summary.to_csv(
        outdir / "gse305146_marker_atac_matched_motif_deviation_summary.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    plot_heatmap(effects, figdir)
    write_note(Path(args.note), len(libraries), summary, effects)
    print(f"Wrote matched motif-deviation summaries for {len(libraries)} libraries")
    print(outdir / "gse305146_marker_atac_matched_motif_deviation_effects.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
