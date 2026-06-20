#!/usr/bin/env python3
"""ChromVAR-style motif deviation in marker-inferred GSE305146 classes."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import subprocess
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analyze_gse305146_marker_atac_motif_instance_peaksets import (
    ANNOTATE_PEAKS,
    FAMILY_ORDER,
    HOMER_BIN,
    MOTIF_FAMILY_LABELS,
    select_motif_records,
)
from infer_gse305146_marker_cell_classes_and_targets import (
    MARKER_SETS,
    infer_classes,
    parse_features as parse_marker_features,
    selected_libraries,
    stream_matrix,
)


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
UNIVERSE_SET_ID = "BACKGROUND_UNIVERSE"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", default="data/geo/GSE305146/metadata/sample_metadata.tsv")
    parser.add_argument("--supplementary", default="data/geo/GSE305146/supplementary")
    parser.add_argument(
        "--background-table",
        default="results/gse305146_scATAC_homer_motif_enrichment/tables/GSE305146_recurrent_peak_background.tsv",
    )
    parser.add_argument(
        "--homer-known",
        default="results/gse305146_scATAC_homer_motif_enrichment/homer_known_motif_summary.tsv",
    )
    parser.add_argument(
        "--homer-out",
        default="results/gse305146_scATAC_homer_motif_enrichment/homer_out",
    )
    parser.add_argument("--outdir", default="results/mechanism_question")
    parser.add_argument("--figdir", default="figures/mechanism_question")
    parser.add_argument("--note", default="notes/gse305146_marker_atac_motif_deviation_first_pass.md")
    parser.add_argument("--q-threshold", type=float, default=0.05)
    parser.add_argument("--region-size", type=int, default=200)
    parser.add_argument("--bin-size", type=int, default=100_000)
    parser.add_argument("--min-total-rna", type=float, default=500.0)
    parser.add_argument("--min-detected-genes", type=int, default=200)
    parser.add_argument("--min-top-score", type=float, default=0.05)
    parser.add_argument("--min-confidence-delta", type=float, default=0.05)
    parser.add_argument("--min-cells-per-donor-class", type=int, default=20)
    parser.add_argument("--limit-libraries", type=int, help="Debug: process first N final-analysis libraries.")
    parser.add_argument("--max-background-regions", type=int, help="Debug: scan only first N recurrent background regions.")
    parser.add_argument("--skip-homer", action="store_true", help="Reuse existing annotatePeaks output.")
    parser.add_argument("--plot-only", action="store_true", help="Regenerate figures from existing output tables.")
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


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", errors="replace") as handle:
        return [{k: (v or "").strip() for k, v in row.items()} for row in csv.DictReader(handle, delimiter="\t")]


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})


def safe_int(value: object) -> int:
    try:
        return int(float(str(value)))
    except ValueError:
        return 0


def read_background_regions(path: Path, region_size: int, max_regions: int | None = None) -> list[dict[str, object]]:
    half = region_size // 2
    regions: list[dict[str, object]] = []
    for rank, row in enumerate(read_tsv(path), start=1):
        chrom = row.get("chrom", "")
        center = safe_int(row.get("center", ""))
        if not chrom.startswith("chr") or center <= 0:
            continue
        start = max(0, center - half)
        end = start + region_size
        peak_id = f"background:{chrom}:{start}-{end}"
        regions.append(
            {
                "peak_id": peak_id,
                "chrom": chrom,
                "start": start,
                "end": end,
                "rank": rank,
                "source_width": row.get("width", ""),
                "source_start": row.get("start", ""),
                "source_end": row.get("end", ""),
                "ds_frequency": row.get("ds_frequency", ""),
                "control_frequency": row.get("control_frequency", ""),
                "delta_frequency_ds_minus_control": row.get("delta_frequency_ds_minus_control", ""),
                "max_frequency": row.get("max_frequency", ""),
            }
        )
        if max_regions and len(regions) >= max_regions:
            break
    return regions


def write_bed(path: Path, regions: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for region in regions:
            handle.write(
                "\t".join(
                    [
                        str(region["chrom"]),
                        str(int(region["start"])),
                        str(int(region["end"])),
                        str(region["peak_id"]),
                        str(int(region["rank"])),
                        ".",
                    ]
                )
                + "\n"
            )


def run_homer_scan(
    bed_path: Path,
    motif_records: list[dict[str, object]],
    output_path: Path,
    log_path: Path,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["PATH"] = f"{HOMER_BIN}:{env.get('PATH', '')}"
    cmd = [
        str(ANNOTATE_PEAKS),
        str(bed_path),
        "hg38",
        "-size",
        "given",
        "-m",
        *[str(record["motif_path"]) for record in motif_records],
        "-nmotifs",
    ]
    with output_path.open("w") as stdout, log_path.open("w") as stderr:
        subprocess.run(cmd, check=True, cwd=Path.cwd(), env=env, stdout=stdout, stderr=stderr)


def parse_homer_scan(
    output_path: Path,
    regions_by_id: dict[str, dict[str, object]],
    motif_records: list[dict[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    motif_by_name = {str(record["motif_name"]): record for record in motif_records}
    table = pd.read_csv(output_path, sep="\t")
    motif_columns: list[tuple[str, str]] = []
    for column in table.columns:
        if "/Homer" not in column:
            continue
        motif_name = column.split(" Distance From Peak", 1)[0]
        motif_columns.append((column, motif_name))

    family_regions: dict[tuple[str, str], dict[str, object]] = {}
    hit_rows: list[dict[str, object]] = []
    peak_column = table.columns[0]
    for row in table.to_dict(orient="records"):
        peak_id = str(row[peak_column])
        region = regions_by_id.get(peak_id)
        if region is None:
            continue
        for column, motif_name in motif_columns:
            count = safe_int(row.get(column, ""))
            if count <= 0:
                continue
            motif = motif_by_name.get(motif_name)
            if motif is None:
                continue
            family = str(motif["motif_family"])
            region_set_id = f"MOTIF_FAMILY|{family}"
            family_regions[(peak_id, family)] = {
                **region,
                "region_set_id": region_set_id,
                "motif_family": family,
                "motif_family_label": MOTIF_FAMILY_LABELS.get(family, family),
                "context_enriched": motif.get("context_enriched", ""),
            }
            hit_rows.append(
                {
                    **region,
                    "region_set_id": region_set_id,
                    "motif_family": family,
                    "motif_family_label": MOTIF_FAMILY_LABELS.get(family, family),
                    "motif_name": motif_name,
                    "motif_count": count,
                    "q_value": motif.get("q_value", ""),
                    "rank": motif.get("rank", ""),
                    "direction": motif.get("direction", ""),
                    "context_enriched": motif.get("context_enriched", ""),
                }
            )
    return list(family_regions.values()), hit_rows


def build_region_set_summary(
    background_regions: list[dict[str, object]],
    family_regions: list[dict[str, object]],
    hit_rows: list[dict[str, object]],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = [
        {
            "region_set_id": UNIVERSE_SET_ID,
            "motif_family": "BACKGROUND",
            "motif_family_label": "Recurrent background",
            "context_enriched": "background",
            "n_regions_in_set": len(background_regions),
            "n_motif_names": 0,
        }
    ]
    family_df = pd.DataFrame(family_regions)
    hit_df = pd.DataFrame(hit_rows)
    if not family_df.empty:
        summary = (
            family_df.groupby(["region_set_id", "motif_family", "motif_family_label"], dropna=False, observed=True)
            .agg(n_regions_in_set=("peak_id", "nunique"))
            .reset_index()
        )
        if not hit_df.empty:
            motif_counts = (
                hit_df.groupby("region_set_id", dropna=False, observed=True)
                .agg(n_motif_names=("motif_name", "nunique"))
                .reset_index()
            )
            enriched = (
                hit_df.groupby("region_set_id", dropna=False, observed=True)
                .agg(context_enriched=("context_enriched", lambda values: "yes" if "yes" in set(values) else "no"))
                .reset_index()
            )
            summary = summary.merge(motif_counts, on="region_set_id", how="left").merge(enriched, on="region_set_id", how="left")
        summary["n_motif_names"] = summary["n_motif_names"].fillna(0).astype(int)
        rows.extend(summary.to_dict(orient="records"))
    out = pd.DataFrame(rows)
    order = {family: index for index, family in enumerate(["BACKGROUND"] + FAMILY_ORDER)}
    out["_order"] = out["motif_family"].map(order).fillna(999).astype(int)
    return out.sort_values(["_order", "region_set_id"]).drop(columns=["_order"])


def build_region_bins(
    background_regions: list[dict[str, object]],
    family_regions: list[dict[str, object]],
    region_set_index: dict[str, int],
    bin_size: int,
) -> dict[str, dict[int, list[tuple[int, int, int]]]]:
    bins: dict[str, dict[int, list[tuple[int, int, int]]]] = defaultdict(lambda: defaultdict(list))
    for region in background_regions:
        chrom = str(region["chrom"])
        start = int(region["start"])
        end = int(region["end"])
        set_index = region_set_index[UNIVERSE_SET_ID]
        for bin_id in range(start // bin_size, end // bin_size + 1):
            bins[chrom][bin_id].append((start, end, set_index))
    for region in family_regions:
        chrom = str(region["chrom"])
        start = int(region["start"])
        end = int(region["end"])
        set_index = region_set_index[str(region["region_set_id"])]
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
    all_peak_counts: np.ndarray | None = None
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
                all_peak_counts = np.zeros(n_cells, dtype=np.float64)
                set_peak_counts = np.zeros((n_cells, n_sets), dtype=np.float32)
                dims_seen = True
                continue
            assert all_peak_counts is not None
            assert set_peak_counts is not None
            row_index = int(fields[0])
            if row_index <= last_gene_expression_row:
                continue
            cell_index = int(fields[1]) - 1
            value = float(fields[2])
            all_peak_counts[cell_index] += value
            for set_index in row_to_sets.get(row_index, []):
                set_peak_counts[cell_index, set_index] += value
    if all_peak_counts is None or set_peak_counts is None:
        raise RuntimeError(f"No Matrix Market dimensions found in {matrix_path}")
    return all_peak_counts, set_peak_counts


def summarize_library(
    library: dict[str, str],
    classes: np.ndarray,
    all_peak_counts: np.ndarray,
    set_peak_counts: np.ndarray,
    region_sets: pd.DataFrame,
    universe_index: int,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    class_values = classes.astype(str)
    for cell_class in sorted(set(class_values)):
        mask = class_values == cell_class
        n_cells = int(mask.sum())
        if n_cells == 0:
            continue
        class_all_peak = float(all_peak_counts[mask].sum())
        universe_count = float(set_peak_counts[mask, universe_index].sum())
        for set_index, region_set in enumerate(region_sets.to_dict(orient="records")):
            set_count = float(set_peak_counts[mask, set_index].sum())
            rows.append(
                {
                    **library,
                    "inferred_cell_class": cell_class,
                    "region_set_id": region_set["region_set_id"],
                    "motif_family": region_set["motif_family"],
                    "motif_family_label": region_set["motif_family_label"],
                    "context_enriched": region_set["context_enriched"],
                    "n_regions_in_set": region_set["n_regions_in_set"],
                    "n_cells": n_cells,
                    "all_peak_counts": class_all_peak,
                    "background_peak_counts": universe_count,
                    "motif_peak_counts": set_count,
                    "fraction_of_background_counts": set_count / universe_count if universe_count > 0 else math.nan,
                    "background_fraction_of_all_peaks": universe_count / class_all_peak if class_all_peak > 0 else math.nan,
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
        "motif_family",
        "motif_family_label",
        "context_enriched",
        "n_regions_in_set",
    ]
    numeric = ["n_cells", "all_peak_counts", "background_peak_counts", "motif_peak_counts"]
    donor = frame.groupby(group_cols, dropna=False, observed=True)[numeric].sum().reset_index()
    donor["fraction_of_background_counts"] = np.where(
        donor["background_peak_counts"] > 0,
        donor["motif_peak_counts"] / donor["background_peak_counts"],
        np.nan,
    )
    donor["background_fraction_of_all_peaks"] = np.where(
        donor["all_peak_counts"] > 0,
        donor["background_peak_counts"] / donor["all_peak_counts"],
        np.nan,
    )
    return donor


def add_deviation_metrics(donor: pd.DataFrame) -> pd.DataFrame:
    donor = donor.copy()
    motif_rows = donor[donor["region_set_id"] != UNIVERSE_SET_ID]
    expected_fraction = (
        motif_rows.groupby("region_set_id", dropna=False, observed=True)
        .agg(
            total_motif_peak_counts=("motif_peak_counts", "sum"),
            total_background_peak_counts=("background_peak_counts", "sum"),
        )
        .reset_index()
    )
    expected_fraction["expected_background_fraction"] = (
        expected_fraction["total_motif_peak_counts"] / expected_fraction["total_background_peak_counts"]
    )
    expected_fraction = expected_fraction[["region_set_id", "expected_background_fraction"]]
    donor = donor.merge(expected_fraction, on="region_set_id", how="left")
    donor.loc[donor["region_set_id"] == UNIVERSE_SET_ID, "expected_background_fraction"] = 1.0
    donor["expected_motif_peak_counts"] = donor["background_peak_counts"] * donor["expected_background_fraction"]
    donor["deviation_raw"] = np.where(
        donor["expected_motif_peak_counts"] > 0,
        (donor["motif_peak_counts"] - donor["expected_motif_peak_counts"]) / donor["expected_motif_peak_counts"],
        np.nan,
    )
    donor["deviation_log2_observed_expected"] = np.log2(
        (donor["motif_peak_counts"].astype(float) + 0.5)
        / (donor["expected_motif_peak_counts"].astype(float) + 0.5)
    )
    variance = donor["expected_motif_peak_counts"] * (1.0 - donor["expected_background_fraction"].clip(upper=0.999999)) + 1.0
    donor["deviation_z"] = (donor["motif_peak_counts"] - donor["expected_motif_peak_counts"]) / np.sqrt(variance)
    return donor


def effect_rows(donor: pd.DataFrame, min_cells: int) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    filtered = donor[(donor["n_cells"] >= min_cells) & (donor["region_set_id"] != UNIVERSE_SET_ID)].copy()
    for (cell_class, set_id), sub in filtered.groupby(["inferred_cell_class", "region_set_id"], dropna=False, observed=True):
        ds = sub[sub["group"] == "DS"]
        control = sub[sub["group"] == "Control"]
        if ds.empty or control.empty:
            continue
        row = {
            "inferred_cell_class": cell_class,
            "region_set_id": set_id,
            "motif_family": sub["motif_family"].iloc[0],
            "motif_family_label": sub["motif_family_label"].iloc[0],
            "context_enriched": sub["context_enriched"].iloc[0],
            "n_regions_in_set": int(sub["n_regions_in_set"].iloc[0]),
            "expected_background_fraction": float(sub["expected_background_fraction"].iloc[0]),
            "n_ds_donors": int(ds["biological_sample_id"].nunique()),
            "n_control_donors": int(control["biological_sample_id"].nunique()),
            "n_ds_cells": int(ds["n_cells"].sum()),
            "n_control_cells": int(control["n_cells"].sum()),
        }
        for metric in [
            "deviation_z",
            "deviation_log2_observed_expected",
            "deviation_raw",
            "fraction_of_background_counts",
            "background_fraction_of_all_peaks",
        ]:
            row[f"mean_ds_{metric}"] = float(ds[metric].mean())
            row[f"mean_control_{metric}"] = float(control[metric].mean())
            row[f"delta_ds_minus_control_{metric}"] = float(ds[metric].mean() - control[metric].mean())
        rows.append(row)
    return pd.DataFrame(rows)


def family_summary(effects: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    subset = effects[effects["inferred_cell_class"].isin(CLASS_ORDER)].copy()
    for set_id, sub in subset.groupby("region_set_id", dropna=False, observed=True):
        z_values = pd.to_numeric(sub["delta_ds_minus_control_deviation_z"], errors="coerce")
        log2_values = pd.to_numeric(sub["delta_ds_minus_control_deviation_log2_observed_expected"], errors="coerce")
        rows.append(
            {
                "region_set_id": set_id,
                "motif_family": sub["motif_family"].iloc[0],
                "motif_family_label": sub["motif_family_label"].iloc[0],
                "context_enriched": sub["context_enriched"].iloc[0],
                "n_regions_in_set": int(sub["n_regions_in_set"].iloc[0]),
                "expected_background_fraction": float(sub["expected_background_fraction"].iloc[0]),
                "n_inferred_classes": int(log2_values.notna().sum()),
                "mean_delta_log2_observed_expected": float(log2_values.mean()),
                "median_delta_log2_observed_expected": float(log2_values.median()),
                "n_classes_log2_positive": int((log2_values > 0).sum()),
                "n_classes_log2_negative": int((log2_values < 0).sum()),
                "max_log2_positive_class": str(sub.loc[log2_values.idxmax(), "inferred_cell_class"]) if len(log2_values) else "",
                "max_delta_log2_observed_expected": float(log2_values.max()) if len(log2_values) else math.nan,
                "max_log2_negative_class": str(sub.loc[log2_values.idxmin(), "inferred_cell_class"]) if len(log2_values) else "",
                "min_delta_log2_observed_expected": float(log2_values.min()) if len(log2_values) else math.nan,
                "mean_delta_deviation_z": float(z_values.mean()),
                "median_delta_deviation_z": float(z_values.median()),
                "n_classes_z_positive": int((z_values > 0).sum()),
                "n_classes_z_negative": int((z_values < 0).sum()),
                "max_z_positive_class": str(sub.loc[z_values.idxmax(), "inferred_cell_class"]) if len(z_values) else "",
                "max_delta_deviation_z": float(z_values.max()) if len(z_values) else math.nan,
                "max_z_negative_class": str(sub.loc[z_values.idxmin(), "inferred_cell_class"]) if len(z_values) else "",
                "min_delta_deviation_z": float(z_values.min()) if len(z_values) else math.nan,
            }
        )
    summary = pd.DataFrame(rows)
    if not summary.empty:
        order = {f"MOTIF_FAMILY|{family}": index for index, family in enumerate(FAMILY_ORDER)}
        summary["_order"] = summary["region_set_id"].map(order).fillna(999).astype(int)
        summary = summary.sort_values(["_order", "region_set_id"]).drop(columns=["_order"])
    return summary


def make_one_heatmap(effects: pd.DataFrame, figdir: Path, value_column: str, suffix: str, label: str, title: str, vmax_cap: float) -> None:
    plot_df = effects[effects["inferred_cell_class"].isin(CLASS_ORDER)].copy()
    row_order = [f"MOTIF_FAMILY|{family}" for family in FAMILY_ORDER if f"MOTIF_FAMILY|{family}" in set(plot_df["region_set_id"])]
    matrix = (
        plot_df.pivot_table(
            index="region_set_id",
            columns="inferred_cell_class",
            values=value_column,
            aggfunc="mean",
        )
        .reindex(index=row_order, columns=CLASS_ORDER)
    )
    labels = (
        plot_df.drop_duplicates("region_set_id")
        .set_index("region_set_id")
        .reindex(row_order)["motif_family_label"]
        .fillna("")
        .tolist()
    )
    values = matrix.to_numpy(dtype=float)
    finite = values[np.isfinite(values)]
    vmax = max(0.025, min(vmax_cap, float(np.nanpercentile(np.abs(finite), 95)) if finite.size else 1.0))
    fig_width = max(11.0, 0.7 * len(CLASS_ORDER) + 4.5)
    fig_height = max(4.8, 0.42 * max(1, len(row_order)) + 2.2)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    image = ax.imshow(matrix, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(CLASS_ORDER)))
    ax.set_xticklabels([CLASS_LABELS.get(item, item) for item in CLASS_ORDER], rotation=45, ha="right")
    ax.set_yticks(range(len(row_order)))
    ax.set_yticklabels(labels)
    ax.set_title(title)
    cbar = fig.colorbar(image, ax=ax, fraction=0.028, pad=0.02)
    cbar.set_label(label)
    ax.grid(which="major", color="white", linewidth=0.8)
    fig.tight_layout()
    figdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(figdir / f"gse305146_marker_atac_motif_deviation_{suffix}_heatmap.png", dpi=300)
    fig.savefig(figdir / f"gse305146_marker_atac_motif_deviation_{suffix}_heatmap.pdf")
    plt.close(fig)


def make_heatmap(effects: pd.DataFrame, figdir: Path) -> None:
    make_one_heatmap(
        effects,
        figdir,
        "delta_ds_minus_control_deviation_log2_observed_expected",
        "log2oe",
        "DS - control delta log2 observed/expected",
        "GSE305146 motif-deviation activity proxy",
        vmax_cap=0.12,
    )
    make_one_heatmap(
        effects,
        figdir,
        "delta_ds_minus_control_deviation_z",
        "z",
        "DS - control delta deviation z",
        "GSE305146 motif-deviation z score proxy",
        vmax_cap=30.0,
    )


def write_note(path: Path, n_libraries: int, region_summary: pd.DataFrame, effects: pd.DataFrame, summary: pd.DataFrame) -> None:
    top_pos = effects.sort_values("delta_ds_minus_control_deviation_log2_observed_expected", ascending=False).head(12)
    top_neg = effects.sort_values("delta_ds_minus_control_deviation_log2_observed_expected", ascending=True).head(12)
    stat = summary[summary["motif_family"] == "STAT_CONTROL"]
    lines = [
        "# GSE305146 Marker-Inferred Motif-Deviation First Pass",
        "",
        "Question: do selected HOMER motif families show DS-control motif-activity shifts within marker-inferred fetal cortical classes?",
        "",
        "## Inputs",
        "",
        f"- Final-analysis libraries processed: {n_libraries}.",
        f"- Recurrent background peaks scanned: {int(region_summary.loc[region_summary['region_set_id'] == UNIVERSE_SET_ID, 'n_regions_in_set'].iloc[0])}.",
        "- Motif universe: selected significant GSE305146 HOMER homeobox/bHLH and RFX/MEF2 motifs plus STAT1/STAT3 non-enriched controls.",
        "- Scoring: motif-bearing recurrent peak counts were compared with expected counts from the full recurrent background peak universe within each donor/class.",
        "- Primary readout: donor/class delta log2 observed/expected motif accessibility. Deviation z is retained as a secondary depth-sensitive diagnostic.",
        "- Interpretation scope: chromVAR-style motif-deviation proxy; not GC-matched chromVAR, not TF binding, and not DNA methylation.",
        "",
        "## Motif Sets",
        "",
    ]
    for _, row in summary.iterrows():
        lines.append(
            f"- {row['motif_family_label']}: {int(row['n_regions_in_set'])} motif-bearing peaks; "
            f"expected background fraction {float(row['expected_background_fraction']):.4g}; "
            f"mean delta log2 observed/expected {float(row['mean_delta_log2_observed_expected']):.3g}; "
            f"{int(row['n_classes_log2_positive'])}/{int(row['n_inferred_classes'])} classes positive; "
            f"mean delta z {float(row['mean_delta_deviation_z']):.3g}; "
            f"context enriched={row['context_enriched']}."
        )
    lines.extend(["", "## Strongest Positive Motif-Deviation Effects", ""])
    for _, row in top_pos.iterrows():
        lines.append(
            f"- {row['motif_family_label']} in {row['inferred_cell_class']}: "
            f"delta log2 observed/expected {float(row['delta_ds_minus_control_deviation_log2_observed_expected']):.3g}; "
            f"delta z {float(row['delta_ds_minus_control_deviation_z']):.3g}; "
            f"DS donors {int(row['n_ds_donors'])}, control donors {int(row['n_control_donors'])}."
        )
    lines.extend(["", "## Strongest Negative Motif-Deviation Effects", ""])
    for _, row in top_neg.iterrows():
        lines.append(
            f"- {row['motif_family_label']} in {row['inferred_cell_class']}: "
            f"delta log2 observed/expected {float(row['delta_ds_minus_control_deviation_log2_observed_expected']):.3g}; "
            f"delta z {float(row['delta_ds_minus_control_deviation_z']):.3g}; "
            f"DS donors {int(row['n_ds_donors'])}, control donors {int(row['n_control_donors'])}."
        )
    if not stat.empty:
        lines.extend(["", "## STAT Control Motifs", ""])
        for _, row in stat.iterrows():
            lines.append(
                f"- {row['motif_family_label']}: mean delta log2 observed/expected {float(row['mean_delta_log2_observed_expected']):.3g}; "
                f"{int(row['n_classes_log2_positive'])}/{int(row['n_inferred_classes'])} classes positive; "
                f"mean delta z {float(row['mean_delta_deviation_z']):.3g}; "
                f"context enriched={row['context_enriched']}."
            )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This pass is closer to motif activity than the earlier HOMER enrichment and motif-instance peak-set projections because it compares motif-bearing peak counts with the full recurrent accessibility background in each donor/class. It still should be described as chromVAR-style rather than full chromVAR because it does not use GC/accessibility-matched background peak sampling.",
            "",
            "The log2 observed/expected readout should be emphasized over pseudobulk z because z is sensitive to donor/class depth and cell number. Large z differences can correspond to small motif-enrichment shifts.",
            "",
            "The strongest robust direction in this first pass is reduced MEF2 motif-deviation in DS across most inferred classes. Developmental homeobox/bHLH motifs show smaller, class-specific positive shifts, especially in cycling/progenitor-like and selected neuron-like contexts, but they are not uniformly positive across all classes. STAT control motifs show a small broad positive deviation, but they remain non-enriched controls and should not be interpreted as a STAT-first direct chromatin mechanism without stronger chromVAR, footprinting, or perturbation support.",
            "",
            "A convincing direct TF-activity claim should require consistency across this deviation proxy, true chromVAR or differential-accessibility modeling, and ideally perturbation evidence. STAT remains a control unless its deviation is both strong and biologically coherent despite lack of HOMER enrichment.",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    figdir = Path(args.figdir)
    if args.plot_only:
        effects = pd.read_csv(outdir / "gse305146_marker_atac_motif_deviation_effects.tsv", sep="\t")
        region_summary = pd.read_csv(outdir / "gse305146_marker_atac_motif_deviation_region_set_summary.tsv", sep="\t")
        library_pb = pd.read_csv(outdir / "gse305146_marker_atac_motif_deviation_library_pseudobulk.tsv", sep="\t")
        n_libraries = int(library_pb["library_prefix"].nunique()) if "library_prefix" in library_pb.columns else 0
        summary = family_summary(effects)
        summary.to_csv(
            outdir / "gse305146_marker_atac_motif_deviation_summary.tsv",
            sep="\t",
            index=False,
            lineterminator="\n",
        )
        make_heatmap(effects, figdir)
        write_note(Path(args.note), n_libraries, region_summary, effects, summary)
        print("Regenerated motif-deviation figures")
        return 0

    scan_dir = outdir / "gse305146_marker_atac_motif_deviation"
    input_dir = scan_dir / "inputs"
    annotate_dir = scan_dir / "annotatePeaks"
    input_dir.mkdir(parents=True, exist_ok=True)
    annotate_dir.mkdir(parents=True, exist_ok=True)

    background_regions = read_background_regions(Path(args.background_table), args.region_size, args.max_background_regions)
    if not background_regions:
        raise SystemExit("No recurrent background regions found")
    motif_records = select_motif_records(Path(args.homer_known), Path(args.homer_out), args.q_threshold)
    write_tsv(
        outdir / "gse305146_marker_atac_motif_deviation_motif_records.tsv",
        motif_records,
        [
            "direction",
            "motif_path",
            "motif_name",
            "motif_family",
            "q_value",
            "rank",
            "target_percent",
            "background_percent",
            "context_enriched",
        ],
    )
    bed_path = input_dir / "gse305146_recurrent_background_hg38_200bp.bed"
    output_path = annotate_dir / "gse305146_recurrent_background_motif_deviation.annotatePeaks.nmotifs.tsv"
    log_path = annotate_dir / "gse305146_recurrent_background_motif_deviation.annotatePeaks.log"
    write_bed(bed_path, background_regions)
    if not args.skip_homer or not output_path.exists():
        run_homer_scan(bed_path, motif_records, output_path, log_path)
    family_regions, hit_rows = parse_homer_scan(output_path, {row["peak_id"]: row for row in background_regions}, motif_records)

    write_tsv(
        outdir / "gse305146_marker_atac_motif_deviation_scan_manifest.tsv",
        [
            {
                "background_table": args.background_table,
                "bed": bed_path,
                "annotatePeaks_output": output_path,
                "annotatePeaks_log": log_path,
                "n_background_regions": len(background_regions),
                "n_motifs_scanned": len(motif_records),
                "n_family_peak_sets": len(family_regions),
                "n_region_motif_hits": len(hit_rows),
            }
        ],
        [
            "background_table",
            "bed",
            "annotatePeaks_output",
            "annotatePeaks_log",
            "n_background_regions",
            "n_motifs_scanned",
            "n_family_peak_sets",
            "n_region_motif_hits",
        ],
    )
    write_tsv(
        outdir / "gse305146_marker_atac_motif_deviation_region_hits.tsv",
        hit_rows,
        [
            "peak_id",
            "chrom",
            "start",
            "end",
            "rank",
            "region_set_id",
            "motif_family",
            "motif_family_label",
            "motif_name",
            "motif_count",
            "q_value",
            "direction",
            "context_enriched",
        ],
    )
    region_summary = build_region_set_summary(background_regions, family_regions, hit_rows)
    region_summary.to_csv(
        outdir / "gse305146_marker_atac_motif_deviation_region_set_summary.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    region_sets = region_summary.to_dict(orient="records")
    region_set_index = {str(row["region_set_id"]): idx for idx, row in enumerate(region_sets)}
    universe_index = region_set_index[UNIVERSE_SET_ID]
    bins = build_region_bins(background_regions, family_regions, region_set_index, args.bin_size)

    libraries = selected_libraries(Path(args.metadata))
    if args.limit_libraries:
        libraries = libraries[: args.limit_libraries]

    supplementary = Path(args.supplementary)
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
        for region_set_id, set_index in region_set_index.items():
            mapping_rows.append(
                {
                    **library,
                    "region_set_id": region_set_id,
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
        all_peak_counts, set_peak_counts = stream_atac_sets(matrix_path, last_gene_row, row_to_sets, len(region_sets))
        library_rows.extend(
            summarize_library(library, classes, all_peak_counts, set_peak_counts, region_summary, universe_index)
        )

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
        "motif_family",
        "motif_family_label",
        "context_enriched",
        "n_regions_in_set",
        "n_cells",
        "all_peak_counts",
        "background_peak_counts",
        "motif_peak_counts",
        "fraction_of_background_counts",
        "background_fraction_of_all_peaks",
    ]
    write_tsv(outdir / "gse305146_marker_atac_motif_deviation_library_pseudobulk.tsv", library_rows, library_fields)
    pd.DataFrame(mapping_rows).to_csv(
        outdir / "gse305146_marker_atac_motif_deviation_mapping_qc.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    donor = add_deviation_metrics(combine_donor_pseudobulk(library_rows))
    donor.to_csv(
        outdir / "gse305146_marker_atac_motif_deviation_donor.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    effects = effect_rows(donor, args.min_cells_per_donor_class)
    effects.to_csv(
        outdir / "gse305146_marker_atac_motif_deviation_effects.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    summary = family_summary(effects)
    summary.to_csv(
        outdir / "gse305146_marker_atac_motif_deviation_summary.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    make_heatmap(effects, figdir)
    write_note(Path(args.note), len(libraries), region_summary, effects, summary)
    print(f"Wrote motif-deviation summaries for {len(libraries)} libraries")
    print(outdir / "gse305146_marker_atac_motif_deviation_effects.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
