#!/usr/bin/env python3
"""Donor-level differential accessibility for recurrent GSE305146 ATAC peaks."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, t as student_t

from analyze_gse305146_marker_atac_motif_deviation import CLASS_LABELS, CLASS_ORDER
from analyze_gse305146_marker_atac_motif_instance_peaksets import FAMILY_ORDER, MOTIF_FAMILY_LABELS
from infer_gse305146_marker_cell_classes_and_targets import (
    MARKER_SETS,
    infer_classes,
    parse_features as parse_marker_features,
    selected_libraries,
    stream_matrix,
)


AUTOSOMES = {f"chr{i}" for i in range(1, 23)}


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
    parser.add_argument("--note", default="notes/gse305146_marker_atac_peak_da_first_pass.md")
    parser.add_argument("--region-size", type=int, default=200)
    parser.add_argument("--bin-size", type=int, default=100_000)
    parser.add_argument("--min-total-rna", type=float, default=500.0)
    parser.add_argument("--min-detected-genes", type=int, default=200)
    parser.add_argument("--min-top-score", type=float, default=0.05)
    parser.add_argument("--min-confidence-delta", type=float, default=0.05)
    parser.add_argument("--min-cells-per-donor-class", type=int, default=20)
    parser.add_argument("--min-donors-per-group", type=int, default=3)
    parser.add_argument("--fdr-threshold", type=float, default=0.05)
    parser.add_argument("--min-abs-beta", type=float, default=0.25)
    parser.add_argument("--include-sex-chromosomes", action="store_true")
    parser.add_argument("--limit-libraries", type=int, help="Debug: process first N final-analysis libraries.")
    parser.add_argument("--max-background-regions", type=int, help="Debug: use first N retained recurrent peaks.")
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


def safe_int(value: object) -> int:
    try:
        return int(float(str(value)))
    except ValueError:
        return 0


def safe_float(value: object) -> float:
    try:
        return float(str(value))
    except ValueError:
        return math.nan


def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, dtype=np.float64)
    q = np.ones_like(p)
    finite = np.isfinite(p)
    if not finite.any():
        return q
    idx = np.where(finite)[0]
    order = idx[np.argsort(p[idx])]
    ranked = p[order]
    n = ranked.size
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    q[order] = np.clip(adjusted, 0.0, 1.0)
    return q


def read_background_regions(args: argparse.Namespace) -> pd.DataFrame:
    half = args.region_size // 2
    rows: list[dict[str, object]] = []
    for rank, row in enumerate(read_tsv(Path(args.background_table)), start=1):
        chrom = row.get("chrom", "")
        if not args.include_sex_chromosomes and chrom not in AUTOSOMES:
            continue
        center = safe_int(row.get("center", ""))
        if not chrom.startswith("chr") or center <= 0:
            continue
        start = max(0, center - half)
        end = start + args.region_size
        rows.append(
            {
                "region_index": len(rows),
                "peak_id": f"background:{chrom}:{start}-{end}",
                "chrom": chrom,
                "start": start,
                "end": end,
                "background_rank": rank,
                "source_width": safe_float(row.get("width", "")),
                "n_source_peaks": safe_float(row.get("n_source_peaks", "")),
                "ds_frequency": safe_float(row.get("ds_frequency", "")),
                "control_frequency": safe_float(row.get("control_frequency", "")),
                "delta_frequency_ds_minus_control": safe_float(row.get("delta_frequency_ds_minus_control", "")),
                "max_frequency": safe_float(row.get("max_frequency", "")),
            }
        )
        if args.max_background_regions and len(rows) >= args.max_background_regions:
            break
    if not rows:
        raise RuntimeError("No recurrent background peaks retained")
    return pd.DataFrame(rows)


def build_region_bins(regions: pd.DataFrame, bin_size: int) -> dict[str, dict[int, list[tuple[int, int, int]]]]:
    bins: dict[str, dict[int, list[tuple[int, int, int]]]] = defaultdict(lambda: defaultdict(list))
    for row in regions.itertuples(index=False):
        chrom = str(row.chrom)
        start = int(row.start)
        end = int(row.end)
        region_index = int(row.region_index)
        for bin_id in range(start // bin_size, end // bin_size + 1):
            bins[chrom][bin_id].append((start, end, region_index))
    return bins


def map_peak_rows_to_regions(
    feature_path: Path,
    bins: dict[str, dict[int, list[tuple[int, int, int]]]],
    bin_size: int,
) -> tuple[dict[int, list[int]], int, int, int, int]:
    row_to_regions: dict[int, list[int]] = {}
    last_gene_expression_row = 0
    n_peak_rows = 0
    mapped_region_ids: set[int] = set()
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
            matches: set[int] = set()
            for bin_id in range(start // bin_size, end // bin_size + 1):
                for region_start, region_end, region_index in bins[chrom].get(bin_id, []):
                    if end >= region_start and start <= region_end:
                        matches.add(region_index)
            if matches:
                ordered = sorted(matches)
                row_to_regions[row_index] = ordered
                mapped_region_ids.update(ordered)
    return row_to_regions, last_gene_expression_row, n_peak_rows, len(row_to_regions), len(mapped_region_ids)


def stream_atac_regions(
    matrix_path: Path,
    classes: np.ndarray,
    last_gene_expression_row: int,
    row_to_regions: dict[int, list[int]],
    n_regions: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int, int]:
    class_to_idx = {name: idx for idx, name in enumerate(CLASS_ORDER)}
    cell_class_idx = np.array([class_to_idx.get(str(cell_class), -1) for cell_class in classes], dtype=np.int16)
    n_classes = len(CLASS_ORDER)
    class_cells = np.array([(cell_class_idx == idx).sum() for idx in range(n_classes)], dtype=np.int32)
    class_peak_depth = np.zeros(n_classes, dtype=np.float64)
    class_region_counts = np.zeros((n_classes, n_regions), dtype=np.float32)
    dims_seen = False
    n_atac_entries = 0
    n_mapped_entries = 0

    with gzip.open(matrix_path, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            fields = line.strip().split()
            if not fields:
                continue
            if not dims_seen:
                dims_seen = True
                continue
            row_index = int(fields[0])
            if row_index <= last_gene_expression_row:
                continue
            cell_index = int(fields[1]) - 1
            class_idx = int(cell_class_idx[cell_index])
            if class_idx < 0:
                continue
            value = float(fields[2])
            n_atac_entries += 1
            class_peak_depth[class_idx] += value
            region_indices = row_to_regions.get(row_index)
            if not region_indices:
                continue
            n_mapped_entries += 1
            for region_index in region_indices:
                class_region_counts[class_idx, region_index] += value
    return class_region_counts, class_peak_depth, class_cells, n_atac_entries, n_mapped_entries


def collect_donor_pseudobulk(
    args: argparse.Namespace,
    regions: pd.DataFrame,
) -> tuple[dict[tuple[str, str], dict[str, object]], pd.DataFrame, pd.DataFrame]:
    supplementary = Path(args.supplementary)
    libraries = selected_libraries(Path(args.metadata))
    if args.limit_libraries:
        libraries = libraries[: args.limit_libraries]
    bins = build_region_bins(regions, args.bin_size)
    n_regions = len(regions)
    marker_names = list(MARKER_SETS)
    donor_store: dict[tuple[str, str], dict[str, object]] = {}
    mapping_rows: list[dict[str, object]] = []
    library_class_rows: list[dict[str, object]] = []

    for idx, library in enumerate(libraries, start=1):
        prefix = library["library_prefix"]
        print(f"[{idx}/{len(libraries)}] {prefix}", flush=True)
        feature_path = supplementary / f"GSE305146_{prefix}_features.tsv.gz"
        matrix_path = supplementary / f"GSE305146_{prefix}_matrix.mtx.gz"
        row_to_marker, row_to_target, last_gene_expression_row, _present, marker_denominators = parse_marker_features(
            feature_path,
            [],
        )
        marker_data = stream_matrix(
            matrix_path,
            row_to_marker,
            row_to_target,
            last_gene_expression_row,
            len(marker_names),
            0,
        )
        classes, _confidence, _top_score, _delta = infer_classes(
            marker_data["total_rna"],
            marker_data["detected_genes"],
            marker_data["marker_counts"],
            marker_denominators,
            args,
        )
        row_to_regions, atac_last_gene_row, n_peak_rows, n_mapped_peak_rows, n_mapped_regions = map_peak_rows_to_regions(
            feature_path,
            bins,
            args.bin_size,
        )
        region_counts, peak_depth, class_cells, n_atac_entries, n_mapped_entries = stream_atac_regions(
            matrix_path,
            classes,
            atac_last_gene_row,
            row_to_regions,
            n_regions,
        )
        mapping_rows.append(
            {
                **library,
                "n_matrix_cells": int(marker_data["n_cells"]),
                "n_peak_feature_rows": n_peak_rows,
                "n_mapped_peak_feature_rows": n_mapped_peak_rows,
                "n_mapped_recurrent_regions": n_mapped_regions,
                "n_atac_entries_in_retained_classes": n_atac_entries,
                "n_mapped_atac_entries_in_retained_classes": n_mapped_entries,
            }
        )
        for class_idx, cell_class in enumerate(CLASS_ORDER):
            n_cells = int(class_cells[class_idx])
            if n_cells <= 0:
                continue
            class_count_sum = float(region_counts[class_idx].sum())
            library_class_rows.append(
                {
                    **library,
                    "inferred_cell_class": cell_class,
                    "n_cells": n_cells,
                    "total_peak_counts": float(peak_depth[class_idx]),
                    "mapped_recurrent_peak_counts": class_count_sum,
                    "mapped_recurrent_fraction": class_count_sum / peak_depth[class_idx] if peak_depth[class_idx] > 0 else math.nan,
                }
            )
            key = (library["biological_sample_id"], cell_class)
            if key not in donor_store:
                donor_store[key] = {
                    **library,
                    "inferred_cell_class": cell_class,
                    "n_cells": 0,
                    "total_peak_counts": 0.0,
                    "counts": np.zeros(n_regions, dtype=np.float32),
                }
            record = donor_store[key]
            record["n_cells"] = int(record["n_cells"]) + n_cells
            record["total_peak_counts"] = float(record["total_peak_counts"]) + float(peak_depth[class_idx])
            record["counts"] += region_counts[class_idx]

    return donor_store, pd.DataFrame(mapping_rows), pd.DataFrame(library_class_rows)


def pcw_adjusted_ols(y: np.ndarray, groups: np.ndarray, pcw: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    valid = np.isfinite(pcw)
    y = y[valid]
    groups = groups[valid].astype(float)
    pcw = pcw[valid].astype(float)
    if y.shape[0] < 5 or len(set(groups.tolist())) < 2:
        n = y.shape[1] if y.ndim == 2 else 0
        return np.full(n, np.nan), np.full(n, np.nan), np.full(n, np.nan)
    pcw_centered = pcw - np.mean(pcw)
    design = np.column_stack([np.ones_like(groups), groups, pcw_centered])
    rank = np.linalg.matrix_rank(design)
    df = y.shape[0] - rank
    if df <= 1:
        return np.full(y.shape[1], np.nan), np.full(y.shape[1], np.nan), np.full(y.shape[1], np.nan)
    xtx_inv = np.linalg.pinv(design.T @ design)
    beta = xtx_inv @ design.T @ y
    residual = y - design @ beta
    sse = np.sum(residual * residual, axis=0)
    sigma2 = sse / df
    se = np.sqrt(np.maximum(sigma2 * xtx_inv[1, 1], 0.0))
    group_beta = beta[1, :]
    t_stat = np.divide(group_beta, se, out=np.zeros_like(group_beta), where=se > 0)
    pvalues = 2.0 * student_t.sf(np.abs(t_stat), df)
    pvalues[~np.isfinite(pvalues)] = 1.0
    return group_beta, t_stat, pvalues


def donor_qc_table(donor_store: dict[tuple[str, str], dict[str, object]]) -> pd.DataFrame:
    rows = []
    for record in donor_store.values():
        rows.append(
            {
                key: value
                for key, value in record.items()
                if key != "counts"
            }
            | {"mapped_recurrent_peak_counts": float(np.asarray(record["counts"]).sum())}
        )
    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    frame["mapped_recurrent_fraction"] = np.where(
        frame["total_peak_counts"] > 0,
        frame["mapped_recurrent_peak_counts"] / frame["total_peak_counts"],
        np.nan,
    )
    return frame


def run_da(
    args: argparse.Namespace,
    regions: pd.DataFrame,
    donor_store: dict[tuple[str, str], dict[str, object]],
    outdir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    all_effects_path = outdir / "gse305146_marker_atac_peak_da_effects.tsv.gz"
    sig_frames: list[pd.DataFrame] = []
    summary_rows: list[dict[str, object]] = []
    top_rows: list[pd.DataFrame] = []
    region_meta = regions.drop(columns=["region_index"]).copy()
    effect_fields = [
        "inferred_cell_class",
        "peak_id",
        "chrom",
        "start",
        "end",
        "background_rank",
        "source_width",
        "n_source_peaks",
        "ds_frequency",
        "control_frequency",
        "delta_frequency_ds_minus_control",
        "max_frequency",
        "n_ds_donors",
        "n_control_donors",
        "n_ds_cells",
        "n_control_cells",
        "mean_ds_log2_peak_cpm",
        "mean_control_log2_peak_cpm",
        "delta_mean_ds_minus_control_log2_peak_cpm",
        "pcw_adjusted_group_beta",
        "pcw_adjusted_group_t",
        "pcw_adjusted_group_p",
        "fdr_bh",
        "significant_da",
        "direction",
        "ds_detected_donor_fraction",
        "control_detected_donor_fraction",
    ]
    with gzip.open(all_effects_path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=effect_fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for cell_class in CLASS_ORDER:
            records = [
                record
                for (_donor, donor_class), record in donor_store.items()
                if donor_class == cell_class
                and int(record["n_cells"]) >= args.min_cells_per_donor_class
                and float(record["total_peak_counts"]) > 0
            ]
            ds_records = [record for record in records if record["group"] == "DS"]
            control_records = [record for record in records if record["group"] == "Control"]
            if len(ds_records) < args.min_donors_per_group or len(control_records) < args.min_donors_per_group:
                summary_rows.append(
                    {
                        "inferred_cell_class": cell_class,
                        "n_tested_peaks": 0,
                        "n_ds_donors": len(ds_records),
                        "n_control_donors": len(control_records),
                        "n_ds_cells": sum(int(record["n_cells"]) for record in ds_records),
                        "n_control_cells": sum(int(record["n_cells"]) for record in control_records),
                        "n_significant_da": 0,
                        "n_ds_open_da": 0,
                        "n_ds_closed_da": 0,
                        "status": "skipped_low_donor_count",
                    }
                )
                continue
            ordered_records = ds_records + control_records
            counts = np.vstack([np.asarray(record["counts"], dtype=np.float64) for record in ordered_records])
            depth = np.array([float(record["total_peak_counts"]) for record in ordered_records], dtype=np.float64)
            groups = np.array([1 if record["group"] == "DS" else 0 for record in ordered_records], dtype=np.int8)
            pcw = np.array([safe_float(record.get("dev_stage_pcw", "")) for record in ordered_records], dtype=np.float64)
            cpm = np.divide(counts * 1_000_000.0, depth[:, None], out=np.zeros_like(counts), where=depth[:, None] > 0)
            log2_cpm = np.log2(cpm + 1.0)
            ds_mask = groups == 1
            control_mask = groups == 0
            mean_ds = log2_cpm[ds_mask].mean(axis=0)
            mean_control = log2_cpm[control_mask].mean(axis=0)
            delta_mean = mean_ds - mean_control
            beta, t_stat, pvalues = pcw_adjusted_ols(log2_cpm, groups, pcw)
            fdr = bh_fdr(pvalues)
            significant = (fdr <= args.fdr_threshold) & (np.abs(beta) >= args.min_abs_beta)
            direction = np.where(beta > 0, "DS_open", "DS_closed")
            direction = np.where(significant, direction, "not_significant")
            ds_detected = (counts[ds_mask] > 0).mean(axis=0)
            control_detected = (counts[control_mask] > 0).mean(axis=0)
            out = region_meta.copy()
            out.insert(0, "inferred_cell_class", cell_class)
            out["n_ds_donors"] = int(ds_mask.sum())
            out["n_control_donors"] = int(control_mask.sum())
            out["n_ds_cells"] = int(sum(int(record["n_cells"]) for record in ds_records))
            out["n_control_cells"] = int(sum(int(record["n_cells"]) for record in control_records))
            out["mean_ds_log2_peak_cpm"] = mean_ds
            out["mean_control_log2_peak_cpm"] = mean_control
            out["delta_mean_ds_minus_control_log2_peak_cpm"] = delta_mean
            out["pcw_adjusted_group_beta"] = beta
            out["pcw_adjusted_group_t"] = t_stat
            out["pcw_adjusted_group_p"] = pvalues
            out["fdr_bh"] = fdr
            out["significant_da"] = significant
            out["direction"] = direction
            out["ds_detected_donor_fraction"] = ds_detected
            out["control_detected_donor_fraction"] = control_detected
            for row in out[effect_fields].to_dict(orient="records"):
                writer.writerow({field: fmt(row.get(field, "")) for field in effect_fields})
            sig = out.loc[significant, effect_fields].copy()
            if not sig.empty:
                sig_frames.append(sig)
                top_open = sig[sig["direction"] == "DS_open"].sort_values("pcw_adjusted_group_beta", ascending=False).head(25)
                top_closed = sig[sig["direction"] == "DS_closed"].sort_values("pcw_adjusted_group_beta", ascending=True).head(25)
                top_rows.extend([top_open, top_closed])
            summary_rows.append(
                {
                    "inferred_cell_class": cell_class,
                    "n_tested_peaks": len(out),
                    "n_ds_donors": int(ds_mask.sum()),
                    "n_control_donors": int(control_mask.sum()),
                    "n_ds_cells": int(sum(int(record["n_cells"]) for record in ds_records)),
                    "n_control_cells": int(sum(int(record["n_cells"]) for record in control_records)),
                    "n_significant_da": int(significant.sum()),
                    "n_ds_open_da": int(((significant) & (beta > 0)).sum()),
                    "n_ds_closed_da": int(((significant) & (beta < 0)).sum()),
                    "median_abs_beta": float(np.nanmedian(np.abs(beta))),
                    "max_positive_beta": float(np.nanmax(beta)),
                    "min_negative_beta": float(np.nanmin(beta)),
                    "status": "tested",
                }
            )
    significant_df = pd.concat(sig_frames, ignore_index=True) if sig_frames else pd.DataFrame(columns=effect_fields)
    top_df = pd.concat(top_rows, ignore_index=True) if top_rows else pd.DataFrame(columns=effect_fields)
    summary_df = pd.DataFrame(summary_rows)
    significant_df.to_csv(
        outdir / "gse305146_marker_atac_peak_da_significant_peaks.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    top_df.to_csv(
        outdir / "gse305146_marker_atac_peak_da_top_peaks.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    summary_df.to_csv(
        outdir / "gse305146_marker_atac_peak_da_class_summary.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    return significant_df, top_df, summary_df


def read_motif_family_masks(motif_hits_path: Path, regions: pd.DataFrame) -> dict[str, np.ndarray]:
    peak_to_index = dict(zip(regions["peak_id"], regions["region_index"]))
    masks = {family: np.zeros(len(regions), dtype=bool) for family in FAMILY_ORDER}
    if not motif_hits_path.exists():
        return masks
    usecols = ["peak_id", "motif_family"]
    for chunk in pd.read_csv(motif_hits_path, sep="\t", usecols=usecols, chunksize=200_000):
        chunk = chunk.dropna()
        for family, sub in chunk.groupby("motif_family", observed=True):
            if family not in masks:
                continue
            indices = [peak_to_index[peak_id] for peak_id in sub["peak_id"].unique() if peak_id in peak_to_index]
            if indices:
                masks[family][indices] = True
    return masks


def motif_enrichment(
    significant_df: pd.DataFrame,
    regions: pd.DataFrame,
    motif_masks: dict[str, np.ndarray],
    args: argparse.Namespace,
) -> pd.DataFrame:
    peak_to_index = dict(zip(regions["peak_id"], regions["region_index"]))
    universe_n = len(regions)
    rows: list[dict[str, object]] = []
    for cell_class in CLASS_ORDER:
        class_sig = significant_df[significant_df["inferred_cell_class"] == cell_class]
        for direction in ["DS_open", "DS_closed"]:
            sig_peaks = class_sig[class_sig["direction"] == direction]["peak_id"].tolist()
            sig_indices = np.array([peak_to_index[peak_id] for peak_id in sig_peaks if peak_id in peak_to_index], dtype=np.int64)
            sig_mask = np.zeros(universe_n, dtype=bool)
            if sig_indices.size:
                sig_mask[sig_indices] = True
            n_sig = int(sig_mask.sum())
            for family in FAMILY_ORDER:
                motif_mask = motif_masks.get(family)
                if motif_mask is None:
                    continue
                motif_n = int(motif_mask.sum())
                a = int(np.logical_and(sig_mask, motif_mask).sum())
                b = n_sig - a
                c = motif_n - a
                d = universe_n - n_sig - c
                if n_sig == 0 or motif_n == 0:
                    odds_ratio = math.nan
                    pvalue = 1.0
                else:
                    odds_ratio, pvalue = fisher_exact([[a, b], [c, d]], alternative="greater")
                rows.append(
                    {
                        "inferred_cell_class": cell_class,
                        "direction": direction,
                        "motif_family": family,
                        "motif_family_label": MOTIF_FAMILY_LABELS.get(family, family),
                        "n_da_peaks": n_sig,
                        "n_motif_hits_in_da": a,
                        "n_motif_hits_in_universe": motif_n,
                        "da_fraction_with_motif": a / n_sig if n_sig else math.nan,
                        "universe_fraction_with_motif": motif_n / universe_n if universe_n else math.nan,
                        "odds_ratio": odds_ratio,
                        "fisher_p": pvalue,
                    }
                )
    out = pd.DataFrame(rows)
    out["fdr_bh"] = bh_fdr(pd.to_numeric(out["fisher_p"], errors="coerce").to_numpy())
    out["significant_enrichment"] = (
        (out["fdr_bh"] <= args.fdr_threshold)
        & (out["odds_ratio"] > 1)
        & (out["n_da_peaks"] > 0)
        & (out["n_motif_hits_in_da"] > 0)
    )
    return out


def plot_da_counts(summary: pd.DataFrame, figdir: Path) -> None:
    figdir.mkdir(parents=True, exist_ok=True)
    frame = summary[summary["inferred_cell_class"].isin(CLASS_ORDER)].copy()
    frame["_order"] = frame["inferred_cell_class"].map({name: idx for idx, name in enumerate(CLASS_ORDER)})
    frame = frame.sort_values("_order")
    labels = [CLASS_LABELS.get(name, name) for name in frame["inferred_cell_class"]]
    x = np.arange(len(frame))
    open_counts = frame["n_ds_open_da"].astype(float).to_numpy()
    closed_counts = -frame["n_ds_closed_da"].astype(float).to_numpy()
    fig, ax = plt.subplots(figsize=(10, 4.8))
    ax.bar(x, open_counts, color="#C34A36", label="DS-open DA peaks")
    ax.bar(x, closed_counts, color="#2F6F9F", label="DS-closed DA peaks")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel("Number of FDR-significant peaks")
    ax.set_title("GSE305146 recurrent-peak differential accessibility")
    ax.legend(frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(figdir / "gse305146_marker_atac_peak_da_counts.png", dpi=300)
    fig.savefig(figdir / "gse305146_marker_atac_peak_da_counts.pdf")
    plt.close(fig)


def plot_motif_heatmap(enrichment: pd.DataFrame, figdir: Path) -> None:
    figdir.mkdir(parents=True, exist_ok=True)
    score = np.zeros((len(FAMILY_ORDER), len(CLASS_ORDER)), dtype=float)
    for family_idx, family in enumerate(FAMILY_ORDER):
        for class_idx, cell_class in enumerate(CLASS_ORDER):
            sub = enrichment[
                (enrichment["motif_family"] == family)
                & (enrichment["inferred_cell_class"] == cell_class)
                & (enrichment["odds_ratio"] > 1)
                & (enrichment["n_motif_hits_in_da"] > 0)
            ].copy()
            if sub.empty:
                continue
            for direction, sign in [("DS_open", 1.0), ("DS_closed", -1.0)]:
                one = sub[sub["direction"] == direction]
                if one.empty:
                    continue
                best = one.sort_values("fdr_bh").iloc[0]
                value = min(-math.log10(max(float(best["fdr_bh"]), 1e-300)), 10.0)
                if abs(value) > abs(score[family_idx, class_idx]):
                    score[family_idx, class_idx] = sign * value
    fig, ax = plt.subplots(figsize=(11, 5.6))
    im = ax.imshow(score, aspect="auto", cmap="RdBu_r", vmin=-10, vmax=10)
    ax.set_xticks(np.arange(len(CLASS_ORDER)))
    ax.set_xticklabels([CLASS_LABELS.get(name, name) for name in CLASS_ORDER], rotation=45, ha="right")
    ax.set_yticks(np.arange(len(FAMILY_ORDER)))
    ax.set_yticklabels([MOTIF_FAMILY_LABELS.get(family, family) for family in FAMILY_ORDER])
    ax.set_title("Motif enrichment among class-resolved DA peaks")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("signed -log10(FDR): red DS-open, blue DS-closed")
    ax.set_xticks(np.arange(-0.5, len(CLASS_ORDER), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(FAMILY_ORDER), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)
    fig.tight_layout()
    fig.savefig(figdir / "gse305146_marker_atac_peak_da_motif_enrichment_heatmap.png", dpi=300)
    fig.savefig(figdir / "gse305146_marker_atac_peak_da_motif_enrichment_heatmap.pdf")
    plt.close(fig)


def write_note(
    path: Path,
    summary: pd.DataFrame,
    enrichment: pd.DataFrame,
    significant: pd.DataFrame,
    args: argparse.Namespace,
) -> None:
    tested = summary[summary["status"] == "tested"].copy()
    tested["_total"] = tested["n_significant_da"].astype(int)
    top_classes = tested.sort_values("_total", ascending=False).head(8)
    enriched = enrichment[
        (enrichment["significant_enrichment"])
        & (enrichment["n_da_peaks"] >= 10)
        & (enrichment["n_motif_hits_in_da"] >= 3)
    ].copy()
    enriched = enriched.sort_values(["fdr_bh", "odds_ratio"], ascending=[True, False]).head(16)
    lines = [
        "# GSE305146 Marker-Inferred Peak-Level Differential Accessibility",
        "",
        "Question: which recurrent fetal-brain ATAC peaks are differentially accessible in DS within marker-inferred cell classes, and do those DA peaks enrich the motif families from the HOMER screen?",
        "",
        "## Methods Snapshot",
        "",
        f"- Retained recurrent peaks: autosomes only unless requested otherwise; FDR <= {args.fdr_threshold:g} and abs(PCW-adjusted group beta) >= {args.min_abs_beta:g}.",
        "- Unit of analysis: biological sample/donor pseudobulk within each marker-inferred cell class.",
        "- Model: peak log2 CPM ~ DS status + post-conception week, fit separately for each inferred class and recurrent peak.",
        "- Motif test: Fisher enrichment of motif-family-bearing peaks among DS-open or DS-closed DA peaks, using the retained recurrent peak background.",
        "",
        "## DA Peak Burden By Inferred Class",
        "",
    ]
    for _, row in top_classes.iterrows():
        lines.append(
            f"- {row['inferred_cell_class']}: {int(row['n_significant_da'])} DA peaks "
            f"({int(row['n_ds_open_da'])} DS-open, {int(row['n_ds_closed_da'])} DS-closed); "
            f"DS donors {int(row['n_ds_donors'])}, controls {int(row['n_control_donors'])}."
        )
    if not significant.empty:
        chrom_summary = (
            significant.groupby(["chrom", "direction"], dropna=False, observed=True)
            .size()
            .reset_index(name="n_significant_peaks")
            .sort_values("n_significant_peaks", ascending=False)
        )
        chrom_text = "; ".join(
            f"{row['chrom']} {row['direction']} {int(row['n_significant_peaks'])}"
            for _, row in chrom_summary.iterrows()
        )
        lines.extend(
            [
                "",
                f"Strict DA chromosome distribution: {chrom_text}.",
            ]
        )
        if set(significant["chrom"].astype(str)) == {"chr21"}:
            lines.append(
                "All strict FDR-significant DA peaks were on chromosome 21, consistent with a strong local trisomy/dosage-linked accessibility component."
            )
    if enriched.empty:
        lines.extend(["", "No motif-family enrichment passed the reporting filter in this first pass."])
    else:
        lines.extend(["", "## Strongest Motif Enrichments In DA Peaks", ""])
        for _, row in enriched.iterrows():
            lines.append(
                f"- {row['motif_family_label']} in {row['inferred_cell_class']} {row['direction']}: "
                f"OR {float(row['odds_ratio']):.3g}, FDR {float(row['fdr_bh']):.3g}, "
                f"{int(row['n_motif_hits_in_da'])}/{int(row['n_da_peaks'])} DA peaks."
            )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This is closer to formal differential accessibility than the motif-deviation proxies because it tests individual recurrent peaks at donor level and adjusts for developmental age. It still uses marker-inferred classes and recurrent-peak overlap from the local 10x features, so it should be interpreted as a first-pass DA screen rather than a fully curated cell-type DA analysis.",
            "",
            "Because strict DA can be sparse and dosage-dominated, motif-family interpretation should be integrated with ranked motif-beta shifts, motif-deviation proxies, HOMER, methylation, and expression layers. A motif family enriched among DS-open or DS-closed DA peaks supports motif-linked regulatory context in that class, but does not prove TF binding or methylation independence.",
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

    regions = read_background_regions(args)
    regions.to_csv(
        outdir / "gse305146_marker_atac_peak_da_tested_regions.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    donor_store, mapping_qc, library_class_qc = collect_donor_pseudobulk(args, regions)
    mapping_qc.to_csv(
        outdir / "gse305146_marker_atac_peak_da_mapping_qc.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    library_class_qc.to_csv(
        outdir / "gse305146_marker_atac_peak_da_library_class_qc.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    donor_qc = donor_qc_table(donor_store)
    donor_qc.to_csv(
        outdir / "gse305146_marker_atac_peak_da_donor_class_qc.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    significant, _top, summary = run_da(args, regions, donor_store, outdir)
    if significant.empty:
        chrom_summary = pd.DataFrame(columns=["chrom", "direction", "n_significant_peaks"])
    else:
        chrom_summary = (
            significant.groupby(["chrom", "direction"], dropna=False, observed=True)
            .size()
            .reset_index(name="n_significant_peaks")
            .sort_values("n_significant_peaks", ascending=False)
        )
    chrom_summary.to_csv(
        outdir / "gse305146_marker_atac_peak_da_significant_chrom_summary.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    motif_masks = read_motif_family_masks(Path(args.motif_hits), regions)
    enrichment = motif_enrichment(significant, regions, motif_masks, args)
    enrichment.to_csv(
        outdir / "gse305146_marker_atac_peak_da_motif_enrichment.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    plot_da_counts(summary, figdir)
    plot_motif_heatmap(enrichment, figdir)
    write_note(Path(args.note), summary, enrichment, significant, args)
    print("Wrote GSE305146 peak-level DA outputs")
    print(outdir / "gse305146_marker_atac_peak_da_class_summary.tsv")
    print(outdir / "gse305146_marker_atac_peak_da_motif_enrichment.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
