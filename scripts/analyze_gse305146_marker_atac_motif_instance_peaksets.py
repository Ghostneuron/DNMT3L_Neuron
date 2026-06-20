#!/usr/bin/env python3
"""Inferred-class accessibility of motif-instance peak subsets in GSE305146."""

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

from infer_gse305146_marker_cell_classes_and_targets import (
    MARKER_SETS,
    infer_classes,
    parse_features as parse_marker_features,
    selected_libraries,
    stream_matrix,
)


HOMER_BIN = Path("/Users/jili/tools/homer/bin")
ANNOTATE_PEAKS = HOMER_BIN / "annotatePeaks.pl"
STAT_CONTROL_MOTIFS = [
    {
        "direction": "DS_open",
        "motif_path": HOMER_BIN.parent / "motifs" / "stat1.motif",
        "motif_name": "STAT1(Stat)/HelaS3-STAT1-ChIP-Seq(GSE12782)/Homer",
        "motif_family": "STAT_CONTROL",
        "q_value": math.nan,
        "rank": "control",
        "context_enriched": "no",
    },
    {
        "direction": "DS_open",
        "motif_path": HOMER_BIN.parent / "motifs" / "stat3.motif",
        "motif_name": "Stat3(Stat)/mES-Stat3-ChIP-Seq(GSE11431)/Homer",
        "motif_family": "STAT_CONTROL",
        "q_value": math.nan,
        "rank": "control",
        "context_enriched": "no",
    },
    {
        "direction": "Control_open",
        "motif_path": HOMER_BIN.parent / "motifs" / "stat1.motif",
        "motif_name": "STAT1(Stat)/HelaS3-STAT1-ChIP-Seq(GSE12782)/Homer",
        "motif_family": "STAT_CONTROL",
        "q_value": math.nan,
        "rank": "control",
        "context_enriched": "no",
    },
    {
        "direction": "Control_open",
        "motif_path": HOMER_BIN.parent / "motifs" / "stat3.motif",
        "motif_name": "Stat3(Stat)/mES-Stat3-ChIP-Seq(GSE11431)/Homer",
        "motif_family": "STAT_CONTROL",
        "q_value": math.nan,
        "rank": "control",
        "context_enriched": "no",
    },
]

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
MOTIF_FAMILY_LABELS = {
    "LHX_HOMEBOX": "LHX homeobox",
    "DLX_HOMEBOX": "DLX homeobox",
    "GSX_HOMEBOX": "GSX homeobox",
    "ISL_HOMEBOX": "ISL homeobox",
    "EMX_HOMEBOX": "EMX homeobox",
    "NKX_HOMEBOX": "NKX homeobox",
    "OLIG_BHLH": "OLIG bHLH",
    "NEURO_BHLH": "Neurogenic bHLH",
    "RFX_XBOX": "RFX/X-box",
    "MEF2": "MEF2",
    "STAT_CONTROL": "STAT control",
}
FAMILY_ORDER = [
    "LHX_HOMEBOX",
    "DLX_HOMEBOX",
    "GSX_HOMEBOX",
    "ISL_HOMEBOX",
    "EMX_HOMEBOX",
    "NKX_HOMEBOX",
    "OLIG_BHLH",
    "NEURO_BHLH",
    "RFX_XBOX",
    "MEF2",
    "STAT_CONTROL",
]
SELECTED_FAMILIES = set(FAMILY_ORDER)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", default="data/geo/GSE305146/metadata/sample_metadata.tsv")
    parser.add_argument("--supplementary", default="data/geo/GSE305146/supplementary")
    parser.add_argument(
        "--ds-open",
        default="results/gse305146_scATAC_homer_motif_enrichment/tables/GSE305146_DS_open_recurrent_peak_calls.tsv",
    )
    parser.add_argument(
        "--control-open",
        default="results/gse305146_scATAC_homer_motif_enrichment/tables/GSE305146_Control_open_recurrent_peak_calls.tsv",
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
    parser.add_argument("--note", default="notes/gse305146_marker_atac_motif_instance_peaksets_first_pass.md")
    parser.add_argument("--q-threshold", type=float, default=0.05)
    parser.add_argument("--bin-size", type=int, default=100_000)
    parser.add_argument("--min-total-rna", type=float, default=500.0)
    parser.add_argument("--min-detected-genes", type=int, default=200)
    parser.add_argument("--min-top-score", type=float, default=0.05)
    parser.add_argument("--min-confidence-delta", type=float, default=0.05)
    parser.add_argument("--min-cells-per-donor-class", type=int, default=20)
    parser.add_argument("--limit-libraries", type=int, help="Debug: process first N final-analysis libraries.")
    parser.add_argument("--skip-homer", action="store_true", help="Reuse existing annotatePeaks outputs.")
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


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(fieldnames=fields, f=handle, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", errors="replace") as handle:
        return [{k: (v or "").strip() for k, v in row.items()} for row in csv.DictReader(handle, delimiter="\t")]


def safe_int(value: str) -> int:
    try:
        return int(float(value))
    except ValueError:
        return 0


def safe_float(value: str) -> float:
    try:
        return float(value)
    except ValueError:
        return math.nan


def motif_family(motif_name: str) -> str:
    lower = motif_name.lower()
    if "stat" in lower:
        return "STAT_CONTROL"
    if "rfx" in lower or "x-box" in lower or "xbox" in lower:
        return "RFX_XBOX"
    if "mef2" in lower:
        return "MEF2"
    if "olig" in lower:
        return "OLIG_BHLH"
    if any(token in lower for token in ["neurod", "atoh", "tcf4", "twist", "ptf1a"]):
        return "NEURO_BHLH"
    if "lhx" in lower:
        return "LHX_HOMEBOX"
    if "dlx" in lower:
        return "DLX_HOMEBOX"
    if "gsx" in lower:
        return "GSX_HOMEBOX"
    if "isl" in lower:
        return "ISL_HOMEBOX"
    if "emx" in lower:
        return "EMX_HOMEBOX"
    if "nkx" in lower:
        return "NKX_HOMEBOX"
    return "OTHER"


def homer_set_dir(homer_out: Path, direction: str) -> Path:
    if direction == "DS_open":
        return homer_out / "GSE305146_fetal_DS_open_recurrent_scATAC_peaks_hg38"
    if direction == "Control_open":
        return homer_out / "GSE305146_fetal_Control_open_recurrent_scATAC_peaks_hg38"
    raise ValueError(direction)


def select_motif_records(homer_known_path: Path, homer_out: Path, q_threshold: float) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for row in read_tsv(homer_known_path):
        direction = row.get("direction", "")
        if direction not in {"DS_open", "Control_open"}:
            continue
        q_value = safe_float(row.get("q_value", ""))
        if not math.isfinite(q_value) or q_value > q_threshold:
            continue
        family = motif_family(row.get("motif_name", ""))
        if family not in SELECTED_FAMILIES or family == "STAT_CONTROL":
            continue
        rank = safe_int(row.get("rank", ""))
        motif_path = homer_set_dir(homer_out, direction) / "knownResults" / f"known{rank}.motif"
        if not motif_path.exists():
            continue
        records.append(
            {
                "direction": direction,
                "motif_path": motif_path,
                "motif_name": row.get("motif_name", ""),
                "motif_family": family,
                "q_value": q_value,
                "rank": rank,
                "target_percent": row.get("target_percent", ""),
                "background_percent": row.get("background_percent", ""),
                "context_enriched": "yes",
            }
        )

    for record in STAT_CONTROL_MOTIFS:
        if Path(record["motif_path"]).exists():
            records.append(record)

    dedup: dict[tuple[str, Path], dict[str, object]] = {}
    for record in records:
        dedup[(str(record["direction"]), Path(record["motif_path"]))] = record
    return list(dedup.values())


def read_peak_regions(path: Path, source_peak_set: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for rank, row in enumerate(read_tsv(path), start=1):
        chrom = row.get("chrom", "")
        start = safe_int(row.get("start", ""))
        end = safe_int(row.get("end", ""))
        if not chrom.startswith("chr") or end <= start:
            continue
        rows.append(
            {
                "peak_id": f"{source_peak_set}:{chrom}:{start}-{end}",
                "source_peak_set": source_peak_set,
                "chrom": chrom,
                "start": start,
                "end": end,
                "rank": rank,
                "ds_frequency": row.get("ds_frequency", ""),
                "control_frequency": row.get("control_frequency", ""),
                "delta_frequency_ds_minus_control": row.get("delta_frequency_ds_minus_control", ""),
            }
        )
    return rows


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
                        str(int(region.get("rank", 0))),
                        "+",
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

    hit_rows: list[dict[str, object]] = []
    region_sets: dict[tuple[str, str, str], dict[str, object]] = {}
    for row in table.to_dict(orient="records"):
        peak_id = str(row[table.columns[0]])
        region = regions_by_id.get(peak_id)
        if region is None:
            continue
        for column, motif_name in motif_columns:
            count = safe_int(str(row.get(column, "")))
            if count <= 0:
                continue
            motif = motif_by_name.get(motif_name)
            if motif is None:
                continue
            source_peak_set = str(region["source_peak_set"])
            family = str(motif["motif_family"])
            set_id = f"{source_peak_set}|{family}"
            set_label = f"{source_peak_set.replace('_', ' ')}: {MOTIF_FAMILY_LABELS.get(family, family)}"
            key = (peak_id, source_peak_set, family)
            region_sets[key] = {
                **region,
                "region_set_id": set_id,
                "region_set_label": set_label,
                "motif_family": family,
                "motif_family_label": MOTIF_FAMILY_LABELS.get(family, family),
                "context_enriched": motif.get("context_enriched", ""),
            }
            hit_rows.append(
                {
                    **region,
                    "region_set_id": set_id,
                    "region_set_label": set_label,
                    "motif_family": family,
                    "motif_family_label": MOTIF_FAMILY_LABELS.get(family, family),
                    "motif_name": motif_name,
                    "motif_count": count,
                    "q_value": motif.get("q_value", ""),
                    "rank": motif.get("rank", ""),
                    "context_enriched": motif.get("context_enriched", ""),
                }
            )
    return list(region_sets.values()), hit_rows


def build_region_bins(
    regions: list[dict[str, object]],
    region_set_index: dict[str, int],
    bin_size: int,
) -> dict[str, dict[int, list[tuple[int, int, int]]]]:
    bins: dict[str, dict[int, list[tuple[int, int, int]]]] = defaultdict(lambda: defaultdict(list))
    for region in regions:
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


def summarize_library(
    library: dict[str, str],
    classes: np.ndarray,
    total_peak_counts: np.ndarray,
    set_peak_counts: np.ndarray,
    region_sets: list[dict[str, object]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for cell_class in sorted(set(classes.astype(str))):
        mask = classes.astype(str) == cell_class
        n_cells = int(mask.sum())
        if n_cells == 0:
            continue
        class_total_peak = float(total_peak_counts[mask].sum())
        for set_index, region_set in enumerate(region_sets):
            values = set_peak_counts[mask, set_index]
            set_count = float(values.sum())
            detected = int(np.count_nonzero(values > 0))
            peak_cpm = set_count / class_total_peak * 1_000_000 if class_total_peak > 0 else math.nan
            rows.append(
                {
                    **library,
                    "inferred_cell_class": cell_class,
                    "region_set_id": region_set["region_set_id"],
                    "region_set_label": region_set["region_set_label"],
                    "source_peak_set": region_set["source_peak_set"],
                    "motif_family": region_set["motif_family"],
                    "motif_family_label": region_set["motif_family_label"],
                    "context_enriched": region_set["context_enriched"],
                    "n_regions_in_set": region_set["n_regions_in_set"],
                    "n_cells": n_cells,
                    "n_cells_with_region_set_signal": detected,
                    "detected_fraction": detected / n_cells if n_cells else math.nan,
                    "total_peak_counts": class_total_peak,
                    "region_set_peak_counts": set_count,
                    "peak_cpm": peak_cpm,
                    "log2_peak_cpm": math.log2(peak_cpm + 1.0) if math.isfinite(peak_cpm) else math.nan,
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
        "region_set_label",
        "source_peak_set",
        "motif_family",
        "motif_family_label",
        "context_enriched",
        "n_regions_in_set",
    ]
    numeric = ["n_cells", "n_cells_with_region_set_signal", "total_peak_counts", "region_set_peak_counts"]
    donor = frame.groupby(group_cols, dropna=False, observed=True)[numeric].sum().reset_index()
    donor["detected_fraction"] = donor["n_cells_with_region_set_signal"] / donor["n_cells"]
    donor["peak_cpm"] = np.where(
        donor["total_peak_counts"] > 0,
        donor["region_set_peak_counts"] / donor["total_peak_counts"] * 1_000_000,
        np.nan,
    )
    donor["log2_peak_cpm"] = np.log2(donor["peak_cpm"] + 1.0)
    return donor


def effect_rows(donor: pd.DataFrame, min_cells: int) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    filtered = donor[donor["n_cells"] >= min_cells].copy()
    for (cell_class, set_id), sub in filtered.groupby(["inferred_cell_class", "region_set_id"], dropna=False, observed=True):
        ds = sub[sub["group"] == "DS"]
        control = sub[sub["group"] == "Control"]
        if ds.empty or control.empty:
            continue
        rows.append(
            {
                "inferred_cell_class": cell_class,
                "region_set_id": set_id,
                "region_set_label": sub["region_set_label"].iloc[0],
                "source_peak_set": sub["source_peak_set"].iloc[0],
                "motif_family": sub["motif_family"].iloc[0],
                "motif_family_label": sub["motif_family_label"].iloc[0],
                "context_enriched": sub["context_enriched"].iloc[0],
                "n_regions_in_set": int(sub["n_regions_in_set"].iloc[0]),
                "n_ds_donors": int(ds["biological_sample_id"].nunique()),
                "n_control_donors": int(control["biological_sample_id"].nunique()),
                "n_ds_cells": int(ds["n_cells"].sum()),
                "n_control_cells": int(control["n_cells"].sum()),
                "mean_ds_log2_peak_cpm": float(ds["log2_peak_cpm"].mean()),
                "mean_control_log2_peak_cpm": float(control["log2_peak_cpm"].mean()),
                "delta_ds_minus_control_log2_peak_cpm": float(ds["log2_peak_cpm"].mean() - control["log2_peak_cpm"].mean()),
                "mean_ds_detected_fraction": float(ds["detected_fraction"].mean()),
                "mean_control_detected_fraction": float(control["detected_fraction"].mean()),
                "delta_detected_fraction": float(ds["detected_fraction"].mean() - control["detected_fraction"].mean()),
                "mean_ds_peak_cpm": float(ds["peak_cpm"].mean()),
                "mean_control_peak_cpm": float(control["peak_cpm"].mean()),
                "delta_ds_minus_control_peak_cpm": float(ds["peak_cpm"].mean() - control["peak_cpm"].mean()),
            }
        )
    return pd.DataFrame(rows)


def make_heatmap(matrix: pd.DataFrame, path_prefix: Path, title: str, vmax: float) -> None:
    matrix = matrix.dropna(how="all").dropna(axis=1, how="all")
    if matrix.empty:
        return
    path_prefix.parent.mkdir(parents=True, exist_ok=True)
    values = matrix.to_numpy(dtype=float)
    fig_w = max(8.0, 0.62 * matrix.shape[1] + 4.0)
    fig_h = max(4.2, 0.34 * matrix.shape[0] + 2.2)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(values, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_title(title, fontsize=13, pad=12)
    ax.set_xticks(range(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(matrix.shape[0]))
    ax.set_yticklabels(matrix.index, fontsize=8)
    ax.set_xticks(np.arange(-0.5, matrix.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, matrix.shape[0], 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.030, pad=0.02)
    cbar.set_label("DS - control delta log2(peak CPM + 1)")
    fig.tight_layout()
    fig.savefig(path_prefix.with_suffix(".png"), dpi=300)
    fig.savefig(path_prefix.with_suffix(".pdf"))
    plt.close(fig)


def write_figures(effects: pd.DataFrame, figdir: Path) -> None:
    class_order = [cell_class for cell_class in CLASS_ORDER if cell_class in set(effects["inferred_cell_class"])]
    row_order: list[str] = []
    for source in ["DS_open_recurrent", "Control_open_recurrent"]:
        for family in FAMILY_ORDER:
            set_id = f"{source}|{family}"
            if set_id in set(effects["region_set_id"]):
                row_order.append(set_id)
    matrix = effects.pivot_table(
        index="region_set_id",
        columns="inferred_cell_class",
        values="delta_ds_minus_control_log2_peak_cpm",
        aggfunc="first",
    ).reindex(index=row_order, columns=class_order)
    labels = effects.drop_duplicates("region_set_id").set_index("region_set_id")["region_set_label"].to_dict()
    matrix.index = [labels.get(idx, idx).replace("_", " ") for idx in matrix.index]
    matrix.columns = [CLASS_LABELS.get(col, col) for col in matrix.columns]
    make_heatmap(
        matrix,
        figdir / "gse305146_marker_atac_motif_instance_peakset_heatmap",
        "GSE305146 motif-instance peak-subset accessibility",
        vmax=1.5,
    )


def write_note(
    path: Path,
    libraries: list[dict[str, str]],
    region_set_summary: pd.DataFrame,
    effects: pd.DataFrame,
) -> None:
    ds = effects[effects["source_peak_set"] == "DS_open_recurrent"].sort_values(
        "delta_ds_minus_control_log2_peak_cpm", ascending=False
    )
    control = effects[effects["source_peak_set"] == "Control_open_recurrent"].sort_values(
        "delta_ds_minus_control_log2_peak_cpm", ascending=True
    )
    stat = effects[effects["motif_family"] == "STAT_CONTROL"].copy()
    stat["abs_delta"] = stat["delta_ds_minus_control_log2_peak_cpm"].abs()
    stat_top = stat.sort_values("abs_delta", ascending=False)

    lines = [
        "# GSE305146 Inferred-Class Motif-Instance Peak-Subset Accessibility",
        "",
        "Question: which HOMER motif families inside DS-open/control-open recurrent peaks carry inferred-cell-class accessibility shifts?",
        "",
        "## Inputs",
        "",
        f"- Final-analysis libraries processed: {len(libraries)}.",
        "- Motif instances: HOMER `annotatePeaks.pl -nmotifs` scans of recurrent DS-open and control-open peak sets.",
        "- Motif families: selected significant homeobox/bHLH families from DS-open peaks, RFX/MEF2 families from control-open peaks, plus STAT1/STAT3 non-enriched control motifs.",
        "- Accessibility: motif-hit peak subsets projected into marker-inferred GSE305146 donor/classes and normalized by total ATAC peak counts.",
        "- Interpretation scope: motif-instance accessibility proxy, not chromVAR, not TF binding, and not DNA methylation.",
        "",
        "## Outputs",
        "",
        "- Motif records: `results/mechanism_question/gse305146_marker_atac_motif_instance_motif_records.tsv`.",
        "- Region motif hits: `results/mechanism_question/gse305146_marker_atac_motif_instance_region_hits.tsv`.",
        "- Region set summary: `results/mechanism_question/gse305146_marker_atac_motif_instance_region_set_summary.tsv`.",
        "- Effects: `results/mechanism_question/gse305146_marker_atac_motif_instance_peakset_effects.tsv`.",
        "- Heatmap: `figures/mechanism_question/gse305146_marker_atac_motif_instance_peakset_heatmap.png`.",
        "",
        "## Motif-Instance Region Sets",
        "",
    ]
    for row in region_set_summary.sort_values(["source_peak_set", "motif_family"]).to_dict(orient="records"):
        lines.append(
            f"- {row['region_set_label']}: {int(row['n_regions_in_set'])} motif-hit regions; "
            f"{int(row['n_motif_names'])} motif(s); context enriched={row['context_enriched']}."
        )

    lines.extend(["", "## Strongest DS-Open Motif-Subset Effects", ""])
    for row in ds.head(10).to_dict(orient="records"):
        lines.append(
            f"- {row['region_set_label']} in {row['inferred_cell_class']}: "
            f"delta {float(row['delta_ds_minus_control_log2_peak_cpm']):.3g}; "
            f"DS donors {int(row['n_ds_donors'])}, control donors {int(row['n_control_donors'])}."
        )

    lines.extend(["", "## Strongest Control-Open Motif-Subset Negative Effects", ""])
    for row in control.head(10).to_dict(orient="records"):
        lines.append(
            f"- {row['region_set_label']} in {row['inferred_cell_class']}: "
            f"delta {float(row['delta_ds_minus_control_log2_peak_cpm']):.3g}; "
            f"DS donors {int(row['n_ds_donors'])}, control donors {int(row['n_control_donors'])}."
        )

    lines.extend(["", "## STAT Control Motifs", ""])
    for row in stat_top.head(8).to_dict(orient="records"):
        lines.append(
            f"- {row['region_set_label']} in {row['inferred_cell_class']}: "
            f"delta {float(row['delta_ds_minus_control_log2_peak_cpm']):.3g}; context enriched={row['context_enriched']}."
        )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This refines the directional peak-set result by asking which motif-instance peak subsets carry the accessibility shift. The positive signal remains concentrated in DS-open developmental homeobox/bHLH motif subsets, while control-open RFX/MEF2 subsets are generally DS-negative.",
            "",
            "STAT1/STAT3 motif-instance subsets are included as non-enriched controls. Any STAT accessibility differences here should be interpreted as ordinary motif-bearing peak accessibility, not evidence that STAT motifs were enriched or that STAT is the dominant direct chromatin mechanism.",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    figdir = Path(args.figdir)
    if args.plot_only:
        effects = pd.read_csv(outdir / "gse305146_marker_atac_motif_instance_peakset_effects.tsv", sep="\t")
        write_figures(effects, figdir)
        print("Regenerated motif-instance peak-subset figures")
        return 0

    supplementary = Path(args.supplementary)
    scan_dir = outdir / "gse305146_marker_atac_motif_instance_peaksets"
    input_dir = scan_dir / "inputs"
    annotate_dir = scan_dir / "annotatePeaks"
    input_dir.mkdir(parents=True, exist_ok=True)
    annotate_dir.mkdir(parents=True, exist_ok=True)

    motif_records = select_motif_records(Path(args.homer_known), Path(args.homer_out), args.q_threshold)
    motif_fields = [
        "direction",
        "motif_path",
        "motif_name",
        "motif_family",
        "q_value",
        "rank",
        "target_percent",
        "background_percent",
        "context_enriched",
    ]
    write_tsv(outdir / "gse305146_marker_atac_motif_instance_motif_records.tsv", motif_records, motif_fields)

    ds_regions = read_peak_regions(Path(args.ds_open), "DS_open_recurrent")
    control_regions = read_peak_regions(Path(args.control_open), "Control_open_recurrent")
    regions_by_direction = {
        "DS_open": ds_regions,
        "Control_open": control_regions,
    }

    all_region_sets: list[dict[str, object]] = []
    all_hit_rows: list[dict[str, object]] = []
    scan_manifest: list[dict[str, object]] = []
    for direction, regions in regions_by_direction.items():
        direction_motifs = [record for record in motif_records if record["direction"] == direction]
        if not direction_motifs:
            continue
        bed_path = input_dir / f"gse305146_{direction}_recurrent_peaks_hg38.bed"
        output_path = annotate_dir / f"gse305146_{direction}_motif_instances.annotatePeaks.nmotifs.tsv"
        log_path = annotate_dir / f"gse305146_{direction}_motif_instances.annotatePeaks.log"
        write_bed(bed_path, regions)
        if not args.skip_homer or not output_path.exists():
            run_homer_scan(bed_path, direction_motifs, output_path, log_path)
        region_sets, hit_rows = parse_homer_scan(output_path, {row["peak_id"]: row for row in regions}, direction_motifs)
        all_region_sets.extend(region_sets)
        all_hit_rows.extend(hit_rows)
        scan_manifest.append(
            {
                "direction": direction,
                "bed": bed_path,
                "annotatePeaks_output": output_path,
                "annotatePeaks_log": log_path,
                "n_regions_scanned": len(regions),
                "n_motifs_scanned": len(direction_motifs),
                "n_region_family_sets": len(region_sets),
                "n_region_motif_hits": len(hit_rows),
            }
        )

    write_tsv(
        outdir / "gse305146_marker_atac_motif_instance_scan_manifest.tsv",
        scan_manifest,
        [
            "direction",
            "bed",
            "annotatePeaks_output",
            "annotatePeaks_log",
            "n_regions_scanned",
            "n_motifs_scanned",
            "n_region_family_sets",
            "n_region_motif_hits",
        ],
    )
    region_hit_fields = [
        "peak_id",
        "source_peak_set",
        "chrom",
        "start",
        "end",
        "region_set_id",
        "region_set_label",
        "motif_family",
        "motif_family_label",
        "motif_name",
        "motif_count",
        "q_value",
        "rank",
        "context_enriched",
    ]
    write_tsv(outdir / "gse305146_marker_atac_motif_instance_region_hits.tsv", all_hit_rows, region_hit_fields)

    region_sets_df = pd.DataFrame(all_region_sets).drop_duplicates(
        ["peak_id", "source_peak_set", "motif_family"]
    )
    region_set_summary = (
        region_sets_df.groupby(
            ["region_set_id", "region_set_label", "source_peak_set", "motif_family", "motif_family_label", "context_enriched"],
            dropna=False,
            observed=True,
        )
        .agg(n_regions_in_set=("peak_id", "nunique"))
        .reset_index()
    )
    hit_frame = pd.DataFrame(all_hit_rows)
    if not hit_frame.empty:
        motif_counts = (
            hit_frame.groupby("region_set_id", dropna=False, observed=True)["motif_name"]
            .nunique()
            .rename("n_motif_names")
            .reset_index()
        )
        region_set_summary = region_set_summary.merge(motif_counts, on="region_set_id", how="left")
    else:
        region_set_summary["n_motif_names"] = 0
    region_set_summary["n_motif_names"] = region_set_summary["n_motif_names"].fillna(0).astype(int)
    region_set_summary.to_csv(
        outdir / "gse305146_marker_atac_motif_instance_region_set_summary.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )

    region_sets = region_set_summary.sort_values(["source_peak_set", "motif_family"]).to_dict(orient="records")
    region_set_index = {str(row["region_set_id"]): idx for idx, row in enumerate(region_sets)}
    bins = build_region_bins(region_sets_df.to_dict(orient="records"), region_set_index, args.bin_size)

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
        total_peak_counts, set_peak_counts = stream_atac_sets(
            matrix_path, last_gene_row, row_to_sets, len(region_sets)
        )
        library_rows.extend(summarize_library(library, classes, total_peak_counts, set_peak_counts, region_sets))

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
        "region_set_label",
        "source_peak_set",
        "motif_family",
        "motif_family_label",
        "context_enriched",
        "n_regions_in_set",
        "n_cells",
        "n_cells_with_region_set_signal",
        "detected_fraction",
        "total_peak_counts",
        "region_set_peak_counts",
        "peak_cpm",
        "log2_peak_cpm",
    ]
    write_tsv(outdir / "gse305146_marker_atac_motif_instance_peakset_library_pseudobulk.tsv", library_rows, library_fields)
    pd.DataFrame(mapping_rows).to_csv(
        outdir / "gse305146_marker_atac_motif_instance_mapping_qc.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    donor = combine_donor_pseudobulk(library_rows)
    donor.to_csv(
        outdir / "gse305146_marker_atac_motif_instance_peakset_donor_pseudobulk.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    effects = effect_rows(donor, args.min_cells_per_donor_class)
    effects.to_csv(
        outdir / "gse305146_marker_atac_motif_instance_peakset_effects.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    write_figures(effects, figdir)
    write_note(Path(args.note), libraries, region_set_summary, effects)
    print(f"Wrote motif-instance peak-subset summaries for {len(libraries)} libraries")
    print(outdir / "gse305146_marker_atac_motif_instance_peakset_effects.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
