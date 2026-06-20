#!/usr/bin/env python3
"""Inferred-cell-class accessibility of GSE305146 HOMER-enriched peak sets."""

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
PEAK_SET_LABELS = {
    "DS_open_recurrent": "DS-open recurrent peaks (homeobox/bHLH enriched)",
    "Control_open_recurrent": "Control-open recurrent peaks (RFX/MEF2 enriched)",
}


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
    parser.add_argument("--outdir", default="results/mechanism_question")
    parser.add_argument("--figdir", default="figures/mechanism_question")
    parser.add_argument("--note", default="notes/gse305146_marker_atac_directional_peaksets_first_pass.md")
    parser.add_argument("--bin-size", type=int, default=100_000)
    parser.add_argument("--min-total-rna", type=float, default=500.0)
    parser.add_argument("--min-detected-genes", type=int, default=200)
    parser.add_argument("--min-top-score", type=float, default=0.05)
    parser.add_argument("--min-confidence-delta", type=float, default=0.05)
    parser.add_argument("--min-cells-per-donor-class", type=int, default=20)
    parser.add_argument("--limit-libraries", type=int, help="Debug: process first N final-analysis libraries.")
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
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
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


def read_peak_regions(path: Path, peak_set: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for index, row in enumerate(read_tsv(path), start=1):
        chrom = row.get("chrom", "")
        start = safe_int(row.get("start", ""))
        end = safe_int(row.get("end", ""))
        if not chrom.startswith("chr") or end <= start:
            continue
        rows.append(
            {
                "region_id": f"{peak_set}:{chrom}:{start}-{end}",
                "peak_set": peak_set,
                "peak_set_label": PEAK_SET_LABELS.get(peak_set, peak_set),
                "chrom": chrom,
                "start": start,
                "end": end,
                "rank": index,
                "ds_frequency": row.get("ds_frequency", ""),
                "control_frequency": row.get("control_frequency", ""),
                "delta_frequency_ds_minus_control": row.get("delta_frequency_ds_minus_control", ""),
                "n_source_peaks": row.get("n_source_peaks", ""),
            }
        )
    return rows


def peak_set_motif_labels(homer_known_path: Path) -> dict[str, str]:
    labels: dict[str, list[str]] = defaultdict(list)
    if not homer_known_path.exists():
        return {}
    for row in read_tsv(homer_known_path):
        direction = row.get("direction", "")
        q_value = float(row.get("q_value", "nan"))
        if not math.isfinite(q_value) or q_value > 0.05:
            continue
        motif = row.get("motif_name", "")
        name = motif.split("(", 1)[0].strip()
        if name and len(labels[direction]) < 8:
            labels[direction].append(name)
    out: dict[str, str] = {}
    if labels.get("DS_open"):
        out["DS_open_recurrent"] = ", ".join(labels["DS_open"])
    if labels.get("Control_open"):
        out["Control_open_recurrent"] = ", ".join(labels["Control_open"])
    return out


def build_region_bins(
    regions: list[dict[str, object]],
    peak_set_index: dict[str, int],
    bin_size: int,
) -> dict[str, dict[int, list[tuple[int, int, int]]]]:
    bins: dict[str, dict[int, list[tuple[int, int, int]]]] = defaultdict(lambda: defaultdict(list))
    for region in regions:
        chrom = str(region["chrom"])
        start = int(region["start"])
        end = int(region["end"])
        set_index = peak_set_index[str(region["peak_set"])]
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
    n_mapped_peak_rows = 0
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
                n_mapped_peak_rows += 1
                for set_index in ordered:
                    mapped_by_set[str(set_index)] += 1

    return row_to_sets, last_gene_expression_row, n_peak_rows, dict(mapped_by_set)


def stream_peak_sets_by_class(
    matrix_path: Path,
    classes: np.ndarray,
    last_gene_expression_row: int,
    row_to_sets: dict[int, list[int]],
    n_sets: int,
) -> tuple[np.ndarray, np.ndarray]:
    n_cells = len(classes)
    total_peak_counts = np.zeros(n_cells, dtype=np.float64)
    set_peak_counts = np.zeros((n_cells, n_sets), dtype=np.float32)
    dims_seen = False

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
            value = float(fields[2])
            total_peak_counts[cell_index] += value
            for set_index in row_to_sets.get(row_index, []):
                set_peak_counts[cell_index, set_index] += value
    return total_peak_counts, set_peak_counts


def summarize_library_peak_sets(
    library: dict[str, str],
    classes: np.ndarray,
    total_peak_counts: np.ndarray,
    set_peak_counts: np.ndarray,
    peak_sets: list[str],
    motif_labels: dict[str, str],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    class_values = sorted(set(classes.astype(str)))
    for cell_class in class_values:
        mask = classes.astype(str) == cell_class
        n_cells = int(mask.sum())
        if n_cells == 0:
            continue
        class_total_peak = float(total_peak_counts[mask].sum())
        for set_index, peak_set in enumerate(peak_sets):
            values = set_peak_counts[mask, set_index]
            set_count = float(values.sum())
            detected = int(np.count_nonzero(values > 0))
            peak_cpm = set_count / class_total_peak * 1_000_000 if class_total_peak > 0 else math.nan
            rows.append(
                {
                    **library,
                    "inferred_cell_class": cell_class,
                    "peak_set": peak_set,
                    "peak_set_label": PEAK_SET_LABELS.get(peak_set, peak_set),
                    "top_known_motifs": motif_labels.get(peak_set, ""),
                    "n_cells": n_cells,
                    "n_cells_with_peak_set_signal": detected,
                    "detected_fraction": detected / n_cells if n_cells else math.nan,
                    "total_peak_counts": class_total_peak,
                    "peak_set_counts": set_count,
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
        "peak_set",
        "peak_set_label",
        "top_known_motifs",
    ]
    numeric = ["n_cells", "n_cells_with_peak_set_signal", "total_peak_counts", "peak_set_counts"]
    donor = frame.groupby(group_cols, dropna=False, observed=True)[numeric].sum().reset_index()
    donor["detected_fraction"] = donor["n_cells_with_peak_set_signal"] / donor["n_cells"]
    donor["peak_cpm"] = np.where(
        donor["total_peak_counts"] > 0,
        donor["peak_set_counts"] / donor["total_peak_counts"] * 1_000_000,
        np.nan,
    )
    donor["log2_peak_cpm"] = np.log2(donor["peak_cpm"] + 1.0)
    return donor


def effect_rows(donor: pd.DataFrame, min_cells: int) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    filtered = donor[donor["n_cells"] >= min_cells].copy()
    for (cell_class, peak_set), sub in filtered.groupby(["inferred_cell_class", "peak_set"], dropna=False, observed=True):
        ds = sub[sub["group"] == "DS"]
        control = sub[sub["group"] == "Control"]
        if ds.empty or control.empty:
            continue
        rows.append(
            {
                "inferred_cell_class": cell_class,
                "peak_set": peak_set,
                "peak_set_label": sub["peak_set_label"].iloc[0],
                "top_known_motifs": sub["top_known_motifs"].iloc[0],
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


def contrast_rows(effects: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for cell_class, sub in effects.groupby("inferred_cell_class", dropna=False, observed=True):
        by_set = {row.peak_set: row for row in sub.itertuples(index=False)}
        if "DS_open_recurrent" not in by_set or "Control_open_recurrent" not in by_set:
            continue
        ds_open = by_set["DS_open_recurrent"]
        control_open = by_set["Control_open_recurrent"]
        rows.append(
            {
                "inferred_cell_class": cell_class,
                "ds_open_delta_log2_peak_cpm": ds_open.delta_ds_minus_control_log2_peak_cpm,
                "control_open_delta_log2_peak_cpm": control_open.delta_ds_minus_control_log2_peak_cpm,
                "ds_open_minus_control_open_delta": ds_open.delta_ds_minus_control_log2_peak_cpm
                - control_open.delta_ds_minus_control_log2_peak_cpm,
                "ds_open_n_ds_donors": ds_open.n_ds_donors,
                "ds_open_n_control_donors": ds_open.n_control_donors,
                "control_open_n_ds_donors": control_open.n_ds_donors,
                "control_open_n_control_donors": control_open.n_control_donors,
                "interpretation": "DS-open peak-set biased toward DS in this class"
                if ds_open.delta_ds_minus_control_log2_peak_cpm - control_open.delta_ds_minus_control_log2_peak_cpm > 0
                else "control-open peak-set relatively stronger in DS-control contrast",
            }
        )
    return pd.DataFrame(rows)


def make_heatmap(matrix: pd.DataFrame, path_prefix: Path, title: str, colorbar_label: str, vmax: float) -> None:
    matrix = matrix.dropna(how="all").dropna(axis=1, how="all")
    if matrix.empty:
        return
    path_prefix.parent.mkdir(parents=True, exist_ok=True)
    values = matrix.to_numpy(dtype=float)
    fig_w = max(7.5, 0.6 * matrix.shape[1] + 3.2)
    fig_h = max(2.8, 0.52 * matrix.shape[0] + 2.1)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(values, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_title(title, fontsize=13, pad=12)
    ax.set_xticks(range(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(matrix.shape[0]))
    ax.set_yticklabels(matrix.index, fontsize=9)
    ax.set_xticks(np.arange(-0.5, matrix.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, matrix.shape[0], 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label(colorbar_label)
    fig.tight_layout()
    fig.savefig(path_prefix.with_suffix(".png"), dpi=300)
    fig.savefig(path_prefix.with_suffix(".pdf"))
    plt.close(fig)


def write_note(
    path: Path,
    libraries: list[dict[str, str]],
    effects: pd.DataFrame,
    contrast: pd.DataFrame,
    region_summary: pd.DataFrame,
) -> None:
    ds_effect = effects[effects["peak_set"] == "DS_open_recurrent"].sort_values(
        "delta_ds_minus_control_log2_peak_cpm", ascending=False
    )
    control_effect = effects[effects["peak_set"] == "Control_open_recurrent"].sort_values(
        "delta_ds_minus_control_log2_peak_cpm", ascending=False
    )
    contrast_top = contrast.sort_values("ds_open_minus_control_open_delta", ascending=False)

    lines = [
        "# GSE305146 Inferred-Class Directional Peak-Set Accessibility",
        "",
        "Question: do inferred fetal cortical classes show DS-control accessibility shifts in the recurrent peak sets that previously drove the GSE305146 HOMER motif screen?",
        "",
        "## Inputs",
        "",
        f"- Final-analysis libraries processed: {len(libraries)}.",
        "- Cell classes: RNA marker-score labels generated per library.",
        "- Peak sets: 5,000 recurrent DS-open peaks and 5,000 recurrent control-open peaks from the previous GSE305146 HOMER input builder.",
        "- Normalization: peak-set counts divided by total ATAC peak counts within each donor/class.",
        "- Interpretation scope: directional peak-set accessibility, not chromVAR, not motif deviation, and not base-resolution methylation.",
        "",
        "## Outputs",
        "",
        "- Library pseudobulk: `results/mechanism_question/gse305146_marker_atac_directional_peakset_library_pseudobulk.tsv`.",
        "- Donor pseudobulk: `results/mechanism_question/gse305146_marker_atac_directional_peakset_donor_pseudobulk.tsv`.",
        "- Effects: `results/mechanism_question/gse305146_marker_atac_directional_peakset_effects.tsv`.",
        "- DS-open versus control-open contrast: `results/mechanism_question/gse305146_marker_atac_directional_peakset_contrast.tsv`.",
        "- Heatmaps: `figures/mechanism_question/gse305146_marker_atac_directional_peakset_*heatmap.png`.",
        "",
        "## Peak-Set Coverage",
        "",
    ]
    for row in region_summary.itertuples(index=False):
        lines.append(
            f"- {row.peak_set}: {int(row.n_input_regions)} input regions; "
            f"mean mapped feature rows per library {row.mean_mapped_feature_rows_per_library:.1f}; "
            f"top motifs: {row.top_known_motifs or 'not available'}."
        )

    lines.extend(["", "## DS-Open HOMER Peak-Set Accessibility", ""])
    for row in ds_effect.head(8).itertuples(index=False):
        lines.append(
            f"- {row.inferred_cell_class}: delta {row.delta_ds_minus_control_log2_peak_cpm:.3g}; "
            f"DS donors {int(row.n_ds_donors)}, control donors {int(row.n_control_donors)}."
        )

    lines.extend(["", "## Control-Open HOMER Peak-Set Accessibility", ""])
    for row in control_effect.head(8).itertuples(index=False):
        lines.append(
            f"- {row.inferred_cell_class}: delta {row.delta_ds_minus_control_log2_peak_cpm:.3g}; "
            f"DS donors {int(row.n_ds_donors)}, control donors {int(row.n_control_donors)}."
        )

    lines.extend(["", "## Directional Contrast", ""])
    lines.append("Positive values mean DS-open peak-set accessibility shifted more toward DS than the control-open peak set in the same inferred class.")
    for row in contrast_top.head(8).itertuples(index=False):
        lines.append(
            f"- {row.inferred_cell_class}: DS-open delta {row.ds_open_delta_log2_peak_cpm:.3g}; "
            f"control-open delta {row.control_open_delta_log2_peak_cpm:.3g}; "
            f"contrast {row.ds_open_minus_control_open_delta:.3g}."
        )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This analysis extends the previous sample-level HOMER result into the recovered multiome matrices. It tests whether the peak sets enriched for DS-open homeobox/bHLH motifs or control-open RFX/MEF2 motifs carry DS-control accessibility shifts within inferred classes.",
            "",
            "Because the peak sets are motif-enriched but individual cells are not scored for motif deviation, this is a bridge analysis rather than true chromVAR. A positive DS-open contrast prioritizes a class for later motif-activity testing; it does not prove direct TF binding or methylation independence.",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def write_figures(effects: pd.DataFrame, contrast: pd.DataFrame, figdir: Path) -> None:
    class_order = [cell_class for cell_class in CLASS_ORDER if cell_class in set(effects["inferred_cell_class"])]
    effect_matrix = effects.pivot_table(
        index="peak_set_label",
        columns="inferred_cell_class",
        values="delta_ds_minus_control_log2_peak_cpm",
        aggfunc="first",
    ).reindex(index=[PEAK_SET_LABELS["DS_open_recurrent"], PEAK_SET_LABELS["Control_open_recurrent"]], columns=class_order)
    effect_matrix.columns = [CLASS_LABELS.get(col, col) for col in effect_matrix.columns]
    make_heatmap(
        effect_matrix,
        figdir / "gse305146_marker_atac_directional_peakset_effect_heatmap",
        "GSE305146 directional recurrent peak-set accessibility",
        "DS - control delta log2(peak CPM + 1)",
        vmax=2.0,
    )

    if contrast.empty:
        return
    contrast_values = contrast.set_index("inferred_cell_class")["ds_open_minus_control_open_delta"].reindex(class_order)
    contrast_matrix = pd.DataFrame(
        [contrast_values.to_numpy()],
        index=["DS-open minus control-open"],
        columns=[CLASS_LABELS.get(col, col) for col in class_order],
    )
    max_abs = float(np.nanmax(np.abs(contrast_values.to_numpy(dtype=float)))) if len(contrast_values) else 1.0
    vmax = max(1.0, min(3.0, math.ceil(max_abs * 10.0) / 10.0))
    make_heatmap(
        contrast_matrix,
        figdir / "gse305146_marker_atac_directional_peakset_contrast_heatmap",
        "GSE305146 DS-open versus control-open peak-set contrast",
        "Delta difference in log2 peak CPM",
        vmax=vmax,
    )


def main() -> int:
    args = parse_args()
    supplementary = Path(args.supplementary)
    outdir = Path(args.outdir)
    figdir = Path(args.figdir)

    if args.plot_only:
        effects = pd.read_csv(outdir / "gse305146_marker_atac_directional_peakset_effects.tsv", sep="\t")
        contrast = pd.read_csv(outdir / "gse305146_marker_atac_directional_peakset_contrast.tsv", sep="\t")
        write_figures(effects, contrast, figdir)
        print("Regenerated directional peak-set figures from existing tables")
        return 0

    peak_sets = ["DS_open_recurrent", "Control_open_recurrent"]
    peak_set_index = {name: idx for idx, name in enumerate(peak_sets)}
    motif_labels = peak_set_motif_labels(Path(args.homer_known))

    regions = read_peak_regions(Path(args.ds_open), "DS_open_recurrent")
    regions.extend(read_peak_regions(Path(args.control_open), "Control_open_recurrent"))
    bins = build_region_bins(regions, peak_set_index, args.bin_size)

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
        for peak_set, set_index in peak_set_index.items():
            mapping_rows.append(
                {
                    **library,
                    "peak_set": peak_set,
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
        total_peak_counts, set_peak_counts = stream_peak_sets_by_class(
            matrix_path,
            classes,
            last_gene_row,
            row_to_sets,
            n_sets=len(peak_sets),
        )
        library_rows.extend(
            summarize_library_peak_sets(library, classes, total_peak_counts, set_peak_counts, peak_sets, motif_labels)
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
        "peak_set",
        "peak_set_label",
        "top_known_motifs",
        "n_cells",
        "n_cells_with_peak_set_signal",
        "detected_fraction",
        "total_peak_counts",
        "peak_set_counts",
        "peak_cpm",
        "log2_peak_cpm",
    ]
    write_tsv(outdir / "gse305146_marker_atac_directional_peakset_library_pseudobulk.tsv", library_rows, library_fields)

    mapping = pd.DataFrame(mapping_rows)
    mapping.to_csv(outdir / "gse305146_marker_atac_directional_peakset_mapping_qc.tsv", sep="\t", index=False)

    donor = combine_donor_pseudobulk(library_rows)
    donor.to_csv(outdir / "gse305146_marker_atac_directional_peakset_donor_pseudobulk.tsv", sep="\t", index=False)
    effects = effect_rows(donor, args.min_cells_per_donor_class)
    effects.to_csv(outdir / "gse305146_marker_atac_directional_peakset_effects.tsv", sep="\t", index=False)
    contrast = contrast_rows(effects)
    contrast.to_csv(outdir / "gse305146_marker_atac_directional_peakset_contrast.tsv", sep="\t", index=False)

    region_counts = pd.DataFrame(regions).groupby("peak_set", dropna=False)["region_id"].nunique()
    region_summary_rows = []
    for peak_set in peak_sets:
        set_mapping = mapping[mapping["peak_set"] == peak_set]
        region_summary_rows.append(
            {
                "peak_set": peak_set,
                "peak_set_label": PEAK_SET_LABELS.get(peak_set, peak_set),
                "n_input_regions": int(region_counts.get(peak_set, 0)),
                "mean_mapped_feature_rows_per_library": float(set_mapping["n_mapped_feature_peak_rows"].mean())
                if not set_mapping.empty
                else math.nan,
                "min_mapped_feature_rows_per_library": int(set_mapping["n_mapped_feature_peak_rows"].min())
                if not set_mapping.empty
                else 0,
                "max_mapped_feature_rows_per_library": int(set_mapping["n_mapped_feature_peak_rows"].max())
                if not set_mapping.empty
                else 0,
                "top_known_motifs": motif_labels.get(peak_set, ""),
            }
        )
    region_summary = pd.DataFrame(region_summary_rows)
    region_summary.to_csv(outdir / "gse305146_marker_atac_directional_peakset_summary.tsv", sep="\t", index=False)

    write_figures(effects, contrast, figdir)

    write_note(Path(args.note), libraries, effects, contrast, region_summary)
    print(f"Wrote directional peak-set summaries for {len(libraries)} libraries")
    print(outdir / "gse305146_marker_atac_directional_peakset_effects.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
