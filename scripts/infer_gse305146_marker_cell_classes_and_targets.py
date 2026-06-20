#!/usr/bin/env python3
"""Marker-based cell-class inference and target summaries for GSE305146 multiome matrices."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


MARKER_SETS: dict[str, list[str]] = {
    "RG": ["SOX2", "PAX6", "VIM", "HES1", "HOPX", "FABP7", "NES", "GLI3"],
    "Cycling_prog": ["MKI67", "TOP2A", "CENPF", "PCLAF", "UBE2C", "BIRC5", "ASPM"],
    "IPC": ["EOMES", "NEUROG2", "NEUROD1", "INSM1", "DLL3", "HES6"],
    "NEU_CUX2": ["CUX1", "CUX2", "SATB2", "POU3F2", "POU3F3"],
    "NEU_TLE4_deep": ["TLE4", "BCL11B", "FEZF2", "TBR1", "FOXP2", "SOX5"],
    "NEU_RORB": ["RORB", "PCP4", "TOX", "FOXP1"],
    "NEU_RELN": ["RELN", "LHX5", "CRABP1", "CALB2"],
    "GABA_SST": ["GAD1", "GAD2", "DLX1", "DLX2", "SST", "CALB2", "RELN"],
    "AST": ["GFAP", "AQP4", "ALDH1L1", "SLC1A3", "SLC1A2", "CLU", "SOX9"],
    "OPC_OL": ["PDGFRA", "CSPG4", "OLIG1", "OLIG2", "SOX10", "BCAN", "MBP", "PLP1", "MOG"],
    "MIC": ["AIF1", "PTPRC", "CX3CR1", "C1QA", "C1QB", "CSF1R", "TYROBP"],
    "VASC": ["PECAM1", "CLDN5", "KDR", "FLT1", "RGS5", "PDGFRB", "COL1A1", "COL3A1"],
}

CORE_GENES = ["DNMT3L", "DNMT3A", "DNMT3B", "STAT1", "STAT3"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", default="data/geo/GSE305146/metadata/sample_metadata.tsv")
    parser.add_argument("--supplementary", default="data/geo/GSE305146/supplementary")
    parser.add_argument(
        "--mechanism-axis",
        default="results/mechanism_question/final_candidate_mechanism_axis_summary.tsv",
    )
    parser.add_argument("--outdir", default="results/mechanism_question")
    parser.add_argument("--figdir", default="figures/mechanism_question")
    parser.add_argument("--note", default="notes/gse305146_marker_cell_class_first_pass.md")
    parser.add_argument("--min-total-rna", type=float, default=500.0)
    parser.add_argument("--min-detected-genes", type=int, default=200)
    parser.add_argument("--min-top-score", type=float, default=0.05)
    parser.add_argument("--min-confidence-delta", type=float, default=0.05)
    parser.add_argument("--min-cells-per-donor-class", type=int, default=20)
    parser.add_argument("--limit-libraries", type=int, help="Debug: process first N final-analysis libraries.")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", errors="replace") as handle:
        return [{k: (v or "").strip() for k, v in row.items()} for row in csv.DictReader(handle, delimiter="\t")]


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})


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


def parse_characteristics(text: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for part in text.split(";"):
        part = part.strip()
        if not part or "=" not in part:
            continue
        key, value = part.split("=", 1)
        values[key.strip().lower()] = value.strip()
    return values


def normalize_group(value: str) -> str:
    lower = value.strip().lower()
    if lower in {"ds", "trisomy 21", "t21", "down syndrome"}:
        return "DS"
    if lower in {"con", "control", "disomic", "euploid"}:
        return "Control"
    return value.strip()


def library_prefix(title: str) -> str:
    match = re.search(r"Library\s+([^,\s]+)_gex", title)
    return match.group(1) if match else ""


def selected_libraries(metadata_path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for row in read_tsv(metadata_path):
        chars = parse_characteristics(row.get("characteristics", ""))
        if chars.get("library type") != "Gene Expression":
            continue
        if chars.get("in final_analysis", "").upper() == "FALSE":
            continue
        prefix = library_prefix(row.get("title", ""))
        if not prefix:
            continue
        rows.append(
            {
                "library_prefix": prefix,
                "geo_sample_accession": row.get("accession", ""),
                "biological_sample_id": chars.get("sample id", ""),
                "sample_name": chars.get("name", ""),
                "group": normalize_group(chars.get("group", chars.get("genotype", ""))),
                "genotype": chars.get("genotype", ""),
                "dev_stage_pcw": chars.get("dev stage_pcw", ""),
                "sex": chars.get("sex", ""),
                "platform_refs": row.get("platform_refs", ""),
            }
        )
    return rows


def load_target_genes(mechanism_axis_path: Path) -> tuple[list[str], dict[str, dict[str, str]]]:
    mechanism_rows = read_tsv(mechanism_axis_path)
    gene_info: dict[str, dict[str, str]] = {}
    for gene in CORE_GENES:
        gene_info[gene] = {
            "gene_symbol": gene,
            "gene_role": "core_target",
            "tier": "",
            "final_mechanism_axis_call": "",
            "mechanism_axis_short": "",
        }
    for row in mechanism_rows:
        gene = row.get("gene_symbol", "").upper()
        if not gene:
            continue
        gene_info[gene] = {
            "gene_symbol": gene,
            "gene_role": "final_candidate",
            "tier": row.get("tier", ""),
            "final_mechanism_axis_call": row.get("final_mechanism_axis_call", ""),
            "mechanism_axis_short": row.get("final_mechanism_axis_call", ""),
        }
    genes = CORE_GENES + [row.get("gene_symbol", "").upper() for row in mechanism_rows if row.get("gene_symbol", "").upper() not in CORE_GENES]
    return genes, gene_info


def parse_features(
    path: Path,
    target_genes: list[str],
) -> tuple[dict[int, list[int]], dict[int, int], int, dict[str, list[str]], dict[str, int]]:
    marker_names = list(MARKER_SETS)
    marker_index = {name: i for i, name in enumerate(marker_names)}
    marker_gene_sets = {name: {gene.upper() for gene in genes} for name, genes in MARKER_SETS.items()}
    target_index = {gene.upper(): i for i, gene in enumerate(target_genes)}
    row_to_marker: dict[int, list[int]] = defaultdict(list)
    row_to_target: dict[int, int] = {}
    present_markers: dict[str, list[str]] = {name: [] for name in marker_names}
    present_targets: dict[str, int] = {}
    last_gene_expression_row = 0

    with gzip.open(path, "rt", errors="replace") as handle:
        for idx, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            ensembl_id, symbol_raw, feature_type = fields[:3]
            if feature_type != "Gene Expression":
                continue
            last_gene_expression_row = idx
            symbol = symbol_raw.upper()
            ensembl = ensembl_id.split(".", 1)[0].upper()
            for name, genes in marker_gene_sets.items():
                if symbol in genes:
                    row_to_marker[idx].append(marker_index[name])
                    present_markers[name].append(symbol)
            if symbol in target_index:
                row_to_target[idx] = target_index[symbol]
                present_targets[symbol] = idx
            elif ensembl in target_index:
                row_to_target[idx] = target_index[ensembl]
                present_targets[ensembl] = idx

    marker_denominators = {name: max(1, len(set(values))) for name, values in present_markers.items()}
    return dict(row_to_marker), row_to_target, last_gene_expression_row, present_markers, marker_denominators


def stream_matrix(
    path: Path,
    row_to_marker: dict[int, list[int]],
    row_to_target: dict[int, int],
    last_gene_expression_row: int,
    n_marker_sets: int,
    n_target_genes: int,
) -> dict[str, np.ndarray | int]:
    total_rna: np.ndarray | None = None
    detected_genes: np.ndarray | None = None
    marker_counts: np.ndarray | None = None
    target_counts: np.ndarray | None = None
    n_features = n_cells = n_entries = 0
    dims_seen = False

    with gzip.open(path, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            fields = line.strip().split()
            if not fields:
                continue
            if not dims_seen:
                n_features, n_cells, n_entries = int(fields[0]), int(fields[1]), int(fields[2])
                total_rna = np.zeros(n_cells, dtype=np.float64)
                detected_genes = np.zeros(n_cells, dtype=np.int32)
                marker_counts = np.zeros((n_cells, n_marker_sets), dtype=np.float32)
                target_counts = np.zeros((n_cells, n_target_genes), dtype=np.float32)
                dims_seen = True
                continue
            row_index = int(fields[0])
            if row_index > last_gene_expression_row:
                continue
            cell_index = int(fields[1]) - 1
            value = float(fields[2])
            assert total_rna is not None
            assert detected_genes is not None
            total_rna[cell_index] += value
            detected_genes[cell_index] += 1
            if row_index in row_to_marker:
                assert marker_counts is not None
                for marker_idx in row_to_marker[row_index]:
                    marker_counts[cell_index, marker_idx] += value
            if row_index in row_to_target:
                assert target_counts is not None
                target_counts[cell_index, row_to_target[row_index]] += value

    if total_rna is None or detected_genes is None or marker_counts is None or target_counts is None:
        raise RuntimeError(f"No Matrix Market dimensions found in {path}")

    return {
        "n_features": n_features,
        "n_cells": n_cells,
        "n_entries": n_entries,
        "total_rna": total_rna,
        "detected_genes": detected_genes,
        "marker_counts": marker_counts,
        "target_counts": target_counts,
    }


def infer_classes(
    total_rna: np.ndarray,
    detected_genes: np.ndarray,
    marker_counts: np.ndarray,
    marker_denominators: dict[str, int],
    args: argparse.Namespace,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    marker_names = list(MARKER_SETS)
    denom = np.array([marker_denominators[name] for name in marker_names], dtype=np.float32)
    denom[denom <= 0] = 1.0
    norm = np.divide(marker_counts, denom[None, :])
    cpm = np.divide(norm * 10_000.0, total_rna[:, None], out=np.zeros_like(norm, dtype=np.float32), where=total_rna[:, None] > 0)
    scores = np.log1p(cpm)
    order = np.argsort(scores, axis=1)
    top_idx = order[:, -1]
    second_idx = order[:, -2] if scores.shape[1] > 1 else top_idx
    top_score = scores[np.arange(scores.shape[0]), top_idx]
    second_score = scores[np.arange(scores.shape[0]), second_idx]
    delta = top_score - second_score
    classes = np.array([marker_names[idx] for idx in top_idx], dtype=object)

    low_quality = (total_rna < args.min_total_rna) | (detected_genes < args.min_detected_genes)
    low_marker = (top_score < args.min_top_score) | (delta < args.min_confidence_delta)
    classes[low_quality] = "low_RNA"
    classes[~low_quality & low_marker] = "ambiguous_marker"
    confidence = np.full(classes.shape, "high", dtype=object)
    confidence[delta < 0.25] = "medium"
    confidence[delta < args.min_confidence_delta] = "low"
    confidence[low_quality] = "low_RNA"
    return classes, confidence, top_score, delta


def summarize_library(
    library: dict[str, str],
    classes: np.ndarray,
    confidence: np.ndarray,
    total_rna: np.ndarray,
    detected_genes: np.ndarray,
    target_counts: np.ndarray,
    target_genes: list[str],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    class_rows: list[dict[str, object]] = []
    pseudobulk_rows: list[dict[str, object]] = []
    for cell_class in sorted(set(classes.astype(str))):
        mask = classes.astype(str) == cell_class
        n_cells = int(mask.sum())
        if n_cells == 0:
            continue
        class_rows.append(
            {
                **library,
                "inferred_cell_class": cell_class,
                "n_cells": n_cells,
                "mean_total_rna": float(total_rna[mask].mean()),
                "median_total_rna": float(np.median(total_rna[mask])),
                "mean_detected_genes": float(detected_genes[mask].mean()),
                "high_confidence_fraction": float(np.mean(confidence[mask] == "high")),
            }
        )
        class_total_rna = float(total_rna[mask].sum())
        for gene_idx, gene in enumerate(target_genes):
            counts = target_counts[mask, gene_idx]
            target_umi = float(counts.sum())
            detected = int(np.count_nonzero(counts > 0))
            cpm = target_umi / class_total_rna * 1_000_000 if class_total_rna > 0 else math.nan
            pseudobulk_rows.append(
                {
                    **library,
                    "inferred_cell_class": cell_class,
                    "gene_symbol": gene,
                    "n_cells": n_cells,
                    "n_detected_cells": detected,
                    "detected_fraction": detected / n_cells if n_cells else math.nan,
                    "total_rna_umi": class_total_rna,
                    "target_umi": target_umi,
                    "cpm": cpm,
                    "log2_cpm": math.log2(cpm + 1.0) if math.isfinite(cpm) else math.nan,
                }
            )
    return class_rows, pseudobulk_rows


def combine_donor_pseudobulk(library_pb: pd.DataFrame, gene_info: dict[str, dict[str, str]]) -> pd.DataFrame:
    group_cols = [
        "biological_sample_id",
        "sample_name",
        "group",
        "genotype",
        "dev_stage_pcw",
        "sex",
        "inferred_cell_class",
        "gene_symbol",
    ]
    numeric = ["n_cells", "n_detected_cells", "total_rna_umi", "target_umi"]
    grouped = library_pb.groupby(group_cols, dropna=False, observed=True)[numeric].sum().reset_index()
    grouped["detected_fraction"] = grouped["n_detected_cells"] / grouped["n_cells"]
    grouped["cpm"] = np.where(grouped["total_rna_umi"] > 0, grouped["target_umi"] / grouped["total_rna_umi"] * 1_000_000, np.nan)
    grouped["log2_cpm"] = np.log2(grouped["cpm"] + 1.0)
    for column in ["gene_role", "tier", "final_mechanism_axis_call"]:
        grouped[column] = grouped["gene_symbol"].map(lambda gene: gene_info.get(str(gene), {}).get(column, ""))
    return grouped


def effect_rows(donor_pb: pd.DataFrame, min_cells: int) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    filtered = donor_pb[donor_pb["n_cells"] >= min_cells].copy()
    for (cell_class, gene), sub in filtered.groupby(["inferred_cell_class", "gene_symbol"], dropna=False, observed=True):
        ds = sub[sub["group"] == "DS"]
        control = sub[sub["group"] == "Control"]
        if ds.empty or control.empty:
            continue
        rows.append(
            {
                "inferred_cell_class": cell_class,
                "gene_symbol": gene,
                "gene_role": sub["gene_role"].iloc[0],
                "tier": sub["tier"].iloc[0],
                "final_mechanism_axis_call": sub["final_mechanism_axis_call"].iloc[0],
                "n_ds_donors": int(ds["biological_sample_id"].nunique()),
                "n_control_donors": int(control["biological_sample_id"].nunique()),
                "n_ds_cells": int(ds["n_cells"].sum()),
                "n_control_cells": int(control["n_cells"].sum()),
                "mean_ds_log2_cpm": float(ds["log2_cpm"].mean()),
                "mean_control_log2_cpm": float(control["log2_cpm"].mean()),
                "delta_ds_minus_control_log2_cpm": float(ds["log2_cpm"].mean() - control["log2_cpm"].mean()),
                "mean_ds_detected_fraction": float(ds["detected_fraction"].mean()),
                "mean_control_detected_fraction": float(control["detected_fraction"].mean()),
                "delta_detected_fraction": float(ds["detected_fraction"].mean() - control["detected_fraction"].mean()),
            }
        )
    return pd.DataFrame(rows)


def mechanism_summary(effects: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    subset = effects[effects["gene_role"] == "final_candidate"].copy()
    subset["delta"] = pd.to_numeric(subset["delta_ds_minus_control_log2_cpm"], errors="coerce")
    for (cell_class, call), sub in subset.groupby(["inferred_cell_class", "final_mechanism_axis_call"], dropna=False, observed=True):
        deltas = sub["delta"].dropna()
        rows.append(
            {
                "inferred_cell_class": cell_class,
                "final_mechanism_axis_call": call,
                "n_genes": int(sub["gene_symbol"].nunique()),
                "mean_delta_log2_cpm": float(deltas.mean()) if len(deltas) else math.nan,
                "median_delta_log2_cpm": float(deltas.median()) if len(deltas) else math.nan,
                "n_up_delta_ge_0_25": int((deltas >= 0.25).sum()),
                "n_down_delta_le_minus_0_25": int((deltas <= -0.25).sum()),
                "genes_up_delta_ge_0_25": ";".join(sub.loc[sub["delta"] >= 0.25, "gene_symbol"].astype(str).tolist()),
                "genes_down_delta_le_minus_0_25": ";".join(sub.loc[sub["delta"] <= -0.25, "gene_symbol"].astype(str).tolist()),
            }
        )
    return pd.DataFrame(rows)


def make_heatmap(matrix: pd.DataFrame, path_prefix: Path, title: str, vmax: float = 1.0) -> None:
    if matrix.empty:
        return
    path_prefix.parent.mkdir(parents=True, exist_ok=True)
    values = matrix.to_numpy(dtype=float)
    fig_w = max(8.0, 0.55 * matrix.shape[1] + 3.0)
    fig_h = max(4.5, 0.32 * matrix.shape[0] + 2.3)
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
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("DS - control delta log2(CPM + 1)")
    fig.tight_layout()
    fig.savefig(path_prefix.with_suffix(".png"), dpi=300)
    fig.savefig(path_prefix.with_suffix(".pdf"))
    plt.close(fig)


def write_note(
    path: Path,
    libraries: list[dict[str, str]],
    class_counts: pd.DataFrame,
    effects: pd.DataFrame,
    mech: pd.DataFrame,
) -> None:
    class_total = class_counts.groupby("inferred_cell_class")["n_cells"].sum().sort_values(ascending=False)
    target_effects = effects[effects["gene_symbol"].isin(CORE_GENES)].copy()
    target_effects["abs_delta"] = target_effects["delta_ds_minus_control_log2_cpm"].abs()
    strongest_targets = target_effects.sort_values("abs_delta", ascending=False).head(10)
    candidate_effects = effects[effects["gene_role"] == "final_candidate"].copy()
    candidate_effects["delta"] = candidate_effects["delta_ds_minus_control_log2_cpm"]
    class_candidate_mean = candidate_effects.groupby("inferred_cell_class")["delta"].mean().sort_values(ascending=False)
    mech_top = mech[mech["n_genes"] >= 4].sort_values("mean_delta_log2_cpm", ascending=False).head(8)

    lines = [
        "# GSE305146 Marker-Based Cell-Class First Pass",
        "",
        "Question: can the now-local GSE305146 multiome matrices be used for cell-context mechanism testing before full chromVAR?",
        "",
        "## Inputs",
        "",
        f"- Final-analysis libraries processed: {len(libraries)}.",
        "- Matrix source: `data/geo/GSE305146/supplementary/*_matrix.mtx.gz`.",
        "- Classifier: marker-score inference from the RNA block of each 10x multiome matrix.",
        "",
        "## Inferred Cell-Class Composition",
        "",
    ]
    for cell_class, n_cells in class_total.head(12).items():
        lines.append(f"- {cell_class}: {int(n_cells)} cells.")
    lines.extend(
        [
            "",
            "## Candidate Expression Pattern",
            "",
            "Gene-weighted final-candidate DS-control means were strongest in:",
        ]
    )
    for cell_class, delta in class_candidate_mean.head(8).items():
        lines.append(f"- {cell_class}: mean delta {delta:.3g}.")
    lines.extend(["", "The weakest candidate means were:"])
    for cell_class, delta in class_candidate_mean.tail(8).sort_values().items():
        lines.append(f"- {cell_class}: mean delta {delta:.3g}.")
    lines.extend(["", "## Strongest Core-Gene Class Effects", ""])
    for row in strongest_targets.itertuples(index=False):
        lines.append(
            f"- {row.gene_symbol} in {row.inferred_cell_class}: delta {row.delta_ds_minus_control_log2_cpm:.3g}; "
            f"DS donors {int(row.n_ds_donors)}, control donors {int(row.n_control_donors)}."
        )
    lines.extend(["", "## Multi-Gene Mechanism-Class Effects", ""])
    for row in mech_top.itertuples(index=False):
        lines.append(
            f"- {row.inferred_cell_class}, {row.final_mechanism_axis_call}: mean delta {row.mean_delta_log2_cpm:.3g}; "
            f"up genes {int(row.n_up_delta_ge_0_25)}, down genes {int(row.n_down_delta_le_minus_0_25)}."
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This is an inferred-label RNA analysis, not chromVAR and not a substitute for curated cell annotations. It establishes a local cell-class scaffold for GSE305146 and identifies which inferred classes carry target/candidate expression shifts. The next step is to connect the same inferred classes to ATAC peak accessibility and motif deviation or motif-set accessibility.",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    supplementary = Path(args.supplementary)
    outdir = Path(args.outdir)
    figdir = Path(args.figdir)

    target_genes, gene_info = load_target_genes(Path(args.mechanism_axis))
    libraries = selected_libraries(Path(args.metadata))
    if args.limit_libraries:
        libraries = libraries[: args.limit_libraries]

    class_rows: list[dict[str, object]] = []
    library_pb_rows: list[dict[str, object]] = []
    marker_rows: list[dict[str, object]] = []

    for index, library in enumerate(libraries, start=1):
        prefix = library["library_prefix"]
        feature_path = supplementary / f"GSE305146_{prefix}_features.tsv.gz"
        matrix_path = supplementary / f"GSE305146_{prefix}_matrix.mtx.gz"
        if not feature_path.exists() or not matrix_path.exists():
            raise FileNotFoundError(f"Missing feature or matrix file for {prefix}")
        print(f"[{index}/{len(libraries)}] {prefix}", flush=True)
        row_to_marker, row_to_target, last_gene_row, present_markers, marker_denominators = parse_features(feature_path, target_genes)
        for marker_class, markers in present_markers.items():
            marker_rows.append(
                {
                    "library_prefix": prefix,
                    "marker_class": marker_class,
                    "n_markers_present": len(set(markers)),
                    "markers_present": ";".join(sorted(set(markers))),
                }
            )
        stats = stream_matrix(
            matrix_path,
            row_to_marker,
            row_to_target,
            last_gene_row,
            n_marker_sets=len(MARKER_SETS),
            n_target_genes=len(target_genes),
        )
        classes, confidence, _top_score, _delta = infer_classes(
            stats["total_rna"],
            stats["detected_genes"],
            stats["marker_counts"],
            marker_denominators,
            args,
        )
        lib_class_rows, lib_pb_rows = summarize_library(
            library,
            classes,
            confidence,
            stats["total_rna"],
            stats["detected_genes"],
            stats["target_counts"],
            target_genes,
        )
        class_rows.extend(lib_class_rows)
        library_pb_rows.extend(lib_pb_rows)

    class_fields = [
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
        "n_cells",
        "mean_total_rna",
        "median_total_rna",
        "mean_detected_genes",
        "high_confidence_fraction",
    ]
    pb_fields = [
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
        "gene_symbol",
        "n_cells",
        "n_detected_cells",
        "detected_fraction",
        "total_rna_umi",
        "target_umi",
        "cpm",
        "log2_cpm",
    ]
    marker_fields = ["library_prefix", "marker_class", "n_markers_present", "markers_present"]
    write_tsv(outdir / "gse305146_marker_cell_class_counts.tsv", class_rows, class_fields)
    write_tsv(outdir / "gse305146_marker_library_pseudobulk.tsv", library_pb_rows, pb_fields)
    write_tsv(outdir / "gse305146_marker_set_recovery.tsv", marker_rows, marker_fields)

    library_pb = pd.DataFrame(library_pb_rows)
    donor_pb = combine_donor_pseudobulk(library_pb, gene_info)
    donor_path = outdir / "gse305146_marker_donor_pseudobulk.tsv"
    donor_pb.to_csv(donor_path, sep="\t", index=False, lineterminator="\n")

    effects = effect_rows(donor_pb, args.min_cells_per_donor_class)
    effects_path = outdir / "gse305146_marker_cell_class_effects.tsv"
    effects.to_csv(effects_path, sep="\t", index=False, lineterminator="\n")

    mech = mechanism_summary(effects)
    mech_path = outdir / "gse305146_marker_mechanism_axis_cell_class_summary.tsv"
    mech.to_csv(mech_path, sep="\t", index=False, lineterminator="\n")

    class_order = (
        pd.DataFrame(class_rows)
        .groupby("inferred_cell_class")["n_cells"]
        .sum()
        .sort_values(ascending=False)
        .index.tolist()
    )
    core = effects[effects["gene_symbol"].isin(CORE_GENES)].copy()
    core_matrix = core.pivot_table(
        index="gene_symbol",
        columns="inferred_cell_class",
        values="delta_ds_minus_control_log2_cpm",
        aggfunc="first",
    ).reindex(index=CORE_GENES, columns=class_order)
    make_heatmap(core_matrix, figdir / "gse305146_marker_core_gene_cell_class_heatmap", "GSE305146 core genes by inferred cell class", vmax=1.0)

    candidate_order = [gene for gene in target_genes if gene not in CORE_GENES]
    cand = effects[effects["gene_role"] == "final_candidate"].copy()
    cand_matrix = cand.pivot_table(
        index="gene_symbol",
        columns="inferred_cell_class",
        values="delta_ds_minus_control_log2_cpm",
        aggfunc="first",
    ).reindex(index=candidate_order, columns=class_order)
    make_heatmap(cand_matrix, figdir / "gse305146_marker_final_candidate_cell_class_heatmap", "GSE305146 final candidates by inferred cell class", vmax=1.2)

    mech_matrix = mech.pivot_table(
        index="final_mechanism_axis_call",
        columns="inferred_cell_class",
        values="mean_delta_log2_cpm",
        aggfunc="first",
    ).reindex(columns=class_order)
    make_heatmap(mech_matrix, figdir / "gse305146_marker_mechanism_axis_cell_class_heatmap", "GSE305146 mechanism classes by inferred cell class", vmax=0.8)

    write_note(Path(args.note), libraries, pd.DataFrame(class_rows), effects, mech)
    print(f"Wrote inferred cell-class summaries for {len(libraries)} libraries")
    print(effects_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
