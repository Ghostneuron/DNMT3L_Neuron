#!/usr/bin/env python3
"""Inferred-cell-class target-window ATAC accessibility for GSE305146."""

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
    CORE_GENES,
    MARKER_SETS,
    infer_classes,
    load_target_genes,
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
MECHANISM_LABELS = {
    "cross_species_dnmt3a_linked_methylation_plus_developmental_cpg_tf_context": "DNMT3A methyl. + dev/CpG",
    "cross_species_dnmt3a_linked_methylation_plus_developmental_tf_context": "DNMT3A methyl. + dev TF",
    "cross_species_dnmt3a_linked_methylation_without_selected_tf_instance": "DNMT3A methyl. only",
    "mouse_dnmt3a_linked_methylation_with_tf_context_human_unresolved": "Mouse methyl. + TF; human unresolved",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", default="data/geo/GSE305146/metadata/sample_metadata.tsv")
    parser.add_argument("--supplementary", default="data/geo/GSE305146/supplementary")
    parser.add_argument("--mechanism-axis", default="results/mechanism_question/final_candidate_mechanism_axis_summary.tsv")
    parser.add_argument("--outdir", default="results/mechanism_question")
    parser.add_argument("--figdir", default="figures/mechanism_question")
    parser.add_argument("--note", default="notes/gse305146_marker_atac_target_window_first_pass.md")
    parser.add_argument("--window-bp", type=int, default=100_000)
    parser.add_argument("--min-total-rna", type=float, default=500.0)
    parser.add_argument("--min-detected-genes", type=int, default=200)
    parser.add_argument("--min-top-score", type=float, default=0.05)
    parser.add_argument("--min-confidence-delta", type=float, default=0.05)
    parser.add_argument("--min-cells-per-donor-class", type=int, default=20)
    parser.add_argument("--limit-libraries", type=int, help="Debug: process first N final-analysis libraries.")
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


def safe_int(value: str) -> int:
    try:
        return int(float(value))
    except ValueError:
        return 0


def parse_window_features(
    feature_path: Path,
    target_genes: list[str],
    window_bp: int,
) -> tuple[dict[int, list[str]], dict[str, dict[str, object]], int]:
    target_set = {gene.upper() for gene in target_genes}
    gene_windows: dict[str, dict[str, object]] = {}
    peak_row_to_genes: dict[int, list[str]] = {}
    last_gene_expression_row = 0

    with gzip.open(feature_path, "rt", errors="replace") as handle:
        for row_index, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            feature_id, symbol_raw, feature_type, chrom, start_raw, end_raw = fields[:6]
            start = safe_int(start_raw)
            end = safe_int(end_raw)
            symbol = symbol_raw.upper()
            if feature_type == "Gene Expression":
                last_gene_expression_row = row_index
                if symbol not in target_set:
                    continue
                if not chrom.startswith("chr") or end <= start:
                    continue
                existing = gene_windows.get(symbol)
                if existing:
                    existing["gene_start"] = min(int(existing["gene_start"]), start)
                    existing["gene_end"] = max(int(existing["gene_end"]), end)
                    existing["window_start"] = max(0, int(existing["gene_start"]) - window_bp)
                    existing["window_end"] = int(existing["gene_end"]) + window_bp
                else:
                    gene_windows[symbol] = {
                        "gene_symbol": symbol,
                        "ensembl_id": feature_id,
                        "chrom": chrom,
                        "gene_start": start,
                        "gene_end": end,
                        "window_start": max(0, start - window_bp),
                        "window_end": end + window_bp,
                    }
                continue
            if feature_type != "Peaks" or end <= start:
                continue
            overlapping: list[str] = []
            for gene, window in gene_windows.items():
                if chrom != window["chrom"]:
                    continue
                if end >= int(window["window_start"]) and start <= int(window["window_end"]):
                    overlapping.append(gene)
            if overlapping:
                peak_row_to_genes[row_index] = overlapping

    for gene, window in gene_windows.items():
        window["n_window_peak_features"] = sum(1 for genes in peak_row_to_genes.values() if gene in genes)
    return peak_row_to_genes, gene_windows, last_gene_expression_row


def stream_atac_by_class(
    matrix_path: Path,
    classes: np.ndarray,
    last_gene_expression_row: int,
    peak_row_to_gene_indices: dict[int, list[int]],
    n_genes: int,
) -> tuple[np.ndarray, np.ndarray]:
    n_cells = classes.shape[0]
    total_peak_counts = np.zeros(n_cells, dtype=np.float64)
    target_peak_counts = np.zeros((n_cells, n_genes), dtype=np.float32)
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
            for gene_index in peak_row_to_gene_indices.get(row_index, []):
                target_peak_counts[cell_index, gene_index] += value

    return total_peak_counts, target_peak_counts


def summarize_library_atac(
    library: dict[str, str],
    classes: np.ndarray,
    total_peak_counts: np.ndarray,
    target_peak_counts: np.ndarray,
    target_genes: list[str],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    class_values = sorted(set(classes.astype(str)))
    for cell_class in class_values:
        mask = classes.astype(str) == cell_class
        n_cells = int(mask.sum())
        if n_cells == 0:
            continue
        class_total_peak = float(total_peak_counts[mask].sum())
        for gene_idx, gene in enumerate(target_genes):
            values = target_peak_counts[mask, gene_idx]
            target_count = float(values.sum())
            detected = int(np.count_nonzero(values > 0))
            peak_cpm = target_count / class_total_peak * 1_000_000 if class_total_peak > 0 else math.nan
            rows.append(
                {
                    **library,
                    "inferred_cell_class": cell_class,
                    "gene_symbol": gene,
                    "n_cells": n_cells,
                    "n_cells_with_window_peak": detected,
                    "detected_fraction": detected / n_cells if n_cells else math.nan,
                    "total_peak_counts": class_total_peak,
                    "target_window_peak_counts": target_count,
                    "peak_cpm": peak_cpm,
                    "log2_peak_cpm": math.log2(peak_cpm + 1.0) if math.isfinite(peak_cpm) else math.nan,
                }
            )
    return rows


def combine_donor_pseudobulk(library_rows: list[dict[str, object]], gene_info: dict[str, dict[str, str]]) -> pd.DataFrame:
    frame = pd.DataFrame(library_rows)
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
    numeric = ["n_cells", "n_cells_with_window_peak", "total_peak_counts", "target_window_peak_counts"]
    donor = frame.groupby(group_cols, dropna=False, observed=True)[numeric].sum().reset_index()
    donor["detected_fraction"] = donor["n_cells_with_window_peak"] / donor["n_cells"]
    donor["peak_cpm"] = np.where(
        donor["total_peak_counts"] > 0,
        donor["target_window_peak_counts"] / donor["total_peak_counts"] * 1_000_000,
        np.nan,
    )
    donor["log2_peak_cpm"] = np.log2(donor["peak_cpm"] + 1.0)
    for column in ["gene_role", "tier", "final_mechanism_axis_call"]:
        donor[column] = donor["gene_symbol"].map(lambda gene: gene_info.get(str(gene), {}).get(column, ""))
    return donor


def effect_rows(donor: pd.DataFrame, min_cells: int) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    filtered = donor[donor["n_cells"] >= min_cells].copy()
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


def mechanism_summary(effects: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    subset = effects[effects["gene_role"] == "final_candidate"].copy()
    subset["delta"] = pd.to_numeric(subset["delta_ds_minus_control_log2_peak_cpm"], errors="coerce")
    for (cell_class, call), sub in subset.groupby(["inferred_cell_class", "final_mechanism_axis_call"], dropna=False, observed=True):
        deltas = sub["delta"].dropna()
        rows.append(
            {
                "inferred_cell_class": cell_class,
                "final_mechanism_axis_call": call,
                "n_genes": int(sub["gene_symbol"].nunique()),
                "mean_delta_log2_peak_cpm": float(deltas.mean()) if len(deltas) else math.nan,
                "median_delta_log2_peak_cpm": float(deltas.median()) if len(deltas) else math.nan,
                "n_up_delta_ge_0_25": int((deltas >= 0.25).sum()),
                "n_down_delta_le_minus_0_25": int((deltas <= -0.25).sum()),
                "genes_up_delta_ge_0_25": ";".join(sub.loc[sub["delta"] >= 0.25, "gene_symbol"].astype(str).tolist()),
                "genes_down_delta_le_minus_0_25": ";".join(sub.loc[sub["delta"] <= -0.25, "gene_symbol"].astype(str).tolist()),
            }
        )
    return pd.DataFrame(rows)


def make_heatmap(matrix: pd.DataFrame, path_prefix: Path, title: str, vmax: float) -> None:
    matrix = matrix.dropna(how="all").dropna(axis=1, how="all")
    if matrix.empty:
        return
    path_prefix.parent.mkdir(parents=True, exist_ok=True)
    values = matrix.to_numpy(dtype=float)
    fig_w = max(7.5, 0.58 * matrix.shape[1] + 3.0)
    fig_h = max(4.2, 0.34 * matrix.shape[0] + 2.2)
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
    cbar = fig.colorbar(im, ax=ax, fraction=0.032, pad=0.02)
    cbar.set_label("DS - control delta log2(peak CPM + 1)")
    fig.tight_layout()
    fig.savefig(path_prefix.with_suffix(".png"), dpi=300)
    fig.savefig(path_prefix.with_suffix(".pdf"))
    plt.close(fig)


def write_note(
    path: Path,
    libraries: list[dict[str, str]],
    effects: pd.DataFrame,
    mech: pd.DataFrame,
    gene_windows: pd.DataFrame,
) -> None:
    core = effects[effects["gene_symbol"].isin(CORE_GENES)].copy()
    core["abs_delta"] = core["delta_ds_minus_control_log2_peak_cpm"].abs()
    strongest_core = core.sort_values("abs_delta", ascending=False).head(10)
    dnmt3l = core[core["gene_symbol"] == "DNMT3L"].sort_values("delta_ds_minus_control_log2_peak_cpm", ascending=False)
    candidate = effects[effects["gene_role"] == "final_candidate"].copy()
    candidate_mean = candidate.groupby("inferred_cell_class")["delta_ds_minus_control_log2_peak_cpm"].mean().sort_values(ascending=False)
    mech_top = mech[mech["n_genes"] >= 4].sort_values("mean_delta_log2_peak_cpm", ascending=False).head(8)
    missing_windows = sorted(set(effects["gene_symbol"]) - set(gene_windows["gene_symbol"])) if not gene_windows.empty else []

    lines = [
        "# GSE305146 Inferred-Class ATAC Target-Window First Pass",
        "",
        "Question: do inferred fetal cortical cell classes show DS-control accessibility shifts around DNMT3L/STAT genes and final candidates?",
        "",
        "## Inputs",
        "",
        f"- Final-analysis libraries processed: {len(libraries)}.",
        "- Inferred classes: RNA marker-score labels generated per library.",
        "- ATAC measurement: peak counts in gene body plus 100 kb flanks, normalized to total peak counts within donor/class.",
        "",
        "## Outputs",
        "",
        "- Library pseudobulk: `results/mechanism_question/gse305146_marker_atac_window_library_pseudobulk.tsv`.",
        "- Donor pseudobulk: `results/mechanism_question/gse305146_marker_atac_window_donor_pseudobulk.tsv`.",
        "- Effects: `results/mechanism_question/gse305146_marker_atac_window_effects.tsv`.",
        "- Mechanism summary: `results/mechanism_question/gse305146_marker_atac_window_mechanism_axis_summary.tsv`.",
        "- Curated heatmaps: `figures/mechanism_question/gse305146_marker_atac_window_*_curated.png`.",
        "",
        "## Target-Window Coverage",
        "",
        f"- Genes with recovered target windows: {gene_windows['gene_symbol'].nunique() if not gene_windows.empty else 0}.",
    ]
    if missing_windows:
        lines.append(f"- Missing windows: {';'.join(missing_windows)}.")
    else:
        lines.append("- Missing windows: none among tested genes.")

    lines.extend(["", "## DNMT3L Accessibility", ""])
    for row in dnmt3l.head(6).itertuples(index=False):
        lines.append(
            f"- {row.inferred_cell_class}: delta {row.delta_ds_minus_control_log2_peak_cpm:.3g}; "
            f"DS donors {int(row.n_ds_donors)}, control donors {int(row.n_control_donors)}."
        )

    lines.extend(["", "## Strongest Core-Gene Accessibility Effects", ""])
    for row in strongest_core.itertuples(index=False):
        lines.append(
            f"- {row.gene_symbol} in {row.inferred_cell_class}: delta {row.delta_ds_minus_control_log2_peak_cpm:.3g}; "
            f"DS donors {int(row.n_ds_donors)}, control donors {int(row.n_control_donors)}."
        )

    lines.extend(["", "## Candidate Accessibility Pattern", ""])
    lines.append("Gene-weighted final-candidate target-window accessibility means were strongest in:")
    for cell_class, delta in candidate_mean.head(8).items():
        lines.append(f"- {cell_class}: mean delta {delta:.3g}.")
    lines.extend(["", "The weakest candidate target-window means were:"])
    for cell_class, delta in candidate_mean.tail(8).sort_values().items():
        lines.append(f"- {cell_class}: mean delta {delta:.3g}.")

    lines.extend(["", "## Multi-Gene Mechanism-Class Accessibility", ""])
    for row in mech_top.itertuples(index=False):
        lines.append(
            f"- {row.inferred_cell_class}, {row.final_mechanism_axis_call}: mean delta {row.mean_delta_log2_peak_cpm:.3g}; "
            f"up genes {int(row.n_up_delta_ge_0_25)}, down genes {int(row.n_down_delta_le_minus_0_25)}."
        )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This is a target-window accessibility proxy, not chromVAR and not base-resolution methylation. It tests whether local candidate/core gene windows show DS-control accessibility shifts within inferred classes. The result should be interpreted as locus-level chromatin context that can guide deeper motif/accessibility analyses.",
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
    gene_index = {gene: idx for idx, gene in enumerate(target_genes)}

    libraries = selected_libraries(Path(args.metadata))
    if args.limit_libraries:
        libraries = libraries[: args.limit_libraries]

    library_rows: list[dict[str, object]] = []
    gene_window_rows: list[dict[str, object]] = []

    for index, library in enumerate(libraries, start=1):
        prefix = library["library_prefix"]
        feature_path = supplementary / f"GSE305146_{prefix}_features.tsv.gz"
        matrix_path = supplementary / f"GSE305146_{prefix}_matrix.mtx.gz"
        if not feature_path.exists() or not matrix_path.exists():
            raise FileNotFoundError(f"Missing feature or matrix for {prefix}")

        print(f"[{index}/{len(libraries)}] {prefix}", flush=True)
        row_to_marker, row_to_target, last_gene_row, _present_markers, marker_denominators = parse_marker_features(
            feature_path, target_genes
        )
        peak_row_to_genes, gene_windows, peak_last_gene_row = parse_window_features(feature_path, target_genes, args.window_bp)
        last_gene_row = max(last_gene_row, peak_last_gene_row)
        peak_row_to_gene_indices = {
            row: [gene_index[gene] for gene in genes if gene in gene_index] for row, genes in peak_row_to_genes.items()
        }
        for row in gene_windows.values():
            gene_window_rows.append({"library_prefix": prefix, **row})

        stats = stream_matrix(
            matrix_path,
            row_to_marker,
            row_to_target,
            last_gene_row,
            n_marker_sets=len(MARKER_SETS),
            n_target_genes=len(target_genes),
        )
        classes, _confidence, _top_score, _delta = infer_classes(
            stats["total_rna"],
            stats["detected_genes"],
            stats["marker_counts"],
            marker_denominators,
            args,
        )
        total_peak_counts, target_peak_counts = stream_atac_by_class(
            matrix_path,
            classes,
            last_gene_row,
            peak_row_to_gene_indices,
            len(target_genes),
        )
        library_rows.extend(summarize_library_atac(library, classes, total_peak_counts, target_peak_counts, target_genes))

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
        "gene_symbol",
        "n_cells",
        "n_cells_with_window_peak",
        "detected_fraction",
        "total_peak_counts",
        "target_window_peak_counts",
        "peak_cpm",
        "log2_peak_cpm",
    ]
    write_tsv(outdir / "gse305146_marker_atac_window_library_pseudobulk.tsv", library_rows, library_fields)

    gene_window_df = pd.DataFrame(gene_window_rows).drop_duplicates()
    gene_window_df.to_csv(outdir / "gse305146_marker_atac_window_gene_windows.tsv", sep="\t", index=False, lineterminator="\n")

    donor = combine_donor_pseudobulk(library_rows, gene_info)
    donor.to_csv(outdir / "gse305146_marker_atac_window_donor_pseudobulk.tsv", sep="\t", index=False, lineterminator="\n")
    effects = effect_rows(donor, args.min_cells_per_donor_class)
    effects.to_csv(outdir / "gse305146_marker_atac_window_effects.tsv", sep="\t", index=False, lineterminator="\n")
    mech = mechanism_summary(effects)
    mech.to_csv(outdir / "gse305146_marker_atac_window_mechanism_axis_summary.tsv", sep="\t", index=False, lineterminator="\n")

    axis = pd.read_csv(args.mechanism_axis, sep="\t")
    class_sizes = axis.groupby("final_mechanism_axis_call")["gene_symbol"].nunique().to_dict()
    class_order = [cell_class for cell_class in CLASS_ORDER if cell_class in set(effects["inferred_cell_class"])]

    core = effects[effects["gene_symbol"].isin(CORE_GENES)].copy()
    core_matrix = core.pivot_table(
        index="gene_symbol",
        columns="inferred_cell_class",
        values="delta_ds_minus_control_log2_peak_cpm",
        aggfunc="first",
    ).reindex(index=CORE_GENES, columns=class_order)
    core_matrix.columns = [CLASS_LABELS.get(col, col) for col in core_matrix.columns]
    make_heatmap(
        core_matrix,
        figdir / "gse305146_marker_atac_window_core_gene_heatmap_curated",
        "GSE305146 target-window accessibility: core genes",
        vmax=1.0,
    )

    candidate_order = axis["gene_symbol"].tolist()
    candidate = effects[effects["gene_role"] == "final_candidate"].copy()
    candidate_matrix = candidate.pivot_table(
        index="gene_symbol",
        columns="inferred_cell_class",
        values="delta_ds_minus_control_log2_peak_cpm",
        aggfunc="first",
    ).reindex(index=candidate_order, columns=class_order)
    candidate_matrix.columns = [CLASS_LABELS.get(col, col) for col in candidate_matrix.columns]
    make_heatmap(
        candidate_matrix,
        figdir / "gse305146_marker_atac_window_candidate_heatmap_curated",
        "GSE305146 target-window accessibility: final candidates",
        vmax=1.2,
    )

    mech_matrix = mech.pivot_table(
        index="final_mechanism_axis_call",
        columns="inferred_cell_class",
        values="mean_delta_log2_peak_cpm",
        aggfunc="first",
    ).reindex(columns=class_order)
    mech_matrix.index = [
        f"{MECHANISM_LABELS.get(idx, idx)} (n={class_sizes.get(idx, 0)})" for idx in mech_matrix.index
    ]
    mech_matrix.columns = [CLASS_LABELS.get(col, col) for col in mech_matrix.columns]
    make_heatmap(
        mech_matrix,
        figdir / "gse305146_marker_atac_window_mechanism_axis_heatmap_curated",
        "GSE305146 target-window accessibility: mechanism classes",
        vmax=0.8,
    )

    write_note(Path(args.note), libraries, effects, mech, gene_window_df)
    print(f"Wrote inferred-class ATAC window summaries for {len(libraries)} libraries")
    print(outdir / "gse305146_marker_atac_window_effects.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
