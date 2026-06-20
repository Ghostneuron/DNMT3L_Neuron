#!/usr/bin/env python3
"""Cell-type-resolved target-gene analysis for the GSE305153 CELLxGENE atlas."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


TARGET_ORDER = ["DNMT3L", "DNMT3A", "DNMT3B", "STAT1", "STAT3"]
DEFAULT_H5AD = "data/cellxgene/GSE305153_complete_fetal_cortex_cellxgene.h5ad"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", default=DEFAULT_H5AD)
    parser.add_argument("--target-genes", default="config/target_genes.tsv")
    parser.add_argument("--outdir", default="results/gse305153_cellxgene_targets")
    parser.add_argument("--figdir", default="figures/gse305153_cellxgene_targets")
    parser.add_argument("--matrix-path", default="raw/X", help="Sparse CSR matrix path inside H5AD.")
    parser.add_argument("--chunk-cells", type=int, default=1500)
    parser.add_argument("--min-cells-per-donor-label", type=int, default=20)
    parser.add_argument("--force-extract", action="store_true")
    return parser.parse_args()


def fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and not math.isfinite(value):
        return ""
    if isinstance(value, float):
        return f"{value:.6g}"
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


def decode_array(values) -> np.ndarray:
    arr = np.asarray(values)
    if arr.dtype.kind in {"S", "O"}:
        return np.array([value.decode() if isinstance(value, bytes) else str(value) for value in arr])
    return arr


def read_h5_column(parent_group, key: str) -> np.ndarray:
    obj = parent_group[key]
    if isinstance(obj, h5py.Group) and {"codes", "categories"}.issubset(obj.keys()):
        codes = np.asarray(obj["codes"])
        categories = decode_array(obj["categories"][:])
        out = np.empty(codes.shape[0], dtype=object)
        valid = codes >= 0
        out[valid] = categories[codes[valid]]
        out[~valid] = ""
        return out
    return decode_array(obj[:])


def read_human_targets(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    by_ensembl: dict[str, str] = {}
    by_symbol: dict[str, str] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["species"] != "human":
                continue
            symbol = row["gene_symbol"]
            ensembl = row.get("ensembl_hint", "").split(".", 1)[0]
            if ensembl:
                by_ensembl[ensembl] = symbol
            by_symbol[symbol.upper()] = symbol
    return by_ensembl, by_symbol


def read_var_targets(f, target_by_ensembl: dict[str, str], target_by_symbol: dict[str, str]) -> list[dict[str, object]]:
    var_group = f["raw/var"] if "raw" in f and "var" in f["raw"] else f["var"]
    ensembl_ids = read_h5_column(var_group, "_index")
    gene_names = read_h5_column(var_group, "gene_name" if "gene_name" in var_group else "feature_name")
    rows: list[dict[str, object]] = []
    seen: set[str] = set()
    for idx, (ensembl, gene_name) in enumerate(zip(ensembl_ids, gene_names)):
        ensembl_clean = str(ensembl).split(".", 1)[0]
        symbol = target_by_ensembl.get(ensembl_clean) or target_by_symbol.get(str(gene_name).upper())
        if symbol and symbol not in seen:
            rows.append(
                {
                    "gene_symbol": symbol,
                    "feature_index": idx,
                    "ensembl_id": ensembl_clean,
                    "gene_name": gene_name,
                }
            )
            seen.add(symbol)
    return sorted(rows, key=lambda row: TARGET_ORDER.index(row["gene_symbol"]) if row["gene_symbol"] in TARGET_ORDER else 99)


def extract_target_counts(
    h5ad_path: Path,
    matrix_path: str,
    target_rows: list[dict[str, object]],
    cache_path: Path,
    chunk_cells: int,
    force: bool,
) -> np.ndarray:
    if cache_path.exists() and not force:
        with np.load(cache_path, allow_pickle=False) as cached:
            return cached["target_counts"]

    target_indices = [int(row["feature_index"]) for row in target_rows]
    with h5py.File(h5ad_path, "r") as f:
        x = f[matrix_path]
        n_cells, _n_features = x.attrs["shape"]
        indptr = x["indptr"][:]
        target_counts = np.zeros((n_cells, len(target_indices)), dtype=np.float32)
        for start in range(0, n_cells, chunk_cells):
            end = min(start + chunk_cells, n_cells)
            data_start = int(indptr[start])
            data_end = int(indptr[end])
            if data_start == data_end:
                continue
            indices = x["indices"][data_start:data_end]
            data = x["data"][data_start:data_end]
            row_offsets = np.diff(indptr[start : end + 1])
            row_ids = np.repeat(np.arange(end - start, dtype=np.int32), row_offsets)
            for gene_pos, feature_index in enumerate(target_indices):
                mask = indices == feature_index
                if np.any(mask):
                    np.add.at(target_counts[start:end, gene_pos], row_ids[mask], data[mask])

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        cache_path,
        target_counts=target_counts,
        gene_symbols=np.array([row["gene_symbol"] for row in target_rows]),
        feature_indices=np.array(target_indices),
    )
    return target_counts


def stage_bin(pcw: object) -> str:
    try:
        value = int(float(pcw))
    except (ValueError, TypeError):
        return "unknown"
    if 10 <= value <= 13:
        return "early_PCW10_13"
    if 16 <= value <= 20:
        return "late_PCW16_20"
    return f"PCW{value}"


def obs_dataframe(f) -> pd.DataFrame:
    columns = [
        "cell_id",
        "group",
        "disease",
        "donor_id",
        "sample_name",
        "dev_PCW",
        "sex",
        "cell_type",
        "author_cell_type",
        "cell_class",
        "cluster_name",
        "cluster_name_subset",
        "cells_in_excitatory_lineage_subset",
        "nCount_RNA",
        "nFeature_RNA",
    ]
    obs = f["obs"]
    data = {key: read_h5_column(obs, key) for key in columns if key in obs}
    df = pd.DataFrame(data)
    df["stage_bin"] = df["dev_PCW"].map(stage_bin)
    df["group"] = df["group"].replace({"CON": "Control"})
    for key in ["cluster_name_subset"]:
        if key in df:
            df[key] = df[key].replace({"": "not_in_excitatory_subset"})
    df["nCount_RNA"] = pd.to_numeric(df["nCount_RNA"], errors="coerce").fillna(0.0)
    return df


def build_cell_values(df: pd.DataFrame, target_counts: np.ndarray, genes: list[str]) -> pd.DataFrame:
    value_df = df[
        [
            "cell_id",
            "group",
            "donor_id",
            "sample_name",
            "dev_PCW",
            "stage_bin",
            "sex",
            "cell_type",
            "author_cell_type",
            "cell_class",
            "cluster_name",
            "cluster_name_subset",
            "nCount_RNA",
        ]
    ].copy()
    for idx, gene in enumerate(genes):
        value_df[gene] = target_counts[:, idx]
    return value_df


def pseudobulk_rows(cell_values: pd.DataFrame, genes: list[str], label_types: list[str]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    base_cols = ["group", "donor_id", "sample_name", "dev_PCW", "stage_bin", "sex"]
    for label_type in label_types:
        group_cols = [label_type, *base_cols]
        grouped = cell_values.groupby(group_cols, dropna=False, observed=True)
        for group_key, sub in grouped:
            key_values = dict(zip(group_cols, group_key))
            total_rna = float(sub["nCount_RNA"].sum())
            n_cells = int(sub.shape[0])
            for gene in genes:
                target_umi = float(sub[gene].sum())
                detected = int((sub[gene] > 0).sum())
                cpm = target_umi / total_rna * 1_000_000 if total_rna > 0 else math.nan
                rows.append(
                    {
                        "label_type": label_type,
                        "label": key_values[label_type],
                        "gene_symbol": gene,
                        "group": key_values["group"],
                        "donor_id": key_values["donor_id"],
                        "sample_name": key_values["sample_name"],
                        "dev_PCW": key_values["dev_PCW"],
                        "stage_bin": key_values["stage_bin"],
                        "sex": key_values["sex"],
                        "n_cells": n_cells,
                        "n_detected_cells": detected,
                        "detected_fraction": detected / n_cells if n_cells else math.nan,
                        "total_rna_umi": total_rna,
                        "target_umi": target_umi,
                        "cpm": cpm,
                        "log2_cpm": math.log2(cpm + 1.0) if math.isfinite(cpm) else math.nan,
                    }
                )
    return rows


def cell_group_rows(cell_values: pd.DataFrame, genes: list[str], label_types: list[str]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for label_type in label_types:
        grouped = cell_values.groupby([label_type, "group"], dropna=False, observed=True)
        for (label, group), sub in grouped:
            total_rna = float(sub["nCount_RNA"].sum())
            n_cells = int(sub.shape[0])
            for gene in genes:
                target_umi = float(sub[gene].sum())
                detected = int((sub[gene] > 0).sum())
                cpm = target_umi / total_rna * 1_000_000 if total_rna > 0 else math.nan
                rows.append(
                    {
                        "label_type": label_type,
                        "label": label,
                        "group": group,
                        "gene_symbol": gene,
                        "n_cells": n_cells,
                        "n_detected_cells": detected,
                        "detected_fraction": detected / n_cells if n_cells else math.nan,
                        "total_rna_umi": total_rna,
                        "target_umi": target_umi,
                        "cpm": cpm,
                        "log2_cpm": math.log2(cpm + 1.0) if math.isfinite(cpm) else math.nan,
                        "mean_count_per_cell": target_umi / n_cells if n_cells else math.nan,
                    }
                )
    return rows


def effect_rows(pseudobulk: pd.DataFrame, min_cells: int) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    pb = pseudobulk[pseudobulk["n_cells"] >= min_cells].copy()
    contexts: list[tuple[str, pd.Series]] = [("all_pcw", pd.Series(True, index=pb.index))]
    for value in ["early_PCW10_13", "late_PCW16_20"]:
        contexts.append((value, pb["stage_bin"] == value))
    for pcw in sorted(pb["dev_PCW"].dropna().unique(), key=lambda x: int(float(x)) if str(x).replace(".", "", 1).isdigit() else 999):
        contexts.append((f"PCW{pcw}", pb["dev_PCW"] == pcw))

    group_cols = ["label_type", "label", "gene_symbol"]
    for context, mask in contexts:
        sub = pb[mask]
        if sub.empty:
            continue
        for key, group_df in sub.groupby(group_cols, dropna=False, observed=True):
            label_type, label, gene = key
            ds = group_df[group_df["group"] == "DS"]
            control = group_df[group_df["group"] == "Control"]
            if ds.empty or control.empty:
                continue
            rows.append(
                {
                    "label_type": label_type,
                    "label": label,
                    "context": context,
                    "gene_symbol": gene,
                    "n_ds_donors": ds["donor_id"].nunique(),
                    "n_control_donors": control["donor_id"].nunique(),
                    "n_ds_cells": int(ds["n_cells"].sum()),
                    "n_control_cells": int(control["n_cells"].sum()),
                    "mean_ds_log2_cpm": float(ds["log2_cpm"].mean()),
                    "mean_control_log2_cpm": float(control["log2_cpm"].mean()),
                    "delta_ds_minus_control_log2_cpm": float(ds["log2_cpm"].mean() - control["log2_cpm"].mean()),
                    "mean_ds_detected_fraction": float(ds["detected_fraction"].mean()),
                    "mean_control_detected_fraction": float(control["detected_fraction"].mean()),
                    "delta_ds_minus_control_detected_fraction": float(ds["detected_fraction"].mean() - control["detected_fraction"].mean()),
                }
            )
    return rows


def plot_outputs(effect_df: pd.DataFrame, cell_group_df: pd.DataFrame, figdir: Path) -> list[Path]:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return []

    figdir.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []

    author_effect = effect_df[
        (effect_df["label_type"] == "author_cell_type")
        & (effect_df["context"] == "all_pcw")
        & (effect_df["gene_symbol"].isin(TARGET_ORDER))
    ].copy()
    if not author_effect.empty:
        labels = list(author_effect["label"].drop_duplicates())
        matrix = []
        for label in labels:
            row = []
            for gene in TARGET_ORDER:
                vals = author_effect[(author_effect["label"] == label) & (author_effect["gene_symbol"] == gene)]
                row.append(float(vals["delta_ds_minus_control_log2_cpm"].iloc[0]) if not vals.empty else math.nan)
            matrix.append(row)
        finite = [value for row in matrix for value in row if math.isfinite(value)]
        vmax = max(0.25, min(3.0, max(abs(value) for value in finite))) if finite else 1.0
        fig, ax = plt.subplots(figsize=(9.2, max(4.8, 0.36 * len(labels) + 1.9)))
        image = ax.imshow(matrix, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax.set_xticks(range(len(TARGET_ORDER)))
        ax.set_xticklabels(TARGET_ORDER, rotation=35, ha="right")
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_title("CELLxGENE GSE305153 DS-control target-gene effects")
        ax.set_xlabel("Target gene")
        ax.set_ylabel("Author cell type")
        for y, row in enumerate(matrix):
            for x, value in enumerate(row):
                if math.isfinite(value):
                    text_color = "white" if abs(value) >= 0.75 * vmax else "black"
                    ax.text(x, y, f"{value:+.2f}", ha="center", va="center", fontsize=9, color=text_color)
        cbar = fig.colorbar(image, ax=ax, shrink=0.75)
        cbar.set_label("Mean donor log2(CPM+1), DS - control")
        fig.tight_layout()
        out = figdir / "target_gene_author_cell_type_effect_heatmap.png"
        fig.savefig(out, dpi=220)
        fig.savefig(out.with_suffix(".pdf"))
        plt.close(fig)
        outputs.extend([out, out.with_suffix(".pdf")])

    dnmt3l = effect_df[
        (effect_df["label_type"] == "cluster_name")
        & (effect_df["context"] == "all_pcw")
        & (effect_df["gene_symbol"] == "DNMT3L")
    ].copy()
    if not dnmt3l.empty:
        dnmt3l = dnmt3l.sort_values("delta_ds_minus_control_log2_cpm")
        fig, ax = plt.subplots(figsize=(8.5, max(4.5, 0.32 * len(dnmt3l) + 1.5)))
        values = dnmt3l["delta_ds_minus_control_log2_cpm"].astype(float).to_numpy()
        colors = ["#2f6f9f" if value >= 0 else "#b65d3b" for value in values]
        ax.barh(range(len(dnmt3l)), values, color=colors)
        ax.axvline(0, color="#333333", linewidth=0.8)
        ax.set_yticks(range(len(dnmt3l)))
        ax.set_yticklabels(dnmt3l["label"], fontsize=8)
        ax.set_xlabel("Mean donor log2(CPM+1), DS - control")
        ax.set_title("DNMT3L effect by fetal cortex cluster")
        fig.tight_layout()
        out = figdir / "dnmt3l_cluster_effects.png"
        fig.savefig(out, dpi=220)
        fig.savefig(out.with_suffix(".pdf"))
        plt.close(fig)
        outputs.extend([out, out.with_suffix(".pdf")])

    det = cell_group_df[
        (cell_group_df["label_type"] == "author_cell_type") & (cell_group_df["gene_symbol"] == "DNMT3L")
    ].copy()
    if not det.empty:
        pivot = det.pivot(index="label", columns="group", values="detected_fraction").fillna(0)
        pivot = pivot.sort_values("DS" if "DS" in pivot.columns else pivot.columns[0])
        fig, ax = plt.subplots(figsize=(8.5, max(4.5, 0.32 * len(pivot) + 1.5)))
        y = np.arange(len(pivot))
        width = 0.38
        ax.barh(y - width / 2, pivot.get("Control", pd.Series(0, index=pivot.index)), height=width, label="Control", color="#7a8c99")
        ax.barh(y + width / 2, pivot.get("DS", pd.Series(0, index=pivot.index)), height=width, label="DS", color="#2f6f9f")
        ax.set_yticks(y)
        ax.set_yticklabels(pivot.index, fontsize=8)
        ax.set_xlabel("DNMT3L detected fraction")
        ax.set_title("DNMT3L detection by author cell type")
        ax.legend(frameon=False)
        fig.tight_layout()
        out = figdir / "dnmt3l_detection_by_author_cell_type.png"
        fig.savefig(out, dpi=220)
        fig.savefig(out.with_suffix(".pdf"))
        plt.close(fig)
        outputs.extend([out, out.with_suffix(".pdf")])

    return outputs


def main() -> int:
    args = parse_args()
    h5ad_path = Path(args.h5ad)
    outdir = Path(args.outdir)
    figdir = Path(args.figdir)
    outdir.mkdir(parents=True, exist_ok=True)

    target_by_ensembl, target_by_symbol = read_human_targets(Path(args.target_genes))
    with h5py.File(h5ad_path, "r") as f:
        target_rows = read_var_targets(f, target_by_ensembl, target_by_symbol)
        obs_df = obs_dataframe(f)

    genes = [row["gene_symbol"] for row in target_rows]
    write_tsv(
        outdir / "target_feature_map.tsv",
        target_rows,
        ["gene_symbol", "feature_index", "ensembl_id", "gene_name"],
    )
    missing = [gene for gene in TARGET_ORDER if gene not in genes]
    if missing:
        raise SystemExit(f"Missing target genes in H5AD var/raw.var: {', '.join(missing)}")

    cache_path = outdir / "target_counts_cache.npz"
    target_counts = extract_target_counts(
        h5ad_path,
        args.matrix_path,
        target_rows,
        cache_path,
        args.chunk_cells,
        args.force_extract,
    )

    cell_values = build_cell_values(obs_df, target_counts, genes)
    label_types = ["cell_type", "author_cell_type", "cell_class", "cluster_name", "cluster_name_subset"]

    pseudobulk = pseudobulk_rows(cell_values, genes, label_types)
    pseudobulk_fields = [
        "label_type",
        "label",
        "gene_symbol",
        "group",
        "donor_id",
        "sample_name",
        "dev_PCW",
        "stage_bin",
        "sex",
        "n_cells",
        "n_detected_cells",
        "detected_fraction",
        "total_rna_umi",
        "target_umi",
        "cpm",
        "log2_cpm",
    ]
    write_tsv(outdir / "target_donor_pseudobulk.tsv", pseudobulk, pseudobulk_fields)

    cell_group = cell_group_rows(cell_values, genes, label_types)
    cell_group_fields = [
        "label_type",
        "label",
        "group",
        "gene_symbol",
        "n_cells",
        "n_detected_cells",
        "detected_fraction",
        "total_rna_umi",
        "target_umi",
        "cpm",
        "log2_cpm",
        "mean_count_per_cell",
    ]
    write_tsv(outdir / "target_cell_group_summary.tsv", cell_group, cell_group_fields)

    effects = effect_rows(pd.DataFrame(pseudobulk), args.min_cells_per_donor_label)
    effect_fields = [
        "label_type",
        "label",
        "context",
        "gene_symbol",
        "n_ds_donors",
        "n_control_donors",
        "n_ds_cells",
        "n_control_cells",
        "mean_ds_log2_cpm",
        "mean_control_log2_cpm",
        "delta_ds_minus_control_log2_cpm",
        "mean_ds_detected_fraction",
        "mean_control_detected_fraction",
        "delta_ds_minus_control_detected_fraction",
    ]
    write_tsv(outdir / "target_celltype_effect_summary.tsv", effects, effect_fields)

    cell_value_fields = [
        "cell_id",
        "group",
        "donor_id",
        "sample_name",
        "dev_PCW",
        "stage_bin",
        "sex",
        "cell_type",
        "author_cell_type",
        "cell_class",
        "cluster_name",
        "cluster_name_subset",
        "nCount_RNA",
        *genes,
    ]
    cell_values[cell_value_fields].to_csv(outdir / "target_cell_values.tsv.gz", sep="\t", index=False, compression="gzip")

    figures = plot_outputs(pd.DataFrame(effects), pd.DataFrame(cell_group), figdir)
    print(
        f"Wrote CELLxGENE target summaries for {len(cell_values)} cells, "
        f"{len(pseudobulk)} donor-label rows, {len(effects)} effects, and {len(figures)} figures."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
