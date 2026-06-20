#!/usr/bin/env python3
"""Integrate mouse methylation-expression candidates with human DS expression."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np


DEFAULT_MOUSE_CANDIDATES = "results/integrated_candidate_genes/mouse_triple_hit_candidate_gene_ranking.tsv"
DEFAULT_OUTDIR = "results/cross_species_candidate_expression"
DEFAULT_GENCODE = "data/reference/gencode/gencode.v39.annotation.gtf.gz"
DEFAULT_GSE301886 = "data/geo/GSE301886/supplementary/GSE301886_T21_raw_counts.txt.gz"
DEFAULT_GSE305139 = "data/geo/GSE305139/supplementary/GSE305139_raw_counts_tissue_bulk.tsv.gz"
DEFAULT_HUMAN_DS_BULK_COUNT_DE = "results/human_ds_bulk_count_de/human_ds_bulk_count_de_all_contexts.tsv"
DEFAULT_H5AD = "data/cellxgene/GSE305153_complete_fetal_cortex_cellxgene.h5ad"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mouse-candidates", default=DEFAULT_MOUSE_CANDIDATES)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--gencode-gtf", default=DEFAULT_GENCODE)
    parser.add_argument("--gse301886-counts", default=DEFAULT_GSE301886)
    parser.add_argument("--gse305139-counts", default=DEFAULT_GSE305139)
    parser.add_argument("--human-ds-bulk-count-de", default=DEFAULT_HUMAN_DS_BULK_COUNT_DE)
    parser.add_argument("--gse305153-h5ad", default=DEFAULT_H5AD)
    parser.add_argument("--matrix-path", default="raw/X")
    parser.add_argument("--chunk-cells", type=int, default=5000)
    parser.add_argument("--effect-threshold", type=float, default=0.25)
    parser.add_argument("--min-donors-per-group", type=int, default=3)
    parser.add_argument("--min-cells-per-donor-label", type=int, default=20)
    parser.add_argument("--max-candidates", type=int, default=0, help="0 means all symbol-mappable candidates.")
    return parser.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="", errors="replace")
    return path.open(newline="", errors="replace")


def fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, (float, np.floating)):
        number = float(value)
        return f"{number:.6g}" if math.isfinite(number) else ""
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    return str(value)


def safe_float(value: object) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return math.nan
    return number if math.isfinite(number) else math.nan


def mean(values: list[float]) -> float:
    clean = [value for value in values if math.isfinite(value)]
    return sum(clean) / len(clean) if clean else math.nan


def direction(delta: float, threshold: float) -> str:
    if not math.isfinite(delta):
        return "not_available"
    if delta >= threshold:
        return "up"
    if delta <= -threshold:
        return "down"
    return "flat"


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})


def parse_characteristics(value: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for part in value.split(";"):
        part = part.strip()
        if not part or "=" not in part:
            continue
        key, val = part.split("=", 1)
        result[key.strip().lower()] = val.strip()
    return result


def load_geo_metadata(accession: str) -> list[dict[str, str]]:
    path = Path("data/geo") / accession / "metadata" / "sample_metadata.tsv"
    rows: list[dict[str, str]] = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.append({**row, **parse_characteristics(row.get("characteristics", ""))})
    return rows


def strip_version(value: str) -> str:
    return value.split(".", 1)[0].strip()


def human_symbol_from_mouse(symbol: str) -> str:
    symbol = symbol.strip()
    if not symbol:
        return ""
    if symbol.startswith("Gm") or symbol.endswith("Rik") or "Rik" in symbol:
        return ""
    if re.match(r"^\d", symbol):
        return ""
    return symbol.upper()


def read_mouse_candidates(path: Path, max_candidates: int) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    seen_human_symbols: set[str] = set()
    with path.open(newline="") as handle:
        for rank, row in enumerate(csv.DictReader(handle, delimiter="\t"), start=1):
            human_symbol = human_symbol_from_mouse(row.get("gene_symbol", ""))
            if not human_symbol or human_symbol in seen_human_symbols:
                continue
            seen_human_symbols.add(human_symbol)
            out = dict(row)
            out["mouse_rank"] = rank
            out["human_symbol"] = human_symbol
            rows.append(out)
            if max_candidates and len(rows) >= max_candidates:
                break
    return rows


def parse_gtf_attributes(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for match in re.finditer(r'([A-Za-z0-9_]+) "([^"]+)"', text):
        attrs[match.group(1)] = match.group(2)
    return attrs


def read_gencode_gene_symbols(path: Path) -> tuple[dict[str, str], set[str]]:
    by_ensembl: dict[str, str] = {}
    symbols: set[str] = set()
    with open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = parse_gtf_attributes(parts[8])
            gene_id = strip_version(attrs.get("gene_id", ""))
            symbol = attrs.get("gene_name", "").upper()
            if gene_id and symbol:
                by_ensembl[gene_id] = symbol
                symbols.add(symbol)
    return by_ensembl, symbols


def metadata_gse301886() -> dict[str, dict[str, object]]:
    by_sample: dict[str, dict[str, object]] = {}
    for row in load_geo_metadata("GSE301886"):
        title = row.get("title", "")
        sample = title.split(",", 1)[1].strip() if "," in title else title.strip()
        karyotype = row.get("karyotype", "")
        group = "DS" if karyotype == "T21" else "Control" if karyotype == "Control" else "Other"
        by_sample[sample] = {
            "sample_label": sample,
            "group": group,
            "context": row.get("tissue", row.get("source_name", "unknown")),
            "pcw": row.get("gestational age", ""),
            "sex": row.get("sex", row.get("Sex", "")),
            "use": True,
        }
    return by_sample


def metadata_gse305139() -> dict[str, dict[str, object]]:
    by_library: dict[str, dict[str, object]] = {}
    for row in load_geo_metadata("GSE305139"):
        library = row.get("library igf_id", "")
        group_raw = row.get("group", "")
        group = "DS" if group_raw == "DS" else "Control" if group_raw in {"CON", "Control"} else "Other"
        by_library[library] = {
            "sample_label": row.get("name", library),
            "group": group,
            "context": "fetal_brain_final",
            "pcw": row.get("dev stage_pcw", ""),
            "sex": row.get("sex", row.get("Sex", "")),
            "use": row.get("in final_analysis", "").upper() == "TRUE",
        }
    return by_library


def read_bulk_counts(
    path: Path,
    dataset: str,
    sample_meta: dict[str, dict[str, object]],
    target_symbols: set[str],
    threshold: float,
    ensembl_to_symbol: dict[str, str] | None = None,
    gene_symbol_col: int | None = None,
    sample_start: int = 1,
) -> tuple[list[dict[str, object]], list[dict[str, object]], set[str]]:
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        samples = header[sample_start:]
        library_sizes = np.zeros(len(samples), dtype=np.float64)
        counts_by_symbol = {symbol: np.zeros(len(samples), dtype=np.float64) for symbol in target_symbols}

        for row in reader:
            if len(row) < sample_start + 1:
                continue
            if gene_symbol_col is not None:
                symbol = row[gene_symbol_col].upper()
            else:
                symbol = (ensembl_to_symbol or {}).get(strip_version(row[0]), "")
            values = np.array([safe_float(value) for value in row[sample_start : sample_start + len(samples)]])
            values[~np.isfinite(values)] = 0.0
            library_sizes += values
            if symbol in counts_by_symbol:
                counts_by_symbol[symbol] += values

    sample_rows: list[dict[str, object]] = []
    present_symbols: set[str] = set()
    for symbol, counts in counts_by_symbol.items():
        if np.sum(counts) <= 0:
            continue
        present_symbols.add(symbol)
        for idx, sample in enumerate(samples):
            meta = sample_meta.get(sample, {})
            if not meta.get("use", False) or meta.get("group") not in {"DS", "Control"}:
                continue
            count = float(counts[idx])
            lib = float(library_sizes[idx])
            cpm = count / lib * 1_000_000 if lib > 0 else math.nan
            log2_cpm = math.log2(cpm + 1.0) if math.isfinite(cpm) else math.nan
            sample_rows.append(
                {
                    "dataset": dataset,
                    "context": meta.get("context", "unknown"),
                    "sample": sample,
                    "sample_label": meta.get("sample_label", sample),
                    "group": meta.get("group", ""),
                    "pcw": meta.get("pcw", ""),
                    "sex": meta.get("sex", ""),
                    "human_symbol": symbol,
                    "count": count,
                    "library_size": lib,
                    "cpm": cpm,
                    "log2_cpm_plus_1": log2_cpm,
                    "detected": "yes" if count > 0 else "no",
                }
            )

    effect_rows = summarize_bulk_effects(sample_rows, threshold)
    return sample_rows, effect_rows, present_symbols


def summarize_bulk_effects(sample_rows: list[dict[str, object]], threshold: float) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str], dict[str, list[dict[str, object]]]] = defaultdict(lambda: defaultdict(list))
    for row in sample_rows:
        grouped[(str(row["dataset"]), str(row["context"]), str(row["human_symbol"]))][str(row["group"])].append(row)

    effects: list[dict[str, object]] = []
    for (dataset, context, symbol), by_group in grouped.items():
        ds = by_group.get("DS", [])
        control = by_group.get("Control", [])
        if not ds or not control:
            continue
        ds_log = [safe_float(row["log2_cpm_plus_1"]) for row in ds]
        control_log = [safe_float(row["log2_cpm_plus_1"]) for row in control]
        ds_cpm = [safe_float(row["cpm"]) for row in ds]
        control_cpm = [safe_float(row["cpm"]) for row in control]
        delta = mean(ds_log) - mean(control_log)
        effects.append(
            {
                "dataset": dataset,
                "context": context,
                "human_symbol": symbol,
                "n_ds": len(ds),
                "n_control": len(control),
                "mean_ds_log2_cpm_plus_1": mean(ds_log),
                "mean_control_log2_cpm_plus_1": mean(control_log),
                "delta_ds_minus_control_log2_cpm_plus_1": delta,
                "mean_ds_cpm": mean(ds_cpm),
                "mean_control_cpm": mean(control_cpm),
                "ds_detected_fraction": sum(1 for row in ds if row["detected"] == "yes") / len(ds),
                "control_detected_fraction": sum(1 for row in control if row["detected"] == "yes") / len(control),
                "direction": direction(delta, threshold),
            }
        )
    return sorted(effects, key=lambda row: (row["dataset"], row["context"], row["human_symbol"]))


def read_bulk_count_de_effects(
    path: Path,
    target_symbols: set[str],
    threshold: float,
) -> tuple[list[dict[str, object]], set[str]]:
    if not path.exists():
        return [], set()
    best_by_key: dict[tuple[str, str, str], dict[str, object]] = {}
    present: set[str] = set()
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            symbol = row.get("gene_symbol", "").upper()
            if symbol not in target_symbols:
                continue
            present.add(symbol)
            dataset = row.get("dataset", "")
            context = row.get("context", "")
            if dataset not in {"GSE301886", "GSE305139"} or not context:
                continue
            effect = safe_float(row.get("effect_log2fc_for_claim"))
            if not math.isfinite(effect):
                continue
            out = {
                "dataset": dataset,
                "context": context,
                "human_symbol": symbol,
                "n_ds": safe_float(row.get("n_ds")),
                "n_control": safe_float(row.get("n_control")),
                "mean_ds_log2_cpm_plus_1": safe_float(row.get("mean_ds_log2_cpm_plus_1")),
                "mean_control_log2_cpm_plus_1": safe_float(row.get("mean_control_log2_cpm_plus_1")),
                "delta_ds_minus_control_log2_cpm_plus_1": effect,
                "descriptive_delta_log2_cpm_plus_1": safe_float(row.get("delta_ds_minus_control_log2_cpm_plus_1")),
                "mean_ds_cpm": safe_float(row.get("mean_ds_cpm")),
                "mean_control_cpm": safe_float(row.get("mean_control_cpm")),
                "ds_detected_fraction": safe_float(row.get("ds_detected_fraction")),
                "control_detected_fraction": safe_float(row.get("control_detected_fraction")),
                "effect_source": row.get("effect_source", ""),
                "edgeR_logFC": safe_float(row.get("logFC")),
                "edgeR_PValue": safe_float(row.get("PValue")),
                "edgeR_FDR": safe_float(row.get("FDR")),
                "tested_by_edgeR": row.get("tested_by_edgeR", ""),
                "direction": direction(effect, threshold),
            }
            key = (dataset, context, symbol)
            current = best_by_key.get(key)
            current_expr = (
                safe_float(current.get("mean_ds_log2_cpm_plus_1")) + safe_float(current.get("mean_control_log2_cpm_plus_1"))
                if current
                else -math.inf
            )
            new_expr = safe_float(out["mean_ds_log2_cpm_plus_1"]) + safe_float(out["mean_control_log2_cpm_plus_1"])
            if current is None or new_expr > current_expr:
                best_by_key[key] = out
    return sorted(best_by_key.values(), key=lambda row: (row["dataset"], row["context"], row["human_symbol"])), present


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


def h5ad_feature_map(h5ad_path: Path, target_symbols: set[str]) -> tuple[dict[int, str], set[str], set[str]]:
    feature_to_symbol: dict[int, str] = {}
    h5_symbols: set[str] = set()
    with h5py.File(h5ad_path, "r") as f:
        var_group = f["raw/var"] if "raw" in f and "var" in f["raw"] else f["var"]
        gene_names = read_h5_column(var_group, "gene_name" if "gene_name" in var_group else "feature_name")
        for idx, gene_name in enumerate(gene_names):
            symbol = str(gene_name).upper()
            h5_symbols.add(symbol)
            if symbol in target_symbols:
                feature_to_symbol[idx] = symbol
    return feature_to_symbol, set(feature_to_symbol.values()), h5_symbols


def build_h5ad_groups(f) -> tuple[list[dict[str, object]], list[np.ndarray], np.ndarray]:
    obs = f["obs"]
    group = read_h5_column(obs, "group")
    group = np.array(["Control" if value == "CON" else str(value) for value in group], dtype=object)
    donor = read_h5_column(obs, "donor_id")
    sample = read_h5_column(obs, "sample_name")
    pcw = read_h5_column(obs, "dev_PCW")
    sex = read_h5_column(obs, "sex") if "sex" in obs else np.array([""] * len(group), dtype=object)
    author_cell_type = read_h5_column(obs, "author_cell_type")
    n_count = np.asarray(read_h5_column(obs, "nCount_RNA"), dtype=np.float64)
    n_cells = len(group)

    group_rows: list[dict[str, object]] = []
    group_ids_by_label_type: list[np.ndarray] = []
    key_to_id: dict[tuple[str, str, str, str], int] = {}

    for label_type, labels in [
        ("all_cells", np.array(["all_cells"] * n_cells, dtype=object)),
        ("author_cell_type", author_cell_type),
    ]:
        group_ids = np.full(n_cells, -1, dtype=np.int32)
        for idx in range(n_cells):
            if group[idx] not in {"DS", "Control"}:
                continue
            label = str(labels[idx]) or "unknown"
            key = (label_type, label, str(donor[idx]), str(group[idx]))
            group_id = key_to_id.get(key)
            if group_id is None:
                group_id = len(group_rows)
                key_to_id[key] = group_id
                group_rows.append(
                    {
                        "label_type": label_type,
                        "label": label,
                        "donor_id": str(donor[idx]),
                        "group": str(group[idx]),
                        "sample_name": str(sample[idx]),
                        "dev_pcw": str(pcw[idx]),
                        "stage_bin": stage_bin(pcw[idx]),
                        "sex": str(sex[idx]),
                        "n_cells": 0,
                        "total_rna_umi": 0.0,
                    }
                )
            group_ids[idx] = group_id
            group_rows[group_id]["n_cells"] = int(group_rows[group_id]["n_cells"]) + 1
            group_rows[group_id]["total_rna_umi"] = float(group_rows[group_id]["total_rna_umi"]) + float(n_count[idx])
        group_ids_by_label_type.append(group_ids)

    return group_rows, group_ids_by_label_type, n_count


def aggregate_h5ad_candidates(
    h5ad_path: Path,
    matrix_path: str,
    target_symbols: set[str],
    chunk_cells: int,
    min_cells_per_donor_label: int,
    threshold: float,
    min_donors: int,
) -> tuple[list[dict[str, object]], list[dict[str, object]], set[str]]:
    feature_to_symbol, present_symbols, _all_h5_symbols = h5ad_feature_map(h5ad_path, target_symbols)
    symbols = sorted(present_symbols)
    symbol_to_pos = {symbol: idx for idx, symbol in enumerate(symbols)}
    feature_indices = np.array(sorted(feature_to_symbol), dtype=np.int64)
    feature_to_pos = np.array([symbol_to_pos[feature_to_symbol[int(idx)]] for idx in feature_indices], dtype=np.int32)

    if not symbols:
        return [], [], set()

    with h5py.File(h5ad_path, "r") as f:
        group_rows, group_ids_by_label_type, _n_count = build_h5ad_groups(f)
        counts = np.zeros((len(group_rows), len(symbols)), dtype=np.float64)
        detected = np.zeros((len(group_rows), len(symbols)), dtype=np.int32)
        x = f[matrix_path]
        shape = x.attrs["shape"]
        n_cells = int(shape[0])
        indptr = x["indptr"][:]

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
            mask = np.isin(indices, feature_indices)
            if not np.any(mask):
                continue
            matched_indices = indices[mask]
            positions = np.searchsorted(feature_indices, matched_indices)
            gene_positions = feature_to_pos[positions]
            local_rows = row_ids[mask]
            values = data[mask]
            for group_ids in group_ids_by_label_type:
                aggregate_group_ids = group_ids[start:end][local_rows]
                valid = aggregate_group_ids >= 0
                if np.any(valid):
                    np.add.at(counts, (aggregate_group_ids[valid], gene_positions[valid]), values[valid])
                    np.add.at(detected, (aggregate_group_ids[valid], gene_positions[valid]), 1)

    donor_rows: list[dict[str, object]] = []
    for group_id, group_row in enumerate(group_rows):
        total_rna = safe_float(group_row["total_rna_umi"])
        n_cells = int(group_row["n_cells"])
        if n_cells < min_cells_per_donor_label:
            continue
        for symbol, pos in symbol_to_pos.items():
            target_umi = float(counts[group_id, pos])
            cpm = target_umi / total_rna * 1_000_000 if total_rna > 0 else math.nan
            donor_rows.append(
                {
                    **group_row,
                    "human_symbol": symbol,
                    "target_umi": target_umi,
                    "cpm": cpm,
                    "log2_cpm_plus_1": math.log2(cpm + 1.0) if math.isfinite(cpm) else math.nan,
                    "n_detected_cells": int(detected[group_id, pos]),
                    "detected_fraction": int(detected[group_id, pos]) / n_cells if n_cells else math.nan,
                }
            )

    effect_rows = summarize_h5ad_effects(donor_rows, threshold, min_donors)
    return donor_rows, effect_rows, present_symbols


def summarize_h5ad_effects(
    donor_rows: list[dict[str, object]],
    threshold: float,
    min_donors: int,
) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str], dict[str, list[dict[str, object]]]] = defaultdict(lambda: defaultdict(list))
    for row in donor_rows:
        grouped[(str(row["label_type"]), str(row["label"]), str(row["human_symbol"]))][str(row["group"])].append(row)

    effects: list[dict[str, object]] = []
    for (label_type, label, symbol), by_group in grouped.items():
        ds = by_group.get("DS", [])
        control = by_group.get("Control", [])
        if len(ds) < min_donors or len(control) < min_donors:
            continue
        ds_log = [safe_float(row["log2_cpm_plus_1"]) for row in ds]
        control_log = [safe_float(row["log2_cpm_plus_1"]) for row in control]
        delta = mean(ds_log) - mean(control_log)
        effects.append(
            {
                "dataset": "GSE305153",
                "label_type": label_type,
                "label": label,
                "context": "all_pcw",
                "human_symbol": symbol,
                "n_ds_donors": len(ds),
                "n_control_donors": len(control),
                "n_ds_cells": sum(int(row["n_cells"]) for row in ds),
                "n_control_cells": sum(int(row["n_cells"]) for row in control),
                "mean_ds_log2_cpm_plus_1": mean(ds_log),
                "mean_control_log2_cpm_plus_1": mean(control_log),
                "delta_ds_minus_control_log2_cpm_plus_1": delta,
                "mean_ds_detected_fraction": mean([safe_float(row["detected_fraction"]) for row in ds]),
                "mean_control_detected_fraction": mean([safe_float(row["detected_fraction"]) for row in control]),
                "delta_detected_fraction": mean([safe_float(row["detected_fraction"]) for row in ds])
                - mean([safe_float(row["detected_fraction"]) for row in control]),
                "direction": direction(delta, threshold),
            }
        )
    return sorted(effects, key=lambda row: (row["label_type"], row["label"], row["human_symbol"]))


def summarize_candidate_support(
    mouse_candidates: list[dict[str, object]],
    bulk_effects: list[dict[str, object]],
    h5ad_effects: list[dict[str, object]],
    bulk_present: set[str],
    h5ad_present: set[str],
    threshold: float,
) -> list[dict[str, object]]:
    bulk_by_symbol: dict[str, list[dict[str, object]]] = defaultdict(list)
    h5_by_symbol: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in bulk_effects:
        bulk_by_symbol[str(row["human_symbol"])].append(row)
    for row in h5ad_effects:
        h5_by_symbol[str(row["human_symbol"])].append(row)

    summary: list[dict[str, object]] = []
    for mouse in mouse_candidates:
        symbol = str(mouse["human_symbol"])
        mouse_direction = str(mouse.get("gse268951_expression_direction", "not_available"))
        bulk_rows = bulk_by_symbol.get(symbol, [])
        h5_rows = h5_by_symbol.get(symbol, [])
        h5_author_rows = [row for row in h5_rows if row.get("label_type") == "author_cell_type"]
        h5_all_rows = [row for row in h5_rows if row.get("label_type") == "all_cells"]

        same_bulk = same_direction_count(bulk_rows, mouse_direction)
        opposite_bulk = opposite_direction_count(bulk_rows, mouse_direction)
        same_celltypes = same_direction_count(h5_author_rows, mouse_direction)
        opposite_celltypes = opposite_direction_count(h5_author_rows, mouse_direction)

        best_bulk_same = best_same_direction_row(bulk_rows, mouse_direction)
        best_h5_same = best_same_direction_row(h5_author_rows, mouse_direction)
        best_h5_abs = max(
            h5_author_rows,
            key=lambda row: abs(safe_float(row.get("delta_ds_minus_control_log2_cpm_plus_1"))),
            default={},
        )
        all_cells = h5_all_rows[0] if h5_all_rows else {}

        human_score = min(2, same_bulk) + min(4, same_celltypes)
        if same_bulk and same_celltypes:
            human_score += 1
        if opposite_bulk > same_bulk and opposite_celltypes > same_celltypes and mouse_direction in {"up", "down"}:
            support_level = "conflicting"
        elif same_bulk and same_celltypes >= 3:
            support_level = "strong"
        elif same_bulk or same_celltypes:
            support_level = "moderate"
        else:
            support_level = "weak_or_none"

        summary.append(
            {
                "mouse_rank": mouse.get("mouse_rank"),
                "mouse_gene_symbol": mouse.get("gene_symbol"),
                "human_symbol": symbol,
                "ortholog_mapping": "symbol_uppercase_first_pass",
                "mouse_evidence_score": safe_float(mouse.get("evidence_score")),
                "mouse_promoter_dmrs": mouse.get("n_promoter_2kb_dmrs"),
                "mouse_gene_body_dmrs": mouse.get("n_gene_body_dmrs"),
                "mouse_gse268951_log2fc": safe_float(mouse.get("gse268951_log2fc_oe_vs_control")),
                "mouse_gse268951_direction": mouse_direction,
                "mouse_dp10_delta": safe_float(mouse.get("dp10_delta_log2_cpm_plus_1")),
                "mouse_dp10_direction": mouse.get("dp10_expression_direction"),
                "present_in_human_bulk_matrices": "yes" if symbol in bulk_present else "no",
                "present_in_gse305153_h5ad": "yes" if symbol in h5ad_present else "no",
                "bulk_contexts_same_direction": same_bulk,
                "bulk_contexts_opposite_direction": opposite_bulk,
                "celltypes_same_direction": same_celltypes,
                "celltypes_opposite_direction": opposite_celltypes,
                "gse305153_all_cells_delta": safe_float(all_cells.get("delta_ds_minus_control_log2_cpm_plus_1")),
                "gse305153_all_cells_direction": all_cells.get("direction", "not_available"),
                "best_bulk_same_direction_context": context_label(best_bulk_same),
                "best_bulk_same_direction_delta": safe_float(best_bulk_same.get("delta_ds_minus_control_log2_cpm_plus_1")),
                "best_celltype_same_direction_label": best_h5_same.get("label", ""),
                "best_celltype_same_direction_delta": safe_float(best_h5_same.get("delta_ds_minus_control_log2_cpm_plus_1")),
                "largest_abs_celltype_label": best_h5_abs.get("label", ""),
                "largest_abs_celltype_delta": safe_float(best_h5_abs.get("delta_ds_minus_control_log2_cpm_plus_1")),
                "human_support_score": human_score,
                "cross_species_priority_score": safe_float(mouse.get("evidence_score")) + human_score,
                "human_support_level": support_level,
                "description": mouse.get("description", ""),
            }
        )

    summary.sort(
        key=lambda row: (
            -safe_float(row["cross_species_priority_score"]),
            str(row["human_support_level"]),
            -safe_float(row["mouse_evidence_score"]),
            str(row["human_symbol"]),
        )
    )
    return summary


def same_direction_count(rows: list[dict[str, object]], reference_direction: str) -> int:
    if reference_direction not in {"up", "down"}:
        return 0
    return sum(1 for row in rows if row.get("direction") == reference_direction)


def opposite_direction_count(rows: list[dict[str, object]], reference_direction: str) -> int:
    if reference_direction == "up":
        return sum(1 for row in rows if row.get("direction") == "down")
    if reference_direction == "down":
        return sum(1 for row in rows if row.get("direction") == "up")
    return 0


def best_same_direction_row(rows: list[dict[str, object]], reference_direction: str) -> dict[str, object]:
    same = [row for row in rows if row.get("direction") == reference_direction]
    if not same:
        return {}
    return max(same, key=lambda row: abs(safe_float(row.get("delta_ds_minus_control_log2_cpm_plus_1"))))


def context_label(row: dict[str, object]) -> str:
    if not row:
        return ""
    return f"{row.get('dataset')}:{row.get('context')}"


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    mouse_candidates = read_mouse_candidates(Path(args.mouse_candidates), args.max_candidates)
    target_symbols = {str(row["human_symbol"]) for row in mouse_candidates}
    ensembl_to_symbol, gencode_symbols = read_gencode_gene_symbols(Path(args.gencode_gtf))

    gse301886_sample_rows, gse301886_effects, gse301886_present = read_bulk_counts(
        Path(args.gse301886_counts),
        "GSE301886",
        metadata_gse301886(),
        target_symbols,
        args.effect_threshold,
        ensembl_to_symbol=ensembl_to_symbol,
        sample_start=1,
    )
    gse305139_sample_rows, gse305139_effects, gse305139_present = read_bulk_counts(
        Path(args.gse305139_counts),
        "GSE305139",
        metadata_gse305139(),
        target_symbols,
        args.effect_threshold,
        gene_symbol_col=1,
        sample_start=2,
    )

    h5ad_donor_rows, h5ad_effects, h5ad_present = aggregate_h5ad_candidates(
        Path(args.gse305153_h5ad),
        args.matrix_path,
        target_symbols,
        args.chunk_cells,
        args.min_cells_per_donor_label,
        args.effect_threshold,
        args.min_donors_per_group,
    )

    bulk_sample_rows = gse301886_sample_rows + gse305139_sample_rows
    count_de_bulk_effects, count_de_present = read_bulk_count_de_effects(
        Path(args.human_ds_bulk_count_de),
        target_symbols,
        args.effect_threshold,
    )
    if count_de_bulk_effects:
        bulk_effects = count_de_bulk_effects
        bulk_present = count_de_present
    else:
        bulk_effects = gse301886_effects + gse305139_effects
        bulk_present = gse301886_present | gse305139_present
    priority_rows = summarize_candidate_support(
        mouse_candidates,
        bulk_effects,
        h5ad_effects,
        bulk_present,
        h5ad_present,
        args.effect_threshold,
    )

    write_tsv(
        outdir / "human_bulk_candidate_sample_values.tsv",
        bulk_sample_rows,
        [
            "dataset",
            "context",
            "sample",
            "sample_label",
            "group",
            "pcw",
            "sex",
            "human_symbol",
            "count",
            "library_size",
            "cpm",
            "log2_cpm_plus_1",
            "detected",
        ],
    )
    write_tsv(
        outdir / "human_bulk_candidate_effects.tsv",
        bulk_effects,
        [
            "dataset",
            "context",
            "human_symbol",
            "n_ds",
            "n_control",
            "mean_ds_log2_cpm_plus_1",
            "mean_control_log2_cpm_plus_1",
            "delta_ds_minus_control_log2_cpm_plus_1",
            "descriptive_delta_log2_cpm_plus_1",
            "mean_ds_cpm",
            "mean_control_cpm",
            "ds_detected_fraction",
            "control_detected_fraction",
            "effect_source",
            "edgeR_logFC",
            "edgeR_PValue",
            "edgeR_FDR",
            "tested_by_edgeR",
            "direction",
        ],
    )
    write_tsv(
        outdir / "gse305153_candidate_donor_pseudobulk.tsv",
        h5ad_donor_rows,
        [
            "label_type",
            "label",
            "donor_id",
            "group",
            "sample_name",
            "dev_pcw",
            "stage_bin",
            "sex",
            "n_cells",
            "total_rna_umi",
            "human_symbol",
            "target_umi",
            "cpm",
            "log2_cpm_plus_1",
            "n_detected_cells",
            "detected_fraction",
        ],
    )
    write_tsv(
        outdir / "gse305153_candidate_celltype_effects.tsv",
        h5ad_effects,
        [
            "dataset",
            "label_type",
            "label",
            "context",
            "human_symbol",
            "n_ds_donors",
            "n_control_donors",
            "n_ds_cells",
            "n_control_cells",
            "mean_ds_log2_cpm_plus_1",
            "mean_control_log2_cpm_plus_1",
            "delta_ds_minus_control_log2_cpm_plus_1",
            "mean_ds_detected_fraction",
            "mean_control_detected_fraction",
            "delta_detected_fraction",
            "direction",
        ],
    )
    write_tsv(
        outdir / "cross_species_candidate_priority.tsv",
        priority_rows,
        [
            "mouse_rank",
            "mouse_gene_symbol",
            "human_symbol",
            "ortholog_mapping",
            "mouse_evidence_score",
            "mouse_promoter_dmrs",
            "mouse_gene_body_dmrs",
            "mouse_gse268951_log2fc",
            "mouse_gse268951_direction",
            "mouse_dp10_delta",
            "mouse_dp10_direction",
            "present_in_human_bulk_matrices",
            "present_in_gse305153_h5ad",
            "bulk_contexts_same_direction",
            "bulk_contexts_opposite_direction",
            "celltypes_same_direction",
            "celltypes_opposite_direction",
            "gse305153_all_cells_delta",
            "gse305153_all_cells_direction",
            "best_bulk_same_direction_context",
            "best_bulk_same_direction_delta",
            "best_celltype_same_direction_label",
            "best_celltype_same_direction_delta",
            "largest_abs_celltype_label",
            "largest_abs_celltype_delta",
            "human_support_score",
            "cross_species_priority_score",
            "human_support_level",
            "description",
        ],
    )

    summary_rows = [
        {"metric": "mouse_symbol_mappable_candidates", "value": len(mouse_candidates)},
        {"metric": "candidate_symbols_in_gencode", "value": len(target_symbols & gencode_symbols)},
        {"metric": "candidate_symbols_in_bulk_matrices", "value": len(bulk_present)},
        {"metric": "candidate_symbols_in_gse305153_h5ad", "value": len(h5ad_present)},
        {"metric": "bulk_effect_rows", "value": len(bulk_effects)},
        {"metric": "gse305153_effect_rows", "value": len(h5ad_effects)},
        {
            "metric": "priority_genes_with_moderate_or_strong_human_support",
            "value": sum(1 for row in priority_rows if row["human_support_level"] in {"moderate", "strong"}),
        },
        {
            "metric": "priority_genes_with_strong_human_support",
            "value": sum(1 for row in priority_rows if row["human_support_level"] == "strong"),
        },
    ]
    write_tsv(outdir / "cross_species_candidate_summary.tsv", summary_rows, ["metric", "value"])


if __name__ == "__main__":
    main()
