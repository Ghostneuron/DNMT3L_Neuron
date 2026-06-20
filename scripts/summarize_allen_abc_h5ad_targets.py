#!/usr/bin/env python3
"""Extract target-gene expression from an Allen ABC H5AD matrix.

This script intentionally uses h5py directly so it can run without anndata or
scanpy. It reads only the requested target gene columns from sparse H5AD
matrices and writes compact per-cell and grouped summaries.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from typing import Any

import h5py
import numpy as np


CELL_VALUE_FIELDS = [
    "dataset",
    "cell_id",
    "gene_symbol",
    "expression",
    "detected",
    "normalization",
]

SUMMARY_FIELDS = [
    "dataset",
    "grouping",
    "group",
    "gene_symbol",
    "normalization",
    "n_cells",
    "detected_n",
    "detected_fraction",
    "mean_expression",
    "median_expression",
    "max_expression",
    "nonzero_mean_expression",
]

QC_FIELDS = [
    "metric",
    "value",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--h5ad",
        required=True,
        help="Allen ABC H5AD expression matrix.",
    )
    parser.add_argument(
        "--dataset",
        default=None,
        help="Dataset label for output rows. Defaults to H5AD filename stem.",
    )
    parser.add_argument(
        "--species",
        required=True,
        choices=["human", "mouse"],
        help="Species to select from config/target_genes.tsv.",
    )
    parser.add_argument(
        "--target-genes",
        default="config/target_genes.tsv",
        help="Project target gene table.",
    )
    parser.add_argument(
        "--outdir",
        default="results/allen_abc/h5ad_targets",
        help="Output directory.",
    )
    parser.add_argument(
        "--group-cols",
        default="anatomical_division_label,library_label",
        help="Comma-separated obs columns for grouped summaries.",
    )
    parser.add_argument(
        "--cell-metadata",
        default=None,
        help="Optional filtered Allen ABC cell metadata CSV to join by cell label.",
    )
    parser.add_argument(
        "--cell-metadata-key",
        default="cell_label",
        help="Join key in --cell-metadata.",
    )
    parser.add_argument(
        "--cell-metadata-columns",
        default=(
            "class,subclass,supertype,cluster,neurotransmitter,"
            "region_of_interest_acronym,dataset_label"
        ),
        help="Comma-separated metadata columns to add to outputs and summaries.",
    )
    parser.add_argument(
        "--chunk-rows",
        type=int,
        default=1000,
        help="CSR row chunk size for target-gene extraction.",
    )
    parser.add_argument(
        "--detection-threshold",
        type=float,
        default=0.0,
        help="Expression value above which a gene is considered detected.",
    )
    parser.add_argument(
        "--skip-cell-values",
        action="store_true",
        help="Do not write long per-cell target-gene values.",
    )
    return parser.parse_args()


def decode_scalar(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.bytes_):
        return value.decode("utf-8")
    return str(value)


def decode_array(dataset: h5py.Dataset) -> list[str]:
    try:
        values = dataset.asstr()[:]
    except (AttributeError, TypeError):
        values = dataset[:]
    return [decode_scalar(value) for value in values]


def decode_dataset_scalar(dataset: h5py.Dataset, default: str = "") -> str:
    try:
        return decode_scalar(dataset[()])
    except Exception:
        return default


def read_obs_column(obs_group: h5py.Group, column: str) -> list[str] | None:
    if column not in obs_group:
        return None
    node = obs_group[column]
    if isinstance(node, h5py.Dataset):
        return decode_array(node)
    if isinstance(node, h5py.Group) and {"codes", "categories"} <= set(node.keys()):
        codes = node["codes"][:]
        categories = decode_array(node["categories"])
        values: list[str] = []
        for code in codes:
            code_int = int(code)
            values.append(categories[code_int] if code_int >= 0 else "")
        return values
    return None


def load_targets(path: Path, species: str) -> list[dict[str, Any]]:
    targets: list[dict[str, Any]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["species"].lower() != species:
                continue
            aliases = [
                alias.upper()
                for alias in row.get("aliases", "").split("|")
                if alias
            ]
            targets.append(
                {
                    "gene_symbol": row["gene_symbol"],
                    "symbols_upper": [row["gene_symbol"].upper(), *aliases],
                    "ensembl_hint": row.get("ensembl_hint", ""),
                }
            )
    return targets


def match_targets(
    var_group: h5py.Group,
    targets: list[dict[str, Any]],
) -> tuple[list[str], list[int], list[dict[str, str]]]:
    gene_symbols = decode_array(var_group["gene_symbol"])
    gene_identifiers = decode_array(var_group["gene_identifier"])
    symbol_to_index = {
        symbol.upper(): idx
        for idx, symbol in enumerate(gene_symbols)
        if symbol
    }
    identifier_to_index = {
        identifier: idx
        for idx, identifier in enumerate(gene_identifiers)
        if identifier
    }

    matched_symbols: list[str] = []
    matched_indices: list[int] = []
    match_rows: list[dict[str, str]] = []
    for target in targets:
        matched_idx: int | None = None
        matched_by = ""
        for symbol_upper in target["symbols_upper"]:
            if symbol_upper in symbol_to_index:
                matched_idx = symbol_to_index[symbol_upper]
                matched_by = "symbol"
                break
        if matched_idx is None and target["ensembl_hint"]:
            matched_idx = identifier_to_index.get(target["ensembl_hint"])
            matched_by = "ensembl_hint" if matched_idx is not None else ""

        if matched_idx is not None:
            matched_symbols.append(gene_symbols[matched_idx])
            matched_indices.append(matched_idx)
            match_rows.append(
                {
                    "target_gene_symbol": target["gene_symbol"],
                    "matched_gene_symbol": gene_symbols[matched_idx],
                    "gene_identifier": gene_identifiers[matched_idx],
                    "gene_index": str(matched_idx),
                    "matched_by": matched_by,
                    "present": "yes",
                }
            )
        else:
            match_rows.append(
                {
                    "target_gene_symbol": target["gene_symbol"],
                    "matched_gene_symbol": "",
                    "gene_identifier": "",
                    "gene_index": "",
                    "matched_by": "",
                    "present": "no",
                }
            )
    return matched_symbols, matched_indices, match_rows


def extract_csr_values(
    x_group: h5py.Group,
    target_indices: list[int],
    chunk_rows: int,
) -> np.ndarray:
    shape = tuple(int(value) for value in x_group.attrs["shape"])
    n_cells = shape[0]
    values = np.zeros((n_cells, len(target_indices)), dtype=np.float32)
    if not target_indices:
        return values

    indptr = x_group["indptr"][:]
    sorted_targets = np.array(sorted(target_indices), dtype=x_group["indices"].dtype)
    sorted_positions = np.array(
        [target_indices.index(int(idx)) for idx in sorted_targets],
        dtype=np.int64,
    )

    for row_start in range(0, n_cells, chunk_rows):
        row_end = min(row_start + chunk_rows, n_cells)
        nz_start = int(indptr[row_start])
        nz_end = int(indptr[row_end])
        if nz_end == nz_start:
            continue
        indices = x_group["indices"][nz_start:nz_end]
        mask = np.isin(indices, sorted_targets)
        if not np.any(mask):
            continue
        data = x_group["data"][nz_start:nz_end]
        row_counts = np.diff(indptr[row_start : row_end + 1])
        row_indices = np.repeat(np.arange(row_start, row_end), row_counts)
        target_cols = indices[mask]
        insertion = np.searchsorted(sorted_targets, target_cols)
        target_positions = sorted_positions[insertion]
        values[row_indices[mask], target_positions] = data[mask]
    return values


def extract_csc_values(
    x_group: h5py.Group,
    target_indices: list[int],
) -> np.ndarray:
    shape = tuple(int(value) for value in x_group.attrs["shape"])
    n_cells = shape[0]
    values = np.zeros((n_cells, len(target_indices)), dtype=np.float32)
    indptr = x_group["indptr"][:]
    for target_pos, gene_idx in enumerate(target_indices):
        nz_start = int(indptr[gene_idx])
        nz_end = int(indptr[gene_idx + 1])
        if nz_end == nz_start:
            continue
        row_indices = x_group["indices"][nz_start:nz_end]
        values[row_indices, target_pos] = x_group["data"][nz_start:nz_end]
    return values


def extract_values(
    x_group: h5py.Group,
    target_indices: list[int],
    chunk_rows: int,
) -> np.ndarray:
    encoding = decode_scalar(x_group.attrs.get("encoding-type", ""))
    if encoding == "csr_matrix":
        return extract_csr_values(x_group, target_indices, chunk_rows)
    if encoding == "csc_matrix":
        return extract_csc_values(x_group, target_indices)
    raise ValueError(f"unsupported X encoding-type: {encoding}")


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def load_external_metadata(
    path: Path,
    key: str,
    columns: list[str],
    cell_ids: list[str],
) -> tuple[dict[str, list[str]], int]:
    wanted = set(cell_ids)
    by_cell: dict[str, dict[str, str]] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or key not in reader.fieldnames:
            raise ValueError(f"{key!r} not found in metadata header")
        available_columns = [
            column
            for column in columns
            if column in reader.fieldnames and column != key
        ]
        for row in reader:
            cell_id = row.get(key, "")
            if cell_id in wanted:
                by_cell[cell_id] = {
                    column: row.get(column, "")
                    for column in available_columns
                }

    values_by_column = {
        column: [
            by_cell.get(cell_id, {}).get(column, "")
            for cell_id in cell_ids
        ]
        for column in available_columns
    }
    return values_by_column, len(by_cell)


def write_cell_values(
    path: Path,
    dataset: str,
    cell_ids: list[str],
    gene_symbols: list[str],
    values: np.ndarray,
    normalization: str,
    detection_threshold: float,
    obs_values: dict[str, list[str]],
) -> None:
    fieldnames = CELL_VALUE_FIELDS + list(obs_values)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for cell_idx, cell_id in enumerate(cell_ids):
            obs_row = {
                column: obs_values[column][cell_idx]
                for column in obs_values
            }
            for gene_idx, gene_symbol in enumerate(gene_symbols):
                expression = float(values[cell_idx, gene_idx])
                writer.writerow(
                    {
                        "dataset": dataset,
                        "cell_id": cell_id,
                        "gene_symbol": gene_symbol,
                        "expression": f"{expression:.6g}",
                        "detected": "1" if expression > detection_threshold else "0",
                        "normalization": normalization,
                        **obs_row,
                    }
                )


def summarize_group(
    dataset: str,
    grouping: str,
    group: str,
    cell_indices: np.ndarray,
    gene_symbols: list[str],
    values: np.ndarray,
    normalization: str,
    detection_threshold: float,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    group_values = values[cell_indices, :]
    for gene_idx, gene_symbol in enumerate(gene_symbols):
        expression = group_values[:, gene_idx]
        detected = expression > detection_threshold
        nonzero = expression[detected]
        rows.append(
            {
                "dataset": dataset,
                "grouping": grouping,
                "group": group,
                "gene_symbol": gene_symbol,
                "normalization": normalization,
                "n_cells": str(expression.size),
                "detected_n": str(int(detected.sum())),
                "detected_fraction": f"{float(detected.mean()):.6g}",
                "mean_expression": f"{float(expression.mean()):.6g}",
                "median_expression": f"{float(np.median(expression)):.6g}",
                "max_expression": f"{float(expression.max()):.6g}",
                "nonzero_mean_expression": f"{float(nonzero.mean()):.6g}"
                if nonzero.size
                else "0",
            }
        )
    return rows


def summarize_values(
    dataset: str,
    gene_symbols: list[str],
    values: np.ndarray,
    normalization: str,
    detection_threshold: float,
    obs_values: dict[str, list[str]],
    group_cols: list[str],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    all_indices = np.arange(values.shape[0])
    rows.extend(
        summarize_group(
            dataset,
            "overall",
            "all",
            all_indices,
            gene_symbols,
            values,
            normalization,
            detection_threshold,
        )
    )

    for group_col in group_cols:
        if group_col not in obs_values:
            continue
        group_to_indices: dict[str, list[int]] = defaultdict(list)
        for idx, group_value in enumerate(obs_values[group_col]):
            group_to_indices[group_value or "NA"].append(idx)
        for group_value in sorted(group_to_indices):
            indices = np.array(group_to_indices[group_value], dtype=np.int64)
            rows.extend(
                summarize_group(
                    dataset,
                    group_col,
                    group_value,
                    indices,
                    gene_symbols,
                    values,
                    normalization,
                    detection_threshold,
                )
            )
    rows.sort(key=lambda row: (row["grouping"], row["group"], row["gene_symbol"]))
    return rows


def main() -> int:
    args = parse_args()
    h5ad_path = Path(args.h5ad)
    dataset = args.dataset or h5ad_path.stem
    outdir = Path(args.outdir) / dataset
    outdir.mkdir(parents=True, exist_ok=True)

    targets = load_targets(Path(args.target_genes), args.species)
    group_cols = [
        column.strip()
        for column in args.group_cols.split(",")
        if column.strip()
    ]
    metadata_columns = [
        column.strip()
        for column in args.cell_metadata_columns.split(",")
        if column.strip()
    ]
    external_metadata_matches = 0

    with h5py.File(h5ad_path, "r") as handle:
        cell_ids = read_obs_column(handle["obs"], "cell_label")
        if cell_ids is None:
            raise ValueError("obs/cell_label is required but was not found")
        obs_values: dict[str, list[str]] = {}
        for column in group_cols:
            values = read_obs_column(handle["obs"], column)
            if values is not None:
                obs_values[column] = values

        if args.cell_metadata:
            external_values, external_metadata_matches = load_external_metadata(
                Path(args.cell_metadata),
                args.cell_metadata_key,
                metadata_columns,
                cell_ids,
            )
            obs_values.update(external_values)

        normalization = (
            decode_dataset_scalar(handle["uns"]["normalization"])
            if "normalization" in handle["uns"]
            else ""
        )
        gene_symbols, target_indices, match_rows = match_targets(handle["var"], targets)
        expression_values = extract_values(handle["X"], target_indices, args.chunk_rows)
        matrix_shape = tuple(int(value) for value in handle["X"].attrs["shape"])
        x_encoding = decode_scalar(handle["X"].attrs.get("encoding-type", ""))

    write_tsv(
        outdir / f"{dataset}_target_gene_matches.tsv",
        [
            "target_gene_symbol",
            "matched_gene_symbol",
            "gene_identifier",
            "gene_index",
            "matched_by",
            "present",
        ],
        match_rows,
    )
    summary_rows = summarize_values(
        dataset,
        gene_symbols,
        expression_values,
        normalization,
        args.detection_threshold,
        obs_values,
        group_cols,
    )
    write_tsv(
        outdir / f"{dataset}_target_gene_summary.tsv",
        SUMMARY_FIELDS,
        summary_rows,
    )
    if not args.skip_cell_values:
        write_cell_values(
            outdir / f"{dataset}_target_gene_cell_values.tsv",
            dataset,
            cell_ids,
            gene_symbols,
            expression_values,
            normalization,
            args.detection_threshold,
            obs_values,
        )

    qc_rows = [
        {"metric": "h5ad_path", "value": str(h5ad_path)},
        {"metric": "dataset", "value": dataset},
        {"metric": "species", "value": args.species},
        {"metric": "x_encoding", "value": x_encoding},
        {"metric": "normalization", "value": normalization},
        {"metric": "n_cells", "value": str(matrix_shape[0])},
        {"metric": "n_genes", "value": str(matrix_shape[1])},
        {"metric": "target_genes_requested", "value": str(len(targets))},
        {"metric": "target_genes_matched", "value": str(len(gene_symbols))},
        {"metric": "group_cols_available", "value": "|".join(obs_values)},
        {"metric": "cell_metadata_path", "value": args.cell_metadata or ""},
        {"metric": "external_metadata_matches", "value": str(external_metadata_matches)},
    ]
    write_tsv(outdir / f"{dataset}_target_gene_qc.tsv", QC_FIELDS, qc_rows)

    print(f"Wrote Allen ABC target-gene summaries to {outdir}")
    print(f"Matrix: {matrix_shape[0]} cells x {matrix_shape[1]} genes ({x_encoding})")
    print(f"Matched target genes: {', '.join(gene_symbols) if gene_symbols else 'none'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
