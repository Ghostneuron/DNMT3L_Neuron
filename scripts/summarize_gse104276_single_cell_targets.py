#!/usr/bin/env python3
"""Summarize GSE104276 fetal PFC single-cell target-gene expression."""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path
from statistics import mean, median

from openpyxl import load_workbook


TARGET_FIELDS = [
    "cell_id",
    "gene_symbol",
    "tpm",
    "detected",
    "cell_type",
    "week",
    "pfc",
    "gene_num",
    "pre_map_reads",
    "aligned_reads",
    "mapping_rate",
]

SUMMARY_FIELDS = [
    "group",
    "gene_symbol",
    "n_cells",
    "detected_n",
    "detected_fraction",
    "mean_tpm",
    "median_tpm",
    "max_tpm",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--matrix",
        default="data/geo/GSE104276/supplementary/GSE104276_all_pfc_2394_UMI_TPM_NOERCC.xls.gz",
        help="GSE104276 all-cell TPM matrix.",
    )
    parser.add_argument(
        "--metadata",
        default="data/geo/GSE104276/supplementary/GSE104276_readme_sample_barcode.xlsx",
        help="GSE104276 sample/barcode workbook.",
    )
    parser.add_argument(
        "--target-genes",
        default="config/target_genes.tsv",
        help="Project target gene manifest.",
    )
    parser.add_argument(
        "--outdir",
        default="results/single_cell/GSE104276",
        help="Output directory.",
    )
    return parser.parse_args()


def load_human_targets(path: Path) -> set[str]:
    targets: set[str] = set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["species"].lower() == "human":
                targets.add(row["gene_symbol"].upper())
                for alias in row.get("aliases", "").split("|"):
                    if alias:
                        targets.add(alias.upper())
    return targets


def load_sample_info(path: Path) -> dict[str, dict[str, str]]:
    workbook = load_workbook(path, read_only=True, data_only=True)
    worksheet = workbook["SampleInfo"]
    metadata: dict[str, dict[str, str]] = {}
    header = [str(value) if value is not None else "cell_id" for value in next(worksheet.iter_rows(values_only=True))]
    header[0] = "cell_id"
    for raw_row in worksheet.iter_rows(values_only=True):
        values = list(raw_row)
        if not values or values[0] is None:
            continue
        row = {key: "" if value is None else str(value) for key, value in zip(header, values)}
        cell_id = row["cell_id"]
        metadata[cell_id] = {
            "cell_type": row.get("cell_types", ""),
            "week": row.get("week", ""),
            "pfc": row.get("pfc", ""),
            "gene_num": row.get("Gene_num", ""),
            "pre_map_reads": row.get("Pre_Map_Reads", ""),
            "aligned_reads": row.get("Aligned_Reads", ""),
            "mapping_rate": row.get("MappingRate", ""),
        }
    return metadata


def extract_target_rows(
    matrix_path: Path,
    target_genes: set[str],
    metadata: dict[str, dict[str, str]],
) -> tuple[list[dict[str, str]], list[str], set[str]]:
    target_rows: list[dict[str, str]] = []
    found_genes: set[str] = set()
    with gzip.open(matrix_path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cells = header
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            gene_symbol = fields[0]
            if gene_symbol.upper() not in target_genes:
                continue
            found_genes.add(gene_symbol.upper())
            values = fields[1:]
            if len(values) != len(cells):
                raise ValueError(
                    f"row for {gene_symbol} has {len(values)} values but {len(cells)} cells"
                )
            for cell_id, value in zip(cells, values):
                meta = metadata.get(cell_id, {})
                tpm = float(value)
                target_rows.append(
                    {
                        "cell_id": cell_id,
                        "gene_symbol": gene_symbol,
                        "tpm": f"{tpm:.6g}",
                        "detected": "1" if tpm > 0 else "0",
                        "cell_type": meta.get("cell_type", ""),
                        "week": meta.get("week", ""),
                        "pfc": meta.get("pfc", ""),
                        "gene_num": meta.get("gene_num", ""),
                        "pre_map_reads": meta.get("pre_map_reads", ""),
                        "aligned_reads": meta.get("aligned_reads", ""),
                        "mapping_rate": meta.get("mapping_rate", ""),
                    }
                )
    return target_rows, cells, found_genes


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def summarize(rows: list[dict[str, str]], group_fields: list[str]) -> list[dict[str, str]]:
    groups: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in rows:
        group = " | ".join(row[field] or "NA" for field in group_fields) if group_fields else "all"
        groups[(group, row["gene_symbol"])].append(float(row["tpm"]))

    output: list[dict[str, str]] = []
    for (group, gene_symbol), values in groups.items():
        detected = [value for value in values if value > 0]
        output.append(
            {
                "group": group,
                "gene_symbol": gene_symbol,
                "n_cells": str(len(values)),
                "detected_n": str(len(detected)),
                "detected_fraction": f"{len(detected) / len(values):.6g}",
                "mean_tpm": f"{mean(values):.6g}",
                "median_tpm": f"{median(values):.6g}",
                "max_tpm": f"{max(values):.6g}",
            }
        )
    output.sort(key=lambda row: (row["group"], row["gene_symbol"]))
    return output


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    target_genes = load_human_targets(Path(args.target_genes))
    metadata = load_sample_info(Path(args.metadata))
    rows, cells, found_genes = extract_target_rows(Path(args.matrix), target_genes, metadata)

    write_tsv(outdir / "GSE104276_target_gene_cell_values.tsv", TARGET_FIELDS, rows)
    write_tsv(
        outdir / "GSE104276_target_gene_overall_summary.tsv",
        SUMMARY_FIELDS,
        summarize(rows, []),
    )
    write_tsv(
        outdir / "GSE104276_target_gene_by_cell_type.tsv",
        SUMMARY_FIELDS,
        summarize(rows, ["cell_type"]),
    )
    write_tsv(
        outdir / "GSE104276_target_gene_by_week.tsv",
        SUMMARY_FIELDS,
        summarize(rows, ["week"]),
    )
    write_tsv(
        outdir / "GSE104276_target_gene_by_cell_type_week.tsv",
        SUMMARY_FIELDS,
        summarize(rows, ["cell_type", "week"]),
    )

    metadata_cells = set(metadata)
    matrix_cells = set(cells)
    qc_rows = [
        {"metric": "matrix_cells", "value": str(len(matrix_cells))},
        {"metric": "metadata_cells", "value": str(len(metadata_cells))},
        {"metric": "cells_with_metadata", "value": str(len(matrix_cells & metadata_cells))},
        {"metric": "matrix_cells_without_metadata", "value": str(len(matrix_cells - metadata_cells))},
        {"metric": "target_genes_found", "value": "|".join(sorted(found_genes))},
    ]
    write_tsv(outdir / "GSE104276_target_gene_qc.tsv", ["metric", "value"], qc_rows)

    print(f"Wrote {len(rows)} target-gene cell rows to {outdir}")
    print(f"Matrix cells: {len(matrix_cells)}; cells with metadata: {len(matrix_cells & metadata_cells)}")
    print(f"Found target genes: {', '.join(sorted(found_genes))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
