#!/usr/bin/env python3
"""Scan processed expression matrices for DNMT3/STAT target genes.

The script is intentionally format-tolerant: point it at a file or directory
of tab/comma-delimited expression matrices and it will search likely
annotation columns for target gene symbols or Ensembl IDs.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import itertools
import math
import re
import sys
from pathlib import Path
from typing import Iterable


TEXT_SUFFIXES = {".txt", ".tsv", ".csv", ".tab", ".gz"}
ANNOTATION_HINTS = (
    "gene",
    "symbol",
    "ensembl",
    "id",
    "name",
    "feature",
    "description",
    "chrom",
    "start",
    "end",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", required=True, help="Dataset label for output files, for example GSE157080.")
    parser.add_argument(
        "--input",
        action="append",
        required=True,
        help="Expression matrix file or directory. Can be supplied multiple times.",
    )
    parser.add_argument("--targets", default="config/target_genes.tsv", help="Target gene TSV.")
    parser.add_argument("--species", choices=["all", "human", "mouse"], default="all", help="Target species filter.")
    parser.add_argument("--outdir", default="results/target_gene_expression", help="Output directory.")
    parser.add_argument(
        "--annotation-columns",
        type=int,
        default=8,
        help="First N columns to search as annotation fields when headers are ambiguous.",
    )
    parser.add_argument(
        "--value-column-regex",
        default=r"(_FPKM$|_TPM$|^TPM$)",
        help="Regex for normalized expression value columns. Falls back to numeric non-annotation columns if no columns match.",
    )
    parser.add_argument("--max-file-mb", type=float, default=750.0, help="Skip input files larger than this.")
    return parser.parse_args()


def normalize_token(value: str) -> str:
    value = value.strip().strip('"').strip("'")
    value = re.sub(r"\.\d+$", "", value)
    return value.upper()


def split_aliases(value: str) -> list[str]:
    aliases = []
    for token in re.split(r"[|;,]\s*|\s+", value.strip()):
        if token:
            aliases.append(token)
    return aliases


def read_targets(path: Path, species: str) -> tuple[dict[str, list[dict[str, str]]], list[dict[str, str]]]:
    alias_map: dict[str, list[dict[str, str]]] = {}
    selected: list[dict[str, str]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if species != "all" and row["species"] != species:
                continue
            selected.append(row)
            aliases = [row["gene_symbol"], *split_aliases(row.get("aliases", "")), *split_aliases(row.get("ensembl_hint", ""))]
            for alias in aliases:
                alias_map.setdefault(normalize_token(alias), []).append(row)
    return alias_map, selected


def iter_input_files(paths: list[str], max_file_mb: float) -> Iterable[Path]:
    max_bytes = max_file_mb * 1024 * 1024
    for value in paths:
        path = Path(value)
        if path.is_dir():
            for item in sorted(path.rglob("*")):
                if item.is_file() and is_text_table(item) and item.stat().st_size <= max_bytes:
                    yield item
        elif path.is_file() and is_text_table(path) and path.stat().st_size <= max_bytes:
            yield path


def is_text_table(path: Path) -> bool:
    if path.suffix.lower() == ".gz":
        return path.with_suffix("").suffix.lower() in TEXT_SUFFIXES
    return path.suffix.lower() in TEXT_SUFFIXES


def open_text(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", newline="", errors="replace")
    return path.open(newline="", errors="replace")


def sniff_delimiter(sample: str, path: Path) -> str:
    if path.suffix.lower() == ".csv" or path.with_suffix("").suffix.lower() == ".csv":
        return ","
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        return dialect.delimiter
    except csv.Error:
        return "\t"


def header_annotation_columns(header: list[str], fallback_n: int) -> set[int]:
    indices = set(range(min(fallback_n, len(header))))
    for idx, name in enumerate(header):
        low = name.lower()
        if any(hint in low for hint in ANNOTATION_HINTS):
            indices.add(idx)
    return indices


def cell_tokens(value: str) -> list[str]:
    pieces = re.split(r"[|;,()\[\]\s]+", value)
    return [normalize_token(piece) for piece in pieces if piece.strip()]


def find_targets(row: list[str], annotation_indices: set[int], alias_map: dict[str, list[dict[str, str]]]) -> list[dict[str, str]]:
    matches = []
    seen = set()
    for idx in sorted(annotation_indices):
        if idx >= len(row):
            continue
        for token in cell_tokens(row[idx]):
            for target in alias_map.get(token, []):
                key = (target["species"], target["gene_symbol"])
                if key in seen:
                    continue
                seen.add(key)
                matches.append(
                    {
                        "species": target["species"],
                        "gene_symbol": target["gene_symbol"],
                        "matched_alias": token,
                        "matched_cell": row[idx],
                        "matched_column_index": str(idx),
                    }
                )
    return matches


def expression_indices(header: list[str], value_column_regex: str) -> set[int]:
    pattern = re.compile(value_column_regex, flags=re.IGNORECASE)
    return {idx for idx, name in enumerate(header) if pattern.search(name)}


def is_number(value: str) -> bool:
    try:
        number = float(value)
    except ValueError:
        return False
    return math.isfinite(number)


def looks_like_no_header_count_row(row: list[str]) -> bool:
    if len(row) < 2:
        return False
    if any(hint in row[0].lower() for hint in ("gene", "symbol", "ensembl", "id")):
        return False
    return is_number(row[1])


def synthetic_count_header(row: list[str]) -> tuple[list[str], set[int], set[int], bool]:
    if len(row) == 4 and row[0].startswith("N_"):
        header = ["gene_id", "unstranded_count", "sense_count", "antisense_count"]
        return header, {0}, {1, 2, 3}, True
    if looks_like_no_header_count_row(row):
        header = ["gene_id", *[f"value_{idx}" for idx in range(1, len(row))]]
        return header, {0}, set(range(1, len(row))), True
    return row, set(), set(), False


def numeric_summary(row: list[str], value_indices: set[int], fallback_excluded_indices: set[int]) -> tuple[int, str, str, str]:
    values = []
    indices = value_indices or set(range(len(row)))
    for idx in sorted(indices):
        if idx >= len(row):
            continue
        if not value_indices and idx in fallback_excluded_indices:
            continue
        value = row[idx]
        try:
            number = float(value)
        except ValueError:
            continue
        if math.isfinite(number):
            values.append(number)
    if not values:
        return 0, "", "", ""
    mean_value = sum(values) / len(values)
    return len(values), f"{mean_value:.6g}", f"{min(values):.6g}", f"{max(values):.6g}"


def safe_name(path: Path) -> str:
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", path.name)
    if len(name) > 140:
        name = name[:140]
    return name


def scan_file(
    path: Path,
    dataset: str,
    alias_map: dict[str, list[dict[str, str]]],
    annotation_columns: int,
    value_column_regex: str,
    outdir: Path,
) -> list[dict[str, str]]:
    with open_text(path) as handle:
        sample = handle.read(8192)
        if not sample:
            return []
        handle.seek(0)
        delimiter = sniff_delimiter(sample, path)
        reader = csv.reader(handle, delimiter=delimiter)
        header = next(reader, None)
        if not header:
            return []
        first_row = header
        header, synthetic_annotation_indices, synthetic_value_indices, header_is_data = synthetic_count_header(first_row)
        annotation_indices = synthetic_annotation_indices or header_annotation_columns(header, annotation_columns)
        value_indices = synthetic_value_indices or expression_indices(header, value_column_regex)
        value_column_names = "|".join(header[idx] for idx in sorted(value_indices))
        target_rows = []
        summary_rows = []
        output_header = [
            "dataset",
            "source_file",
            "row_number",
            "matched_species",
            "matched_gene",
            "matched_alias",
            "matched_cell",
            *header,
        ]

        row_iter = itertools.chain([first_row], reader) if header_is_data else reader
        first_row_number = 1 if header_is_data else 2
        for row_number, row in enumerate(row_iter, start=first_row_number):
            matches = find_targets(row, annotation_indices, alias_map)
            if not matches:
                continue
            numeric_n, numeric_mean, numeric_min, numeric_max = numeric_summary(row, value_indices, annotation_indices)
            for match in matches:
                target_rows.append(
                    {
                        "dataset": dataset,
                        "source_file": str(path),
                        "row_number": str(row_number),
                        "matched_species": match["species"],
                        "matched_gene": match["gene_symbol"],
                        "matched_alias": match["matched_alias"],
                        "matched_cell": match["matched_cell"],
                        "_row_values": row,
                    }
                )
                summary_rows.append(
                    {
                        "dataset": dataset,
                        "source_file": str(path),
                        "row_number": str(row_number),
                        "species": match["species"],
                        "gene_symbol": match["gene_symbol"],
                        "matched_alias": match["matched_alias"],
                        "numeric_value_count": str(numeric_n),
                        "numeric_mean": numeric_mean,
                        "numeric_min": numeric_min,
                        "numeric_max": numeric_max,
                        "value_columns": value_column_names,
                    }
                )

    if target_rows:
        out_path = outdir / f"{dataset}_{safe_name(path)}_target_rows.tsv"
        with out_path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(output_header)
            for item in target_rows:
                writer.writerow(
                    [
                        item["dataset"],
                        item["source_file"],
                        item["row_number"],
                        item["matched_species"],
                        item["matched_gene"],
                        item["matched_alias"],
                        item["matched_cell"],
                        *item["_row_values"],
                    ]
                )
    return summary_rows


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def best_rows_by_gene(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    best: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        key = (row["species"], row["gene_symbol"])
        try:
            value = float(row["numeric_mean"]) if row["numeric_mean"] else float("-inf")
        except ValueError:
            value = float("-inf")
        current = best.get(key)
        if current is None:
            best[key] = row
            continue
        try:
            current_value = float(current["numeric_mean"]) if current["numeric_mean"] else float("-inf")
        except ValueError:
            current_value = float("-inf")
        if value > current_value:
            best[key] = row
    return [best[key] for key in sorted(best)]


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    alias_map, targets = read_targets(Path(args.targets), args.species)

    all_summary_rows: list[dict[str, str]] = []
    files = list(iter_input_files(args.input, args.max_file_mb))
    if not files:
        print("No input expression tables found.", file=sys.stderr)
        return 1

    for path in files:
        print(f"Scanning {path}", file=sys.stderr)
        try:
            all_summary_rows.extend(
                scan_file(path, args.dataset, alias_map, args.annotation_columns, args.value_column_regex, outdir)
            )
        except UnicodeDecodeError:
            print(f"Skipping non-text file: {path}", file=sys.stderr)
        except csv.Error as exc:
            print(f"Skipping malformed table {path}: {exc}", file=sys.stderr)

    summary_path = outdir / f"{args.dataset}_target_expression_summary.tsv"
    fields = [
        "dataset",
        "source_file",
        "row_number",
        "species",
        "gene_symbol",
        "matched_alias",
        "numeric_value_count",
        "numeric_mean",
        "numeric_min",
        "numeric_max",
        "value_columns",
    ]
    write_tsv(summary_path, all_summary_rows, fields)
    write_tsv(outdir / f"{args.dataset}_target_gene_best_rows.tsv", best_rows_by_gene(all_summary_rows), fields)

    found = {(row["species"], row["gene_symbol"]) for row in all_summary_rows}
    presence_rows = [
        {
            "dataset": args.dataset,
            "species": target["species"],
            "gene_symbol": target["gene_symbol"],
            "found": "yes" if (target["species"], target["gene_symbol"]) in found else "no",
        }
        for target in targets
    ]
    write_tsv(outdir / f"{args.dataset}_target_gene_presence.tsv", presence_rows, ["dataset", "species", "gene_symbol", "found"])
    print(f"Wrote {len(all_summary_rows)} target rows to {summary_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
