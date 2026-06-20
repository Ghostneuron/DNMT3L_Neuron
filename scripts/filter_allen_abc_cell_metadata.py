#!/usr/bin/env python3
"""Filter large Allen ABC cell metadata tables to the cells in one H5AD file."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

import h5py
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--h5ad",
        required=True,
        help="Allen ABC H5AD file containing obs/cell_label.",
    )
    parser.add_argument(
        "--metadata",
        required=True,
        help="Large Allen ABC cell metadata CSV.",
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Filtered metadata CSV/TSV output.",
    )
    parser.add_argument(
        "--key",
        default="cell_label",
        help="Metadata key column to match against H5AD obs/cell_label.",
    )
    parser.add_argument(
        "--qc-out",
        default=None,
        help="Optional TSV with filtering metrics.",
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


def load_h5ad_cell_labels(path: Path) -> list[str]:
    with h5py.File(path, "r") as handle:
        return decode_array(handle["obs"]["cell_label"])


def write_qc(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "value"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    h5ad_path = Path(args.h5ad)
    metadata_path = Path(args.metadata)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    cell_labels = load_h5ad_cell_labels(h5ad_path)
    wanted = set(cell_labels)
    matched: set[str] = set()
    n_metadata_rows = 0

    with metadata_path.open(newline="") as in_handle, out_path.open("w", newline="") as out_handle:
        reader = csv.DictReader(in_handle)
        if reader.fieldnames is None or args.key not in reader.fieldnames:
            raise ValueError(f"{args.key!r} not found in metadata header")
        writer = csv.DictWriter(out_handle, fieldnames=reader.fieldnames)
        writer.writeheader()
        for row in reader:
            n_metadata_rows += 1
            cell_label = row.get(args.key, "")
            if cell_label in wanted:
                writer.writerow(row)
                matched.add(cell_label)

    qc_rows = [
        {"metric": "h5ad_path", "value": str(h5ad_path)},
        {"metric": "metadata_path", "value": str(metadata_path)},
        {"metric": "out_path", "value": str(out_path)},
        {"metric": "h5ad_cells", "value": str(len(cell_labels))},
        {"metric": "metadata_rows_scanned", "value": str(n_metadata_rows)},
        {"metric": "matched_cells", "value": str(len(matched))},
        {"metric": "h5ad_cells_without_metadata", "value": str(len(wanted - matched))},
    ]
    if args.qc_out:
        write_qc(Path(args.qc_out), qc_rows)

    print(f"Filtered {len(matched)} metadata rows for {len(cell_labels)} H5AD cells.")
    print(f"Output: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
