#!/usr/bin/env python3
"""Extract the GSE157080 author count-based DE table from Table S3."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        default="data/geo/GSE157080/supplementary/1-s2.0-S0301008221000423-mmc5.xlsx",
    )
    parser.add_argument("--sheet", default="CON-VS-DNMT3L total")
    parser.add_argument(
        "--output",
        default="data/geo/GSE157080/supplementary/GSE157080_Table_S3_author_count_de.tsv",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    source = Path(args.input)
    output = Path(args.output)
    if not source.exists():
        raise FileNotFoundError(source)

    table = pd.read_excel(source, sheet_name=args.sheet, dtype=str)
    table = table.fillna("")
    table["source_workbook"] = source.name
    table["source_sheet"] = args.sheet
    table["analysis_note"] = (
        "Author supplied supplementary Table S3 CON vs DNMT3L differential-expression table; "
        "the accompanying report describes read-count input for differential analysis."
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL, lineterminator="\n")
    print(f"Wrote {len(table):,} rows to {output}")


if __name__ == "__main__":
    main()
