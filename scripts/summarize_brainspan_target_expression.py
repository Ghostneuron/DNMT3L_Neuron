#!/usr/bin/env python3
"""Summarize BrainSpan target-gene expression by age and brain region."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean, median


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        default="data/allen/brainspan/RNASeq_Gencode_v10/brainspan_target_gene_rpkm.tsv",
        help="Target-gene BrainSpan RPKM table from download_allen_brainspan_targets.py.",
    )
    parser.add_argument(
        "--outdir",
        default="results/allen_brainspan",
        help="Output directory for summary tables.",
    )
    return parser.parse_args()


def as_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def summarize(values: list[float]) -> dict[str, str]:
    detected = [value for value in values if value > 0]
    return {
        "n": str(len(values)),
        "mean_rpkm": f"{mean(values):.6g}",
        "median_rpkm": f"{median(values):.6g}",
        "min_rpkm": f"{min(values):.6g}",
        "max_rpkm": f"{max(values):.6g}",
        "detected_n": str(len(detected)),
        "detected_fraction": f"{len(detected) / len(values):.6g}",
    }


def life_stage(row: dict[str, str]) -> str:
    tags = row["age_tags"].lower()
    if row["age_embryonic"] == "True":
        return "prenatal"
    if "adult" in tags:
        return "adult"
    if "adolescent" in tags or "adolesence" in tags:
        return "adolescent"
    if "childhood" in tags:
        return "childhood"
    if "infant" in tags:
        return "infancy"
    return "postnatal_other"


def write_group_summary(
    rows: list[dict[str, str]],
    group_fields: list[str],
    path: Path,
    sort_fields: list[str] | None = None,
) -> None:
    groups: dict[tuple[str, ...], list[float]] = defaultdict(list)
    metadata: dict[tuple[str, ...], dict[str, str]] = {}
    for row in rows:
        key = tuple(row[field] for field in group_fields)
        groups[key].append(as_float(row["rpkm"]))
        metadata.setdefault(key, {field: row[field] for field in group_fields})

    output_rows: list[dict[str, str]] = []
    for key, values in groups.items():
        out = dict(metadata[key])
        out.update(summarize(values))
        output_rows.append(out)

    if sort_fields:
        output_rows.sort(key=lambda row: tuple(row[field] for field in sort_fields))
    else:
        output_rows.sort(key=lambda row: tuple(row[field] for field in group_fields))

    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=group_fields
            + [
                "n",
                "mean_rpkm",
                "median_rpkm",
                "min_rpkm",
                "max_rpkm",
                "detected_n",
                "detected_fraction",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(output_rows)


def write_prenatal_age_adjusted_region_summary(
    rows: list[dict[str, str]],
    path: Path,
    pseudocount: float = 0.1,
) -> None:
    prenatal_rows = [row for row in rows if row["life_stage"] == "prenatal"]
    output_rows: list[dict[str, str]] = []
    for gene_symbol in sorted({row["gene_symbol"] for row in prenatal_rows}):
        gene_rows = [row for row in prenatal_rows if row["gene_symbol"] == gene_symbol]
        xs = [as_float(row["age_days"]) for row in gene_rows]
        ys = [math.log2(as_float(row["rpkm"]) + pseudocount) for row in gene_rows]
        mean_x = mean(xs)
        mean_y = mean(ys)
        denominator = sum((x - mean_x) ** 2 for x in xs)
        slope = (
            sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys)) / denominator
            if denominator
            else 0.0
        )
        intercept = mean_y - slope * mean_x

        by_region: dict[tuple[str, str], list[tuple[float, float]]] = defaultdict(list)
        for row, x, y in zip(gene_rows, xs, ys):
            expected = intercept + slope * x
            key = (row["structure_acronym"], row["structure_name"])
            by_region[key].append((y - expected, as_float(row["rpkm"])))

        for (structure_acronym, structure_name), values in by_region.items():
            residuals = [value[0] for value in values]
            raw_values = [value[1] for value in values]
            output_rows.append(
                {
                    "gene_symbol": gene_symbol,
                    "structure_acronym": structure_acronym,
                    "structure_name": structure_name,
                    "n": str(len(values)),
                    "median_age_adjusted_log2_residual": f"{median(residuals):.6g}",
                    "mean_age_adjusted_log2_residual": f"{mean(residuals):.6g}",
                    "median_raw_rpkm": f"{median(raw_values):.6g}",
                    "mean_raw_rpkm": f"{mean(raw_values):.6g}",
                    "slope_log2_per_week": f"{slope * 7:.6g}",
                }
            )

    output_rows.sort(
        key=lambda row: (
            row["gene_symbol"],
            -float(row["median_age_adjusted_log2_residual"]),
            row["structure_acronym"],
        )
    )
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "gene_symbol",
                "structure_acronym",
                "structure_name",
                "n",
                "median_age_adjusted_log2_residual",
                "mean_age_adjusted_log2_residual",
                "median_raw_rpkm",
                "mean_raw_rpkm",
                "slope_log2_per_week",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(output_rows)


def main() -> int:
    args = parse_args()
    input_path = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with input_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    for row in rows:
        row["developmental_period"] = (
            "prenatal" if row["age_embryonic"] == "True" else "postnatal"
        )
        row["life_stage"] = life_stage(row)
        row["age_sort_days"] = f"{as_float(row['age_days']):012.3f}"

    write_group_summary(
        rows,
        ["gene_symbol"],
        outdir / "brainspan_target_gene_overall_summary.tsv",
    )
    write_group_summary(
        rows,
        ["gene_symbol", "developmental_period"],
        outdir / "brainspan_target_gene_period_summary.tsv",
    )
    write_group_summary(
        rows,
        ["gene_symbol", "life_stage"],
        outdir / "brainspan_target_gene_life_stage_summary.tsv",
    )
    write_group_summary(
        rows,
        ["gene_symbol", "age_sort_days", "age_name", "developmental_period"],
        outdir / "brainspan_target_gene_age_summary.tsv",
        sort_fields=["gene_symbol", "age_sort_days"],
    )
    write_group_summary(
        rows,
        ["gene_symbol", "structure_acronym", "structure_name"],
        outdir / "brainspan_target_gene_region_summary.tsv",
    )

    adult_rows = [row for row in rows if row["life_stage"] == "adult"]
    if adult_rows:
        write_group_summary(
            adult_rows,
            ["gene_symbol"],
            outdir / "brainspan_adult_target_gene_overall_summary.tsv",
        )
        write_group_summary(
            adult_rows,
            ["gene_symbol", "structure_acronym", "structure_name"],
            outdir / "brainspan_adult_target_gene_region_summary.tsv",
        )

    prenatal_rows = [row for row in rows if row["life_stage"] == "prenatal"]
    if prenatal_rows:
        write_group_summary(
            prenatal_rows,
            ["gene_symbol"],
            outdir / "brainspan_prenatal_target_gene_overall_summary.tsv",
        )
        write_group_summary(
            prenatal_rows,
            ["gene_symbol", "structure_acronym", "structure_name"],
            outdir / "brainspan_prenatal_target_gene_region_summary.tsv",
        )
        write_prenatal_age_adjusted_region_summary(
            rows,
            outdir / "brainspan_prenatal_age_adjusted_region_summary.tsv",
        )

    print(f"Read {len(rows)} target-gene rows from {input_path}")
    print(f"Wrote summaries to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
