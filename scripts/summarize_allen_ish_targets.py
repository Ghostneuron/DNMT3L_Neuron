#!/usr/bin/env python3
"""Summarize Allen ISH target-gene metadata and structure-level expression."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from statistics import mean, median


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--indir", default="data/allen/ish")
    parser.add_argument("--outdir", default="results/allen_ish")
    parser.add_argument("--top-n", type=int, default=10)
    return parser.parse_args()


def as_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def summarize_availability(indir: Path, outdir: Path) -> None:
    rows = read_tsv(indir / "allen_ish_target_availability.tsv")
    write_tsv(
        outdir / "allen_ish_target_availability_summary.tsv",
        rows[0].keys() if rows else [],
        rows,
    )


def summarize_experiments(indir: Path, outdir: Path) -> None:
    rows = read_tsv(indir / "allen_ish_section_datasets.tsv")
    groups: dict[tuple[str, str, str, str], list[dict[str, str]]] = {}
    for row in rows:
        key = (
            row["product_name"],
            row["species"],
            row["gene_symbol"],
            row["age_name"],
        )
        groups.setdefault(key, []).append(row)

    output: list[dict[str, str]] = []
    for (product_name, species, gene_symbol, age_name), values in groups.items():
        delegates = [row for row in values if row["delegate"].lower() == "true"]
        output.append(
            {
                "product_name": product_name,
                "species": species,
                "gene_symbol": gene_symbol,
                "age_name": age_name,
                "experiments": str(len(values)),
                "delegate_experiments": str(len(delegates)),
                "planes": "|".join(sorted({row["plane_of_section"] for row in values})),
            }
        )
    output.sort(key=lambda row: (row["product_name"], row["gene_symbol"], row["age_name"]))
    write_tsv(
        outdir / "allen_ish_experiment_age_summary.tsv",
        [
            "product_name",
            "species",
            "gene_symbol",
            "age_name",
            "experiments",
            "delegate_experiments",
            "planes",
        ],
        output,
    )


def summarize_top_structures(indir: Path, outdir: Path, top_n: int) -> None:
    rows = read_tsv(indir / "allen_ish_structure_unionizes.tsv")
    groups: dict[tuple[str, str, str, str, str], list[dict[str, str]]] = {}
    for row in rows:
        if row["delegate"].lower() != "true":
            continue
        key = (
            row["product_name"],
            row["gene_symbol"],
            row["age_name"],
            row["section_data_set_id"],
            row["plane_of_section"],
        )
        groups.setdefault(key, []).append(row)

    output: list[dict[str, str]] = []
    for key, values in groups.items():
        values.sort(key=lambda row: as_float(row["expression_energy"]), reverse=True)
        for rank, row in enumerate(values[:top_n], start=1):
            output.append(
                {
                    "product_name": row["product_name"],
                    "species": row["species"],
                    "gene_symbol": row["gene_symbol"],
                    "age_name": row["age_name"],
                    "age_days": row["age_days"],
                    "age_embryonic": row["age_embryonic"],
                    "section_data_set_id": row["section_data_set_id"],
                    "plane_of_section": row["plane_of_section"],
                    "rank": str(rank),
                    "structure_acronym": row["structure_acronym"],
                    "structure_name": row["structure_name"],
                    "expression_energy": row["expression_energy"],
                    "expression_density": row["expression_density"],
                    "voxel_energy_mean": row["voxel_energy_mean"],
                }
            )
    output.sort(
        key=lambda row: (
            row["product_name"],
            row["gene_symbol"],
            as_float(row["age_days"]),
            row["section_data_set_id"],
            int(row["rank"]),
        )
    )
    write_tsv(
        outdir / "allen_ish_delegate_top_structures.tsv",
        [
            "product_name",
            "species",
            "gene_symbol",
            "age_name",
            "age_days",
            "age_embryonic",
            "section_data_set_id",
            "plane_of_section",
            "rank",
            "structure_acronym",
            "structure_name",
            "expression_energy",
            "expression_density",
            "voxel_energy_mean",
        ],
        output,
    )


def summarize_dataset_energy(indir: Path, outdir: Path) -> None:
    rows = read_tsv(indir / "allen_ish_structure_unionizes.tsv")
    groups: dict[tuple[str, str, str, str, str], list[float]] = {}
    metadata: dict[tuple[str, str, str, str, str], dict[str, str]] = {}
    for row in rows:
        if row["delegate"].lower() != "true":
            continue
        key = (
            row["product_name"],
            row["gene_symbol"],
            row["age_name"],
            row["age_days"],
            row["section_data_set_id"],
        )
        groups.setdefault(key, []).append(as_float(row["expression_energy"]))
        metadata.setdefault(key, row)

    output: list[dict[str, str]] = []
    for key, values in groups.items():
        row = metadata[key]
        detected = [value for value in values if value > 0]
        output.append(
            {
                "product_name": row["product_name"],
                "species": row["species"],
                "gene_symbol": row["gene_symbol"],
                "age_name": row["age_name"],
                "age_days": row["age_days"],
                "age_embryonic": row["age_embryonic"],
                "section_data_set_id": row["section_data_set_id"],
                "structures": str(len(values)),
                "detected_structures": str(len(detected)),
                "median_expression_energy": f"{median(values):.6g}",
                "mean_expression_energy": f"{mean(values):.6g}",
                "max_expression_energy": f"{max(values):.6g}",
            }
        )
    output.sort(
        key=lambda row: (row["product_name"], row["gene_symbol"], as_float(row["age_days"]))
    )
    write_tsv(
        outdir / "allen_ish_delegate_dataset_energy_summary.tsv",
        [
            "product_name",
            "species",
            "gene_symbol",
            "age_name",
            "age_days",
            "age_embryonic",
            "section_data_set_id",
            "structures",
            "detected_structures",
            "median_expression_energy",
            "mean_expression_energy",
            "max_expression_energy",
        ],
        output,
    )


def main() -> int:
    args = parse_args()
    indir = Path(args.indir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    summarize_availability(indir, outdir)
    summarize_experiments(indir, outdir)
    summarize_top_structures(indir, outdir, args.top_n)
    summarize_dataset_energy(indir, outdir)
    print(f"Wrote Allen ISH summaries to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
