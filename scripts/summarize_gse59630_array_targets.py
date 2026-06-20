#!/usr/bin/env python3
"""Summarize target genes in GSE59630 developmental Down syndrome brain array data."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
import statistics
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-genes", default="config/target_genes.tsv")
    parser.add_argument("--metadata", default="data/geo/GSE59630/metadata/sample_metadata.tsv")
    parser.add_argument("--series-matrix", default="data/geo/GSE59630/series_matrix/GSE59630_series_matrix.txt.gz")
    parser.add_argument("--platform-table", default="data/geo/GPL5175/annotation/GPL5175_table.txt")
    parser.add_argument("--outdir", default="results/ds_fetal_expression")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def parse_characteristics(text: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for part in text.split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
            values[key.strip().lower()] = value.strip()
    return values


def load_human_targets(path: Path) -> dict[str, dict[str, str]]:
    targets: dict[str, dict[str, str]] = {}
    for row in read_tsv(path):
        if row.get("species") == "human":
            targets[row["gene_symbol"]] = row
    return targets


def load_platform_targets(path: Path, targets: dict[str, dict[str, str]]) -> dict[str, dict[str, str]]:
    symbols = set(targets)
    rows: dict[str, dict[str, str]] = {}
    with path.open(newline="") as handle:
        filtered = (line for line in handle if not line.startswith("#"))
        reader = csv.DictReader(filtered, delimiter="\t")
        for row in reader:
            assignment = row.get("gene_assignment", "")
            mrna_assignment = row.get("mrna_assignment", "")
            for symbol in symbols:
                ensembl = targets[symbol].get("ensembl_hint", "")
                if f"// {symbol} //" in assignment or f"gene:{ensembl}" in mrna_assignment:
                    rows[row["ID"]] = {
                        "probe_id": row["ID"],
                        "gene_symbol": symbol,
                        "ensembl_hint": ensembl,
                        "gene_assignment": assignment,
                    }
    return rows


def load_metadata(path: Path) -> dict[str, dict[str, str]]:
    metadata: dict[str, dict[str, str]] = {}
    for row in read_tsv(path):
        chars = parse_characteristics(row.get("characteristics", ""))
        accession = row.get("accession", "")
        region = chars.get("region", "")
        age_label = chars.get("age", "")
        metadata[accession] = {
            "geo_sample_accession": accession,
            "source_name": row.get("source_name", ""),
            "region": region,
            "analysis_region": normalize_region(region),
            "group": normalize_group(chars.get("disease status", "")),
            "sex": chars.get("sex", ""),
            "age_label": age_label,
            "life_stage": life_stage(age_label),
            "age_sort": f"{age_sort_value(age_label):.6g}" if age_sort_value(age_label) is not None else "",
            "stage": chars.get("stage", ""),
            "rin": chars.get("rna integrity number", ""),
        }
    return metadata


def normalize_region(region: str) -> str:
    if region in {"FC", "DFC"}:
        return "frontal_cortex"
    if region == "CBC":
        return "cerebellar_cortex"
    return region


def normalize_group(group: str) -> str:
    group_norm = group.strip().lower()
    if group_norm in {"ds", "down syndrome", "t21", "trisomy 21"}:
        return "DS"
    if group_norm in {"ctl", "control", "con", "euploid"}:
        return "Control"
    return group.strip()


def age_sort_value(age: str) -> float | None:
    match = re.fullmatch(r"(\d+(?:\.\d+)?)(wg|mo|yr)", age.strip(), flags=re.IGNORECASE)
    if not match:
        return None
    value = float(match.group(1))
    unit = match.group(2).lower()
    if unit == "wg":
        return -40.0 + value
    if unit == "mo":
        return value / 12.0
    if unit == "yr":
        return value
    return None


def life_stage(age: str) -> str:
    match = re.fullmatch(r"(\d+(?:\.\d+)?)(wg|mo|yr)", age.strip(), flags=re.IGNORECASE)
    if not match:
        return "unknown"
    value = float(match.group(1))
    unit = match.group(2).lower()
    if unit == "wg":
        return "prenatal"
    if unit == "mo":
        return "infant"
    if value < 18:
        return "child_adolescent"
    return "adult"


def load_series_matrix(path: Path, wanted_ids: set[str]) -> tuple[list[str], dict[str, list[float]]]:
    sample_ids: list[str] = []
    values: dict[str, list[float]] = {}
    in_table = False
    with gzip.open(path, "rt", newline="") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line == "!series_matrix_table_begin":
                in_table = True
                continue
            if line == "!series_matrix_table_end":
                break
            if not in_table:
                continue
            fields = [field.strip('"') for field in line.split("\t")]
            if fields and fields[0] == "ID_REF":
                sample_ids = fields[1:]
                continue
            probe_id = fields[0]
            if probe_id in wanted_ids:
                values[probe_id] = [float(value) for value in fields[1:]]
    return sample_ids, values


def fmt(value: float | int | None) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def mean(values: list[float]) -> float:
    return statistics.fmean(values) if values else math.nan


def median(values: list[float]) -> float:
    return statistics.median(values) if values else math.nan


def summarize_by_keys(rows: list[dict[str, object]], keys: list[str]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, ...], list[dict[str, object]]] = {}
    for row in rows:
        key = tuple(str(row.get(field, "")) for field in keys)
        grouped.setdefault(key, []).append(row)
    summaries: list[dict[str, object]] = []
    for key, group_rows in sorted(grouped.items()):
        expression = [float(row["expression_value"]) for row in group_rows]
        summary = dict(zip(keys, key, strict=False))
        summary.update(
            {
                "n": len(group_rows),
                "mean_expression": fmt(mean(expression)),
                "median_expression": fmt(median(expression)),
                "min_expression": fmt(min(expression)),
                "max_expression": fmt(max(expression)),
            }
        )
        summaries.append(summary)
    return summaries


def summarize_effects(rows: list[dict[str, object]], keys: list[str]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, ...], dict[str, list[float]]] = {}
    for row in rows:
        key = tuple(str(row.get(field, "")) for field in keys)
        group = str(row.get("group", ""))
        grouped.setdefault(key, {"DS": [], "Control": []})
        if group in {"DS", "Control"}:
            grouped[key][group].append(float(row["expression_value"]))
    summaries: list[dict[str, object]] = []
    for key, values in sorted(grouped.items()):
        if not values["DS"] or not values["Control"]:
            continue
        summary = dict(zip(keys, key, strict=False))
        summary.update(
            {
                "n_ds": len(values["DS"]),
                "n_control": len(values["Control"]),
                "mean_expression_ds": fmt(mean(values["DS"])),
                "mean_expression_control": fmt(mean(values["Control"])),
                "delta_ds_minus_control": fmt(mean(values["DS"]) - mean(values["Control"])),
            }
        )
        summaries.append(summary)
    return summaries


def main() -> int:
    args = parse_args()
    targets = load_human_targets(Path(args.target_genes))
    probe_targets = load_platform_targets(Path(args.platform_table), targets)
    metadata = load_metadata(Path(args.metadata))
    sample_ids, matrix_values = load_series_matrix(Path(args.series_matrix), set(probe_targets))

    sample_rows: list[dict[str, object]] = []
    for probe_id, target in sorted(probe_targets.items(), key=lambda item: item[1]["gene_symbol"]):
        values = matrix_values.get(probe_id)
        if values is None:
            continue
        for sample_id, value in zip(sample_ids, values, strict=False):
            meta = metadata.get(sample_id, {})
            sample_rows.append(
                {
                    "dataset": "GSE59630",
                    "gene_symbol": target["gene_symbol"],
                    "probe_id": probe_id,
                    "ensembl_hint": target["ensembl_hint"],
                    "geo_sample_accession": sample_id,
                    "source_name": meta.get("source_name", ""),
                    "region": meta.get("region", ""),
                    "analysis_region": meta.get("analysis_region", ""),
                    "group": meta.get("group", ""),
                    "sex": meta.get("sex", ""),
                    "age_label": meta.get("age_label", ""),
                    "life_stage": meta.get("life_stage", ""),
                    "age_sort": meta.get("age_sort", ""),
                    "stage": meta.get("stage", ""),
                    "rin": meta.get("rin", ""),
                    "expression_value": fmt(value),
                }
            )

    outdir = Path(args.outdir)
    write_tsv(
        outdir / "GSE59630_target_probe_map.tsv",
        list(probe_targets.values()),
        ["probe_id", "gene_symbol", "ensembl_hint", "gene_assignment"],
    )
    write_tsv(
        outdir / "GSE59630_target_sample_values.tsv",
        sample_rows,
        [
            "dataset",
            "gene_symbol",
            "probe_id",
            "ensembl_hint",
            "geo_sample_accession",
            "source_name",
            "region",
            "analysis_region",
            "group",
            "sex",
            "age_label",
            "life_stage",
            "age_sort",
            "stage",
            "rin",
            "expression_value",
        ],
    )
    group_summary = summarize_by_keys(sample_rows, ["dataset", "gene_symbol", "life_stage", "analysis_region", "group"])
    effect_summary = summarize_effects(sample_rows, ["dataset", "gene_symbol", "life_stage", "analysis_region"])
    prenatal_effect_summary = summarize_effects(
        [row for row in sample_rows if row.get("life_stage") == "prenatal"],
        ["dataset", "gene_symbol", "analysis_region"],
    )
    write_tsv(
        outdir / "GSE59630_target_group_summary.tsv",
        group_summary,
        [
            "dataset",
            "gene_symbol",
            "life_stage",
            "analysis_region",
            "group",
            "n",
            "mean_expression",
            "median_expression",
            "min_expression",
            "max_expression",
        ],
    )
    write_tsv(
        outdir / "GSE59630_target_effect_summary.tsv",
        effect_summary,
        [
            "dataset",
            "gene_symbol",
            "life_stage",
            "analysis_region",
            "n_ds",
            "n_control",
            "mean_expression_ds",
            "mean_expression_control",
            "delta_ds_minus_control",
        ],
    )
    write_tsv(
        outdir / "GSE59630_target_prenatal_effect_summary.tsv",
        prenatal_effect_summary,
        [
            "dataset",
            "gene_symbol",
            "analysis_region",
            "n_ds",
            "n_control",
            "mean_expression_ds",
            "mean_expression_control",
            "delta_ds_minus_control",
        ],
    )
    print(f"Wrote GSE59630 target summaries to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
