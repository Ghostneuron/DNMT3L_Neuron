#!/usr/bin/env python3
"""Stream target-gene rows from GSE305146 fetal cortex 10x multiome matrices."""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import math
import re
import statistics
import urllib.request
from pathlib import Path


BASE_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE305nnn/GSE305146/suppl"
USER_AGENT = "DNMT3L-Neuron-GSE305146-targets/0.1"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-genes", default="config/target_genes.tsv")
    parser.add_argument("--metadata", default="data/geo/GSE305146/metadata/sample_metadata.tsv")
    parser.add_argument("--local-supplementary", default="data/geo/GSE305146/supplementary")
    parser.add_argument("--outdir", default="results/ds_fetal_expression")
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--limit-libraries", type=int, help="Debug: process only the first N selected libraries.")
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


def normalize_group(group: str) -> str:
    group_norm = group.strip().lower()
    if group_norm in {"ds", "trisomy 21", "t21", "down syndrome"}:
        return "DS"
    if group_norm in {"con", "control", "disomic", "euploid"}:
        return "Control"
    return group.strip()


def load_human_targets(path: Path) -> dict[str, str]:
    targets: dict[str, str] = {}
    for row in read_tsv(path):
        if row.get("species") == "human":
            targets[row["ensembl_hint"]] = row["gene_symbol"]
    return targets


def selected_libraries(metadata_path: Path) -> list[dict[str, str]]:
    libraries: list[dict[str, str]] = []
    for row in read_tsv(metadata_path):
        chars = parse_characteristics(row.get("characteristics", ""))
        if chars.get("library type") != "Gene Expression":
            continue
        if chars.get("in final_analysis", "").upper() == "FALSE":
            continue
        prefix = library_prefix(row.get("title", ""))
        if not prefix:
            continue
        libraries.append(
            {
                "library_prefix": prefix,
                "geo_sample_accession": row.get("accession", ""),
                "biological_sample_id": chars.get("sample id", ""),
                "sample_name": chars.get("name", ""),
                "group": normalize_group(chars.get("group", chars.get("genotype", ""))),
                "genotype": chars.get("genotype", ""),
                "dev_stage_pcw": chars.get("dev stage_pcw", ""),
                "sex": chars.get("sex", ""),
                "platform_refs": row.get("platform_refs", ""),
            }
        )
    return libraries


def library_prefix(title: str) -> str:
    match = re.search(r"Library\s+([^,\s]+)_gex", title)
    return match.group(1) if match else ""


def open_text_gzip(path_or_url: str, timeout: int) -> io.TextIOBase:
    if path_or_url.startswith("http"):
        request = urllib.request.Request(path_or_url, headers={"User-Agent": USER_AGENT})
        response = urllib.request.urlopen(request, timeout=timeout)
        return io.TextIOWrapper(gzip.GzipFile(fileobj=response))
    return gzip.open(path_or_url, "rt")


def feature_url(prefix: str) -> str:
    return f"{BASE_URL}/GSE305146_{prefix}_features.tsv.gz"


def matrix_url(prefix: str) -> str:
    return f"{BASE_URL}/GSE305146_{prefix}_matrix.mtx.gz"


def local_or_url(local_root: Path, prefix: str, kind: str) -> str:
    local_path = local_root / f"GSE305146_{prefix}_{kind}"
    if local_path.exists():
        return str(local_path)
    if kind == "features.tsv.gz":
        return feature_url(prefix)
    if kind == "matrix.mtx.gz":
        return matrix_url(prefix)
    raise ValueError(kind)


def parse_features(path_or_url: str, targets: dict[str, str], timeout: int) -> tuple[dict[int, str], int, list[dict[str, object]]]:
    row_to_symbol: dict[int, str] = {}
    target_feature_rows: list[dict[str, object]] = []
    last_gene_expression_row = 0
    gene_expression_after_peak = False
    seen_peak = False

    with open_text_gzip(path_or_url, timeout) as handle:
        for index, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            ensembl, symbol, feature_type = fields[:3]
            if feature_type == "Gene Expression":
                if seen_peak:
                    gene_expression_after_peak = True
                last_gene_expression_row = index
            else:
                seen_peak = True
            if feature_type == "Gene Expression" and (ensembl in targets or symbol in targets.values()):
                target_symbol = targets.get(ensembl, symbol)
                row_to_symbol[index] = target_symbol
                target_feature_rows.append(
                    {
                        "feature_row": index,
                        "ensembl_id": ensembl,
                        "gene_symbol": target_symbol,
                        "feature_type": feature_type,
                    }
                )
    if gene_expression_after_peak:
        raise RuntimeError(f"Unexpected feature ordering in {path_or_url}: Gene Expression rows after Peaks")
    return row_to_symbol, last_gene_expression_row, target_feature_rows


def stream_matrix_targets(
    path_or_url: str,
    row_to_symbol: dict[int, str],
    last_gene_expression_row: int,
    timeout: int,
) -> dict[str, object]:
    target_umi = {symbol: 0.0 for symbol in row_to_symbol.values()}
    target_detected_cells = {symbol: 0 for symbol in row_to_symbol.values()}
    total_gene_expression_umi = 0.0
    n_features = n_cells = n_entries = 0
    dims_seen = False

    with open_text_gzip(path_or_url, timeout) as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            fields = line.strip().split()
            if not fields:
                continue
            if not dims_seen:
                n_features, n_cells, n_entries = (int(fields[0]), int(fields[1]), int(fields[2]))
                dims_seen = True
                continue
            row_index = int(fields[0])
            value = float(fields[2])
            if row_index <= last_gene_expression_row:
                total_gene_expression_umi += value
            symbol = row_to_symbol.get(row_index)
            if symbol:
                target_umi[symbol] += value
                target_detected_cells[symbol] += 1

    return {
        "n_features": n_features,
        "n_cells": n_cells,
        "n_entries": n_entries,
        "total_gene_expression_umi": total_gene_expression_umi,
        "target_umi": target_umi,
        "target_detected_cells": target_detected_cells,
    }


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


def build_library_rows(args: argparse.Namespace) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    targets = load_human_targets(Path(args.target_genes))
    libraries = selected_libraries(Path(args.metadata))
    if args.limit_libraries:
        libraries = libraries[: args.limit_libraries]

    library_rows: list[dict[str, object]] = []
    readiness_rows: list[dict[str, object]] = []
    feature_rows: list[dict[str, object]] = []
    local_root = Path(args.local_supplementary)

    for idx, library in enumerate(libraries, start=1):
        prefix = library["library_prefix"]
        features_path = local_or_url(local_root, prefix, "features.tsv.gz")
        matrix_path = local_or_url(local_root, prefix, "matrix.mtx.gz")
        print(f"[{idx}/{len(libraries)}] streaming {prefix}", flush=True)

        row_to_symbol, last_gene_row, target_feature_rows = parse_features(features_path, targets, args.timeout)
        for feature in target_feature_rows:
            feature_rows.append({"library_prefix": prefix, **feature})

        matrix_stats = stream_matrix_targets(matrix_path, row_to_symbol, last_gene_row, args.timeout)
        total_umi = float(matrix_stats["total_gene_expression_umi"])
        n_cells = int(matrix_stats["n_cells"])
        target_umi = matrix_stats["target_umi"]
        target_detected = matrix_stats["target_detected_cells"]
        assert isinstance(target_umi, dict)
        assert isinstance(target_detected, dict)

        readiness_rows.append(
            {
                "library_prefix": prefix,
                "biological_sample_id": library["biological_sample_id"],
                "group": library["group"],
                "dev_stage_pcw": library["dev_stage_pcw"],
                "n_features": matrix_stats["n_features"],
                "n_cells": n_cells,
                "n_entries": matrix_stats["n_entries"],
                "last_gene_expression_row": last_gene_row,
                "target_features_found": len(row_to_symbol),
                "total_gene_expression_umi": fmt(total_umi),
            }
        )

        for ensembl, symbol in targets.items():
            umi = float(target_umi.get(symbol, 0.0))
            detected_cells = int(target_detected.get(symbol, 0))
            cpm = umi / total_umi * 1_000_000 if total_umi else math.nan
            library_rows.append(
                {
                    "dataset": "GSE305146",
                    "library_prefix": prefix,
                    "geo_sample_accession": library["geo_sample_accession"],
                    "biological_sample_id": library["biological_sample_id"],
                    "sample_name": library["sample_name"],
                    "group": library["group"],
                    "genotype": library["genotype"],
                    "dev_stage_pcw": library["dev_stage_pcw"],
                    "sex": library["sex"],
                    "platform_refs": library["platform_refs"],
                    "gene_symbol": symbol,
                    "ensembl_id": ensembl,
                    "n_cells": n_cells,
                    "total_gene_expression_umi": fmt(total_umi),
                    "target_umi": fmt(umi),
                    "target_cpm": fmt(cpm),
                    "target_log2_cpm": fmt(math.log2(cpm + 1.0) if not math.isnan(cpm) else math.nan),
                    "detected_cells": detected_cells,
                    "detected_fraction": fmt(detected_cells / n_cells if n_cells else math.nan),
                }
            )
    return library_rows, readiness_rows, feature_rows


def aggregate_biological_samples(library_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str], list[dict[str, object]]] = {}
    for row in library_rows:
        key = (
            str(row.get("biological_sample_id", "")),
            str(row.get("gene_symbol", "")),
            str(row.get("group", "")),
        )
        grouped.setdefault(key, []).append(row)

    sample_rows: list[dict[str, object]] = []
    for (sample_id, gene, group), rows in sorted(grouped.items()):
        total_umi = sum(float(row["total_gene_expression_umi"]) for row in rows)
        target_umi = sum(float(row["target_umi"]) for row in rows)
        n_cells = sum(int(row["n_cells"]) for row in rows)
        detected_cells = sum(int(row["detected_cells"]) for row in rows)
        cpm = target_umi / total_umi * 1_000_000 if total_umi else math.nan
        first = rows[0]
        sample_rows.append(
            {
                "dataset": "GSE305146",
                "biological_sample_id": sample_id,
                "group": group,
                "genotype": first.get("genotype", ""),
                "dev_stage_pcw": first.get("dev_stage_pcw", ""),
                "sex": first.get("sex", ""),
                "gene_symbol": gene,
                "library_count": len(rows),
                "n_cells": n_cells,
                "total_gene_expression_umi": fmt(total_umi),
                "target_umi": fmt(target_umi),
                "target_cpm": fmt(cpm),
                "target_log2_cpm": fmt(math.log2(cpm + 1.0) if not math.isnan(cpm) else math.nan),
                "detected_cells": detected_cells,
                "detected_fraction": fmt(detected_cells / n_cells if n_cells else math.nan),
            }
        )
    return sample_rows


def summarize_by_keys(rows: list[dict[str, object]], keys: list[str]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, ...], list[dict[str, object]]] = {}
    for row in rows:
        key = tuple(str(row.get(field, "")) for field in keys)
        grouped.setdefault(key, []).append(row)
    summaries: list[dict[str, object]] = []
    for key, group_rows in sorted(grouped.items()):
        log_values = [float(row["target_log2_cpm"]) for row in group_rows if row.get("target_log2_cpm") != ""]
        cpm_values = [float(row["target_cpm"]) for row in group_rows if row.get("target_cpm") != ""]
        detected_fraction = [float(row["detected_fraction"]) for row in group_rows if row.get("detected_fraction") != ""]
        summary = dict(zip(keys, key, strict=False))
        summary.update(
            {
                "n": len(group_rows),
                "mean_cpm": fmt(mean(cpm_values)),
                "median_cpm": fmt(median(cpm_values)),
                "mean_log2_cpm": fmt(mean(log_values)),
                "median_log2_cpm": fmt(median(log_values)),
                "mean_detected_fraction": fmt(mean(detected_fraction)),
                "median_detected_fraction": fmt(median(detected_fraction)),
            }
        )
        summaries.append(summary)
    return summaries


def summarize_effects(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[str, dict[str, list[dict[str, object]]]] = {}
    for row in rows:
        gene = str(row.get("gene_symbol", ""))
        group = str(row.get("group", ""))
        grouped.setdefault(gene, {"DS": [], "Control": []})
        if group in {"DS", "Control"}:
            grouped[gene][group].append(row)
    summaries: list[dict[str, object]] = []
    for gene, values in sorted(grouped.items()):
        if not values["DS"] or not values["Control"]:
            continue
        ds_log = [float(row["target_log2_cpm"]) for row in values["DS"]]
        control_log = [float(row["target_log2_cpm"]) for row in values["Control"]]
        ds_det = [float(row["detected_fraction"]) for row in values["DS"]]
        control_det = [float(row["detected_fraction"]) for row in values["Control"]]
        summaries.append(
            {
                "dataset": "GSE305146",
                "gene_symbol": gene,
                "n_ds": len(values["DS"]),
                "n_control": len(values["Control"]),
                "mean_log2_cpm_ds": fmt(mean(ds_log)),
                "mean_log2_cpm_control": fmt(mean(control_log)),
                "delta_log2_cpm_ds_minus_control": fmt(mean(ds_log) - mean(control_log)),
                "mean_detected_fraction_ds": fmt(mean(ds_det)),
                "mean_detected_fraction_control": fmt(mean(control_det)),
                "delta_detected_fraction_ds_minus_control": fmt(mean(ds_det) - mean(control_det)),
            }
        )
    return summaries


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    library_rows, readiness_rows, feature_rows = build_library_rows(args)
    sample_rows = aggregate_biological_samples(library_rows)
    group_summary = summarize_by_keys(sample_rows, ["dataset", "gene_symbol", "group"])
    pcw_summary = summarize_by_keys(sample_rows, ["dataset", "gene_symbol", "dev_stage_pcw", "group"])
    effect_summary = summarize_effects(sample_rows)

    write_tsv(
        outdir / "GSE305146_target_library_readiness.tsv",
        readiness_rows,
        [
            "library_prefix",
            "biological_sample_id",
            "group",
            "dev_stage_pcw",
            "n_features",
            "n_cells",
            "n_entries",
            "last_gene_expression_row",
            "target_features_found",
            "total_gene_expression_umi",
        ],
    )
    write_tsv(
        outdir / "GSE305146_target_feature_rows.tsv",
        feature_rows,
        ["library_prefix", "feature_row", "ensembl_id", "gene_symbol", "feature_type"],
    )
    write_tsv(
        outdir / "GSE305146_target_library_values.tsv",
        library_rows,
        [
            "dataset",
            "library_prefix",
            "geo_sample_accession",
            "biological_sample_id",
            "sample_name",
            "group",
            "genotype",
            "dev_stage_pcw",
            "sex",
            "platform_refs",
            "gene_symbol",
            "ensembl_id",
            "n_cells",
            "total_gene_expression_umi",
            "target_umi",
            "target_cpm",
            "target_log2_cpm",
            "detected_cells",
            "detected_fraction",
        ],
    )
    write_tsv(
        outdir / "GSE305146_target_biological_sample_values.tsv",
        sample_rows,
        [
            "dataset",
            "biological_sample_id",
            "group",
            "genotype",
            "dev_stage_pcw",
            "sex",
            "gene_symbol",
            "library_count",
            "n_cells",
            "total_gene_expression_umi",
            "target_umi",
            "target_cpm",
            "target_log2_cpm",
            "detected_cells",
            "detected_fraction",
        ],
    )
    write_tsv(
        outdir / "GSE305146_target_group_summary.tsv",
        group_summary,
        [
            "dataset",
            "gene_symbol",
            "group",
            "n",
            "mean_cpm",
            "median_cpm",
            "mean_log2_cpm",
            "median_log2_cpm",
            "mean_detected_fraction",
            "median_detected_fraction",
        ],
    )
    write_tsv(
        outdir / "GSE305146_target_pcw_summary.tsv",
        pcw_summary,
        [
            "dataset",
            "gene_symbol",
            "dev_stage_pcw",
            "group",
            "n",
            "mean_cpm",
            "median_cpm",
            "mean_log2_cpm",
            "median_log2_cpm",
            "mean_detected_fraction",
            "median_detected_fraction",
        ],
    )
    write_tsv(
        outdir / "GSE305146_target_effect_summary.tsv",
        effect_summary,
        [
            "dataset",
            "gene_symbol",
            "n_ds",
            "n_control",
            "mean_log2_cpm_ds",
            "mean_log2_cpm_control",
            "delta_log2_cpm_ds_minus_control",
            "mean_detected_fraction_ds",
            "mean_detected_fraction_control",
            "delta_detected_fraction_ds_minus_control",
        ],
    )
    print(f"Wrote GSE305146 target summaries to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
