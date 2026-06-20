#!/usr/bin/env python3
"""Summarize target-gene expression in fetal Down syndrome brain datasets.

The first-pass DS question is whether DNMT3L and related genes are higher in
trisomy 21 fetal brain, and whether that differs by age or region. This script
uses the small processed bulk RNA-seq count tables that are public in GEO.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import statistics
from pathlib import Path


DEFAULT_DATASETS = {
    "GSE301886": {
        "counts": "data/geo/GSE301886/supplementary/GSE301886_T21_raw_counts.txt.gz",
        "metadata": "data/geo/GSE301886/metadata/sample_metadata.tsv",
        "id_column": "gene.id",
        "gene_name_column": "",
    },
    "GSE305139": {
        "counts": "data/geo/GSE305139/supplementary/GSE305139_raw_counts_tissue_bulk.tsv.gz",
        "metadata": "data/geo/GSE305139/metadata/sample_metadata.tsv",
        "id_column": "gene_id",
        "gene_name_column": "gene_name",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-genes", default="config/target_genes.tsv")
    parser.add_argument("--outdir", default="results/ds_fetal_expression")
    parser.add_argument(
        "--human-ds-bulk-count-de",
        default="results/human_ds_bulk_count_de/human_ds_bulk_count_de_all_contexts.tsv",
    )
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


def clean_ensembl(gene_id: str) -> str:
    return gene_id.split(".", 1)[0]


def as_float(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def mean(values: list[float]) -> float:
    return statistics.fmean(values) if values else math.nan


def median(values: list[float]) -> float:
    return statistics.median(values) if values else math.nan


def fmt(value: float | int | None) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def load_human_targets(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    by_ensembl: dict[str, str] = {}
    by_symbol: dict[str, str] = {}
    for row in read_tsv(path):
        if row.get("species") != "human":
            continue
        symbol = row["gene_symbol"]
        ensembl = row["ensembl_hint"]
        by_ensembl[ensembl] = symbol
        by_symbol[symbol.upper()] = symbol
    return by_ensembl, by_symbol


def load_metadata(dataset: str, path: Path) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    for row in read_tsv(path):
        chars = parse_characteristics(row.get("characteristics", ""))
        if dataset == "GSE301886":
            sample_label = row.get("title", "").split(",")[-1].strip()
            group = chars.get("karyotype", "")
            age = chars.get("gestational age", "")
            final_analysis = "TRUE"
        elif dataset == "GSE305139":
            sample_label = chars.get("library igf_id", "")
            group = chars.get("group", chars.get("genotype", ""))
            age = chars.get("dev stage_pcw", "")
            final_analysis = chars.get("in final_analysis", "")
        else:
            continue

        rows[sample_label] = {
            "geo_sample_accession": row.get("accession", ""),
            "region": chars.get("tissue", row.get("source_name", "")),
            "group": normalize_group(group),
            "gestational_age_pcw": age,
            "sex": chars.get("sex", chars.get("Sex", "")),
            "final_analysis": final_analysis,
        }
    return rows


def normalize_group(group: str) -> str:
    group_norm = group.strip().lower()
    if group_norm in {"t21", "trisomy 21", "ds", "down syndrome"}:
        return "DS"
    if group_norm in {"control", "con", "disomic", "euploid"}:
        return "Control"
    return group.strip()


def included_for_summary(meta: dict[str, str]) -> bool:
    if meta.get("group") not in {"DS", "Control"}:
        return False
    final_analysis = meta.get("final_analysis", "")
    return final_analysis.upper() != "FALSE"


def read_count_targets(
    dataset: str,
    path: Path,
    target_by_ensembl: dict[str, str],
    target_by_symbol: dict[str, str],
) -> tuple[list[str], dict[str, float], list[dict[str, object]]]:
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        id_index = 0
        gene_name_index = header.index("gene_name") if "gene_name" in header else None
        first_sample_index = 2 if gene_name_index is not None else 1
        sample_labels = header[first_sample_index:]
        library_sizes = {sample: 0.0 for sample in sample_labels}
        target_rows: list[dict[str, object]] = []

        for row in reader:
            if not row:
                continue
            gene_id_raw = row[id_index]
            ensembl = clean_ensembl(gene_id_raw)
            gene_name = row[gene_name_index] if gene_name_index is not None else ""
            symbol = target_by_ensembl.get(ensembl) or target_by_symbol.get(gene_name.upper())
            counts = [as_float(value) or 0.0 for value in row[first_sample_index:]]
            for sample, count in zip(sample_labels, counts, strict=False):
                library_sizes[sample] += count
            if symbol:
                target_rows.append(
                    {
                        "dataset": dataset,
                        "gene_symbol": symbol,
                        "ensembl_id": ensembl,
                        "gene_id_raw": gene_id_raw,
                        "counts": dict(zip(sample_labels, counts, strict=False)),
                    }
                )
    return sample_labels, library_sizes, target_rows


def build_sample_values(args: argparse.Namespace) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    target_by_ensembl, target_by_symbol = load_human_targets(Path(args.target_genes))
    sample_values: list[dict[str, object]] = []
    dataset_rows: list[dict[str, object]] = []

    for dataset, spec in DEFAULT_DATASETS.items():
        counts_path = Path(spec["counts"])
        metadata_path = Path(spec["metadata"])
        if not counts_path.exists() or not metadata_path.exists():
            dataset_rows.append(
                {
                    "dataset": dataset,
                    "counts_path": str(counts_path),
                    "metadata_path": str(metadata_path),
                    "status": "missing_input",
                    "samples_in_counts": 0,
                    "samples_with_metadata": 0,
                    "target_gene_rows": 0,
                }
            )
            continue
        metadata = load_metadata(dataset, metadata_path)
        sample_labels, library_sizes, target_rows = read_count_targets(
            dataset, counts_path, target_by_ensembl, target_by_symbol
        )
        dataset_rows.append(
            {
                "dataset": dataset,
                "counts_path": str(counts_path),
                "metadata_path": str(metadata_path),
                "status": "ready",
                "samples_in_counts": len(sample_labels),
                "samples_with_metadata": sum(1 for sample in sample_labels if sample in metadata),
                "target_gene_rows": len(target_rows),
            }
        )

        for target in target_rows:
            counts_by_sample = target["counts"]
            assert isinstance(counts_by_sample, dict)
            for sample in sample_labels:
                meta = metadata.get(sample, {})
                count = float(counts_by_sample.get(sample, 0.0))
                library_size = library_sizes[sample]
                cpm = count / library_size * 1_000_000 if library_size else math.nan
                log2_cpm = math.log2(cpm + 1.0) if not math.isnan(cpm) else math.nan
                sample_values.append(
                    {
                        "dataset": dataset,
                        "gene_symbol": target["gene_symbol"],
                        "ensembl_id": target["ensembl_id"],
                        "gene_id_raw": target["gene_id_raw"],
                        "sample_label": sample,
                        "geo_sample_accession": meta.get("geo_sample_accession", ""),
                        "region": meta.get("region", ""),
                        "group": meta.get("group", ""),
                        "gestational_age_pcw": meta.get("gestational_age_pcw", ""),
                        "age_week_pcw": rounded_age(meta.get("gestational_age_pcw", "")),
                        "sex": meta.get("sex", ""),
                        "final_analysis": meta.get("final_analysis", ""),
                        "included_for_summary": "TRUE" if included_for_summary(meta) else "FALSE",
                        "library_size": fmt(library_size),
                        "count": fmt(count),
                        "cpm": fmt(cpm),
                        "log2_cpm": fmt(log2_cpm),
                        "detected": "TRUE" if count > 0 else "FALSE",
                    }
                )
    return sample_values, dataset_rows


def rounded_age(value: str) -> str:
    age = as_float(value)
    if age is None:
        return ""
    return str(int(round(age)))


def summarize_groups(
    sample_values: list[dict[str, object]],
    count_de: dict[tuple[str, str, str], dict[str, str]] | None = None,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    included = [row for row in sample_values if row.get("included_for_summary") == "TRUE"]
    group_keys = ["dataset", "gene_symbol", "region", "group"]
    age_keys = ["dataset", "gene_symbol", "region", "age_week_pcw", "group"]
    group_summary = summarize_by_keys(included, group_keys)
    age_summary = summarize_by_keys([row for row in included if row.get("age_week_pcw")], age_keys)
    effect_summary = summarize_effects(included, count_de)
    return group_summary, age_summary, effect_summary


def summarize_by_keys(rows: list[dict[str, object]], keys: list[str]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, ...], list[dict[str, object]]] = {}
    for row in rows:
        key = tuple(str(row.get(field, "")) for field in keys)
        grouped.setdefault(key, []).append(row)

    summaries: list[dict[str, object]] = []
    for key, group_rows in sorted(grouped.items()):
        log_values = [float(row["log2_cpm"]) for row in group_rows if row.get("log2_cpm") != ""]
        cpm_values = [float(row["cpm"]) for row in group_rows if row.get("cpm") != ""]
        age_values = [float(row["gestational_age_pcw"]) for row in group_rows if as_float(str(row.get("gestational_age_pcw", ""))) is not None]
        detected = sum(1 for row in group_rows if row.get("detected") == "TRUE")
        summary = dict(zip(keys, key, strict=False))
        summary.update(
            {
                "n": len(group_rows),
                "age_min_pcw": fmt(min(age_values) if age_values else math.nan),
                "age_max_pcw": fmt(max(age_values) if age_values else math.nan),
                "mean_cpm": fmt(mean(cpm_values)),
                "median_cpm": fmt(median(cpm_values)),
                "mean_log2_cpm": fmt(mean(log_values)),
                "median_log2_cpm": fmt(median(log_values)),
                "detected_n": detected,
                "detected_fraction": fmt(detected / len(group_rows) if group_rows else math.nan),
            }
        )
        summaries.append(summary)
    return summaries


def de_context_key(dataset: str, region: str) -> tuple[str, str]:
    aliases = {
        ("GSE305139", "fetal brain"): ("GSE305139", "fetal_brain_final"),
    }
    return aliases.get((dataset, region), (dataset, region))


def load_count_de(path: Path) -> dict[tuple[str, str, str], dict[str, str]]:
    out: dict[tuple[str, str, str], dict[str, str]] = {}
    if not path.exists():
        return out
    for row in read_tsv(path):
        gene = row.get("gene_symbol", "")
        if not gene:
            continue
        key = (row.get("dataset", ""), row.get("context", ""), gene)
        current = out.get(key)
        current_mean = (as_float(current.get("mean_ds_cpm", "")) if current else None) or 0
        current_mean += (as_float(current.get("mean_control_cpm", "")) if current else None) or 0
        row_mean = (as_float(row.get("mean_ds_cpm", "")) or 0) + (as_float(row.get("mean_control_cpm", "")) or 0)
        if current is None or row_mean > current_mean:
            out[key] = row
    return out


def summarize_effects(
    rows: list[dict[str, object]],
    count_de: dict[tuple[str, str, str], dict[str, str]] | None = None,
) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str], dict[str, list[dict[str, object]]]] = {}
    for row in rows:
        key = (str(row.get("dataset", "")), str(row.get("gene_symbol", "")), str(row.get("region", "")))
        grouped.setdefault(key, {"DS": [], "Control": []})
        group = str(row.get("group", ""))
        if group in {"DS", "Control"}:
            grouped[key][group].append(row)

    summaries: list[dict[str, object]] = []
    for (dataset, gene, region), groups in sorted(grouped.items()):
        ds_rows = groups["DS"]
        control_rows = groups["Control"]
        if not ds_rows or not control_rows:
            continue
        ds_log = [float(row["log2_cpm"]) for row in ds_rows if row.get("log2_cpm") != ""]
        control_log = [float(row["log2_cpm"]) for row in control_rows if row.get("log2_cpm") != ""]
        ds_cpm = [float(row["cpm"]) for row in ds_rows if row.get("cpm") != ""]
        control_cpm = [float(row["cpm"]) for row in control_rows if row.get("cpm") != ""]
        descriptive_delta = mean(ds_log) - mean(control_log)
        de_row = (count_de or {}).get((*de_context_key(dataset, region), gene), {})
        edgeR_effect = as_float(de_row.get("effect_log2fc_for_claim", ""))
        summaries.append(
            {
                "dataset": dataset,
                "gene_symbol": gene,
                "region": region,
                "n_ds": len(ds_rows),
                "n_control": len(control_rows),
                "mean_log2_cpm_ds": fmt(mean(ds_log)),
                "mean_log2_cpm_control": fmt(mean(control_log)),
                "delta_log2_cpm_ds_minus_control": fmt(descriptive_delta),
                "effect_log2fc_for_claim": fmt(edgeR_effect if edgeR_effect is not None else descriptive_delta),
                "effect_source": de_row.get("effect_source", "mean_log2CPM_descriptive"),
                "edgeR_logFC": fmt(as_float(de_row.get("logFC", ""))),
                "edgeR_PValue": fmt(as_float(de_row.get("PValue", ""))),
                "edgeR_FDR": fmt(as_float(de_row.get("FDR", ""))),
                "tested_by_edgeR": de_row.get("tested_by_edgeR", ""),
                "mean_cpm_ds": fmt(mean(ds_cpm)),
                "mean_cpm_control": fmt(mean(control_cpm)),
                "delta_cpm_ds_minus_control": fmt(mean(ds_cpm) - mean(control_cpm)),
                "detected_fraction_ds": fmt(sum(row.get("detected") == "TRUE" for row in ds_rows) / len(ds_rows)),
                "detected_fraction_control": fmt(sum(row.get("detected") == "TRUE" for row in control_rows) / len(control_rows)),
            }
        )
    return summaries


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    sample_values, dataset_rows = build_sample_values(args)
    count_de = load_count_de(Path(args.human_ds_bulk_count_de))
    group_summary, age_summary, effect_summary = summarize_groups(sample_values, count_de)

    write_tsv(
        outdir / "ds_bulk_dataset_readiness.tsv",
        dataset_rows,
        [
            "dataset",
            "counts_path",
            "metadata_path",
            "status",
            "samples_in_counts",
            "samples_with_metadata",
            "target_gene_rows",
        ],
    )
    write_tsv(
        outdir / "ds_target_sample_values.tsv",
        sample_values,
        [
            "dataset",
            "gene_symbol",
            "ensembl_id",
            "gene_id_raw",
            "sample_label",
            "geo_sample_accession",
            "region",
            "group",
            "gestational_age_pcw",
            "age_week_pcw",
            "sex",
            "final_analysis",
            "included_for_summary",
            "library_size",
            "count",
            "cpm",
            "log2_cpm",
            "detected",
        ],
    )
    write_tsv(
        outdir / "ds_target_group_summary.tsv",
        group_summary,
        [
            "dataset",
            "gene_symbol",
            "region",
            "group",
            "n",
            "age_min_pcw",
            "age_max_pcw",
            "mean_cpm",
            "median_cpm",
            "mean_log2_cpm",
            "median_log2_cpm",
            "detected_n",
            "detected_fraction",
        ],
    )
    write_tsv(
        outdir / "ds_target_age_week_summary.tsv",
        age_summary,
        [
            "dataset",
            "gene_symbol",
            "region",
            "age_week_pcw",
            "group",
            "n",
            "age_min_pcw",
            "age_max_pcw",
            "mean_cpm",
            "median_cpm",
            "mean_log2_cpm",
            "median_log2_cpm",
            "detected_n",
            "detected_fraction",
        ],
    )
    write_tsv(
        outdir / "ds_target_region_effect_summary.tsv",
        effect_summary,
        [
            "dataset",
            "gene_symbol",
            "region",
            "n_ds",
            "n_control",
            "mean_log2_cpm_ds",
            "mean_log2_cpm_control",
            "delta_log2_cpm_ds_minus_control",
            "effect_log2fc_for_claim",
            "effect_source",
            "edgeR_logFC",
            "edgeR_PValue",
            "edgeR_FDR",
            "tested_by_edgeR",
            "mean_cpm_ds",
            "mean_cpm_control",
            "delta_cpm_ds_minus_control",
            "detected_fraction_ds",
            "detected_fraction_control",
        ],
    )

    print(f"Wrote DS fetal expression summaries to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
