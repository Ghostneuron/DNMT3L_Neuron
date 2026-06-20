#!/usr/bin/env python3
"""Summarize global CpG methylation in GSE168276 SH-SY5Y WGBS samples."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


DEFAULT_METADATA = "data/geo/GSE168276/metadata/sample_metadata.tsv"
DEFAULT_SUPPLEMENTARY_DIR = "data/geo/GSE168276/supplementary"
DEFAULT_OUTDIR = "results/shsy5y_gse168276"
PHASE_ORDER = {"Growth": 0, "Phase 1": 1, "Phase 3": 2}
GENOTYPE_ORDER = {"GFP": 0, "DNMT3L+GFP": 1}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", default=DEFAULT_METADATA)
    parser.add_argument("--supplementary-dir", default=DEFAULT_SUPPLEMENTARY_DIR)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--min-coverage", type=int, default=5)
    parser.add_argument("--workers", type=int, default=4)
    return parser.parse_args()


def parse_characteristics(value: str) -> dict[str, str]:
    fields = {}
    for token in value.split(";"):
        if "=" not in token:
            continue
        key, raw = token.split("=", 1)
        fields[key.strip().lower()] = raw.strip()
    return fields


def read_samples(metadata_path: Path, supplementary_dir: Path) -> list[dict[str, str]]:
    samples = []
    with metadata_path.open(newline="", errors="replace") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            chars = parse_characteristics(row["characteristics"])
            url = row["supplementary_data"]
            filename = Path(url).name
            path = supplementary_dir / filename
            samples.append(
                {
                    "accession": row["accession"],
                    "sample": row["title"],
                    "genotype": chars.get("genotype", ""),
                    "phase": chars.get("differentiation phase", ""),
                    "replicate": chars.get("cell line replicate", ""),
                    "path": str(path),
                }
            )
    return samples


def summarize_sample(sample: dict[str, str], min_coverage: int) -> dict[str, str]:
    path = Path(sample["path"])
    total_cpg_lines = 0
    covered_cpg_lines = 0
    methylated_reads = 0
    unmethylated_reads = 0
    methylated_reads_all = 0
    unmethylated_reads_all = 0
    beta_sum = 0.0
    coverage_sum = 0

    with gzip.open(path, "rt", errors="replace") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7 or parts[5] != "CG":
                continue
            total_cpg_lines += 1
            methylated = int(parts[3])
            unmethylated = int(parts[4])
            coverage = methylated + unmethylated
            methylated_reads_all += methylated
            unmethylated_reads_all += unmethylated
            if coverage < min_coverage:
                continue
            covered_cpg_lines += 1
            methylated_reads += methylated
            unmethylated_reads += unmethylated
            coverage_sum += coverage
            beta_sum += methylated / coverage

    weighted_beta = methylated_reads / (methylated_reads + unmethylated_reads) if methylated_reads + unmethylated_reads else float("nan")
    weighted_beta_all = (
        methylated_reads_all / (methylated_reads_all + unmethylated_reads_all)
        if methylated_reads_all + unmethylated_reads_all
        else float("nan")
    )
    mean_beta = beta_sum / covered_cpg_lines if covered_cpg_lines else float("nan")
    mean_coverage = coverage_sum / covered_cpg_lines if covered_cpg_lines else float("nan")

    return {
        **{key: sample[key] for key in ["accession", "sample", "genotype", "phase", "replicate"]},
        "file": str(path),
        "min_coverage": str(min_coverage),
        "total_cpg_lines": str(total_cpg_lines),
        "covered_cpg_lines": str(covered_cpg_lines),
        "covered_cpg_fraction": fmt(covered_cpg_lines / total_cpg_lines if total_cpg_lines else float("nan")),
        "methylated_reads_cov_ge_min": str(methylated_reads),
        "unmethylated_reads_cov_ge_min": str(unmethylated_reads),
        "weighted_methylation_cov_ge_min": fmt(weighted_beta),
        "mean_cpg_methylation_cov_ge_min": fmt(mean_beta),
        "mean_coverage_cov_ge_min": fmt(mean_coverage),
        "weighted_methylation_all_cpg_lines": fmt(weighted_beta_all),
    }


def safe_float(value: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return number if math.isfinite(number) else float("nan")


def mean(values: list[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    return sum(finite) / len(finite) if finite else float("nan")


def stdev(values: list[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    if len(finite) < 2:
        return float("nan")
    avg = mean(finite)
    return math.sqrt(sum((value - avg) ** 2 for value in finite) / (len(finite) - 1))


def fmt(value: float) -> str:
    return f"{value:.6g}" if math.isfinite(value) else ""


def phase_summary(sample_rows: list[dict[str, str]]) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    group_rows = []
    comparison_rows = []
    grouped: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in sample_rows:
        grouped.setdefault((row["phase"], row["genotype"]), []).append(row)

    for (phase, genotype), rows in sorted(grouped.items(), key=lambda item: (PHASE_ORDER.get(item[0][0], 99), GENOTYPE_ORDER.get(item[0][1], 99))):
        weighted = [safe_float(row["weighted_methylation_cov_ge_min"]) for row in rows]
        unweighted = [safe_float(row["mean_cpg_methylation_cov_ge_min"]) for row in rows]
        coverage = [safe_float(row["mean_coverage_cov_ge_min"]) for row in rows]
        covered = [safe_float(row["covered_cpg_lines"]) for row in rows]
        group_rows.append(
            {
                "phase": phase,
                "genotype": genotype,
                "n_samples": str(len(rows)),
                "mean_weighted_methylation": fmt(mean(weighted)),
                "sd_weighted_methylation": fmt(stdev(weighted)),
                "mean_unweighted_cpg_methylation": fmt(mean(unweighted)),
                "sd_unweighted_cpg_methylation": fmt(stdev(unweighted)),
                "mean_covered_cpg_lines": fmt(mean(covered)),
                "mean_coverage": fmt(mean(coverage)),
                "samples": ";".join(row["sample"] for row in sorted(rows, key=lambda r: r["replicate"])),
            }
        )

    for phase in sorted({row["phase"] for row in sample_rows}, key=lambda value: PHASE_ORDER.get(value, 99)):
        dnmt3l = [row for row in sample_rows if row["phase"] == phase and row["genotype"] == "DNMT3L+GFP"]
        gfp = [row for row in sample_rows if row["phase"] == phase and row["genotype"] == "GFP"]
        dnmt3l_weighted = [safe_float(row["weighted_methylation_cov_ge_min"]) for row in dnmt3l]
        gfp_weighted = [safe_float(row["weighted_methylation_cov_ge_min"]) for row in gfp]
        dnmt3l_unweighted = [safe_float(row["mean_cpg_methylation_cov_ge_min"]) for row in dnmt3l]
        gfp_unweighted = [safe_float(row["mean_cpg_methylation_cov_ge_min"]) for row in gfp]
        comparison_rows.append(
            {
                "phase": phase,
                "n_dnmt3l_gfp": str(len(dnmt3l)),
                "n_gfp": str(len(gfp)),
                "dnmt3l_gfp_mean_weighted_methylation": fmt(mean(dnmt3l_weighted)),
                "gfp_mean_weighted_methylation": fmt(mean(gfp_weighted)),
                "delta_weighted_methylation_dnmt3l_minus_gfp": fmt(mean(dnmt3l_weighted) - mean(gfp_weighted)),
                "dnmt3l_gfp_mean_unweighted_cpg_methylation": fmt(mean(dnmt3l_unweighted)),
                "gfp_mean_unweighted_cpg_methylation": fmt(mean(gfp_unweighted)),
                "delta_unweighted_cpg_methylation_dnmt3l_minus_gfp": fmt(mean(dnmt3l_unweighted) - mean(gfp_unweighted)),
            }
        )
    return group_rows, comparison_rows


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    samples = read_samples(Path(args.metadata), Path(args.supplementary_dir))
    missing = [sample["path"] for sample in samples if not Path(sample["path"]).exists()]
    if missing:
        raise FileNotFoundError(f"Missing {len(missing)} GSE168276 CpG reports: {missing[:3]}")

    sample_rows = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(summarize_sample, sample, args.min_coverage): sample for sample in samples}
        for future in as_completed(futures):
            row = future.result()
            print(f"Summarized {row['sample']}", flush=True)
            sample_rows.append(row)

    sample_rows.sort(key=lambda row: (PHASE_ORDER.get(row["phase"], 99), GENOTYPE_ORDER.get(row["genotype"], 99), int(row["replicate"])))
    group_rows, comparison_rows = phase_summary(sample_rows)

    outdir = Path(args.outdir)
    sample_fields = [
        "accession",
        "sample",
        "genotype",
        "phase",
        "replicate",
        "file",
        "min_coverage",
        "total_cpg_lines",
        "covered_cpg_lines",
        "covered_cpg_fraction",
        "methylated_reads_cov_ge_min",
        "unmethylated_reads_cov_ge_min",
        "weighted_methylation_cov_ge_min",
        "mean_cpg_methylation_cov_ge_min",
        "mean_coverage_cov_ge_min",
        "weighted_methylation_all_cpg_lines",
    ]
    group_fields = [
        "phase",
        "genotype",
        "n_samples",
        "mean_weighted_methylation",
        "sd_weighted_methylation",
        "mean_unweighted_cpg_methylation",
        "sd_unweighted_cpg_methylation",
        "mean_covered_cpg_lines",
        "mean_coverage",
        "samples",
    ]
    comparison_fields = [
        "phase",
        "n_dnmt3l_gfp",
        "n_gfp",
        "dnmt3l_gfp_mean_weighted_methylation",
        "gfp_mean_weighted_methylation",
        "delta_weighted_methylation_dnmt3l_minus_gfp",
        "dnmt3l_gfp_mean_unweighted_cpg_methylation",
        "gfp_mean_unweighted_cpg_methylation",
        "delta_unweighted_cpg_methylation_dnmt3l_minus_gfp",
    ]
    write_tsv(outdir / "gse168276_sample_global_methylation.tsv", sample_rows, sample_fields)
    write_tsv(outdir / "gse168276_phase_genotype_global_methylation.tsv", group_rows, group_fields)
    write_tsv(outdir / "gse168276_phase_dnmt3l_vs_gfp_global_delta.tsv", comparison_rows, comparison_fields)
    print(f"Wrote {len(sample_rows)} sample rows")
    print(f"Wrote {len(group_rows)} phase/genotype rows")
    print(f"Wrote {len(comparison_rows)} phase comparison rows")


if __name__ == "__main__":
    main()
