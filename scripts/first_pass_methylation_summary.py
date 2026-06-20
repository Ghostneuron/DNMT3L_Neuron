#!/usr/bin/env python3
"""First-pass methylation delta summaries for direct DNMT3L perturbation datasets."""

from __future__ import annotations

import argparse
import csv
import gzip
import heapq
import math
from pathlib import Path


HUMAN_SAMPLES = [
    {"sample": "Z1", "condition": "control", "meth": "SAMPLE 1 Methylated signal", "unmeth": "SAMPLE 1 Unmethylated Signal", "pval": "SAMPLE 1 Detection Pval"},
    {"sample": "Z2", "condition": "control", "meth": "SAMPLE 2 Methylated signal", "unmeth": "SAMPLE 2 Unmethylated Signal", "pval": "SAMPLE 2 Detection Pval"},
    {"sample": "Z3", "condition": "control", "meth": "SAMPLE 3 Methylated signal", "unmeth": "SAMPLE 3 Unmethylated Signal", "pval": "SAMPLE 3 Detection Pval"},
    {"sample": "D1", "condition": "DNMT3L_overexpression", "meth": "SAMPLE 4 Methylated signal", "unmeth": "SAMPLE 4 Unmethylated Signal", "pval": "SAMPLE 4 Detection Pval"},
    {"sample": "D2", "condition": "DNMT3L_overexpression", "meth": "SAMPLE 5 Methylated signal", "unmeth": "SAMPLE 5 Unmethylated Signal", "pval": "SAMPLE 5 Detection Pval"},
    {"sample": "D3", "condition": "DNMT3L_overexpression", "meth": "SAMPLE 6 Methylated signal", "unmeth": "SAMPLE 6 Unmethylated Signal", "pval": "SAMPLE 6 Detection Pval"},
]

MOUSE_SAMPLES = [
    {"sample": "B-3", "condition": "control", "sex": "female"},
    {"sample": "B-4", "condition": "control", "sex": "female"},
    {"sample": "B-5", "condition": "control", "sex": "female"},
    {"sample": "CC1rest", "condition": "control", "sex": "male"},
    {"sample": "CC3rest", "condition": "control", "sex": "male"},
    {"sample": "CF2rest", "condition": "control", "sex": "male"},
    {"sample": "CC2rest", "condition": "Dnmt3l_overexpression", "sex": "male"},
    {"sample": "CD5rest", "condition": "Dnmt3l_overexpression", "sex": "male"},
    {"sample": "CF1rest", "condition": "Dnmt3l_overexpression", "sex": "male"},
    {"sample": "D-1", "condition": "Dnmt3l_overexpression", "sex": "female"},
    {"sample": "CD6rest", "condition": "Dnmt3l_overexpression", "sex": "female"},
    {"sample": "CD16rest", "condition": "Dnmt3l_overexpression", "sex": "female"},
]

COMPARISONS = {
    "GSE157078": [
        {"comparison": "human_all_control_vs_DNMT3L_overexpression", "control": ["Z1", "Z2", "Z3"], "case": ["D1", "D2", "D3"]},
    ],
    "GSE271949": [
        {
            "comparison": "mouse_all_control_vs_Dnmt3l_overexpression",
            "control": ["B-3", "B-4", "B-5", "CC1rest", "CC3rest", "CF2rest"],
            "case": ["CC2rest", "CD5rest", "CF1rest", "D-1", "CD6rest", "CD16rest"],
        },
        {
            "comparison": "mouse_male_control_vs_Dnmt3l_overexpression",
            "control": ["CC1rest", "CC3rest", "CF2rest"],
            "case": ["CC2rest", "CD5rest", "CF1rest"],
        },
        {
            "comparison": "mouse_female_control_vs_Dnmt3l_overexpression",
            "control": ["B-3", "B-4", "B-5"],
            "case": ["D-1", "CD6rest", "CD16rest"],
        },
    ],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--human-input", default="data/geo/GSE157078/supplementary/GSE157078_signals.txt.gz")
    parser.add_argument("--mouse-input", default="data/geo/GSE271949/supplementary/GSE271949_all_sample.CpG.txt.gz")
    parser.add_argument("--outdir", default="results/first_pass_methylation")
    parser.add_argument("--top-n", type=int, default=100)
    parser.add_argument("--detection-p", type=float, default=0.01)
    parser.add_argument("--beta-offset", type=float, default=100.0)
    return parser.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="", errors="replace")
    return path.open(newline="", errors="replace")


def safe_float(value: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return number if math.isfinite(number) else float("nan")


def mean(values: list[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    return sum(finite) / len(finite) if finite else float("nan")


def fmt(value: float) -> str:
    return f"{value:.6g}" if math.isfinite(value) else ""


def split_chrbase(value: str) -> tuple[str, str]:
    if "." not in value:
        return "", ""
    chrom, pos = value.rsplit(".", 1)
    return chrom, pos


def update_summary(summary: dict[str, float], delta: float) -> None:
    summary["n_loci"] += 1
    summary["sum_delta_beta"] += delta
    for threshold in (0.05, 0.10, 0.20):
        if delta >= threshold:
            summary[f"hyper_delta_ge_{threshold:.2f}"] += 1
        if delta <= -threshold:
            summary[f"hypo_delta_le_neg_{threshold:.2f}"] += 1


def push_top(heap: list[tuple[float, int, dict[str, str]]], row: dict[str, str], top_n: int, counter: int) -> None:
    score = safe_float(row["abs_delta_beta"])
    item = (score, counter, row)
    if len(heap) < top_n:
        heapq.heappush(heap, item)
    elif score > heap[0][0]:
        heapq.heapreplace(heap, item)


def make_delta_row(dataset: str, comparison: dict[str, list[str]], locus_id: str, values: dict[str, float]) -> dict[str, str] | None:
    control_values = [values[sample] for sample in comparison["control"]]
    case_values = [values[sample] for sample in comparison["case"]]
    if any(not math.isfinite(value) for value in control_values + case_values):
        return None
    control_mean = mean(control_values)
    case_mean = mean(case_values)
    delta = case_mean - control_mean
    chrom, pos = split_chrbase(locus_id) if dataset == "GSE271949" else ("", "")
    return {
        "dataset": dataset,
        "comparison": comparison["comparison"],
        "locus_id": locus_id,
        "chrom": chrom,
        "position": pos,
        "mean_control_beta": fmt(control_mean),
        "mean_case_beta": fmt(case_mean),
        "delta_beta": fmt(delta),
        "abs_delta_beta": fmt(abs(delta)),
        "control_values_beta": ",".join(fmt(value) for value in control_values),
        "case_values_beta": ",".join(fmt(value) for value in case_values),
    }


def human_beta_values(row: dict[str, str], detection_p: float, beta_offset: float) -> dict[str, float] | None:
    values = {}
    for sample in HUMAN_SAMPLES:
        pvalue = safe_float(row.get(sample["pval"], ""))
        if not math.isfinite(pvalue) or pvalue > detection_p:
            return None
        methylated = safe_float(row.get(sample["meth"], ""))
        unmethylated = safe_float(row.get(sample["unmeth"], ""))
        denominator = methylated + unmethylated + beta_offset
        if not math.isfinite(denominator) or denominator <= 0:
            return None
        values[sample["sample"]] = methylated / denominator
    return values


def mouse_beta_values(row: dict[str, str]) -> dict[str, float] | None:
    values = {}
    for sample in MOUSE_SAMPLES:
        value = safe_float(row.get(sample["sample"], ""))
        if not math.isfinite(value):
            return None
        values[sample["sample"]] = value / 100.0
    return values


def analyze_table(
    dataset: str,
    path: Path,
    outdir: Path,
    top_n: int,
    value_parser,
    locus_column: str,
    detection_p: float,
    beta_offset: float,
) -> list[dict[str, str]]:
    summaries = {
        comparison["comparison"]: {
            "dataset": dataset,
            "comparison": comparison["comparison"],
            "n_loci": 0,
            "sum_delta_beta": 0.0,
            "hyper_delta_ge_0.05": 0,
            "hypo_delta_le_neg_0.05": 0,
            "hyper_delta_ge_0.10": 0,
            "hypo_delta_le_neg_0.10": 0,
            "hyper_delta_ge_0.20": 0,
            "hypo_delta_le_neg_0.20": 0,
        }
        for comparison in COMPARISONS[dataset]
    }
    top_heaps = {comparison["comparison"]: [] for comparison in COMPARISONS[dataset]}
    counter = 0

    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            values = value_parser(row, detection_p, beta_offset) if dataset == "GSE157078" else value_parser(row)
            if values is None:
                continue
            locus_id = row[locus_column]
            for comparison in COMPARISONS[dataset]:
                delta_row = make_delta_row(dataset, comparison, locus_id, values)
                if delta_row is None:
                    continue
                delta = safe_float(delta_row["delta_beta"])
                update_summary(summaries[comparison["comparison"]], delta)
                push_top(top_heaps[comparison["comparison"]], delta_row, top_n=top_n, counter=counter)
                counter += 1

    for comparison_name, heap in top_heaps.items():
        rows = [item[2] for item in sorted(heap, reverse=True)]
        write_tsv(
            outdir / f"{dataset}_{comparison_name}_top_abs_delta.tsv",
            rows,
            [
                "dataset",
                "comparison",
                "locus_id",
                "chrom",
                "position",
                "mean_control_beta",
                "mean_case_beta",
                "delta_beta",
                "abs_delta_beta",
                "control_values_beta",
                "case_values_beta",
            ],
        )

    summary_rows = []
    for summary in summaries.values():
        n_loci = int(summary["n_loci"])
        mean_delta = summary["sum_delta_beta"] / n_loci if n_loci else float("nan")
        summary_rows.append(
            {
                **{key: str(value) for key, value in summary.items() if key != "sum_delta_beta"},
                "mean_delta_beta": fmt(mean_delta),
            }
        )
    return summary_rows


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    all_summary_rows = []
    all_summary_rows.extend(
        analyze_table(
            dataset="GSE157078",
            path=Path(args.human_input),
            outdir=outdir,
            top_n=args.top_n,
            value_parser=human_beta_values,
            locus_column="ID_REF",
            detection_p=args.detection_p,
            beta_offset=args.beta_offset,
        )
    )
    all_summary_rows.extend(
        analyze_table(
            dataset="GSE271949",
            path=Path(args.mouse_input),
            outdir=outdir,
            top_n=args.top_n,
            value_parser=mouse_beta_values,
            locus_column="chrBase",
            detection_p=args.detection_p,
            beta_offset=args.beta_offset,
        )
    )
    write_tsv(
        outdir / "methylation_delta_summary.tsv",
        all_summary_rows,
        [
            "dataset",
            "comparison",
            "n_loci",
            "mean_delta_beta",
            "hyper_delta_ge_0.05",
            "hypo_delta_le_neg_0.05",
            "hyper_delta_ge_0.10",
            "hypo_delta_le_neg_0.10",
            "hyper_delta_ge_0.20",
            "hypo_delta_le_neg_0.20",
        ],
    )
    print(f"Wrote methylation summaries for {len(all_summary_rows)} comparisons to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
