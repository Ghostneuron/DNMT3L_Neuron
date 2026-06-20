#!/usr/bin/env python3
"""Compare Dp10/Dp16 mouse WGBS methylation shifts with Dnmt3l-OE CpGs."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path


DEFAULT_DNMT3L_DMCS = "results/mouse_candidate_dmrs_sex_consistent/GSE271949_candidate_dmcs.tsv"
DEFAULT_OUTDIR = "results/dp10_methylation_dnmt3l_concordance"


TRACKS = [
    {
        "model": "Dp10",
        "path": "data/geo/GSE74505/supplementary/raw_wgbs/GSM1921944_dp10_vs_control_differential_methylation.bedgraph.gz",
    },
    {
        "model": "Dp16",
        "path": "data/geo/GSE74505/supplementary/raw_wgbs/GSM1921945_dp16_vs_control_differential_methylation.bedgraph.gz",
    },
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dnmt3l-dmcs", default=DEFAULT_DNMT3L_DMCS)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
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


def median(values: list[float]) -> float:
    finite = sorted(value for value in values if math.isfinite(value))
    if not finite:
        return float("nan")
    midpoint = len(finite) // 2
    if len(finite) % 2:
        return finite[midpoint]
    return (finite[midpoint - 1] + finite[midpoint]) / 2


def fmt(value: float) -> str:
    return f"{value:.6g}" if math.isfinite(value) else ""


def read_dnmt3l_dmcs(path: Path) -> dict[tuple[str, int, int], dict[str, str]]:
    dmcs = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            key = (row["chrom"], int(row["start"]), int(row["end"]))
            dmcs[key] = row
    return dmcs


def read_track_overlaps(model: str, path: Path, dmcs: dict[tuple[str, int, int], dict[str, str]]) -> list[dict[str, str]]:
    rows = []
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, value = line.rstrip("\n").split("\t")[:4]
            key = (chrom, int(start), int(end))
            dmc = dmcs.get(key)
            if not dmc:
                continue
            model_delta = safe_float(value)
            expected_sign = 1 if dmc["direction"] == "hyper" else -1
            signed_model_delta = model_delta * expected_sign
            rows.append(
                {
                    "model": model,
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "locus_id": dmc["locus_id"],
                    "dnmt3l_direction": dmc["direction"],
                    "dnmt3l_delta_beta": dmc["delta_beta"],
                    "dnmt3l_male_delta_beta": dmc["male_delta_beta"],
                    "dnmt3l_female_delta_beta": dmc["female_delta_beta"],
                    "model_delta_beta": fmt(model_delta),
                    "signed_model_delta_by_dnmt3l_direction": fmt(signed_model_delta),
                    "concordant": "yes" if signed_model_delta > 0 else "no",
                }
            )
    return rows


def summarize(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    by_model: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_model[row["model"]].append(row)
    summaries = []
    for model, model_rows in sorted(by_model.items()):
        dnmt3l_deltas = [safe_float(row["dnmt3l_delta_beta"]) for row in model_rows]
        model_deltas = [safe_float(row["model_delta_beta"]) for row in model_rows]
        signed = [safe_float(row["signed_model_delta_by_dnmt3l_direction"]) for row in model_rows]
        hyper_rows = [row for row in model_rows if row["dnmt3l_direction"] == "hyper"]
        hypo_rows = [row for row in model_rows if row["dnmt3l_direction"] == "hypo"]
        concordant = [value for value in signed if value > 0]
        discordant = [value for value in signed if value < 0]
        summaries.append(
            {
                "model": model,
                "n_overlap_cpgs": str(len(model_rows)),
                "n_dnmt3l_hyper_cpgs": str(len(hyper_rows)),
                "n_dnmt3l_hypo_cpgs": str(len(hypo_rows)),
                "mean_dnmt3l_delta_beta": fmt(mean(dnmt3l_deltas)),
                "mean_model_delta_beta": fmt(mean(model_deltas)),
                "median_model_delta_beta": fmt(median(model_deltas)),
                "mean_signed_model_delta_by_dnmt3l_direction": fmt(mean(signed)),
                "median_signed_model_delta_by_dnmt3l_direction": fmt(median(signed)),
                "concordant_n": str(len(concordant)),
                "discordant_n": str(len(discordant)),
                "concordant_fraction": fmt(len(concordant) / len(model_rows)) if model_rows else "",
                "same_direction_ge_0.02": str(sum(1 for value in signed if value >= 0.02)),
                "same_direction_ge_0.05": str(sum(1 for value in signed if value >= 0.05)),
                "opposite_direction_ge_0.02": str(sum(1 for value in signed if value <= -0.02)),
                "opposite_direction_ge_0.05": str(sum(1 for value in signed if value <= -0.05)),
                "mean_abs_model_delta_beta": fmt(mean([abs(value) for value in model_deltas])),
            }
        )
    return summaries


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    dmcs = read_dnmt3l_dmcs(Path(args.dnmt3l_dmcs))
    detail_rows = []
    for track in TRACKS:
        path = Path(track["path"])
        if path.exists():
            detail_rows.extend(read_track_overlaps(track["model"], path, dmcs))

    summary_rows = summarize(detail_rows)
    summary_fields = [
        "model",
        "n_overlap_cpgs",
        "n_dnmt3l_hyper_cpgs",
        "n_dnmt3l_hypo_cpgs",
        "mean_dnmt3l_delta_beta",
        "mean_model_delta_beta",
        "median_model_delta_beta",
        "mean_signed_model_delta_by_dnmt3l_direction",
        "median_signed_model_delta_by_dnmt3l_direction",
        "concordant_n",
        "discordant_n",
        "concordant_fraction",
        "same_direction_ge_0.02",
        "same_direction_ge_0.05",
        "opposite_direction_ge_0.02",
        "opposite_direction_ge_0.05",
        "mean_abs_model_delta_beta",
    ]
    detail_fields = [
        "model",
        "chrom",
        "start",
        "end",
        "locus_id",
        "dnmt3l_direction",
        "dnmt3l_delta_beta",
        "dnmt3l_male_delta_beta",
        "dnmt3l_female_delta_beta",
        "model_delta_beta",
        "signed_model_delta_by_dnmt3l_direction",
        "concordant",
    ]
    write_tsv(outdir / "dp10_dp16_dnmt3l_cpg_concordance_summary.tsv", summary_rows, summary_fields)
    write_tsv(outdir / "dp10_dp16_dnmt3l_cpg_values.tsv", detail_rows, detail_fields)
    print(f"Wrote {len(summary_rows)} model summaries")
    print(f"Wrote {len(detail_rows)} CpG overlap rows")


if __name__ == "__main__":
    main()
