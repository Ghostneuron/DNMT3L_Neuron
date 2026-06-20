#!/usr/bin/env python3
"""Call first-pass candidate DMCs and DMRs from mouse GSE271949 CpG methylation."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
from pathlib import Path


CONTROL_MALE = ["CC1rest", "CC3rest", "CF2rest"]
CONTROL_FEMALE = ["B-3", "B-4", "B-5"]
CASE_MALE = ["CC2rest", "CD5rest", "CF1rest"]
CASE_FEMALE = ["D-1", "CD6rest", "CD16rest"]
CONTROL_ALL = CONTROL_MALE + CONTROL_FEMALE
CASE_ALL = CASE_MALE + CASE_FEMALE


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="data/geo/GSE271949/supplementary/GSE271949_all_sample.CpG.txt.gz")
    parser.add_argument("--outdir", default="results/mouse_candidate_dmrs")
    parser.add_argument("--min-delta", type=float, default=0.10, help="Minimum absolute delta beta for candidate CpGs.")
    parser.add_argument("--max-gap", type=int, default=500, help="Maximum gap between candidate CpGs in one DMR.")
    parser.add_argument("--min-cpgs", type=int, default=3, help="Minimum candidate CpGs in one DMR.")
    parser.add_argument(
        "--require-sex-consistent-cpgs",
        action="store_true",
        help="Only group CpGs where male and female deltas have the same direction as the all-sample delta.",
    )
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


def split_chrbase(value: str) -> tuple[str, int]:
    chrom, pos = value.rsplit(".", 1)
    return chrom, int(pos)


def chrom_sort_key(chrom: str) -> tuple[int, int | str]:
    name = chrom.removeprefix("chr")
    if name.isdigit():
        return (0, int(name))
    special = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    return (0, special[name]) if name in special else (1, name)


def beta_values(row: dict[str, str], samples: list[str]) -> list[float]:
    return [safe_float(row.get(sample, "")) / 100.0 for sample in samples]


def direction(delta: float, min_delta: float) -> str:
    if delta >= min_delta:
        return "hyper"
    if delta <= -min_delta:
        return "hypo"
    return "none"


def direction_matches(delta: float, expected: str) -> bool:
    if expected == "hyper":
        return delta > 0
    if expected == "hypo":
        return delta < 0
    return False


def read_candidate_dmcs(path: Path, min_delta: float) -> list[dict[str, object]]:
    candidates: list[dict[str, object]] = []
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            locus = row["chrBase"]
            chrom, pos = split_chrbase(locus)

            control_all = beta_values(row, CONTROL_ALL)
            case_all = beta_values(row, CASE_ALL)
            control_male = beta_values(row, CONTROL_MALE)
            case_male = beta_values(row, CASE_MALE)
            control_female = beta_values(row, CONTROL_FEMALE)
            case_female = beta_values(row, CASE_FEMALE)

            mean_control = mean(control_all)
            mean_case = mean(case_all)
            delta = mean_case - mean_control
            call = direction(delta, min_delta)
            if call == "none":
                continue

            male_delta = mean(case_male) - mean(control_male)
            female_delta = mean(case_female) - mean(control_female)
            sex_consistent = direction_matches(male_delta, call) and direction_matches(female_delta, call)
            candidates.append(
                {
                    "chrom": chrom,
                    "position": pos,
                    "start": pos - 1,
                    "end": pos,
                    "locus_id": locus,
                    "direction": call,
                    "mean_control_beta": mean_control,
                    "mean_case_beta": mean_case,
                    "delta_beta": delta,
                    "mean_control_male_beta": mean(control_male),
                    "mean_case_male_beta": mean(case_male),
                    "male_delta_beta": male_delta,
                    "mean_control_female_beta": mean(control_female),
                    "mean_case_female_beta": mean(case_female),
                    "female_delta_beta": female_delta,
                    "sex_consistent": sex_consistent,
                }
            )
    candidates.sort(key=lambda row: (chrom_sort_key(str(row["chrom"])), int(row["position"])))
    return candidates


def cpg_for_output(row: dict[str, object]) -> dict[str, str]:
    return {
        "chrom": str(row["chrom"]),
        "start": str(row["start"]),
        "end": str(row["end"]),
        "locus_id": str(row["locus_id"]),
        "direction": str(row["direction"]),
        "mean_control_beta": fmt(float(row["mean_control_beta"])),
        "mean_case_beta": fmt(float(row["mean_case_beta"])),
        "delta_beta": fmt(float(row["delta_beta"])),
        "mean_control_male_beta": fmt(float(row["mean_control_male_beta"])),
        "mean_case_male_beta": fmt(float(row["mean_case_male_beta"])),
        "male_delta_beta": fmt(float(row["male_delta_beta"])),
        "mean_control_female_beta": fmt(float(row["mean_control_female_beta"])),
        "mean_case_female_beta": fmt(float(row["mean_case_female_beta"])),
        "female_delta_beta": fmt(float(row["female_delta_beta"])),
        "sex_consistent": "yes" if row["sex_consistent"] else "no",
    }


def should_group(prev: dict[str, object], row: dict[str, object], max_gap: int) -> bool:
    return (
        prev["chrom"] == row["chrom"]
        and prev["direction"] == row["direction"]
        and int(row["position"]) - int(prev["position"]) <= max_gap
    )


def summarize_region(rows: list[dict[str, object]], region_id: int) -> dict[str, str]:
    chrom = str(rows[0]["chrom"])
    start = min(int(row["start"]) for row in rows)
    end = max(int(row["end"]) for row in rows)
    deltas = [float(row["delta_beta"]) for row in rows]
    male_deltas = [float(row["male_delta_beta"]) for row in rows]
    female_deltas = [float(row["female_delta_beta"]) for row in rows]
    control = [float(row["mean_control_beta"]) for row in rows]
    case = [float(row["mean_case_beta"]) for row in rows]
    direction_value = str(rows[0]["direction"])
    sex_consistent_count = sum(1 for row in rows if row["sex_consistent"])
    return {
        "dmr_id": f"GSE271949_DMR_{region_id:06d}",
        "chrom": chrom,
        "start": str(start),
        "end": str(end),
        "direction": direction_value,
        "n_cpgs": str(len(rows)),
        "width_bp": str(end - start),
        "mean_control_beta": fmt(mean(control)),
        "mean_case_beta": fmt(mean(case)),
        "mean_delta_beta": fmt(mean(deltas)),
        "max_abs_delta_beta": fmt(max(abs(delta) for delta in deltas)),
        "mean_male_delta_beta": fmt(mean(male_deltas)),
        "mean_female_delta_beta": fmt(mean(female_deltas)),
        "sex_consistent_cpgs": str(sex_consistent_count),
        "sex_consistent_fraction": fmt(sex_consistent_count / len(rows)),
        "first_locus_id": str(rows[0]["locus_id"]),
        "last_locus_id": str(rows[-1]["locus_id"]),
    }


def call_dmrs(
    candidates: list[dict[str, object]],
    max_gap: int,
    min_cpgs: int,
    require_sex_consistent_cpgs: bool,
) -> list[dict[str, str]]:
    regions: list[dict[str, str]] = []
    current: list[dict[str, object]] = []
    region_id = 1

    grouped_candidates = [row for row in candidates if row["sex_consistent"]] if require_sex_consistent_cpgs else candidates
    for row in grouped_candidates:
        if not current:
            current = [row]
            continue
        if should_group(current[-1], row, max_gap):
            current.append(row)
            continue
        if len(current) >= min_cpgs:
            regions.append(summarize_region(current, region_id))
            region_id += 1
        current = [row]

    if len(current) >= min_cpgs:
        regions.append(summarize_region(current, region_id))
    return regions


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_bed(path: Path, regions: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        for region in regions:
            score = min(1000, int(abs(float(region["mean_delta_beta"])) * 1000))
            strand = "."
            handle.write(
                "\t".join(
                    [
                        region["chrom"],
                        region["start"],
                        region["end"],
                        region["dmr_id"],
                        str(score),
                        strand,
                        region["direction"],
                        region["n_cpgs"],
                        region["mean_delta_beta"],
                    ]
                )
                + "\n"
            )


def write_summary(path: Path, candidates: list[dict[str, object]], regions: list[dict[str, str]], args: argparse.Namespace) -> None:
    hyper_cpgs = sum(1 for row in candidates if row["direction"] == "hyper")
    hypo_cpgs = sum(1 for row in candidates if row["direction"] == "hypo")
    sex_consistent = sum(1 for row in candidates if row["sex_consistent"])
    hyper_regions = sum(1 for row in regions if row["direction"] == "hyper")
    hypo_regions = sum(1 for row in regions if row["direction"] == "hypo")
    rows = [
        {
            "dataset": "GSE271949",
            "min_delta": fmt(args.min_delta),
            "max_gap": str(args.max_gap),
            "min_cpgs": str(args.min_cpgs),
            "require_sex_consistent_cpgs": "yes" if args.require_sex_consistent_cpgs else "no",
            "candidate_cpgs": str(len(candidates)),
            "hyper_candidate_cpgs": str(hyper_cpgs),
            "hypo_candidate_cpgs": str(hypo_cpgs),
            "sex_consistent_candidate_cpgs": str(sex_consistent),
            "candidate_dmrs": str(len(regions)),
            "hyper_candidate_dmrs": str(hyper_regions),
            "hypo_candidate_dmrs": str(hypo_regions),
        }
    ]
    write_tsv(
        path,
        rows,
        [
            "dataset",
            "min_delta",
            "max_gap",
            "min_cpgs",
            "require_sex_consistent_cpgs",
            "candidate_cpgs",
            "hyper_candidate_cpgs",
            "hypo_candidate_cpgs",
            "sex_consistent_candidate_cpgs",
            "candidate_dmrs",
            "hyper_candidate_dmrs",
            "hypo_candidate_dmrs",
        ],
    )


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    candidates = read_candidate_dmcs(Path(args.input), min_delta=args.min_delta)
    regions = call_dmrs(
        candidates,
        max_gap=args.max_gap,
        min_cpgs=args.min_cpgs,
        require_sex_consistent_cpgs=args.require_sex_consistent_cpgs,
    )

    cpg_fields = [
        "chrom",
        "start",
        "end",
        "locus_id",
        "direction",
        "mean_control_beta",
        "mean_case_beta",
        "delta_beta",
        "mean_control_male_beta",
        "mean_case_male_beta",
        "male_delta_beta",
        "mean_control_female_beta",
        "mean_case_female_beta",
        "female_delta_beta",
        "sex_consistent",
    ]
    region_fields = [
        "dmr_id",
        "chrom",
        "start",
        "end",
        "direction",
        "n_cpgs",
        "width_bp",
        "mean_control_beta",
        "mean_case_beta",
        "mean_delta_beta",
        "max_abs_delta_beta",
        "mean_male_delta_beta",
        "mean_female_delta_beta",
        "sex_consistent_cpgs",
        "sex_consistent_fraction",
        "first_locus_id",
        "last_locus_id",
    ]
    write_tsv(outdir / "GSE271949_candidate_dmcs.tsv", [cpg_for_output(row) for row in candidates], cpg_fields)
    write_tsv(outdir / "GSE271949_candidate_dmrs.tsv", regions, region_fields)
    write_bed(outdir / "GSE271949_candidate_dmrs.bed", regions)
    write_summary(outdir / "GSE271949_candidate_dmr_summary.tsv", candidates, regions, args)

    print(
        f"Wrote {len(candidates)} candidate CpGs and {len(regions)} candidate DMRs to {outdir}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
