#!/usr/bin/env python3
"""Summarize candidate DMR overlaps with chromatin annotation BED files."""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import defaultdict
from glob import glob
from pathlib import Path


DMR_SETS = [
    {
        "dataset": "GSE157078",
        "species": "human",
        "genome_build": "hg19",
        "label": "human_EPIC_candidate_DMRs",
        "path": "results/human_candidate_dmrs/GSE157078_candidate_dmrs_hg19.bed",
    },
    {
        "dataset": "GSE271949",
        "species": "mouse",
        "genome_build": "mm10",
        "label": "mouse_capture_candidate_DMRs_sex_consistent",
        "path": "results/mouse_candidate_dmrs_sex_consistent/GSE271949_candidate_dmrs.bed",
    },
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources", default="config/chromatin_overlap_sources.tsv")
    parser.add_argument("--outdir", default="results/chromatin_overlaps")
    return parser.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="", errors="replace")
    return path.open(newline="", errors="replace")


def chrom_sort_key(chrom: str) -> tuple[int, int | str]:
    name = chrom.removeprefix("chr")
    if name.isdigit():
        return (0, int(name))
    special = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    return (0, special[name]) if name in special else (1, name)


def parse_bed_line(line: str) -> tuple[str, int, int, list[str]] | None:
    if not line.strip() or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
        return None
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 3:
        return None
    try:
        start = int(float(parts[1]))
        end = int(float(parts[2]))
    except ValueError:
        return None
    if end <= start:
        return None
    return parts[0], start, end, parts


def load_intervals(path: Path) -> dict[str, list[tuple[int, int]]]:
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with open_text(path) as handle:
        for line in handle:
            parsed = parse_bed_line(line)
            if parsed is None:
                continue
            chrom, start, end, _parts = parsed
            intervals[chrom].append((start, end))
    for chrom in intervals:
        intervals[chrom].sort()
    return intervals


def load_dmrs(path: Path) -> list[dict[str, str | int]]:
    dmrs = []
    with path.open(newline="") as handle:
        for line in handle:
            parsed = parse_bed_line(line)
            if parsed is None:
                continue
            chrom, start, end, parts = parsed
            dmrs.append(
                {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "dmr_id": parts[3] if len(parts) > 3 else f"{chrom}:{start}-{end}",
                    "direction": parts[6] if len(parts) > 6 else "unknown",
                }
            )
    dmrs.sort(key=lambda row: (chrom_sort_key(str(row["chrom"])), int(row["start"])))
    return dmrs


def overlaps(intervals: list[tuple[int, int]], start: int, end: int, pointer: int) -> tuple[bool, int]:
    while pointer < len(intervals) and intervals[pointer][1] <= start:
        pointer += 1
    idx = pointer
    while idx < len(intervals) and intervals[idx][0] < end:
        if intervals[idx][1] > start:
            return True, pointer
        idx += 1
    return False, pointer


def summarize_overlap(dmrs: list[dict[str, str | int]], intervals_by_chrom: dict[str, list[tuple[int, int]]]) -> dict[str, int]:
    pointers = {chrom: 0 for chrom in intervals_by_chrom}
    summary = {
        "total_dmrs": len(dmrs),
        "overlapping_dmrs": 0,
        "hyper_total": 0,
        "hyper_overlapping": 0,
        "hypo_total": 0,
        "hypo_overlapping": 0,
    }
    for dmr in dmrs:
        direction = str(dmr["direction"])
        if direction == "hyper":
            summary["hyper_total"] += 1
        elif direction == "hypo":
            summary["hypo_total"] += 1

        chrom = str(dmr["chrom"])
        intervals = intervals_by_chrom.get(chrom)
        has_overlap = False
        if intervals:
            has_overlap, pointers[chrom] = overlaps(intervals, int(dmr["start"]), int(dmr["end"]), pointers.get(chrom, 0))
        if not has_overlap:
            continue
        summary["overlapping_dmrs"] += 1
        if direction == "hyper":
            summary["hyper_overlapping"] += 1
        elif direction == "hypo":
            summary["hypo_overlapping"] += 1
    return summary


def read_sources(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def label_for_file(source: dict[str, str], path: Path) -> str:
    matches = glob(source["path_glob"])
    if len(matches) == 1:
        return source["label"]
    name = path.name.removesuffix(".gz").removesuffix(".bed")
    return f"{source['label']}::{name}"


def fraction(numerator: int, denominator: int) -> str:
    if denominator == 0:
        return ""
    return f"{numerator / denominator:.6g}"


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    sources = read_sources(Path(args.sources))
    dmr_sets = [
        {**dmr_set, "dmrs": load_dmrs(Path(dmr_set["path"]))}
        for dmr_set in DMR_SETS
        if Path(dmr_set["path"]).exists()
    ]
    rows = []
    for source in sources:
        source_paths = [Path(path) for path in sorted(glob(source["path_glob"]))]
        if not source_paths:
            rows.append(
                {
                    "source_label": source["label"],
                    "annotation_file": source["path_glob"],
                    "status": "missing",
                }
            )
            continue
        for source_path in source_paths:
            intervals = load_intervals(source_path)
            interval_count = sum(len(values) for values in intervals.values())
            for dmr_set in dmr_sets:
                if dmr_set["species"] != source["species"] or dmr_set["genome_build"] != source["genome_build"]:
                    continue
                summary = summarize_overlap(dmr_set["dmrs"], intervals)
                rows.append(
                    {
                        "dmr_dataset": dmr_set["dataset"],
                        "dmr_label": dmr_set["label"],
                        "species": source["species"],
                        "genome_build": source["genome_build"],
                        "source_type": source["source_type"],
                        "source_label": label_for_file(source, source_path),
                        "annotation_file": str(source_path),
                        "annotation_intervals": str(interval_count),
                        "total_dmrs": str(summary["total_dmrs"]),
                        "overlapping_dmrs": str(summary["overlapping_dmrs"]),
                        "overlap_fraction": fraction(summary["overlapping_dmrs"], summary["total_dmrs"]),
                        "hyper_total": str(summary["hyper_total"]),
                        "hyper_overlapping": str(summary["hyper_overlapping"]),
                        "hyper_overlap_fraction": fraction(summary["hyper_overlapping"], summary["hyper_total"]),
                        "hypo_total": str(summary["hypo_total"]),
                        "hypo_overlapping": str(summary["hypo_overlapping"]),
                        "hypo_overlap_fraction": fraction(summary["hypo_overlapping"], summary["hypo_total"]),
                        "status": "ok",
                    }
                )
    outdir = Path(args.outdir)
    write_tsv(
        outdir / "dmr_chromatin_overlap_summary.tsv",
        rows,
        [
            "dmr_dataset",
            "dmr_label",
            "species",
            "genome_build",
            "source_type",
            "source_label",
            "annotation_file",
            "annotation_intervals",
            "total_dmrs",
            "overlapping_dmrs",
            "overlap_fraction",
            "hyper_total",
            "hyper_overlapping",
            "hyper_overlap_fraction",
            "hypo_total",
            "hypo_overlapping",
            "hypo_overlap_fraction",
            "status",
        ],
    )
    print(f"Wrote {len(rows)} overlap summaries to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
