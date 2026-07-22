#!/usr/bin/env python3
"""Prepare GSE305146 recurrent scATAC peak-call sets for HOMER motif enrichment."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import random
import re
from collections import defaultdict
from pathlib import Path


HOMER_BIN = os.environ.get("HOMER_FIND_MOTIFS", "findMotifsGenome.pl")
STANDARD_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
CHROM_RANK = {chrom: rank for rank, chrom in enumerate(STANDARD_CHROMS)}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", default="data/geo/GSE305146/metadata/sample_metadata.tsv")
    parser.add_argument("--supplementary", default="data/geo/GSE305146/supplementary")
    parser.add_argument("--outdir", default="results/gse305146_scATAC_homer_motif_enrichment")
    parser.add_argument("--genome", default="hg38")
    parser.add_argument("--merge-gap", type=int, default=250)
    parser.add_argument("--region-size", type=int, default=200)
    parser.add_argument("--top-n", type=int, default=5000)
    parser.add_argument("--background-max-regions", type=int, default=50000)
    parser.add_argument("--min-group-frequency", type=float, default=0.20)
    parser.add_argument("--min-delta-frequency", type=float, default=0.15)
    parser.add_argument("--min-background-frequency", type=float, default=0.10)
    parser.add_argument("--seed", type=int, default=305146)
    parser.add_argument("--homer-bin", default=HOMER_BIN)
    parser.add_argument("--homer-cpus", type=int, default=8)
    parser.add_argument("--homer-len", default="8,10,12")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", errors="replace") as handle:
        return [{k: (v or "").strip() for k, v in row.items()} for row in csv.DictReader(handle, delimiter="\t")]


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def parse_characteristics(text: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for part in text.split(";"):
        part = part.strip()
        if not part or "=" not in part:
            continue
        key, value = part.split("=", 1)
        values[key.strip().lower()] = value.strip()
    return values


def normalize_group(group: str) -> str:
    value = group.strip().lower()
    if value in {"ds", "trisomy 21", "t21", "down syndrome"}:
        return "DS"
    if value in {"con", "control", "disomic", "euploid"}:
        return "Control"
    return group.strip()


def selected_atac_libraries(metadata_path: Path) -> list[dict[str, str]]:
    libraries: list[dict[str, str]] = []
    for row in read_tsv(metadata_path):
        chars = parse_characteristics(row.get("characteristics", ""))
        if chars.get("library type") != "Chromatin Accessibility":
            continue
        if chars.get("in final_analysis", "").upper() == "FALSE":
            continue
        match = re.search(r"Library\s+([^,\s]+)_atac", row.get("title", ""))
        if not match:
            continue
        group = normalize_group(chars.get("group", chars.get("genotype", "")))
        if group not in {"DS", "Control"}:
            continue
        libraries.append(
            {
                "library_prefix": match.group(1),
                "geo_sample_accession": row.get("accession", ""),
                "biological_sample_id": chars.get("sample id", ""),
                "sample_name": chars.get("name", ""),
                "group": group,
                "genotype": chars.get("genotype", ""),
                "dev_stage_pcw": chars.get("dev stage_pcw", ""),
                "sex": chars.get("sex", ""),
            }
        )
    return libraries


def safe_int(value: str) -> int:
    try:
        return int(float(value))
    except ValueError:
        return 0


def read_peak_intervals(
    feature_path: Path,
    group: str,
    bit_index: int,
    intervals_by_chrom: dict[str, list[tuple[int, int, int, int]]],
) -> int:
    ds_bit = (1 << bit_index) if group == "DS" else 0
    control_bit = (1 << bit_index) if group == "Control" else 0
    n_peaks = 0
    with gzip.open(feature_path, "rt", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6 or fields[2] != "Peaks":
                continue
            chrom = fields[3]
            if chrom not in CHROM_RANK:
                continue
            start = safe_int(fields[4])
            end = safe_int(fields[5])
            if end <= start:
                continue
            intervals_by_chrom[chrom].append((start, end, ds_bit, control_bit))
            n_peaks += 1
    return n_peaks


def merge_intervals(
    intervals_by_chrom: dict[str, list[tuple[int, int, int, int]]],
    merge_gap: int,
    n_ds: int,
    n_control: int,
) -> list[dict[str, object]]:
    merged: list[dict[str, object]] = []
    for chrom in STANDARD_CHROMS:
        intervals = intervals_by_chrom.get(chrom, [])
        if not intervals:
            continue
        intervals.sort(key=lambda row: (row[0], row[1]))
        current_start, current_end, current_ds_mask, current_control_mask = intervals[0]
        n_source_peaks = 1
        for start, end, ds_mask, control_mask in intervals[1:]:
            if start <= current_end + merge_gap:
                current_end = max(current_end, end)
                current_ds_mask |= ds_mask
                current_control_mask |= control_mask
                n_source_peaks += 1
                continue
            merged.append(region_row(chrom, current_start, current_end, current_ds_mask, current_control_mask, n_source_peaks, n_ds, n_control))
            current_start, current_end, current_ds_mask, current_control_mask = start, end, ds_mask, control_mask
            n_source_peaks = 1
        merged.append(region_row(chrom, current_start, current_end, current_ds_mask, current_control_mask, n_source_peaks, n_ds, n_control))
    return merged


def region_row(
    chrom: str,
    start: int,
    end: int,
    ds_mask: int,
    control_mask: int,
    n_source_peaks: int,
    n_ds: int,
    n_control: int,
) -> dict[str, object]:
    ds_count = int(ds_mask).bit_count()
    control_count = int(control_mask).bit_count()
    ds_frequency = ds_count / n_ds if n_ds else 0.0
    control_frequency = control_count / n_control if n_control else 0.0
    return {
        "chrom": chrom,
        "start": start,
        "end": end,
        "width": end - start,
        "center": (start + end) // 2,
        "n_source_peaks": n_source_peaks,
        "ds_library_count": ds_count,
        "control_library_count": control_count,
        "ds_frequency": ds_frequency,
        "control_frequency": control_frequency,
        "delta_frequency_ds_minus_control": ds_frequency - control_frequency,
        "max_frequency": max(ds_frequency, control_frequency),
    }


def sort_regions(regions: list[dict[str, object]]) -> list[dict[str, object]]:
    return sorted(regions, key=lambda row: (CHROM_RANK[str(row["chrom"])], int(row["start"]), int(row["end"])))


def select_directional_regions(
    merged: list[dict[str, object]],
    direction: str,
    min_group_frequency: float,
    min_delta_frequency: float,
    top_n: int,
) -> list[dict[str, object]]:
    if direction == "DS_open":
        selected = [
            row
            for row in merged
            if float(row["ds_frequency"]) >= min_group_frequency
            and float(row["delta_frequency_ds_minus_control"]) >= min_delta_frequency
        ]
        selected.sort(
            key=lambda row: (
                -float(row["delta_frequency_ds_minus_control"]),
                -float(row["ds_frequency"]),
                float(row["control_frequency"]),
                CHROM_RANK[str(row["chrom"])],
                int(row["start"]),
            )
        )
    elif direction == "Control_open":
        selected = [
            row
            for row in merged
            if float(row["control_frequency"]) >= min_group_frequency
            and -float(row["delta_frequency_ds_minus_control"]) >= min_delta_frequency
        ]
        selected.sort(
            key=lambda row: (
                float(row["delta_frequency_ds_minus_control"]),
                -float(row["control_frequency"]),
                float(row["ds_frequency"]),
                CHROM_RANK[str(row["chrom"])],
                int(row["start"]),
            )
        )
    else:
        raise ValueError(f"Unknown direction: {direction}")

    if selected:
        return selected[:top_n]

    if direction == "DS_open":
        fallback = sorted(
            (row for row in merged if float(row["delta_frequency_ds_minus_control"]) > 0),
            key=lambda row: (-float(row["delta_frequency_ds_minus_control"]), -float(row["ds_frequency"])),
        )
    else:
        fallback = sorted(
            (row for row in merged if float(row["delta_frequency_ds_minus_control"]) < 0),
            key=lambda row: (float(row["delta_frequency_ds_minus_control"]), -float(row["control_frequency"])),
        )
    return fallback[:top_n]


def select_background(
    merged: list[dict[str, object]],
    min_background_frequency: float,
    max_regions: int,
    seed: int,
) -> list[dict[str, object]]:
    background = [row for row in merged if float(row["max_frequency"]) >= min_background_frequency]
    background = sort_regions(background)
    if len(background) <= max_regions:
        return background
    rng = random.Random(seed)
    sampled = rng.sample(background, max_regions)
    return sort_regions(sampled)


def write_bed(path: Path, rows: list[dict[str, object]], region_size: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    half = region_size // 2
    with path.open("w", newline="") as handle:
        for idx, row in enumerate(sort_regions(rows), start=1):
            center = int(row["center"])
            start = max(0, center - half)
            end = start + region_size
            name = (
                f"{row['chrom']}:{row['start']}-{row['end']}"
                f"|ds={float(row['ds_frequency']):.3f}"
                f"|control={float(row['control_frequency']):.3f}"
            )
            handle.write(f"{row['chrom']}\t{start}\t{end}\t{name}\t0\t.\n")


def write_region_table(path: Path, rows: list[dict[str, object]]) -> None:
    fields = [
        "chrom",
        "start",
        "end",
        "width",
        "center",
        "n_source_peaks",
        "ds_library_count",
        "control_library_count",
        "ds_frequency",
        "control_frequency",
        "delta_frequency_ds_minus_control",
        "max_frequency",
    ]
    formatted = []
    for row in sort_regions(rows):
        row_out = dict(row)
        for key in ["ds_frequency", "control_frequency", "delta_frequency_ds_minus_control", "max_frequency"]:
            row_out[key] = f"{float(row[key]):.6g}"
        formatted.append(row_out)
    write_tsv(path, formatted, fields)


def homer_command(row: dict[str, object], homer_bin: str, cpus: int, motif_len: str) -> str:
    return (
        f"{homer_bin} {row['input_bed']} {row['genome']} {row['homer_outdir']}"
        f" -size {row['window_size']} -mask -p {cpus} -len {motif_len} -bg {row['background_bed']}"
    )


def main() -> None:
    args = parse_args()
    metadata_path = Path(args.metadata)
    supplementary = Path(args.supplementary)
    outdir = Path(args.outdir)
    input_dir = outdir / "inputs"
    background_dir = outdir / "backgrounds"
    table_dir = outdir / "tables"
    homer_outdir = outdir / "homer_out"

    libraries = selected_atac_libraries(metadata_path)
    if not libraries:
        raise SystemExit(f"No final-analysis ATAC libraries found in {metadata_path}")
    n_ds = sum(1 for row in libraries if row["group"] == "DS")
    n_control = sum(1 for row in libraries if row["group"] == "Control")

    intervals_by_chrom: dict[str, list[tuple[int, int, int, int]]] = defaultdict(list)
    library_rows: list[dict[str, object]] = []
    for bit_index, library in enumerate(libraries):
        prefix = library["library_prefix"]
        feature_path = supplementary / f"GSE305146_{prefix}_features.tsv.gz"
        if not feature_path.exists():
            raise SystemExit(f"Missing feature file for {prefix}: {feature_path}")
        n_peaks = read_peak_intervals(feature_path, library["group"], bit_index, intervals_by_chrom)
        library_rows.append({**library, "feature_file": str(feature_path), "standard_chrom_peak_features": n_peaks})

    merged = merge_intervals(intervals_by_chrom, args.merge_gap, n_ds, n_control)
    ds_open = select_directional_regions(merged, "DS_open", args.min_group_frequency, args.min_delta_frequency, args.top_n)
    control_open = select_directional_regions(merged, "Control_open", args.min_group_frequency, args.min_delta_frequency, args.top_n)
    background = select_background(merged, args.min_background_frequency, args.background_max_regions, args.seed)

    ds_bed = input_dir / "GSE305146_DS_open_recurrent_peaks_hg38.bed"
    control_bed = input_dir / "GSE305146_Control_open_recurrent_peaks_hg38.bed"
    background_bed = background_dir / "GSE305146_recurrent_peak_background_hg38.bed"
    write_bed(ds_bed, ds_open, args.region_size)
    write_bed(control_bed, control_open, args.region_size)
    write_bed(background_bed, background, args.region_size)

    write_region_table(table_dir / "GSE305146_all_merged_recurrent_peak_calls.tsv", merged)
    write_region_table(table_dir / "GSE305146_DS_open_recurrent_peak_calls.tsv", ds_open)
    write_region_table(table_dir / "GSE305146_Control_open_recurrent_peak_calls.tsv", control_open)
    write_region_table(table_dir / "GSE305146_recurrent_peak_background.tsv", background)

    write_tsv(
        outdir / "library_manifest.tsv",
        library_rows,
        [
            "library_prefix",
            "geo_sample_accession",
            "biological_sample_id",
            "sample_name",
            "group",
            "genotype",
            "dev_stage_pcw",
            "sex",
            "standard_chrom_peak_features",
            "feature_file",
        ],
    )

    set_rows: list[dict[str, object]] = [
        {
            "set_id": "GSE305146_fetal_DS_open_recurrent_scATAC_peaks_hg38",
            "dataset": "GSE305146_fetal_DS_multiome",
            "species": "Homo_sapiens",
            "genome": args.genome,
            "level": "recurrent_scATAC_peak_call",
            "direction": "DS_open",
            "n_regions": len(ds_open),
            "window_size": args.region_size,
            "input_bed": str(ds_bed),
            "background_bed": str(background_bed),
            "homer_outdir": str(homer_outdir / "GSE305146_fetal_DS_open_recurrent_scATAC_peaks_hg38"),
        },
        {
            "set_id": "GSE305146_fetal_Control_open_recurrent_scATAC_peaks_hg38",
            "dataset": "GSE305146_fetal_DS_multiome",
            "species": "Homo_sapiens",
            "genome": args.genome,
            "level": "recurrent_scATAC_peak_call",
            "direction": "Control_open",
            "n_regions": len(control_open),
            "window_size": args.region_size,
            "input_bed": str(control_bed),
            "background_bed": str(background_bed),
            "homer_outdir": str(homer_outdir / "GSE305146_fetal_Control_open_recurrent_scATAC_peaks_hg38"),
        },
    ]
    for row in set_rows:
        row["command"] = homer_command(row, args.homer_bin, args.homer_cpus, args.homer_len)

    write_tsv(
        outdir / "homer_command_manifest.tsv",
        set_rows,
        [
            "set_id",
            "dataset",
            "species",
            "genome",
            "level",
            "direction",
            "n_regions",
            "window_size",
            "input_bed",
            "background_bed",
            "homer_outdir",
            "command",
        ],
    )

    summary_rows = [
        {"metric": "final_analysis_atac_libraries", "value": len(libraries)},
        {"metric": "ds_libraries", "value": n_ds},
        {"metric": "control_libraries", "value": n_control},
        {"metric": "merged_standard_chrom_peak_regions", "value": len(merged)},
        {"metric": "ds_open_regions", "value": len(ds_open)},
        {"metric": "control_open_regions", "value": len(control_open)},
        {"metric": "background_regions", "value": len(background)},
        {"metric": "merge_gap_bp", "value": args.merge_gap},
        {"metric": "region_size_bp", "value": args.region_size},
        {"metric": "min_group_frequency", "value": args.min_group_frequency},
        {"metric": "min_delta_frequency", "value": args.min_delta_frequency},
        {"metric": "min_background_frequency", "value": args.min_background_frequency},
        {"metric": "background_seed", "value": args.seed},
    ]
    write_tsv(outdir / "peak_set_summary.tsv", summary_rows, ["metric", "value"])

    print(f"Wrote HOMER inputs to {outdir}")
    print(f"DS-open regions: {len(ds_open)}")
    print(f"Control-open regions: {len(control_open)}")
    print(f"Background regions: {len(background)}")


if __name__ == "__main__":
    main()
