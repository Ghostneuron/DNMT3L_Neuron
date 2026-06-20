#!/usr/bin/env python3
"""Prepare HOMER inputs for DS/control brain methylation contrasts.

This is a first-pass, array-based DS methylation motif workflow.  It computes
context-specific DS-control delta beta values, writes effect-size DMC tables,
clusters nearby same-direction DMCs into DMR-like regions, and emits HOMER
foreground/background BEDs plus a command manifest.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
import shlex
from collections import defaultdict
from pathlib import Path

from compare_ds_methylation_to_dnmt3l_signature import (
    build_context_specs,
    eligible_contexts,
    gse63347_samples,
    gse73747_samples,
    gse74486_xml_samples,
    mean_beta_for_processed_row,
    mean_beta_for_signal_row,
    open_text,
    processed_matrix_sample_specs,
    safe_float,
    signal_matrix_sample_specs,
    slug,
)
from prepare_homer_motif_inputs import chrom_sort_key, read_human_epic_coords, write_bed


DEFAULT_CONTEXTS = [
    "GSE73747_fetal_all",
    "GSE73747_fetal_frontal",
    "GSE73747_fetal_temporal",
    "GSE74486_fetus_cerebrum_whole_cerebrum",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="results/ds_homer_motif_enrichment")
    parser.add_argument("--human-manifest", default="data/geo/GPL21145/GPL21145_MethylationEPIC_15073387_v-1-0.csv.gz")
    parser.add_argument("--genome", default="hg19")
    parser.add_argument("--window-size", type=int, default=200)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--motif-lengths", default="8,10")
    parser.add_argument("--min-group-samples", type=int, default=2)
    parser.add_argument("--min-detected-fraction", type=float, default=0.7)
    parser.add_argument("--detection-p", type=float, default=0.01)
    parser.add_argument("--gse74486-detection-p", type=float, default=0.05)
    parser.add_argument("--beta-offset", type=float, default=100.0)
    parser.add_argument("--min-abs-delta", type=float, default=0.05)
    parser.add_argument("--cluster-gap-bp", type=int, default=500)
    parser.add_argument("--min-dmr-probes", type=int, default=2)
    parser.add_argument(
        "--context",
        action="append",
        default=[],
        help="Context label to include. Can be repeated. Defaults to fetal-priority contexts.",
    )
    parser.add_argument(
        "--include-dmc-commands",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Also include DMC-level HOMER commands in the manifest.",
    )
    return parser.parse_args()


def fmt(value: float) -> str:
    return f"{value:.6g}" if math.isfinite(value) else ""


def mean(values: list[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    return sum(finite) / len(finite) if finite else float("nan")


def sample_threshold(n_samples: int, fraction: float, minimum: int) -> int:
    return max(minimum, math.ceil(n_samples * fraction))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def dataset_specs(args: argparse.Namespace) -> list[dict[str, object]]:
    return [
        {
            "dataset": "GSE73747",
            "kind": "signal",
            "matrix_path": Path("data/geo/GSE73747/supplementary/GSE73747_unmethylated_methylated.txt.gz"),
            "samples": gse73747_samples(Path("data/geo/GSE73747/metadata/sample_metadata.tsv")),
            "detection_p": args.detection_p,
            "compact_headers": False,
        },
        {
            "dataset": "GSE63347",
            "kind": "signal",
            "matrix_path": Path("data/geo/GSE63347/supplementary/GSE63347_datSignal.csv.gz"),
            "samples": gse63347_samples(Path("data/geo/GSE63347/metadata/sample_metadata.tsv")),
            "detection_p": args.detection_p,
            "compact_headers": True,
        },
        {
            "dataset": "GSE74486",
            "kind": "processed",
            "matrix_path": Path(
                "data/geo/GSE74486/supplementary/"
                "GSE74486_DS_processed_norm_batch_corrected_filtered_matrix_2.txt.gz"
            ),
            "samples": gse74486_xml_samples(Path("data/geo/GSE74486/metadata/GSE74486_family.xml")),
            "detection_p": args.gse74486_detection_p,
        },
    ]


def context_sample_counts(contexts: dict[str, dict[str, list[dict[str, str]]]]) -> dict[str, dict[str, int]]:
    return {
        context: {"Control": len(groups["Control"]), "DS": len(groups["DS"])}
        for context, groups in contexts.items()
    }


def context_label_for_set_id(context: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", context).strip("_")


def make_dmc_row(
    dataset: str,
    context: str,
    probe_id: str,
    coord: tuple[str, int, int],
    control_mean: float,
    ds_mean: float,
    control_n: int,
    ds_n: int,
) -> dict[str, str]:
    delta = ds_mean - control_mean
    chrom, start, end = coord
    direction = "hyper" if delta > 0 else "hypo"
    return {
        "dataset": dataset,
        "context": context,
        "probe_id": probe_id,
        "chrom": chrom,
        "start": str(start),
        "end": str(end),
        "direction": direction,
        "delta_beta": fmt(delta),
        "abs_delta_beta": fmt(abs(delta)),
        "control_mean_beta": fmt(control_mean),
        "ds_mean_beta": fmt(ds_mean),
        "control_n_detected": str(control_n),
        "ds_n_detected": str(ds_n),
    }


def scan_signal_dataset(
    spec: dict[str, object],
    coords: dict[str, tuple[str, int, int]],
    selected_contexts: set[str],
    args: argparse.Namespace,
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    dataset = str(spec["dataset"])
    samples = list(spec["samples"])  # type: ignore[arg-type]
    contexts = eligible_contexts(samples, args.min_group_samples)
    contexts = {context: groups for context, groups in contexts.items() if context in selected_contexts}
    counts = context_sample_counts(contexts)
    background_rows: dict[str, list[tuple[str, int, int, str, int, str]]] = defaultdict(list)
    dmc_rows: list[dict[str, str]] = []

    with open_text(Path(spec["matrix_path"])) as handle:  # type: ignore[arg-type]
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        sample_specs = signal_matrix_sample_specs(header, samples, compact=bool(spec["compact_headers"]))
        context_specs = build_context_specs(contexts, sample_specs)
        for row in reader:
            probe_id = row[0]
            coord = coords.get(probe_id)
            if not coord:
                continue
            for context, groups in context_specs.items():
                control_mean, control_n = mean_beta_for_signal_row(
                    row,
                    groups["Control"],
                    float(spec["detection_p"]),
                    args.beta_offset,
                )
                ds_mean, ds_n = mean_beta_for_signal_row(
                    row,
                    groups["DS"],
                    float(spec["detection_p"]),
                    args.beta_offset,
                )
                control_min = sample_threshold(counts[context]["Control"], args.min_detected_fraction, args.min_group_samples)
                ds_min = sample_threshold(counts[context]["DS"], args.min_detected_fraction, args.min_group_samples)
                if control_n < control_min or ds_n < ds_min:
                    continue
                background_rows[context].append((*coord, probe_id, 0, "."))
                delta = ds_mean - control_mean
                if math.isfinite(delta) and abs(delta) >= args.min_abs_delta:
                    dmc_rows.append(make_dmc_row(dataset, context, probe_id, coord, control_mean, ds_mean, control_n, ds_n))
    return dmc_rows, background_summary_rows(dataset, background_rows, counts), writeable_background_rows(background_rows)


def scan_processed_dataset(
    spec: dict[str, object],
    coords: dict[str, tuple[str, int, int]],
    selected_contexts: set[str],
    args: argparse.Namespace,
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    dataset = str(spec["dataset"])
    samples = list(spec["samples"])  # type: ignore[arg-type]
    contexts = eligible_contexts(samples, args.min_group_samples)
    contexts = {context: groups for context, groups in contexts.items() if context in selected_contexts}
    counts = context_sample_counts(contexts)
    background_rows: dict[str, list[tuple[str, int, int, str, int, str]]] = defaultdict(list)
    dmc_rows: list[dict[str, str]] = []

    with open_text(Path(spec["matrix_path"])) as handle:  # type: ignore[arg-type]
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        sample_specs = processed_matrix_sample_specs(header, samples)
        context_specs = build_context_specs(contexts, sample_specs)
        for row in reader:
            probe_id = row[0]
            coord = coords.get(probe_id)
            if not coord:
                continue
            for context, groups in context_specs.items():
                control_mean, control_n = mean_beta_for_processed_row(row, groups["Control"], float(spec["detection_p"]))
                ds_mean, ds_n = mean_beta_for_processed_row(row, groups["DS"], float(spec["detection_p"]))
                control_min = sample_threshold(counts[context]["Control"], args.min_detected_fraction, args.min_group_samples)
                ds_min = sample_threshold(counts[context]["DS"], args.min_detected_fraction, args.min_group_samples)
                if control_n < control_min or ds_n < ds_min:
                    continue
                background_rows[context].append((*coord, probe_id, 0, "."))
                delta = ds_mean - control_mean
                if math.isfinite(delta) and abs(delta) >= args.min_abs_delta:
                    dmc_rows.append(make_dmc_row(dataset, context, probe_id, coord, control_mean, ds_mean, control_n, ds_n))
    return dmc_rows, background_summary_rows(dataset, background_rows, counts), writeable_background_rows(background_rows)


def background_summary_rows(
    dataset: str,
    background_rows: dict[str, list[tuple[str, int, int, str, int, str]]],
    counts: dict[str, dict[str, int]],
) -> list[dict[str, str]]:
    rows = []
    for context, beds in sorted(background_rows.items()):
        rows.append(
            {
                "dataset": dataset,
                "context": context,
                "n_control": str(counts.get(context, {}).get("Control", "")),
                "n_ds": str(counts.get(context, {}).get("DS", "")),
                "n_background_probes": str(len(beds)),
            }
        )
    return rows


def writeable_background_rows(background_rows: dict[str, list[tuple[str, int, int, str, int, str]]]) -> list[dict[str, str]]:
    rows = []
    for context, beds in background_rows.items():
        for chrom, start, end, name, score, strand in beds:
            rows.append(
                {
                    "context": context,
                    "chrom": chrom,
                    "start": str(start),
                    "end": str(end),
                    "name": name,
                    "score": str(score),
                    "strand": strand,
                }
            )
    return rows


def cluster_dmrs(dmc_rows: list[dict[str, str]], max_gap: int, min_probes: int) -> list[dict[str, str]]:
    by_context_direction: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in dmc_rows:
        by_context_direction[(row["dataset"], row["context"], row["direction"])].append(row)

    dmr_rows: list[dict[str, str]] = []
    for (dataset, context, direction), rows in sorted(by_context_direction.items()):
        rows = sorted(rows, key=lambda row: (chrom_sort_key(row["chrom"]), int(row["start"]), int(row["end"])))
        clusters: list[list[dict[str, str]]] = []
        current: list[dict[str, str]] = []
        for row in rows:
            if (
                current
                and row["chrom"] == current[-1]["chrom"]
                and int(row["start"]) - int(current[-1]["end"]) <= max_gap
            ):
                current.append(row)
            else:
                if current:
                    clusters.append(current)
                current = [row]
        if current:
            clusters.append(current)

        kept_index = 0
        for cluster in clusters:
            if len(cluster) < min_probes:
                continue
            kept_index += 1
            starts = [int(row["start"]) for row in cluster]
            ends = [int(row["end"]) for row in cluster]
            deltas = [safe_float(row["delta_beta"]) for row in cluster]
            abs_deltas = [abs(value) for value in deltas if math.isfinite(value)]
            chrom = cluster[0]["chrom"]
            dmr_id = f"{context_label_for_set_id(context)}_{direction}_dmr_{kept_index:05d}"
            dmr_rows.append(
                {
                    "dataset": dataset,
                    "context": context,
                    "dmr_id": dmr_id,
                    "chrom": chrom,
                    "start": str(min(starts)),
                    "end": str(max(ends)),
                    "direction": direction,
                    "n_probes": str(len(cluster)),
                    "mean_delta_beta": fmt(mean(deltas)),
                    "mean_abs_delta_beta": fmt(mean(abs_deltas)),
                    "max_abs_delta_beta": fmt(max(abs_deltas) if abs_deltas else float("nan")),
                    "probe_ids": ";".join(row["probe_id"] for row in cluster),
                }
            )
    return dmr_rows


def bed_from_dmc(row: dict[str, str]) -> tuple[str, int, int, str, int, str]:
    score = min(1000, int(safe_float(row["abs_delta_beta"]) * 1000))
    return row["chrom"], int(row["start"]), int(row["end"]), row["probe_id"], score, "."


def bed_from_dmr(row: dict[str, str]) -> tuple[str, int, int, str, int, str]:
    score = min(1000, int(safe_float(row["mean_abs_delta_beta"]) * 1000))
    return row["chrom"], int(row["start"]), max(int(row["start"]) + 1, int(row["end"])), row["dmr_id"], score, "."


def split_rows(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    return {
        "all": rows,
        "hyper": [row for row in rows if row["direction"] == "hyper"],
        "hypo": [row for row in rows if row["direction"] == "hypo"],
    }


def write_foregrounds_and_commands(
    outdir: Path,
    rows: list[dict[str, str]],
    level: str,
    genome: str,
    background_paths: dict[str, Path],
    commands: list[dict[str, str]],
    window_size: int,
    threads: int,
    motif_lengths: str,
) -> None:
    input_dir = outdir / "inputs"
    homer_out_base = outdir / "homer_out"
    by_context: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_context[(row["dataset"], row["context"])].append(row)

    for (dataset, context), context_rows in sorted(by_context.items()):
        for direction, direction_rows in split_rows(context_rows).items():
            if not direction_rows:
                continue
            set_id = f"{dataset}_{context_label_for_set_id(context)}_{level}_{direction}_{genome}"
            bed_path = input_dir / f"{set_id}.bed"
            if level == "DMC":
                bed_rows = [bed_from_dmc(row) for row in direction_rows]
            else:
                bed_rows = [bed_from_dmr(row) for row in direction_rows]
            write_bed(bed_path, bed_rows)
            homer_outdir = homer_out_base / set_id
            background = background_paths[context]
            command_parts = [
                "findMotifsGenome.pl",
                str(bed_path),
                genome,
                str(homer_outdir),
                "-size",
                str(window_size),
                "-mask",
                "-p",
                str(threads),
                "-len",
                motif_lengths,
                "-bg",
                str(background),
            ]
            commands.append(
                {
                    "set_id": set_id,
                    "dataset": f"{dataset}_{context}",
                    "species": "Homo_sapiens",
                    "genome": genome,
                    "level": level,
                    "direction": direction,
                    "n_regions": str(len(bed_rows)),
                    "window_size": str(window_size),
                    "input_bed": str(bed_path),
                    "background_bed": str(background),
                    "homer_outdir": str(homer_outdir),
                    "command": " ".join(shlex.quote(part) for part in command_parts),
                }
            )


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    fields = [
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
    ]
    write_tsv(path, rows, fields)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    selected_contexts = set(args.context or DEFAULT_CONTEXTS)
    coords = read_human_epic_coords(Path(args.human_manifest))

    all_dmc_rows: list[dict[str, str]] = []
    all_background_summary_rows: list[dict[str, str]] = []
    background_bed_rows: dict[str, list[tuple[str, int, int, str, int, str]]] = defaultdict(list)

    for spec in dataset_specs(args):
        if not Path(spec["matrix_path"]).exists():  # type: ignore[arg-type]
            continue
        if spec["kind"] == "signal":
            dmc_rows, summary_rows, raw_bg_rows = scan_signal_dataset(spec, coords, selected_contexts, args)
        else:
            dmc_rows, summary_rows, raw_bg_rows = scan_processed_dataset(spec, coords, selected_contexts, args)
        all_dmc_rows.extend(dmc_rows)
        all_background_summary_rows.extend(summary_rows)
        for row in raw_bg_rows:
            background_bed_rows[row["context"]].append(
                (row["chrom"], int(row["start"]), int(row["end"]), row["name"], int(row["score"]), row["strand"])
            )

    dmr_rows = cluster_dmrs(all_dmc_rows, args.cluster_gap_bp, args.min_dmr_probes)

    background_dir = outdir / "backgrounds"
    background_paths: dict[str, Path] = {}
    for context, bed_rows in sorted(background_bed_rows.items()):
        path = background_dir / f"{context_label_for_set_id(context)}_background_{args.genome}.bed"
        write_bed(path, bed_rows)
        background_paths[context] = path

    commands: list[dict[str, str]] = []
    write_foregrounds_and_commands(
        outdir=outdir,
        rows=dmr_rows,
        level="DMR",
        genome=args.genome,
        background_paths=background_paths,
        commands=commands,
        window_size=args.window_size,
        threads=args.threads,
        motif_lengths=args.motif_lengths,
    )
    if args.include_dmc_commands:
        write_foregrounds_and_commands(
            outdir=outdir,
            rows=all_dmc_rows,
            level="DMC",
            genome=args.genome,
            background_paths=background_paths,
            commands=commands,
            window_size=args.window_size,
            threads=args.threads,
            motif_lengths=args.motif_lengths,
        )

    write_tsv(
        outdir / "ds_dmc_candidates.tsv",
        all_dmc_rows,
        [
            "dataset",
            "context",
            "probe_id",
            "chrom",
            "start",
            "end",
            "direction",
            "delta_beta",
            "abs_delta_beta",
            "control_mean_beta",
            "ds_mean_beta",
            "control_n_detected",
            "ds_n_detected",
        ],
    )
    write_tsv(
        outdir / "ds_dmr_like_regions.tsv",
        dmr_rows,
        [
            "dataset",
            "context",
            "dmr_id",
            "chrom",
            "start",
            "end",
            "direction",
            "n_probes",
            "mean_delta_beta",
            "mean_abs_delta_beta",
            "max_abs_delta_beta",
            "probe_ids",
        ],
    )
    write_tsv(
        outdir / "ds_context_background_summary.tsv",
        all_background_summary_rows,
        ["dataset", "context", "n_control", "n_ds", "n_background_probes"],
    )
    write_manifest(outdir / "homer_command_manifest.tsv", commands)

    print(f"Wrote {len(all_dmc_rows)} DS DMC candidate rows")
    print(f"Wrote {len(dmr_rows)} DS DMR-like region rows")
    print(f"Wrote {len(commands)} HOMER command rows to {outdir / 'homer_command_manifest.tsv'}")


if __name__ == "__main__":
    main()
