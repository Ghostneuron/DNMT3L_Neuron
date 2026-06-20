#!/usr/bin/env python3
"""Prepare a tractable first-pass HOMER motif set for mouse Dnmt3l-OE DMRs."""

from __future__ import annotations

import argparse
import csv
import random
import shlex
from pathlib import Path

from prepare_homer_motif_inputs import chrom_sort_key, open_text, safe_float, write_bed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="results/mouse_homer_motif_enrichment")
    parser.add_argument(
        "--dmrs",
        default="results/mouse_candidate_dmrs_sex_consistent/GSE271949_candidate_dmrs.tsv",
    )
    parser.add_argument(
        "--background",
        default="results/homer_motif_enrichment/backgrounds/GSE271949_capture_all_cpg_background_mm10.bed",
    )
    parser.add_argument("--genome", default="mm10")
    parser.add_argument("--background-size", type=int, default=400_000)
    parser.add_argument("--seed", type=int, default=271949)
    parser.add_argument("--top-hyper", type=int, default=5000)
    parser.add_argument("--window-size", type=int, default=200)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--motif-lengths", default="8,10")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_bed(path: Path) -> list[tuple[str, int, int, str, int, str]]:
    rows = []
    with path.open(newline="") as handle:
        for line in handle:
            if not line.strip():
                continue
            chrom, start, end, name, score, strand = line.rstrip("\n").split("\t")[:6]
            rows.append((chrom, int(start), int(end), name, int(score), strand))
    return rows


def sample_background(rows: list[tuple[str, int, int, str, int, str]], n_rows: int, seed: int) -> list[tuple[str, int, int, str, int, str]]:
    if n_rows <= 0 or n_rows >= len(rows):
        sampled = rows
    else:
        rng = random.Random(seed)
        sampled = rng.sample(rows, n_rows)
    return sorted(sampled, key=lambda row: (chrom_sort_key(row[0]), row[1], row[2], row[3]))


def dmr_score(row: dict[str, str]) -> int:
    return min(1000, int(abs(safe_float(row.get("mean_delta_beta", ""))) * 1000))


def dmr_effect_key(row: dict[str, str]) -> tuple[float, float, int]:
    return (
        abs(safe_float(row.get("max_abs_delta_beta", ""))),
        abs(safe_float(row.get("mean_delta_beta", ""))),
        int(float(row.get("n_cpgs", "0") or 0)),
    )


def dmr_to_bed(row: dict[str, str]) -> tuple[str, int, int, str, int, str]:
    start = int(float(row["start"]))
    end = int(float(row["end"]))
    return row["chrom"], start, max(start + 1, end), row["dmr_id"], dmr_score(row), "."


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
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def add_command(
    commands: list[dict[str, str]],
    outdir: Path,
    genome: str,
    set_label: str,
    direction: str,
    foreground_rows: list[dict[str, str]],
    background_bed: Path,
    window_size: int,
    threads: int,
    motif_lengths: str,
) -> None:
    set_id = f"GSE271949_mouse_Dnmt3l_OE_{set_label}_{genome}"
    bed_path = outdir / "inputs" / f"{set_id}.bed"
    write_bed(bed_path, [dmr_to_bed(row) for row in foreground_rows])
    homer_outdir = outdir / "homer_out" / set_id
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
        str(background_bed),
    ]
    commands.append(
        {
            "set_id": set_id,
            "dataset": f"GSE271949_mouse_Dnmt3l_OE_{set_label}",
            "species": "Mus_musculus",
            "genome": genome,
            "level": "DMR",
            "direction": direction,
            "n_regions": str(len(foreground_rows)),
            "window_size": str(window_size),
            "input_bed": str(bed_path),
            "background_bed": str(background_bed),
            "homer_outdir": str(homer_outdir),
            "command": " ".join(shlex.quote(part) for part in command_parts),
        }
    )


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    dmrs = read_tsv(Path(args.dmrs))
    hyper = [row for row in dmrs if row.get("direction") == "hyper"]
    hypo = [row for row in dmrs if row.get("direction") == "hypo"]
    hyper_top = sorted(hyper, key=dmr_effect_key, reverse=True)[: args.top_hyper]

    background_rows = read_bed(Path(args.background))
    sampled_background = sample_background(background_rows, args.background_size, args.seed)
    background_bed = outdir / "backgrounds" / f"GSE271949_capture_sample{len(sampled_background)}_seed{args.seed}_{args.genome}.bed"
    write_bed(background_bed, sampled_background)

    commands: list[dict[str, str]] = []
    add_command(
        commands,
        outdir,
        args.genome,
        f"DMR_hyper_top{len(hyper_top)}",
        "hyper_top",
        hyper_top,
        background_bed,
        args.window_size,
        args.threads,
        args.motif_lengths,
    )
    add_command(
        commands,
        outdir,
        args.genome,
        "DMR_hypo_all",
        "hypo",
        hypo,
        background_bed,
        args.window_size,
        args.threads,
        args.motif_lengths,
    )
    add_command(
        commands,
        outdir,
        args.genome,
        "DMR_hyper_all",
        "hyper",
        hyper,
        background_bed,
        args.window_size,
        args.threads,
        args.motif_lengths,
    )

    write_manifest(outdir / "homer_command_manifest.tsv", commands)
    print(f"Wrote sampled background with {len(sampled_background)} regions")
    print(f"Wrote {len(commands)} HOMER command rows to {outdir / 'homer_command_manifest.tsv'}")
    print(f"Mouse DMR counts: hyper={len(hyper)}, hypo={len(hypo)}, hyper_top={len(hyper_top)}")


if __name__ == "__main__":
    main()
