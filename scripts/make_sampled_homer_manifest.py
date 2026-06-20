#!/usr/bin/env python3
"""Create a HOMER manifest with deterministic sampled backgrounds.

This keeps large array-derived motif tests tractable while retaining a
context-matched detected-probe background for each contrast.
"""

from __future__ import annotations

import argparse
import csv
import random
import shlex
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--background-size", type=int, default=100000)
    parser.add_argument("--seed", type=int, default=74486)
    parser.add_argument(
        "--de-novo-count",
        type=int,
        default=8,
        help="Pass HOMER -S with this many de novo motifs per run.",
    )
    parser.add_argument(
        "--direction",
        action="append",
        default=[],
        help="Direction to keep. Can be repeated. Defaults to hyper and hypo.",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sample_background(path: Path, out_path: Path, size: int, seed: int) -> tuple[int, int]:
    rows = path.read_text(errors="replace").splitlines()
    rng = random.Random(seed)
    sampled = rows if len(rows) <= size else rng.sample(rows, size)
    sampled.sort(key=lambda line: line.split("\t")[:3])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(sampled) + ("\n" if sampled else ""))
    return len(rows), len(sampled)


def rebuild_command(row: dict[str, str], de_novo_count: int) -> str:
    parts = [
        "findMotifsGenome.pl",
        row["input_bed"],
        row["genome"],
        row["homer_outdir"],
        "-size",
        row["window_size"],
        "-mask",
        "-p",
        "8",
        "-len",
        "8,10",
        "-S",
        str(de_novo_count),
        "-bg",
        row["background_bed"],
    ]
    return " ".join(shlex.quote(part) for part in parts)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    directions = set(args.direction or ["hyper", "hypo"])
    rows = [row for row in read_tsv(Path(args.manifest)) if row["direction"] in directions]

    sampled_paths: dict[str, Path] = {}
    background_rows: list[dict[str, str]] = []
    for row in rows:
        original = Path(row["background_bed"])
        if str(original) not in sampled_paths:
            sampled = outdir / "backgrounds" / f"{original.stem}_sample{args.background_size}_seed{args.seed}.bed"
            total, kept = sample_background(original, sampled, args.background_size, args.seed + len(sampled_paths))
            sampled_paths[str(original)] = sampled
            background_rows.append(
                {
                    "original_background": str(original),
                    "sampled_background": str(sampled),
                    "total_rows": str(total),
                    "sampled_rows": str(kept),
                }
            )

    new_rows: list[dict[str, str]] = []
    for row in rows:
        new_row = dict(row)
        new_row["background_bed"] = str(sampled_paths[row["background_bed"]])
        new_row["homer_outdir"] = str(outdir / "homer_out" / row["set_id"])
        new_row["command"] = rebuild_command(new_row, args.de_novo_count)
        new_rows.append(new_row)

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
    write_tsv(outdir / "homer_command_manifest.tsv", new_rows, fields)
    write_tsv(
        outdir / "sampled_background_summary.tsv",
        background_rows,
        ["original_background", "sampled_background", "total_rows", "sampled_rows"],
    )
    print(f"Wrote {len(new_rows)} HOMER command rows to {outdir / 'homer_command_manifest.tsv'}")
    print(f"Wrote {len(background_rows)} sampled backgrounds")


if __name__ == "__main__":
    main()
