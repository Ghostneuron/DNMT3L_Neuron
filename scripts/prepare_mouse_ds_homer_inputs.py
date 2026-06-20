#!/usr/bin/env python3
"""Prepare first-pass HOMER inputs for mouse DS-model WGBS tracks."""

from __future__ import annotations

import argparse
import csv
import gzip
from dataclasses import dataclass
from pathlib import Path


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


@dataclass(frozen=True)
class CpG:
    chrom: str
    start: int
    end: int
    delta: float
    model: str

    @property
    def key(self) -> tuple[str, int, int]:
        return (self.chrom, self.start, self.end)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="results/mouse_ds_homer_motif_enrichment")
    parser.add_argument("--genome", default="mm10")
    parser.add_argument("--top-n", type=int, default=5000)
    parser.add_argument("--size", default="200")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--motif-lengths", default="8,10")
    return parser.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="", errors="replace")
    return path.open(newline="", errors="replace")


def read_track(model: str, path: Path) -> list[CpG]:
    rows = []
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, value = line.rstrip("\n").split("\t")[:4]
            rows.append(CpG(chrom=chrom, start=int(start), end=int(end), delta=float(value), model=model))
    return rows


def bed_row(cpg: CpG, name: str) -> tuple[str, int, int, str, int, str]:
    return (cpg.chrom, cpg.start, cpg.end, name, 0, "+")


def write_bed(path: Path, rows: list[tuple[str, int, int, str, int, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        for row in rows:
            handle.write("\t".join(str(value) for value in row) + "\n")


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "set_id",
        "dataset",
        "genome",
        "level",
        "direction",
        "n_regions",
        "input_bed",
        "background_bed",
        "homer_outdir",
        "command",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def make_command(input_bed: Path, genome: str, outdir: Path, background_bed: Path, args: argparse.Namespace) -> str:
    return " ".join(
        [
            "findMotifsGenome.pl",
            str(input_bed),
            genome,
            str(outdir),
            "-size",
            args.size,
            "-mask",
            "-p",
            str(args.threads),
            "-len",
            args.motif_lengths,
            "-bg",
            str(background_bed),
        ]
    )


def add_manifest_row(
    manifest_rows: list[dict[str, str]],
    *,
    model: str,
    label: str,
    direction: str,
    selected: list[CpG],
    input_bed: Path,
    background_bed: Path,
    genome: str,
    outdir: Path,
    args: argparse.Namespace,
) -> None:
    set_id = f"GSE74505_{model}_DS_model_{label}_{genome}"
    homer_outdir = outdir / "homer_out" / set_id
    manifest_rows.append(
        {
            "set_id": set_id,
            "dataset": f"GSE74505_{model}_DS_model",
            "genome": genome,
            "level": "DMC",
            "direction": direction,
            "n_regions": str(len(selected)),
            "input_bed": str(input_bed),
            "background_bed": str(background_bed),
            "homer_outdir": str(homer_outdir),
            "command": make_command(input_bed, genome, homer_outdir, background_bed, args),
        }
    )


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    input_dir = outdir / "inputs"
    background_dir = outdir / "backgrounds"
    manifest_rows: list[dict[str, str]] = []

    by_model = {track["model"]: read_track(track["model"], Path(track["path"])) for track in TRACKS}

    background_by_key: dict[tuple[str, int, int], CpG] = {}
    for rows in by_model.values():
        for cpg in rows:
            background_by_key.setdefault(cpg.key, cpg)
    background = sorted(background_by_key.values(), key=lambda cpg: (cpg.chrom, cpg.start, cpg.end))
    background_bed = background_dir / f"GSE74505_Dp10_Dp16_union_differential_CpG_background_{args.genome}.bed"
    write_bed(
        background_bed,
        [bed_row(cpg, f"bg_{index}") for index, cpg in enumerate(background, start=1)],
    )

    for model, rows in by_model.items():
        hyper = sorted((cpg for cpg in rows if cpg.delta > 0), key=lambda cpg: abs(cpg.delta), reverse=True)
        hypo = sorted((cpg for cpg in rows if cpg.delta < 0), key=lambda cpg: abs(cpg.delta), reverse=True)
        for direction, selected_all in [("hyper", hyper), ("hypo", hypo)]:
            selected = selected_all[: args.top_n]
            label = f"DMC_{direction}_top{len(selected)}"
            input_bed = input_dir / f"GSE74505_{model}_DS_model_{label}_{args.genome}.bed"
            write_bed(
                input_bed,
                [bed_row(cpg, f"{model}_{direction}_{index}") for index, cpg in enumerate(selected, start=1)],
            )
            add_manifest_row(
                manifest_rows,
                model=model,
                label=label,
                direction=f"{direction}_top",
                selected=selected,
                input_bed=input_bed,
                background_bed=background_bed,
                genome=args.genome,
                outdir=outdir,
                args=args,
            )

    write_manifest(outdir / "homer_command_manifest.tsv", manifest_rows)
    print(f"Wrote common background with {len(background)} CpGs")
    print(f"Wrote {len(manifest_rows)} HOMER command rows to {outdir / 'homer_command_manifest.tsv'}")
    for model, rows in by_model.items():
        hyper_n = sum(1 for cpg in rows if cpg.delta > 0)
        hypo_n = sum(1 for cpg in rows if cpg.delta < 0)
        print(f"{model}: total={len(rows)} hyper={hyper_n} hypo={hypo_n}")


if __name__ == "__main__":
    main()
