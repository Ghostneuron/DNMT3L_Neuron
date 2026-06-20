#!/usr/bin/env python3
"""Prepare HOMER motif-enrichment inputs for DNMT3L-sensitive methylation loci.

The script makes foreground BEDs for differential CpGs and DMRs, split by
methylation direction, and matched assay-background BEDs from all measured CpGs.
The resulting command manifest can be run with HOMER's findMotifsGenome.pl to
identify known TFBS motifs and de novo motifs.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import shlex
from pathlib import Path


HUMAN_DETECTION_COLUMNS = [
    "SAMPLE 1 Detection Pval",
    "SAMPLE 2 Detection Pval",
    "SAMPLE 3 Detection Pval",
    "SAMPLE 4 Detection Pval",
    "SAMPLE 5 Detection Pval",
    "SAMPLE 6 Detection Pval",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="results/homer_motif_enrichment")
    parser.add_argument("--window-size", type=int, default=200)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--motif-lengths", default="8,10,12")
    parser.add_argument("--human-genome", default="hg19")
    parser.add_argument("--mouse-genome", default="mm10")
    parser.add_argument("--human-dmrs", default="results/human_candidate_dmrs/GSE157078_candidate_dmrs.tsv")
    parser.add_argument("--human-dmcs", default="results/human_candidate_dmrs/GSE157078_candidate_dmcs.tsv")
    parser.add_argument("--human-signals", default="data/geo/GSE157078/supplementary/GSE157078_signals.txt.gz")
    parser.add_argument(
        "--human-manifest",
        default="data/geo/GPL21145/GPL21145_MethylationEPIC_15073387_v-1-0.csv.gz",
    )
    parser.add_argument(
        "--human-detection-p",
        type=float,
        default=0.01,
        help="Detection-p cutoff for human EPIC background probes.",
    )
    parser.add_argument(
        "--mouse-dmrs",
        default="results/mouse_candidate_dmrs_sex_consistent/GSE271949_candidate_dmrs.tsv",
    )
    parser.add_argument(
        "--mouse-dmcs",
        default="results/mouse_candidate_dmrs_sex_consistent/GSE271949_candidate_dmcs.tsv",
    )
    parser.add_argument(
        "--mouse-cpg-table",
        default="data/geo/GSE271949/supplementary/GSE271949_all_sample.CpG.txt.gz",
    )
    parser.add_argument(
        "--mouse-require-sex-consistent",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use only sex-consistent mouse candidate CpGs for the foreground DMC set.",
    )
    parser.add_argument(
        "--include-all-direction",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Also emit combined hyper+hypo foreground sets.",
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


def chrom_name(value: str) -> str:
    value = (value or "").strip()
    if not value:
        return ""
    return value if value.startswith("chr") else f"chr{value}"


def chrom_sort_key(chrom: str) -> tuple[int, int | str]:
    name = chrom.removeprefix("chr")
    if name.isdigit():
        return (0, int(name))
    special = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    return (0, special[name]) if name in special else (1, name)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_bed(path: Path, rows: list[tuple[str, int, int, str, int, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = sorted(rows, key=lambda row: (chrom_sort_key(row[0]), row[1], row[2], row[3]))
    with path.open("w", newline="") as handle:
        for chrom, start, end, name, score, strand in rows:
            if chrom and end > start:
                handle.write(f"{chrom}\t{max(0, start)}\t{end}\t{name}\t{score}\t{strand}\n")


def region_to_bed(row: dict[str, str], name: str) -> tuple[str, int, int, str, int, str]:
    chrom = row["chrom"]
    start = int(float(row["start"]))
    end = int(float(row["end"]))
    score = min(1000, int(abs(safe_float(row.get("mean_delta_beta", row.get("delta_beta", "0")))) * 1000))
    return chrom, max(0, start), max(start + 1, end), name, score, "."


def split_candidate_rows(rows: list[dict[str, str]], include_all: bool) -> dict[str, list[dict[str, str]]]:
    split = {
        "hyper": [row for row in rows if row.get("direction") == "hyper"],
        "hypo": [row for row in rows if row.get("direction") == "hypo"],
    }
    if include_all:
        split = {"all": rows, **split}
    return split


def read_human_epic_coords(manifest_path: Path) -> dict[str, tuple[str, int, int]]:
    coords: dict[str, tuple[str, int, int]] = {}
    with open_text(manifest_path) as handle:
        for line in handle:
            if line.startswith("IlmnID,"):
                header = next(csv.reader([line]))
                break
        else:
            raise RuntimeError(f"Could not find IlmnID header in {manifest_path}")

        reader = csv.DictReader(handle, fieldnames=header)
        for row in reader:
            probe_id = row.get("Name") or row.get("IlmnID")
            chrom = chrom_name(row.get("CHR", ""))
            pos = (row.get("MAPINFO") or "").strip()
            if not probe_id or not chrom or not pos.isdigit():
                continue
            start = int(pos) - 1
            coords[probe_id] = (chrom, start, start + 1)
    return coords


def human_detection_passed(signal_path: Path, detection_p: float) -> set[str]:
    passed: set[str] = set()
    with open_text(signal_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            probe_id = row.get("ID_REF", "")
            if not probe_id:
                continue
            pvals = [safe_float(row.get(column, "")) for column in HUMAN_DETECTION_COLUMNS]
            if all(math.isfinite(pval) and pval <= detection_p for pval in pvals):
                passed.add(probe_id)
    return passed


def build_human_background(manifest_path: Path, signal_path: Path, detection_p: float) -> list[tuple[str, int, int, str, int, str]]:
    coords = read_human_epic_coords(manifest_path)
    passed = human_detection_passed(signal_path, detection_p)
    rows = []
    for probe_id in sorted(passed):
        coord = coords.get(probe_id)
        if coord:
            rows.append((*coord, probe_id, 0, "."))
    return rows


def split_chrbase(value: str) -> tuple[str, int]:
    chrom, pos = value.rsplit(".", 1)
    return chrom, int(pos)


def build_mouse_background(cpg_table: Path) -> list[tuple[str, int, int, str, int, str]]:
    rows = []
    with open_text(cpg_table) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            locus = row.get("chrBase", "")
            if not locus or "." not in locus:
                continue
            chrom, pos = split_chrbase(locus)
            rows.append((chrom, pos - 1, pos, locus, 0, "."))
    return rows


def foreground_from_dmcs(
    rows: list[dict[str, str]],
    id_column: str,
    require_sex_consistent: bool = False,
) -> list[dict[str, str]]:
    if not require_sex_consistent:
        return rows
    return [row for row in rows if row.get("sex_consistent", "").lower() == "yes" and row.get(id_column)]


def write_foreground_sets(
    input_dir: Path,
    dataset_label: str,
    species: str,
    genome: str,
    level: str,
    rows: list[dict[str, str]],
    id_column: str,
    background_bed: Path,
    include_all: bool,
    commands: list[dict[str, str]],
    window_size: int,
    threads: int,
    motif_lengths: str,
) -> None:
    for direction, direction_rows in split_candidate_rows(rows, include_all).items():
        if not direction_rows:
            continue
        set_id = f"{dataset_label}_{level}_{direction}_{genome}"
        bed_path = input_dir / f"{set_id}.bed"
        bed_rows = [
            region_to_bed(row, row.get(id_column, f"{set_id}_{idx}"))
            for idx, row in enumerate(direction_rows, start=1)
        ]
        write_bed(bed_path, bed_rows)
        outdir = Path("results/homer_motif_enrichment") / "homer_out" / set_id
        command_parts = [
            "findMotifsGenome.pl",
            str(bed_path),
            genome,
            str(outdir),
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
                "dataset": dataset_label,
                "species": species,
                "genome": genome,
                "level": level,
                "direction": direction,
                "n_regions": str(len(bed_rows)),
                "window_size": str(window_size),
                "input_bed": str(bed_path),
                "background_bed": str(background_bed),
                "homer_outdir": str(outdir),
                "command": " ".join(shlex.quote(part) for part in command_parts),
            }
        )


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
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
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_shell_script(path: Path, commands: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        handle.write("#!/usr/bin/env bash\n")
        handle.write("set -euo pipefail\n\n")
        for row in commands:
            handle.write(f"echo 'Running {row['set_id']}'\n")
            handle.write(f"{row['command']}\n\n")


def main() -> None:
    args = parse_args()
    if args.window_size <= 0 or args.window_size % 2:
        raise SystemExit("--window-size must be a positive even integer")

    outdir = Path(args.outdir)
    input_dir = outdir / "inputs"
    background_dir = outdir / "backgrounds"
    commands: list[dict[str, str]] = []

    human_background = background_dir / f"GSE157078_EPIC_detected_background_{args.human_genome}.bed"
    write_bed(
        human_background,
        build_human_background(Path(args.human_manifest), Path(args.human_signals), args.human_detection_p),
    )

    mouse_background = background_dir / f"GSE271949_capture_all_cpg_background_{args.mouse_genome}.bed"
    write_bed(mouse_background, build_mouse_background(Path(args.mouse_cpg_table)))

    human_dmcs = read_tsv(Path(args.human_dmcs))
    human_dmrs = read_tsv(Path(args.human_dmrs))
    mouse_dmcs = foreground_from_dmcs(
        read_tsv(Path(args.mouse_dmcs)),
        id_column="locus_id",
        require_sex_consistent=args.mouse_require_sex_consistent,
    )
    mouse_dmrs = read_tsv(Path(args.mouse_dmrs))

    write_foreground_sets(
        input_dir=input_dir,
        dataset_label="GSE157078_human_DNMT3L_OE",
        species="Homo_sapiens",
        genome=args.human_genome,
        level="DMC",
        rows=human_dmcs,
        id_column="probe_id",
        background_bed=human_background,
        include_all=args.include_all_direction,
        commands=commands,
        window_size=args.window_size,
        threads=args.threads,
        motif_lengths=args.motif_lengths,
    )
    write_foreground_sets(
        input_dir=input_dir,
        dataset_label="GSE157078_human_DNMT3L_OE",
        species="Homo_sapiens",
        genome=args.human_genome,
        level="DMR",
        rows=human_dmrs,
        id_column="dmr_id",
        background_bed=human_background,
        include_all=args.include_all_direction,
        commands=commands,
        window_size=args.window_size,
        threads=args.threads,
        motif_lengths=args.motif_lengths,
    )
    write_foreground_sets(
        input_dir=input_dir,
        dataset_label="GSE271949_mouse_Dnmt3l_OE_sex_consistent",
        species="Mus_musculus",
        genome=args.mouse_genome,
        level="DMC",
        rows=mouse_dmcs,
        id_column="locus_id",
        background_bed=mouse_background,
        include_all=args.include_all_direction,
        commands=commands,
        window_size=args.window_size,
        threads=args.threads,
        motif_lengths=args.motif_lengths,
    )
    write_foreground_sets(
        input_dir=input_dir,
        dataset_label="GSE271949_mouse_Dnmt3l_OE_sex_consistent",
        species="Mus_musculus",
        genome=args.mouse_genome,
        level="DMR",
        rows=mouse_dmrs,
        id_column="dmr_id",
        background_bed=mouse_background,
        include_all=args.include_all_direction,
        commands=commands,
        window_size=args.window_size,
        threads=args.threads,
        motif_lengths=args.motif_lengths,
    )

    write_manifest(outdir / "homer_command_manifest.tsv", commands)
    write_shell_script(outdir / "run_homer_motif_enrichment.sh", commands)

    print(f"Wrote {len(commands)} HOMER foreground command rows")
    print(f"Wrote inputs to {input_dir}")
    print(f"Wrote backgrounds to {background_dir}")
    print(f"Wrote command manifest to {outdir / 'homer_command_manifest.tsv'}")


if __name__ == "__main__":
    main()
