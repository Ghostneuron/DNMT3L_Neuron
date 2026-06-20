#!/usr/bin/env python3
"""Summarize HOMER known and de novo motif results into compact TSV tables."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", default="results/homer_motif_enrichment/homer_command_manifest.tsv")
    parser.add_argument("--outdir", default="results/homer_motif_enrichment")
    parser.add_argument("--top-n", type=int, default=25)
    return parser.parse_args()


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_known_results(path: Path, set_row: dict[str, str], top_n: int) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = []
        for rank, row in enumerate(reader, start=1):
            if rank > top_n:
                break
            rows.append(
                {
                    "set_id": set_row["set_id"],
                    "dataset": set_row["dataset"],
                    "genome": set_row["genome"],
                    "level": set_row["level"],
                    "direction": set_row["direction"],
                    "rank": str(rank),
                    "motif_name": row.get("Motif Name", ""),
                    "consensus": row.get("Consensus", ""),
                    "p_value": row.get("P-value", ""),
                    "log_p_value": row.get("Log P-value", ""),
                    "q_value": row.get("q-value (Benjamini)", ""),
                    "target_percent": row.get("% of Target Sequences with Motif", ""),
                    "background_percent": row.get("% of Background Sequences with Motif", ""),
                }
            )
    return rows


def parse_denovo_motif(path: Path) -> dict[str, str]:
    motif = {
        "motif_file": str(path),
        "motif_name": "",
        "consensus": "",
        "log_p_value": "",
        "target_percent": "",
        "background_percent": "",
        "best_match": "",
    }
    if not path.exists():
        return motif
    first_line = path.read_text(errors="replace").splitlines()[0] if path.stat().st_size else ""
    if not first_line.startswith(">"):
        return motif
    parts = first_line[1:].split("\t")
    motif["motif_name"] = parts[0] if parts else ""
    motif["consensus"] = parts[0].split(",")[0] if parts else ""
    for part in parts:
        if part.startswith("log P-value:"):
            motif["log_p_value"] = part.replace("log P-value:", "").strip()
        elif "% of Targets" in part:
            motif["target_percent"] = part.strip()
        elif "% of Background" in part:
            motif["background_percent"] = part.strip()
        elif "BestGuess:" in part:
            motif["best_match"] = re.sub(r".*BestGuess:", "", part).strip()
    return motif


def parse_denovo_results(homer_outdir: Path, set_row: dict[str, str], top_n: int) -> list[dict[str, str]]:
    motif_dir = homer_outdir / "homerResults"
    if not motif_dir.exists():
        return []
    rows = []
    for rank, path in enumerate(sorted(motif_dir.glob("motif*.motif")), start=1):
        if rank > top_n:
            break
        motif = parse_denovo_motif(path)
        rows.append(
            {
                "set_id": set_row["set_id"],
                "dataset": set_row["dataset"],
                "genome": set_row["genome"],
                "level": set_row["level"],
                "direction": set_row["direction"],
                "rank": str(rank),
                **motif,
            }
        )
    return rows


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    known_rows: list[dict[str, str]] = []
    denovo_rows: list[dict[str, str]] = []
    for row in read_manifest(Path(args.manifest)):
        homer_outdir = Path(row["homer_outdir"])
        known_rows.extend(parse_known_results(homer_outdir / "knownResults.txt", row, args.top_n))
        denovo_rows.extend(parse_denovo_results(homer_outdir, row, args.top_n))

    outdir = Path(args.outdir)
    write_tsv(
        outdir / "homer_known_motif_summary.tsv",
        known_rows,
        [
            "set_id",
            "dataset",
            "genome",
            "level",
            "direction",
            "rank",
            "motif_name",
            "consensus",
            "p_value",
            "log_p_value",
            "q_value",
            "target_percent",
            "background_percent",
        ],
    )
    write_tsv(
        outdir / "homer_denovo_motif_summary.tsv",
        denovo_rows,
        [
            "set_id",
            "dataset",
            "genome",
            "level",
            "direction",
            "rank",
            "motif_file",
            "motif_name",
            "consensus",
            "log_p_value",
            "target_percent",
            "background_percent",
            "best_match",
        ],
    )
    print(f"Wrote {len(known_rows)} known motif rows")
    print(f"Wrote {len(denovo_rows)} de novo motif rows")


if __name__ == "__main__":
    main()
