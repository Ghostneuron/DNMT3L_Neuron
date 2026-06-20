#!/usr/bin/env python3
"""Run HOMER motif-enrichment commands prepared by prepare_homer_motif_inputs.py."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", default="results/homer_motif_enrichment/homer_command_manifest.tsv")
    parser.add_argument(
        "--set-id",
        action="append",
        default=[],
        help="Specific HOMER set_id to run. Can be passed multiple times. Defaults to all.",
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--continue-on-error", action="store_true")
    return parser.parse_args()


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_status(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    with path.open(newline="") as handle:
        return {
            row["set_id"]: row
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("set_id")
        }


def write_status(path: Path, rows: list[dict[str, str]], status_by_set: dict[str, dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["set_id", "status", "returncode", "homer_outdir", "command"]
    ordered_rows = [status_by_set[row["set_id"]] for row in rows if row["set_id"] in status_by_set]
    manifest_sets = {row["set_id"] for row in rows}
    ordered_rows.extend(row for set_id, row in sorted(status_by_set.items()) if set_id not in manifest_sets)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(ordered_rows)


def main() -> None:
    args = parse_args()
    manifest = Path(args.manifest)
    rows = read_manifest(manifest)
    if args.set_id:
        selected = set(args.set_id)
        rows = [row for row in rows if row["set_id"] in selected]
    if not rows:
        raise SystemExit("No HOMER command rows selected")

    status_path = manifest.parent / "homer_run_status.tsv"
    status_by_set = read_status(status_path)

    if shutil.which("findMotifsGenome.pl") is None:
        for row in rows:
            status_by_set[row["set_id"]] = {
                "set_id": row["set_id"],
                "status": "not_run_homer_not_on_path",
                "returncode": "",
                "homer_outdir": row["homer_outdir"],
                "command": row["command"],
            }
        write_status(status_path, read_manifest(manifest), status_by_set)
        print("findMotifsGenome.pl is not on PATH; wrote not-run status.", file=sys.stderr)
        raise SystemExit(127)

    for row in rows:
        print(f"Running {row['set_id']}")
        if args.dry_run:
            status_by_set[row["set_id"]] = {
                "set_id": row["set_id"],
                "status": "dry_run",
                "returncode": "0",
                "homer_outdir": row["homer_outdir"],
                "command": row["command"],
            }
            continue
        completed = subprocess.run(row["command"], shell=True, check=False)
        status = "completed" if completed.returncode == 0 else "failed"
        status_by_set[row["set_id"]] = {
            "set_id": row["set_id"],
            "status": status,
            "returncode": str(completed.returncode),
            "homer_outdir": row["homer_outdir"],
            "command": row["command"],
        }
        write_status(status_path, read_manifest(manifest), status_by_set)
        if completed.returncode != 0 and not args.continue_on_error:
            raise SystemExit(completed.returncode)

    write_status(status_path, read_manifest(manifest), status_by_set)
    print(f"Wrote run status to {status_path}")


if __name__ == "__main__":
    main()
