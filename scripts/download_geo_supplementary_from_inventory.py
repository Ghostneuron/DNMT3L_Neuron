#!/usr/bin/env python3
"""Download selected GEO supplementary files from an inventory TSV."""

from __future__ import annotations

import argparse
import csv
import shutil
import urllib.parse
import urllib.request
from pathlib import Path


USER_AGENT = "DNMT3L-Neuron-GEO-inventory-download/0.1"
DEFAULT_STRATEGIES = ["download_now_processed", "download_now_matrix_index"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", default="results/ds_ncbi_geo/high_priority_supplementary_inventory.tsv")
    parser.add_argument("--outroot", default="data/geo")
    parser.add_argument("--strategies", nargs="*", default=DEFAULT_STRATEGIES)
    parser.add_argument("--accessions", nargs="*", help="Optional subset of accessions.")
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--limit", type=int, help="Optional maximum number of files to download.")
    return parser.parse_args()


def filename_from_url(url: str) -> str:
    return Path(urllib.parse.urlparse(url).path).name


def download(url: str, dest: Path, timeout: int, force: bool) -> tuple[str, str]:
    if not url.startswith("http"):
        return "skipped_non_http", ""
    if dest.exists() and not force:
        return "already_exists", str(dest.stat().st_size)
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as response, tmp.open("wb") as handle:
            shutil.copyfileobj(response, handle)
        tmp.replace(dest)
        return "downloaded", str(dest.stat().st_size)
    except Exception as exc:
        if tmp.exists():
            tmp.unlink()
        return f"error:{type(exc).__name__}:{exc}", ""


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    wanted_strategies = set(args.strategies)
    wanted_accessions = {acc.upper() for acc in args.accessions or []}
    outroot = Path(args.outroot)

    selected: list[dict[str, str]] = []
    with Path(args.inventory).open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            accession = row.get("accession", "").upper()
            if wanted_accessions and accession not in wanted_accessions:
                continue
            if row.get("download_strategy") not in wanted_strategies:
                continue
            url = row.get("url", "")
            filename = filename_from_url(url)
            if not filename:
                continue
            selected.append({**row, "filename": filename})

    if args.limit:
        selected = selected[: args.limit]

    log_rows: list[dict[str, str]] = []
    for index, row in enumerate(selected, start=1):
        accession = row["accession"].upper()
        filename = row["filename"]
        dest = outroot / accession / "supplementary" / filename
        print(f"[{index}/{len(selected)}] {accession} {filename}", flush=True)
        status, bytes_written = download(row["url"], dest, args.timeout, args.force)
        log_rows.append(
            {
                "accession": accession,
                "filename": filename,
                "download_strategy": row.get("download_strategy", ""),
                "url": row.get("url", ""),
                "status": status,
                "path": str(dest) if status not in {"skipped_non_http"} else "",
                "bytes": bytes_written,
            }
        )

    log_path = outroot / "ds_ncbi_geo_inventory_download_log.tsv"
    fields = ["accession", "filename", "download_strategy", "url", "status", "path", "bytes"]
    write_tsv(log_path, log_rows, fields)
    print(f"Wrote download log for {len(log_rows)} files to {log_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
