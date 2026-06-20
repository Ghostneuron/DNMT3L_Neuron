#!/usr/bin/env python3
"""Download missing GSE305146 Matrix Market files listed in local GEO metadata."""

from __future__ import annotations

import argparse
import csv
import shutil
import urllib.parse
import urllib.request
from pathlib import Path


USER_AGENT = "DNMT3L-Neuron-GSE305146-matrix-download/0.1"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", default="data/geo/GSE305146/metadata/supplementary_files.tsv")
    parser.add_argument("--outdir", default="data/geo/GSE305146/supplementary")
    parser.add_argument("--log", default="data/geo/GSE305146/metadata/matrix_download_log.tsv")
    parser.add_argument("--timeout", type=int, default=240)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--limit", type=int, help="Optional debug limit.")
    return parser.parse_args()


def filename_from_url(url: str) -> str:
    return Path(urllib.parse.urlparse(url).path).name


def read_matrix_urls(path: Path) -> list[str]:
    urls: list[str] = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            url = row.get("url", "")
            if url.endswith("_matrix.mtx.gz"):
                urls.append(url)
    return urls


def download(url: str, dest: Path, timeout: int, force: bool) -> dict[str, object]:
    filename = filename_from_url(url)
    if dest.exists() and not force and dest.stat().st_size > 0:
        return {
            "filename": filename,
            "url": url,
            "status": "already_exists",
            "bytes": dest.stat().st_size,
            "expected_bytes": "",
            "path": str(dest),
        }

    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".tmp")
    if tmp.exists():
        tmp.unlink()

    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as response:
            expected = response.headers.get("Content-Length", "")
            with tmp.open("wb") as handle:
                shutil.copyfileobj(response, handle, length=1024 * 1024)
        tmp.replace(dest)
        status = "downloaded"
        bytes_written = dest.stat().st_size
        if expected and bytes_written != int(expected):
            status = "downloaded_size_mismatch"
        return {
            "filename": filename,
            "url": url,
            "status": status,
            "bytes": bytes_written,
            "expected_bytes": expected,
            "path": str(dest),
        }
    except Exception as exc:
        if tmp.exists():
            tmp.unlink()
        return {
            "filename": filename,
            "url": url,
            "status": f"error:{type(exc).__name__}:{exc}",
            "bytes": "",
            "expected_bytes": "",
            "path": str(dest),
        }


def write_log(path: Path, rows: list[dict[str, object]]) -> None:
    fields = ["filename", "url", "status", "bytes", "expected_bytes", "path"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    urls = read_matrix_urls(Path(args.metadata))
    if args.limit:
        urls = urls[: args.limit]

    outdir = Path(args.outdir)
    rows: list[dict[str, object]] = []
    for index, url in enumerate(urls, start=1):
        filename = filename_from_url(url)
        dest = outdir / filename
        print(f"[{index}/{len(urls)}] {filename}", flush=True)
        row = download(url, dest, args.timeout, args.force)
        rows.append(row)
        print(f"  {row['status']} {row['bytes']}", flush=True)
        write_log(Path(args.log), rows)

    n_downloaded = sum(1 for row in rows if str(row["status"]).startswith("downloaded"))
    n_existing = sum(1 for row in rows if row["status"] == "already_exists")
    n_errors = sum(1 for row in rows if str(row["status"]).startswith("error:"))
    print(f"Downloaded={n_downloaded}; already_exists={n_existing}; errors={n_errors}", flush=True)
    print(args.log, flush=True)
    return 1 if n_errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
