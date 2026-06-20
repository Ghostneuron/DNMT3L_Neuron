#!/usr/bin/env python3
"""Inventory GEO supplementary files and assign a first-pass download strategy."""

from __future__ import annotations

import argparse
import csv
import re
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path


USER_AGENT = "DNMT3L-Neuron-GEO-supplementary-inventory/0.1"

RAW_OR_LARGE_RE = re.compile(
    r"(raw\.tar|\.tar$|\.bw$|\.bigwig$|\.hic$|\.bam$|\.cram$|\.fastq|\.fq|"
    r"\.sra$|\.cel(\.gz)?$|\.idat(\.gz)?$|\.rda(\.gz)?$|\.rds$|\.h5$|"
    r"\.h5ad$|\.loom$|\.vcf(\.gz)?$|cpg_report|bismark|\.cov(\.gz)?)",
    flags=re.IGNORECASE,
)
PROCESSED_RE = re.compile(
    r"(counts?|rawcounts?|readspergene|fpkm|tpm|rpkm|deg|deseq|edgeR|limma|"
    r"expression|normalized|non.normalized|matrix|beta|methylation|sample|"
    r"metadata|annotation|result|\.txt(\.gz)?$|\.tsv(\.gz)?$|\.csv(\.gz)?$|"
    r"\.tab(\.gz)?$|\.xls(x)?$)",
    flags=re.IGNORECASE,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accession-table", help="TSV containing an accession column.")
    parser.add_argument("--accession-column", default="accession", help="Column name for accessions.")
    parser.add_argument("--accessions", nargs="*", help="Accessions to inventory.")
    parser.add_argument("--metadata-root", default="data/geo", help="Root containing data/geo/<GSE>/metadata.")
    parser.add_argument("--out", default="results/ds_ncbi_geo/high_priority_supplementary_inventory.tsv")
    parser.add_argument("--max-download-mb", type=float, default=250.0)
    parser.add_argument("--timeout", type=int, default=30)
    parser.add_argument(
        "--skip-size-probe",
        action="store_true",
        help="Do not issue remote HEAD requests; classify by filename only.",
    )
    return parser.parse_args()


def read_accessions(args: argparse.Namespace) -> list[str]:
    accessions: list[str] = []
    if args.accession_table:
        with Path(args.accession_table).open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                accession = row.get(args.accession_column, "").strip().upper()
                if accession:
                    accessions.append(accession)
    for accession in args.accessions or []:
        accessions.append(accession.strip().upper())
    return sorted(set(accessions))


def remote_size(url: str, timeout: int) -> tuple[int | None, str]:
    if not url.startswith("http"):
        return None, "non_http"
    try:
        req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(req, timeout=timeout) as response:
            length = response.headers.get("Content-Length")
            return (int(length) if length else None), "ok"
    except Exception as exc:
        return None, f"head_error:{type(exc).__name__}"


def filename_from_url(url: str) -> str:
    if not url.startswith("http"):
        return url
    return Path(urllib.parse.urlparse(url).path).name


def strategy(filename: str, size: int | None, max_download_mb: float) -> str:
    lower = filename.lower()
    if filename in {"", "NONE"}:
        return "no_file"
    if "matrix.mtx" in lower:
        if size is not None and size / 1024 / 1024 <= max_download_mb:
            return "stream_or_download_matrix"
        return "stream_matrix_skip_large"
    if re.search(r"(features|barcodes).*\.tsv", lower):
        return "download_now_matrix_index"
    if RAW_OR_LARGE_RE.search(lower):
        return "skip_raw_or_large_object"
    if size is not None and size / 1024 / 1024 > max_download_mb:
        return "skip_over_size_cap"
    if PROCESSED_RE.search(lower):
        return "download_now_processed"
    return "inspect"


def should_probe_size(filename: str) -> bool:
    lower = filename.lower()
    if filename in {"", "NONE"}:
        return False
    if RAW_OR_LARGE_RE.search(lower) and "matrix.mtx" not in lower:
        return False
    return True


def read_series_title(metadata_root: Path, accession: str) -> str:
    path = metadata_root / accession / "metadata" / "series_metadata.tsv"
    if not path.exists():
        return ""
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return rows[0].get("title", "") if rows else ""


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    metadata_root = Path(args.metadata_root)
    rows: list[dict[str, str]] = []

    for accession in read_accessions(args):
        supp_path = metadata_root / accession / "metadata" / "supplementary_files.tsv"
        title = read_series_title(metadata_root, accession)
        if not supp_path.exists():
            rows.append(
                {
                    "accession": accession,
                    "series_title": title,
                    "context": "",
                    "sample_or_series": "",
                    "filename": "",
                    "remote_size_bytes": "",
                    "remote_size_mb": "",
                    "head_status": "missing_supplementary_table",
                    "download_strategy": "metadata_missing",
                    "url": "",
                }
            )
            continue
        with supp_path.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                url = row.get("url", "")
                filename = filename_from_url(url)
                if args.skip_size_probe:
                    size, head_status = None, "size_probe_skipped_by_user"
                elif should_probe_size(filename):
                    size, head_status = remote_size(url, args.timeout)
                else:
                    size, head_status = None, "size_probe_skipped_by_pattern"
                rows.append(
                    {
                        "accession": accession,
                        "series_title": title,
                        "context": row.get("context", ""),
                        "sample_or_series": row.get("accession", ""),
                        "filename": filename,
                        "remote_size_bytes": "" if size is None else str(size),
                        "remote_size_mb": "" if size is None else f"{size / 1024 / 1024:.2f}",
                        "head_status": head_status,
                        "download_strategy": strategy(filename, size, args.max_download_mb),
                        "url": url,
                    }
                )

    fields = [
        "accession",
        "series_title",
        "context",
        "sample_or_series",
        "filename",
        "remote_size_bytes",
        "remote_size_mb",
        "head_status",
        "download_strategy",
        "url",
    ]
    write_tsv(Path(args.out), rows, fields)

    print(f"Wrote {len(rows)} supplementary rows to {args.out}")
    for strategy_name, count in Counter(row["download_strategy"] for row in rows).most_common():
        print(f"{strategy_name}\t{count}")
    by_accession: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if row["download_strategy"] in {
            "download_now_processed",
            "download_now_matrix_index",
            "stream_or_download_matrix",
        }:
            by_accession[row["accession"]].append(row)
    print("first-pass-downloadable-accessions")
    for accession in sorted(by_accession):
        total = sum(int(row["remote_size_bytes"] or 0) for row in by_accession[accession])
        print(f"{accession}\t{len(by_accession[accession])}\t{total / 1024 / 1024:.2f} MB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
