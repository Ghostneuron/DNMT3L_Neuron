#!/usr/bin/env python3
"""Download the public Zenodo code package for Science aea1549 postnatal PFC."""

from __future__ import annotations

import argparse
import csv
import json
import re
import urllib.request
import zipfile
from pathlib import Path


RECORD_URL = "https://zenodo.org/api/records/17527842"
ACCESSION = "ZENODO_17527842"
DEFAULT_OUTDIR = Path("data/zenodo") / ACCESSION
USER_AGENT = "DNMT3L-Neuron-Zenodo17527842/0.1"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument(
        "--include-large-motif",
        action="store_true",
        help="Also download Motif_analysis.zip (~710 MB).",
    )
    parser.add_argument(
        "--no-extract",
        action="store_true",
        help="Download archives without extracting the small code zips.",
    )
    return parser.parse_args()


def request_bytes(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req) as response:
        return response.read()


def safe_local_name(file_key: str) -> str:
    return file_key.replace(" ", "_")


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def maybe_extract(zip_path: Path, extract_dir: Path) -> None:
    if zip_path.suffix.lower() != ".zip":
        return
    with zipfile.ZipFile(zip_path) as archive:
        archive.extractall(extract_dir)


def audit_sample_ids(preprocessing_code: Path, outdir: Path) -> None:
    if not preprocessing_code.exists():
        return
    code = preprocessing_code.read_text()
    pattern = re.compile(r"[DH]23(?:0215|0223|0307|0313|0315|0330|0406|0629)E\d|H240110E[2-6]")
    seen: list[str] = []
    for match in pattern.finditer(code):
        sample_id = match.group(0)
        if sample_id not in seen:
            seen.append(sample_id)
    rows = []
    for sample_id in seen:
        rows.append(
            {
                "sample_id_in_code": sample_id,
                "condition_in_code": "Ts21" if sample_id.startswith("D") else "Control",
                "donor_prefix_in_code": sample_id.split("E")[0],
            }
        )
    write_tsv(
        outdir / "results/postnatal_pfc_science_aea1549/code_sample_id_audit.tsv",
        rows,
        ["sample_id_in_code", "condition_in_code", "donor_prefix_in_code"],
    )


def main() -> None:
    args = parse_args()
    metadata_dir = args.outdir / "metadata"
    supp_dir = args.outdir / "supplementary"
    extract_dir = supp_dir / "extracted"
    metadata_dir.mkdir(parents=True, exist_ok=True)
    supp_dir.mkdir(parents=True, exist_ok=True)

    record = json.loads(request_bytes(RECORD_URL).decode("utf-8"))
    (metadata_dir / "record.json").write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")

    manifest_rows: list[dict[str, str]] = []
    for file_info in record.get("files", []):
        key = file_info.get("key", "")
        url = file_info.get("links", {}).get("self", "")
        is_large_motif = key == "Motif_analysis.zip"
        should_download = key == "readme.txt" or (not is_large_motif or args.include_large_motif)
        local_path = ""
        status = "available_not_downloaded"
        if should_download:
            target = metadata_dir / key if key == "readme.txt" else supp_dir / safe_local_name(key)
            target.write_bytes(request_bytes(url))
            local_path = str(target)
            status = "downloaded"
            if key != "readme.txt" and not args.no_extract and not is_large_motif:
                maybe_extract(target, extract_dir)
        manifest_rows.append(
            {
                "record_id": str(record.get("id", "")),
                "doi": record.get("doi", ""),
                "file_key": key,
                "size_bytes": str(file_info.get("size", "")),
                "checksum": file_info.get("checksum", ""),
                "url": url,
                "local_path": local_path,
                "download_status": status,
                "notes": "large motif archive; download only if motif-output internals are needed"
                if is_large_motif
                else "public code/readme artifact",
            }
        )

    write_tsv(
        metadata_dir / "file_manifest.tsv",
        manifest_rows,
        [
            "record_id",
            "doi",
            "file_key",
            "size_bytes",
            "checksum",
            "url",
            "local_path",
            "download_status",
            "notes",
        ],
    )

    audit_sample_ids(
        extract_dir / "RNA preprocessing" / "step2_QC.R",
        Path("."),
    )
    sample_audit = Path("results/postnatal_pfc_science_aea1549/code_sample_id_audit.tsv")
    sample_rows = list(csv.DictReader(sample_audit.open(), delimiter="\t")) if sample_audit.exists() else []
    summary_rows = [
        {"metric": "zenodo_record_id", "value": str(record.get("id", ""))},
        {"metric": "zenodo_doi", "value": record.get("doi", "")},
        {"metric": "publication_title", "value": record.get("metadata", {}).get("title", "")},
        {"metric": "downloaded_file_count", "value": str(sum(r["download_status"] == "downloaded" for r in manifest_rows))},
        {
            "metric": "available_not_downloaded_file_count",
            "value": str(sum(r["download_status"] == "available_not_downloaded" for r in manifest_rows)),
        },
        {"metric": "code_sample_ids_detected", "value": str(len(sample_rows))},
        {"metric": "ts21_code_sample_ids_detected", "value": str(sum(r["condition_in_code"] == "Ts21" for r in sample_rows))},
        {
            "metric": "control_code_sample_ids_detected",
            "value": str(sum(r["condition_in_code"] == "Control" for r in sample_rows)),
        },
    ]
    write_tsv(
        Path("results/postnatal_pfc_science_aea1549/dataset_audit_summary.tsv"),
        summary_rows,
        ["metric", "value"],
    )


if __name__ == "__main__":
    main()
