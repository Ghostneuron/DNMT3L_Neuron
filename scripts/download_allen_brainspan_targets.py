#!/usr/bin/env python3
"""Download BrainSpan metadata and extract target-gene RNA-seq RPKM values.

The Allen BrainSpan expression service can be fragile, but each specimen has a
well-known gene-level RNA-seq file. This script streams those files and keeps
only the target genes from config/target_genes.tsv.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import sys
import time
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


API_BASE = "https://api.brain-map.org"
QUERY_URL = f"{API_BASE}/api/v2/data/query.json"
PRODUCT_ID = 24
GENE_FILE_TYPE_ID = 121134332


METADATA_FIELDS = [
    "specimen_id",
    "specimen_name",
    "external_specimen_name",
    "donor_id",
    "donor_name",
    "external_donor_name",
    "sex",
    "age_name",
    "age_days",
    "age_embryonic",
    "age_tags",
    "structure_id",
    "structure_acronym",
    "structure_name",
    "hemisphere",
    "rna_integrity_number",
    "tissue_ph",
    "gene_file_download_link",
    "gene_file_path",
]

EXPRESSION_FIELDS = METADATA_FIELDS + [
    "gene_symbol",
    "ensembl_id",
    "entrez_id",
    "rpkm",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--target-genes",
        default="config/target_genes.tsv",
        help="Target gene manifest with species, gene_symbol, and ensembl_hint columns.",
    )
    parser.add_argument(
        "--outdir",
        default="data/allen/brainspan/RNASeq_Gencode_v10",
        help="Output directory for BrainSpan metadata and target-gene table.",
    )
    parser.add_argument(
        "--page-size",
        type=int,
        default=200,
        help="Allen API page size for specimen metadata.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional specimen limit for smoke tests.",
    )
    parser.add_argument(
        "--sleep",
        type=float,
        default=0.05,
        help="Delay between specimen file requests.",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=60.0,
        help="Network timeout in seconds.",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=3,
        help="Retries for each network request.",
    )
    parser.add_argument(
        "--metadata-only",
        action="store_true",
        help="Only fetch specimen metadata; skip expression file streaming.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite target-gene output instead of resuming from existing rows.",
    )
    return parser.parse_args()


def fetch_url(url: str, timeout: float, retries: int) -> bytes:
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            request = Request(url, headers={"User-Agent": "DNMT3L_Neuron/0.1"})
            with urlopen(request, timeout=timeout) as response:
                return response.read()
        except (HTTPError, URLError, TimeoutError) as exc:
            last_error = exc
            if attempt == retries:
                break
            time.sleep(min(2.0 * attempt, 10.0))
    raise RuntimeError(f"failed to fetch {url}: {last_error}")


def fetch_json(criteria: str, timeout: float, retries: int) -> dict[str, Any]:
    url = f"{QUERY_URL}?{urlencode({'criteria': criteria})}"
    payload = fetch_url(url, timeout=timeout, retries=retries)
    data = json.loads(payload.decode("utf-8"))
    if not data.get("success"):
        raise RuntimeError(f"Allen API query failed: {data.get('msg')}")
    return data


def load_human_targets(path: Path) -> tuple[set[str], set[str]]:
    symbols: set[str] = set()
    ensembl_ids: set[str] = set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["species"].lower() != "human":
                continue
            symbols.add(row["gene_symbol"].upper())
            if row.get("ensembl_hint"):
                ensembl_ids.add(row["ensembl_hint"])
            for alias in row.get("aliases", "").split("|"):
                if alias:
                    symbols.add(alias.upper())
    return symbols, ensembl_ids


def get_gene_file(specimen: dict[str, Any]) -> dict[str, Any] | None:
    for file_info in specimen.get("well_known_files", []):
        path = file_info.get("path", "")
        if (
            file_info.get("well_known_file_type_id") == GENE_FILE_TYPE_ID
            or path.endswith("_RNASeq_Gene.txt")
        ):
            return file_info
    return None


def flatten_specimen(specimen: dict[str, Any]) -> dict[str, str]:
    donor = specimen.get("donor") or {}
    age = donor.get("age") or {}
    structure = specimen.get("structure") or {}
    gene_file = get_gene_file(specimen) or {}
    download_link = gene_file.get("download_link", "")
    if download_link and download_link.startswith("/"):
        download_link = f"{API_BASE}{download_link}"
    return {
        "specimen_id": str(specimen.get("id", "")),
        "specimen_name": str(specimen.get("name", "")),
        "external_specimen_name": str(specimen.get("external_specimen_name", "")),
        "donor_id": str(specimen.get("donor_id", "")),
        "donor_name": str(donor.get("name", "")),
        "external_donor_name": str(donor.get("external_donor_name", "")),
        "sex": str(donor.get("sex", "")),
        "age_name": str(age.get("name", "")),
        "age_days": str(age.get("days", "")),
        "age_embryonic": str(age.get("embryonic", "")),
        "age_tags": str(age.get("tags", "")),
        "structure_id": str(specimen.get("structure_id", "")),
        "structure_acronym": str(structure.get("acronym", "")),
        "structure_name": str(structure.get("name", "")),
        "hemisphere": str(specimen.get("hemisphere", "")),
        "rna_integrity_number": str(specimen.get("rna_integrity_number", "")),
        "tissue_ph": str(specimen.get("tissue_ph", "")),
        "gene_file_download_link": download_link,
        "gene_file_path": str(gene_file.get("path", "")),
    }


def fetch_specimens(page_size: int, timeout: float, retries: int) -> list[dict[str, Any]]:
    specimens: list[dict[str, Any]] = []
    start = 0
    total: int | None = None
    while total is None or start < total:
        criteria = (
            "model::Specimen,"
            "rma::criteria,"
            f"donor(products[id$eq{PRODUCT_ID}]),"
            "well_known_files(well_known_file_type[name$eqRNASeqSummarizedToGenes]),"
            "rma::include,well_known_files,donor(age),structure,"
            f"rma::options[start_row$eq{start}][num_rows$eq{page_size}]"
        )
        data = fetch_json(criteria, timeout=timeout, retries=retries)
        total = int(data["total_rows"])
        page = data.get("msg", [])
        specimens.extend(page)
        if not page:
            break
        start += len(page)
    return specimens


def write_metadata(rows: list[dict[str, str]], path: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=METADATA_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def completed_specimens(path: Path) -> set[str]:
    if not path.exists():
        return set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["specimen_id"] for row in reader if row.get("specimen_id")}


def iter_target_rows(
    metadata: dict[str, str],
    target_symbols: set[str],
    target_ensembl_ids: set[str],
    timeout: float,
    retries: int,
) -> list[dict[str, str]]:
    url = metadata["gene_file_download_link"]
    if not url:
        raise RuntimeError(f"missing gene file link for specimen {metadata['specimen_id']}")

    request = Request(url, headers={"User-Agent": "DNMT3L_Neuron/0.1"})
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            with urlopen(request, timeout=timeout) as response:
                text = io.TextIOWrapper(response, encoding="utf-8", errors="replace")
                header = next(text).rstrip("\n").split("\t")
                if header[:3] != ["EntrezID", "GencodeID", "RPKM"]:
                    raise RuntimeError(f"unexpected header for {url}: {header}")
                rows: list[dict[str, str]] = []
                for line in text:
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 3:
                        continue
                    entrez_id, gencode_id, rpkm = fields[:3]
                    if "|" in gencode_id:
                        ensembl_id, gene_symbol = gencode_id.split("|", 1)
                    else:
                        ensembl_id, gene_symbol = gencode_id, ""
                    if (
                        gene_symbol.upper() not in target_symbols
                        and ensembl_id not in target_ensembl_ids
                    ):
                        continue
                    row = dict(metadata)
                    row.update(
                        {
                            "gene_symbol": gene_symbol,
                            "ensembl_id": ensembl_id,
                            "entrez_id": entrez_id,
                            "rpkm": rpkm,
                        }
                    )
                    rows.append(row)
                return rows
        except (HTTPError, URLError, TimeoutError, RuntimeError) as exc:
            last_error = exc
            if attempt == retries:
                break
            time.sleep(min(2.0 * attempt, 10.0))
    raise RuntimeError(f"failed to stream {url}: {last_error}")


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    target_symbols, target_ensembl_ids = load_human_targets(Path(args.target_genes))
    specimens = fetch_specimens(args.page_size, args.timeout, args.retries)
    if args.limit is not None:
        specimens = specimens[: args.limit]

    metadata_rows = [flatten_specimen(specimen) for specimen in specimens]
    metadata_path = outdir / "brainspan_rnaseq_specimen_metadata.tsv"
    write_metadata(metadata_rows, metadata_path)
    print(f"Wrote {len(metadata_rows)} specimen metadata rows to {metadata_path}")

    if args.metadata_only:
        return 0

    expression_path = outdir / "brainspan_target_gene_rpkm.tsv"
    failures_path = outdir / "brainspan_target_gene_download_failures.tsv"
    if args.force and expression_path.exists():
        expression_path.unlink()
    done = completed_specimens(expression_path)
    mode = "a" if expression_path.exists() else "w"

    total_target_rows = 0
    failures: list[dict[str, str]] = []
    with expression_path.open(mode, newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=EXPRESSION_FIELDS, delimiter="\t")
        if mode == "w":
            writer.writeheader()
        for index, metadata in enumerate(metadata_rows, start=1):
            specimen_id = metadata["specimen_id"]
            if specimen_id in done:
                continue
            try:
                rows = iter_target_rows(
                    metadata,
                    target_symbols,
                    target_ensembl_ids,
                    timeout=args.timeout,
                    retries=args.retries,
                )
                writer.writerows(rows)
                handle.flush()
                total_target_rows += len(rows)
            except RuntimeError as exc:
                failures.append(
                    {
                        "specimen_id": specimen_id,
                        "specimen_name": metadata["specimen_name"],
                        "gene_file_download_link": metadata["gene_file_download_link"],
                        "error": str(exc),
                    }
                )
            if index % 25 == 0:
                print(f"Processed {index}/{len(metadata_rows)} specimens", file=sys.stderr)
            if args.sleep:
                time.sleep(args.sleep)

    with failures_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "specimen_id",
                "specimen_name",
                "gene_file_download_link",
                "error",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(failures)

    existing = sum(1 for _ in expression_path.open()) - 1
    print(f"Wrote/retained {existing} target-gene rows to {expression_path}")
    print(f"Wrote {len(failures)} failures to {failures_path}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
