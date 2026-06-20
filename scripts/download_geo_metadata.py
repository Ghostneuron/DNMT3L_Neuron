#!/usr/bin/env python3
"""Download and summarize GEO Series MINiML metadata.

This script reads `config/datasets.tsv`, downloads GEO Series MINiML files,
parses sample-level metadata, and writes small TSV summaries under
`data/geo/<GSE>/metadata/`.

It intentionally does not download large supplementary payloads unless
`--download-supplementary` is supplied.
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import sys
import tarfile
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Iterable


USER_AGENT = "DNMT3L-Neuron-metadata/0.1"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", default="config/datasets.tsv", help="Dataset manifest TSV.")
    parser.add_argument("--outdir", default="data/geo", help="Output root for GEO downloads.")
    parser.add_argument(
        "--accessions",
        nargs="*",
        help="Optional GEO accessions to process. Defaults to all GEO GSE rows in the manifest.",
    )
    parser.add_argument("--force", action="store_true", help="Redownload existing MINiML archives.")
    parser.add_argument("--timeout", type=int, default=60, help="Network timeout in seconds.")
    parser.add_argument(
        "--download-supplementary",
        action="store_true",
        help="Download supplementary files listed in MINiML metadata.",
    )
    parser.add_argument(
        "--max-file-mb",
        type=float,
        default=1024.0,
        help="Skip supplementary downloads above this size when Content-Length is available.",
    )
    parser.add_argument(
        "--include-supplementary-regex",
        help="Only download supplementary URLs whose URL or filename matches this regex.",
    )
    parser.add_argument(
        "--exclude-supplementary-regex",
        help="Skip supplementary URLs whose URL or filename matches this regex.",
    )
    return parser.parse_args()


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [row for row in reader]


def geo_bucket(accession: str) -> str:
    match = re.fullmatch(r"(GSE)(\d+)", accession.strip(), flags=re.IGNORECASE)
    if not match:
        raise ValueError(f"Not a GEO Series accession: {accession}")
    digits = match.group(2)
    if len(digits) <= 3:
        return "GSEnnn"
    return f"GSE{digits[:-3]}nnn"


def miniml_url(accession: str) -> str:
    acc = accession.upper()
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{geo_bucket(acc)}/{acc}/miniml/{acc}_family.xml.tgz"


def geo_url(accession: str) -> str:
    return f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession.upper()}"


def request(url: str, timeout: int, method: str = "GET") -> urllib.request.Request:
    return urllib.request.Request(url, method=method, headers={"User-Agent": USER_AGENT})


def download(url: str, dest: Path, timeout: int, force: bool = False) -> None:
    if dest.exists() and not force:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    with urllib.request.urlopen(request(url, timeout), timeout=timeout) as response, tmp.open("wb") as handle:
        shutil.copyfileobj(response, handle)
    tmp.replace(dest)


def extract_xml(archive: Path, dest: Path, force: bool = False) -> None:
    if dest.exists() and not force:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "r:gz") as tar:
        members = [m for m in tar.getmembers() if m.name.endswith(".xml")]
        if not members:
            raise RuntimeError(f"No XML file found in {archive}")
        member = members[0]
        source = tar.extractfile(member)
        if source is None:
            raise RuntimeError(f"Could not extract {member.name} from {archive}")
        with dest.open("wb") as handle:
            shutil.copyfileobj(source, handle)


def local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def children(elem: ET.Element, name: str) -> Iterable[ET.Element]:
    for child in elem:
        if local_name(child.tag) == name:
            yield child


def first_text(elem: ET.Element, name: str) -> str:
    for child in children(elem, name):
        if child.text:
            return " ".join(child.text.split())
    return ""


def all_text(elem: ET.Element, name: str) -> list[str]:
    values: list[str] = []
    for child in elem.iter():
        if local_name(child.tag) == name and child.text:
            values.append(" ".join(child.text.split()))
    return values


def collect_platform_refs(sample: ET.Element) -> str:
    refs = []
    for child in sample.iter():
        if local_name(child.tag) == "Platform-Ref":
            ref = child.attrib.get("ref", "")
            if ref:
                refs.append(ref)
    return "|".join(sorted(set(refs)))


def collect_characteristics(sample: ET.Element) -> str:
    parts = []
    for child in sample.iter():
        if local_name(child.tag) == "Characteristics":
            tag = child.attrib.get("tag", "").strip()
            value = " ".join((child.text or "").split())
            if tag and value:
                parts.append(f"{tag}={value}")
            elif value:
                parts.append(value)
    return "; ".join(parts)


def collect_supplementary(elem: ET.Element) -> list[str]:
    urls = []
    for child in elem.iter():
        if local_name(child.tag) == "Supplementary-Data" and child.text:
            urls.append(child.text.strip())
    return sorted(set(urls))


def parse_miniml(xml_path: Path) -> tuple[dict[str, str], list[dict[str, str]], list[dict[str, str]]]:
    root = ET.parse(xml_path).getroot()
    series = next((elem for elem in root.iter() if local_name(elem.tag) == "Series"), None)
    if series is None:
        raise RuntimeError(f"No Series element found in {xml_path}")

    series_acc = series.attrib.get("iid", "")
    series_summary = {
        "accession": series_acc,
        "title": first_text(series, "Title"),
        "summary": " ".join(all_text(series, "Summary")),
        "overall_design": " ".join(all_text(series, "Overall-Design")),
        "geo_url": geo_url(series_acc) if series_acc else "",
    }

    samples: list[dict[str, str]] = []
    supplementary: list[dict[str, str]] = []

    for url in collect_supplementary(series):
        supplementary.append({"context": "series", "accession": series_acc, "url": normalize_geo_ftp(url)})

    for sample in root.iter():
        if local_name(sample.tag) != "Sample":
            continue
        sample_acc = sample.attrib.get("iid", "")
        channel = next(children(sample, "Channel"), None)
        channel = channel if channel is not None else sample
        row = {
            "accession": sample_acc,
            "title": first_text(sample, "Title"),
            "organism": first_text(channel, "Organism"),
            "source_name": first_text(channel, "Source"),
            "molecule": first_text(channel, "Molecule"),
            "platform_refs": collect_platform_refs(sample),
            "characteristics": collect_characteristics(sample),
            "supplementary_data": "|".join(normalize_geo_ftp(url) for url in collect_supplementary(sample)),
        }
        samples.append(row)
        for url in collect_supplementary(sample):
            supplementary.append({"context": "sample", "accession": sample_acc, "url": normalize_geo_ftp(url)})

    series_summary["sample_count"] = str(len(samples))
    return series_summary, samples, supplementary


def normalize_geo_ftp(url: str) -> str:
    if url.startswith("ftp://ftp.ncbi.nlm.nih.gov/"):
        return "https://ftp.ncbi.nlm.nih.gov/" + url.removeprefix("ftp://ftp.ncbi.nlm.nih.gov/")
    return url


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def remote_size(url: str, timeout: int) -> int | None:
    try:
        with urllib.request.urlopen(request(url, timeout, method="HEAD"), timeout=timeout) as response:
            length = response.headers.get("Content-Length")
            return int(length) if length else None
    except Exception:
        return None


def download_supplementary(rows: list[dict[str, str]], dest_dir: Path, timeout: int, max_file_mb: float) -> list[dict[str, str]]:
    log_rows = []
    dest_dir.mkdir(parents=True, exist_ok=True)
    max_bytes = int(max_file_mb * 1024 * 1024)
    for row in rows:
        if row.get("status", "").startswith("skipped_"):
            log_rows.append(row)
            continue
        url = row["url"]
        if not url.startswith("http"):
            log_rows.append({**row, "status": "skipped_non_http", "path": "", "bytes": ""})
            continue
        size = remote_size(url, timeout)
        if size is not None and size > max_bytes:
            log_rows.append({**row, "status": "skipped_too_large", "path": "", "bytes": str(size)})
            continue
        filename = Path(urllib.parse.urlparse(url).path).name or f"{row['accession']}_supplementary"
        dest = dest_dir / filename
        try:
            download(url, dest, timeout=timeout, force=False)
            log_rows.append({**row, "status": "downloaded", "path": str(dest), "bytes": str(dest.stat().st_size)})
        except Exception as exc:
            log_rows.append({**row, "status": f"error: {exc}", "path": str(dest), "bytes": str(size or "")})
    return log_rows


def filter_supplementary_rows(
    rows: list[dict[str, str]],
    include_regex: str | None,
    exclude_regex: str | None,
) -> list[dict[str, str]]:
    include = re.compile(include_regex, flags=re.IGNORECASE) if include_regex else None
    exclude = re.compile(exclude_regex, flags=re.IGNORECASE) if exclude_regex else None
    filtered = []
    for row in rows:
        url = row.get("url", "")
        filename = Path(urllib.parse.urlparse(url).path).name
        target = f"{url} {filename}"
        if exclude and exclude.search(target):
            filtered.append({**row, "status": "skipped_excluded", "path": "", "bytes": ""})
            continue
        if include and not include.search(target):
            filtered.append({**row, "status": "skipped_not_included", "path": "", "bytes": ""})
            continue
        filtered.append(row)
    return filtered


def main() -> int:
    args = parse_args()
    manifest_path = Path(args.manifest)
    outdir = Path(args.outdir)
    requested = {acc.upper() for acc in args.accessions or []}

    rows = read_manifest(manifest_path)
    known_accessions = {row.get("accession", "").upper() for row in rows}
    for accession in sorted(requested - known_accessions):
        if accession.startswith("GSE"):
            rows.append(
                {
                    "accession": accession,
                    "priority": "external_requested",
                    "species": "",
                    "modality": "",
                }
            )
    geo_rows = [
        row
        for row in rows
        if row.get("accession", "").upper().startswith("GSE")
        and (not requested or row.get("accession", "").upper() in requested)
    ]
    if not geo_rows:
        print("No GEO GSE accessions selected.", file=sys.stderr)
        return 1

    summary_rows = []
    for row in geo_rows:
        acc = row["accession"].upper()
        base = outdir / acc
        archive = base / "raw" / f"{acc}_family.xml.tgz"
        xml_path = base / "metadata" / f"{acc}_family.xml"
        print(f"[{acc}] downloading metadata", file=sys.stderr)
        try:
            download(miniml_url(acc), archive, timeout=args.timeout, force=args.force)
            extract_xml(archive, xml_path, force=args.force)
            series, samples, supplementary = parse_miniml(xml_path)
        except urllib.error.HTTPError as exc:
            print(f"[{acc}] HTTP error: {exc}", file=sys.stderr)
            continue
        except Exception as exc:
            print(f"[{acc}] error: {exc}", file=sys.stderr)
            continue

        write_tsv(
            base / "metadata" / "series_metadata.tsv",
            [series],
            ["accession", "title", "summary", "overall_design", "sample_count", "geo_url"],
        )
        write_tsv(
            base / "metadata" / "sample_metadata.tsv",
            samples,
            [
                "accession",
                "title",
                "organism",
                "source_name",
                "molecule",
                "platform_refs",
                "characteristics",
                "supplementary_data",
            ],
        )
        write_tsv(base / "metadata" / "supplementary_files.tsv", supplementary, ["context", "accession", "url"])

        if args.download_supplementary:
            filtered_supplementary = filter_supplementary_rows(
                supplementary,
                include_regex=args.include_supplementary_regex,
                exclude_regex=args.exclude_supplementary_regex,
            )
            download_log = download_supplementary(
                filtered_supplementary,
                base / "supplementary",
                timeout=args.timeout,
                max_file_mb=args.max_file_mb,
            )
            write_tsv(
                base / "metadata" / "supplementary_download_log.tsv",
                download_log,
                ["context", "accession", "url", "status", "path", "bytes"],
            )

        summary_rows.append(
            {
                "accession": acc,
                "title": series.get("title", ""),
                "sample_count": series.get("sample_count", ""),
                "priority": row.get("priority", ""),
                "species": row.get("species", ""),
                "modality": row.get("modality", ""),
                "metadata_dir": str(base / "metadata"),
            }
        )

    write_tsv(
        outdir / "geo_metadata_download_summary.tsv",
        summary_rows,
        ["accession", "title", "sample_count", "priority", "species", "modality", "metadata_dir"],
    )
    print(f"Wrote {len(summary_rows)} dataset summaries to {outdir}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
