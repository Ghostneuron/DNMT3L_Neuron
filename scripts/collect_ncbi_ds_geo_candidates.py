#!/usr/bin/env python3
"""Collect Down syndrome/trisomy 21 brain-neural GEO candidates from NCBI.

The goal is broad recall, not final curation. The script queries NCBI GDS with
multiple overlapping terms, de-duplicates GEO Series records, scores each record
for brain/neural relevance, and writes compact TSV inventories for follow-up.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Iterable


USER_AGENT = "DNMT3L-Neuron-NCBI-GEO-search/0.1"

DEFAULT_QUERIES = [
    "trisomy brain human AND gse[ETYP]",
    '"Down syndrome" brain human AND gse[ETYP]',
    '"trisomy 21" brain human AND gse[ETYP]',
    '"Down syndrome" cortex human AND gse[ETYP]',
    '"trisomy 21" cortex human AND gse[ETYP]',
    '"Down syndrome" "fetal brain" human AND gse[ETYP]',
    '"trisomy 21" "fetal brain" human AND gse[ETYP]',
    '"Down syndrome" neuron human AND gse[ETYP]',
    '"trisomy 21" neuron human AND gse[ETYP]',
    '"Down syndrome" neural human AND gse[ETYP]',
    '"trisomy 21" neural human AND gse[ETYP]',
    '"Down syndrome" iPSC neural human AND gse[ETYP]',
    '"trisomy 21" iPSC neural human AND gse[ETYP]',
    '"Down syndrome" organoid human AND gse[ETYP]',
    '"trisomy 21" organoid human AND gse[ETYP]',
    '"Down syndrome" methylation brain human AND gse[ETYP]',
    '"trisomy 21" methylation brain human AND gse[ETYP]',
    '"Down syndrome" ATAC brain human AND gse[ETYP]',
    '"trisomy 21" ATAC brain human AND gse[ETYP]',
]

DIRECT_BRAIN_RE = re.compile(
    r"\b(brain|cortex|cortical|cerebral|cerebrum|cerebellum|hippocamp|prefrontal|"
    r"pfc|forebrain|telencephalon|fetal brain|foetal brain)\b",
    flags=re.IGNORECASE,
)
NEURAL_MODEL_RE = re.compile(
    r"\b(neural|neuron|neuronal|neuro|npc|nsc|progenitor|organoid|iPSC|stem cell|"
    r"astrocyte|glia|glial|xenograft)\b",
    flags=re.IGNORECASE,
)
DS_RE = re.compile(r"\b(down syndrome|trisomy 21|t21|hsa21)\b", flags=re.IGNORECASE)
CHR21_RE = re.compile(r"\b(chromosome 21|chr21|human chromosome 21|HSA21)\b", flags=re.IGNORECASE)
METHYLATION_RE = re.compile(r"\b(methylation|bisulfite|methylome|450k|epic|wgbs)\b", flags=re.IGNORECASE)
SINGLE_CELL_RE = re.compile(r"\b(single.cell|single nucleus|snrna|scrna|multiome|10x|mtx)\b", flags=re.IGNORECASE)
REGULATORY_RE = re.compile(r"\b(atac|chip|cut&tag|cut&run|hic|hi-c|chromatin|occupancy)\b", flags=re.IGNORECASE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="results/ds_ncbi_geo", help="Output directory.")
    parser.add_argument("--retmax", type=int, default=500, help="Maximum NCBI hits per query.")
    parser.add_argument("--sleep", type=float, default=0.34, help="Delay between NCBI requests.")
    parser.add_argument("--timeout", type=int, default=60, help="Network timeout in seconds.")
    parser.add_argument("--query", action="append", help="Additional query term. Can be repeated.")
    return parser.parse_args()


def ncbi_json(url: str, timeout: int) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8"))


def esearch(term: str, retmax: int, timeout: int) -> tuple[list[str], int]:
    params = urllib.parse.urlencode(
        {
            "db": "gds",
            "term": term,
            "retmode": "json",
            "retmax": str(retmax),
        }
    )
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?{params}"
    data = ncbi_json(url, timeout)
    result = data.get("esearchresult", {})
    return result.get("idlist", []), int(result.get("count", 0))


def chunks(items: list[str], size: int) -> Iterable[list[str]]:
    for start in range(0, len(items), size):
        yield items[start : start + size]


def esummary(uids: list[str], timeout: int, sleep: float) -> dict[str, dict]:
    records: dict[str, dict] = {}
    for group in chunks(uids, 100):
        params = urllib.parse.urlencode({"db": "gds", "id": ",".join(group), "retmode": "json"})
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?{params}"
        data = ncbi_json(url, timeout).get("result", {})
        for uid in data.get("uids", []):
            item = data.get(uid, {})
            if item:
                records[uid] = item
        time.sleep(sleep)
    return records


def text_blob(record: dict) -> str:
    sample_titles = " ".join(sample.get("title", "") for sample in record.get("samples", []))
    return " ".join(
        str(record.get(key, ""))
        for key in ["title", "summary", "gdstype", "suppfile", "taxon"]
    ) + " " + sample_titles


def score_record(record: dict) -> dict[str, str]:
    blob = text_blob(record)
    direct_brain = bool(DIRECT_BRAIN_RE.search(blob))
    neural_model = bool(NEURAL_MODEL_RE.search(blob))
    ds_related = bool(DS_RE.search(blob))
    chromosome21_context = bool(CHR21_RE.search(blob))
    methylation = bool(METHYLATION_RE.search(blob))
    single_cell = bool(SINGLE_CELL_RE.search(blob))
    regulatory = bool(REGULATORY_RE.search(blob))

    score = 0
    score += 5 if direct_brain else 0
    score += 3 if neural_model else 0
    score += 3 if ds_related else 0
    score += 1 if chromosome21_context and not ds_related else 0
    score += 2 if methylation else 0
    score += 2 if single_cell else 0
    score += 1 if regulatory else 0
    if record.get("taxon") == "Homo sapiens":
        score += 1

    if direct_brain and single_cell:
        category = "direct_brain_single_cell"
    elif direct_brain and methylation:
        category = "direct_brain_methylation"
    elif direct_brain:
        category = "direct_brain_expression"
    elif neural_model and single_cell:
        category = "neural_model_single_cell"
    elif neural_model and methylation:
        category = "neural_model_methylation"
    elif neural_model:
        category = "neural_model_expression"
    elif methylation:
        category = "non_neural_methylation_or_context"
    else:
        category = "low_relevance_or_context"

    return {
        "relevance_score": str(score),
        "category": category,
        "ds_related": "yes" if ds_related else "no",
        "chromosome21_context": "yes" if chromosome21_context else "no",
        "direct_brain": "yes" if direct_brain else "no",
        "neural_model": "yes" if neural_model else "no",
        "methylation_related": "yes" if methylation else "no",
        "single_cell_related": "yes" if single_cell else "no",
        "regulatory_related": "yes" if regulatory else "no",
    }


def geo_url(accession: str) -> str:
    return f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}"


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    queries = DEFAULT_QUERIES + (args.query or [])

    uid_to_queries: dict[str, list[str]] = {}
    query_rows: list[dict[str, str]] = []
    for query in queries:
        ids, count = esearch(query, args.retmax, args.timeout)
        for uid in ids:
            uid_to_queries.setdefault(uid, []).append(query)
        query_rows.append({"query": query, "ncbi_count": str(count), "uids_returned": str(len(ids))})
        time.sleep(args.sleep)

    records = esummary(sorted(uid_to_queries), args.timeout, args.sleep)
    rows: list[dict[str, str]] = []
    for uid, record in records.items():
        if record.get("entrytype") != "GSE":
            continue
        accession = record.get("accession", "")
        sample_titles = "; ".join(sample.get("title", "") for sample in record.get("samples", [])[:8])
        row = {
            "uid": uid,
            "accession": accession,
            "title": record.get("title", ""),
            "summary": " ".join(str(record.get("summary", "")).split()),
            "taxon": record.get("taxon", ""),
            "gdstype": record.get("gdstype", ""),
            "sample_count": str(record.get("n_samples", "")),
            "suppfile": record.get("suppfile", ""),
            "publication_date": record.get("pdat", ""),
            "pubmed_ids": "|".join(record.get("pubmedids", [])),
            "gpl": record.get("gpl", ""),
            "ftplink": record.get("ftplink", ""),
            "geo_url": geo_url(accession),
            "matched_queries": " | ".join(uid_to_queries.get(uid, [])),
            "first_sample_titles": sample_titles,
        }
        row.update(score_record(record))
        rows.append(row)

    rows.sort(key=lambda r: (-int(r["relevance_score"]), r["accession"]))
    high_priority = [
        row
        for row in rows
        if row["taxon"] == "Homo sapiens"
        and row["ds_related"] == "yes"
        and (row["direct_brain"] == "yes" or row["neural_model"] == "yes")
    ]
    direct_brain = [row for row in high_priority if row["direct_brain"] == "yes"]

    fields = [
        "uid",
        "accession",
        "relevance_score",
        "category",
        "ds_related",
        "chromosome21_context",
        "direct_brain",
        "neural_model",
        "methylation_related",
        "single_cell_related",
        "regulatory_related",
        "title",
        "taxon",
        "gdstype",
        "sample_count",
        "suppfile",
        "publication_date",
        "pubmed_ids",
        "gpl",
        "ftplink",
        "geo_url",
        "matched_queries",
        "first_sample_titles",
        "summary",
    ]
    write_tsv(outdir / "ncbi_geo_search_queries.tsv", query_rows, ["query", "ncbi_count", "uids_returned"])
    write_tsv(outdir / "ncbi_geo_all_gse_candidates.tsv", rows, fields)
    write_tsv(outdir / "ncbi_geo_high_priority_brain_neural_candidates.tsv", high_priority, fields)
    write_tsv(outdir / "ncbi_geo_direct_brain_candidates.tsv", direct_brain, fields)

    print(f"Wrote {len(rows)} GSE candidates to {outdir / 'ncbi_geo_all_gse_candidates.tsv'}")
    print(f"High-priority brain/neural DS candidates: {len(high_priority)}")
    print(f"Direct brain candidates: {len(direct_brain)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
