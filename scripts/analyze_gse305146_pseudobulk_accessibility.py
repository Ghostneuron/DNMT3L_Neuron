#!/usr/bin/env python3
"""Pseudo-bulk target-region accessibility from GSE305146 fetal brain multiome matrices."""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import math
import re
import statistics
import urllib.request
from collections import defaultdict
from pathlib import Path

from scipy import stats


BASE_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE305nnn/GSE305146/suppl"
USER_AGENT = "DNMT3L-Neuron-GSE305146-accessibility/0.1"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-genes", default="config/target_genes.tsv")
    parser.add_argument("--final-candidates", default="results/manuscript_summary/final_candidate_manuscript_summary.tsv")
    parser.add_argument("--metadata", default="data/geo/GSE305146/metadata/sample_metadata.tsv")
    parser.add_argument("--local-supplementary", default="data/geo/GSE305146/supplementary")
    parser.add_argument("--outdir", default="results/mechanism_question")
    parser.add_argument("--note", default="notes/gse305146_scATAC_accessibility_proxy_first_pass.md")
    parser.add_argument("--window-bp", type=int, default=100_000)
    parser.add_argument("--timeout", type=int, default=240)
    parser.add_argument("--limit-libraries", type=int, help="Debug: process only the first N selected libraries.")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", errors="replace") as handle:
        return [{k: (v or "").strip() for k, v in row.items()} for row in csv.DictReader(handle, delimiter="\t")]


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fieldnames})


def parse_characteristics(text: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for part in text.split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
            values[key.strip().lower()] = value.strip()
    return values


def normalize_group(group: str) -> str:
    group_norm = group.strip().lower()
    if group_norm in {"ds", "trisomy 21", "t21", "down syndrome"}:
        return "DS"
    if group_norm in {"con", "control", "disomic", "euploid"}:
        return "Control"
    return group.strip()


def library_prefix(title: str) -> str:
    match = re.search(r"Library\s+([^,\s]+)_gex", title)
    return match.group(1) if match else ""


def selected_libraries(metadata_path: Path) -> list[dict[str, str]]:
    libraries: list[dict[str, str]] = []
    for row in read_tsv(metadata_path):
        chars = parse_characteristics(row.get("characteristics", ""))
        if chars.get("library type") != "Gene Expression":
            continue
        if chars.get("in final_analysis", "").upper() == "FALSE":
            continue
        prefix = library_prefix(row.get("title", ""))
        if not prefix:
            continue
        libraries.append(
            {
                "library_prefix": prefix,
                "geo_sample_accession": row.get("accession", ""),
                "biological_sample_id": chars.get("sample id", ""),
                "sample_name": chars.get("name", ""),
                "group": normalize_group(chars.get("group", chars.get("genotype", ""))),
                "genotype": chars.get("genotype", ""),
                "dev_stage_pcw": chars.get("dev stage_pcw", ""),
                "sex": chars.get("sex", ""),
                "platform_refs": row.get("platform_refs", ""),
            }
        )
    return libraries


def load_targets(target_gene_path: Path, final_candidate_path: Path) -> dict[str, dict[str, str]]:
    targets: dict[str, dict[str, str]] = {}
    for row in read_tsv(target_gene_path):
        if row.get("species") != "human":
            continue
        symbol = row["gene_symbol"].upper()
        targets[symbol] = {
            "gene_symbol": symbol,
            "ensembl_id": row.get("ensembl_hint", ""),
            "target_class": "mechanism_core",
            "tier": "",
            "candidate_group": "",
            "mechanism_class": "",
        }

    for row in read_tsv(final_candidate_path):
        symbol = row.get("gene_symbol", "").upper()
        if not symbol:
            continue
        targets.setdefault(
            symbol,
            {
                "gene_symbol": symbol,
                "ensembl_id": "",
                "target_class": "final_candidate",
                "tier": row.get("tier", ""),
                "candidate_group": row.get("candidate_group", ""),
                "mechanism_class": row.get("mechanism_class", ""),
            },
        )
        if targets[symbol]["target_class"] == "mechanism_core":
            targets[symbol]["target_class"] = "mechanism_core_and_final_candidate"
        targets[symbol]["tier"] = row.get("tier", targets[symbol].get("tier", ""))
        targets[symbol]["candidate_group"] = row.get("candidate_group", targets[symbol].get("candidate_group", ""))
        targets[symbol]["mechanism_class"] = row.get("mechanism_class", targets[symbol].get("mechanism_class", ""))
    return targets


def open_text_gzip(path_or_url: str, timeout: int) -> io.TextIOBase:
    if path_or_url.startswith("http"):
        request = urllib.request.Request(path_or_url, headers={"User-Agent": USER_AGENT})
        response = urllib.request.urlopen(request, timeout=timeout)
        return io.TextIOWrapper(gzip.GzipFile(fileobj=response))
    return gzip.open(path_or_url, "rt")


def local_or_url(local_root: Path, prefix: str, kind: str) -> str:
    local_path = local_root / f"GSE305146_{prefix}_{kind}"
    if local_path.exists():
        return str(local_path)
    return f"{BASE_URL}/GSE305146_{prefix}_{kind}"


def safe_int(value: str) -> int:
    try:
        return int(float(value))
    except ValueError:
        return 0


def parse_features(
    path_or_url: str,
    targets: dict[str, dict[str, str]],
    window_bp: int,
    timeout: int,
) -> tuple[dict[int, list[str]], dict[str, dict[str, object]], int, int]:
    by_ensembl = {row["ensembl_id"]: symbol for symbol, row in targets.items() if row.get("ensembl_id")}
    target_symbols = set(targets)
    gene_windows: dict[str, dict[str, object]] = {}
    peak_row_to_symbols: dict[int, list[str]] = {}
    last_gene_expression_row = 0
    n_peak_features = 0

    with open_text_gzip(path_or_url, timeout) as handle:
        for index, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            feature_id, symbol_raw, feature_type, chrom, start_raw, end_raw = fields[:6]
            symbol = symbol_raw.upper()
            start = safe_int(start_raw)
            end = safe_int(end_raw)
            if feature_type == "Gene Expression":
                last_gene_expression_row = index
                target_symbol = by_ensembl.get(feature_id) or (symbol if symbol in target_symbols else "")
                if target_symbol:
                    gene_windows[target_symbol] = {
                        "gene_symbol": target_symbol,
                        "ensembl_id": feature_id,
                        "chrom": chrom,
                        "gene_start": start,
                        "gene_end": end,
                        "window_start": max(0, start - window_bp),
                        "window_end": end + window_bp,
                    }
                continue

            if feature_type != "Peaks":
                continue
            n_peak_features += 1
            overlapping = []
            for target_symbol, window in gene_windows.items():
                if chrom != window["chrom"]:
                    continue
                if end >= int(window["window_start"]) and start <= int(window["window_end"]):
                    overlapping.append(target_symbol)
            if overlapping:
                peak_row_to_symbols[index] = overlapping

    for target_symbol, window in gene_windows.items():
        window["n_window_peak_features"] = sum(1 for symbols in peak_row_to_symbols.values() if target_symbol in symbols)
    return peak_row_to_symbols, gene_windows, last_gene_expression_row, n_peak_features


def stream_matrix_accessibility(
    path_or_url: str,
    peak_row_to_symbols: dict[int, list[str]],
    last_gene_expression_row: int,
    target_symbols: set[str],
    timeout: int,
) -> dict[str, object]:
    peak_counts = {symbol: 0.0 for symbol in target_symbols}
    peak_nonzero_entries = {symbol: 0 for symbol in target_symbols}
    total_peak_counts = 0.0
    n_features = n_cells = n_entries = 0
    dims_seen = False

    with open_text_gzip(path_or_url, timeout) as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            fields = line.strip().split()
            if not fields:
                continue
            if not dims_seen:
                n_features, n_cells, n_entries = (int(fields[0]), int(fields[1]), int(fields[2]))
                dims_seen = True
                continue
            row_index = int(fields[0])
            value = float(fields[2])
            if row_index > last_gene_expression_row:
                total_peak_counts += value
            for symbol in peak_row_to_symbols.get(row_index, []):
                peak_counts[symbol] += value
                peak_nonzero_entries[symbol] += 1

    return {
        "n_features": n_features,
        "n_cells": n_cells,
        "n_entries": n_entries,
        "total_peak_counts": total_peak_counts,
        "peak_counts": peak_counts,
        "peak_nonzero_entries": peak_nonzero_entries,
    }


def fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}" if math.isfinite(value) else ""
    return str(value)


def mean(values: list[float]) -> float:
    return statistics.fmean(values) if values else math.nan


def median(values: list[float]) -> float:
    return statistics.median(values) if values else math.nan


def safe_float(value: object) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return math.nan
    return number if math.isfinite(number) else math.nan


def safe_ttest(a: list[float], b: list[float]) -> float:
    a = [value for value in a if math.isfinite(value)]
    b = [value for value in b if math.isfinite(value)]
    if len(a) < 2 or len(b) < 2:
        return math.nan
    _stat, pvalue = stats.ttest_ind(a, b, equal_var=False)
    return pvalue if math.isfinite(pvalue) else math.nan


def bh_fdr(pvalues: list[float]) -> list[float]:
    indexed = [(idx, p) for idx, p in enumerate(pvalues) if math.isfinite(p)]
    fdr = [math.nan] * len(pvalues)
    if not indexed:
        return fdr
    indexed.sort(key=lambda item: item[1])
    n = len(indexed)
    running = 1.0
    for rank, (idx, pvalue) in reversed(list(enumerate(indexed, start=1))):
        running = min(running, pvalue * n / rank)
        fdr[idx] = running
    return fdr


def build_library_rows(args: argparse.Namespace) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    targets = load_targets(Path(args.target_genes), Path(args.final_candidates))
    libraries = selected_libraries(Path(args.metadata))
    if args.limit_libraries:
        libraries = libraries[: args.limit_libraries]

    local_root = Path(args.local_supplementary)
    library_rows: list[dict[str, object]] = []
    readiness_rows: list[dict[str, object]] = []
    window_rows: list[dict[str, object]] = []

    for idx, library in enumerate(libraries, start=1):
        prefix = library["library_prefix"]
        print(f"[{idx}/{len(libraries)}] streaming accessibility {prefix}", flush=True)
        features_path = local_or_url(local_root, prefix, "features.tsv.gz")
        matrix_path = local_or_url(local_root, prefix, "matrix.mtx.gz")
        peak_row_to_symbols, gene_windows, last_gene_row, n_peak_features = parse_features(
            features_path, targets, args.window_bp, args.timeout
        )
        matrix_stats = stream_matrix_accessibility(
            matrix_path, peak_row_to_symbols, last_gene_row, set(gene_windows), args.timeout
        )
        total_peak_counts = float(matrix_stats["total_peak_counts"])
        peak_counts = matrix_stats["peak_counts"]
        peak_nonzero = matrix_stats["peak_nonzero_entries"]
        assert isinstance(peak_counts, dict)
        assert isinstance(peak_nonzero, dict)

        readiness_rows.append(
            {
                "library_prefix": prefix,
                "biological_sample_id": library["biological_sample_id"],
                "group": library["group"],
                "dev_stage_pcw": library["dev_stage_pcw"],
                "n_features": matrix_stats["n_features"],
                "n_cells": matrix_stats["n_cells"],
                "n_entries": matrix_stats["n_entries"],
                "last_gene_expression_row": last_gene_row,
                "n_peak_features": n_peak_features,
                "target_genes_with_coordinates": len(gene_windows),
                "target_peak_rows": len(peak_row_to_symbols),
                "total_peak_counts": total_peak_counts,
            }
        )

        for symbol, window in sorted(gene_windows.items()):
            target = targets.get(symbol, {})
            count = float(peak_counts.get(symbol, 0.0))
            cpm = count / total_peak_counts * 1_000_000 if total_peak_counts else math.nan
            row = {
                "dataset": "GSE305146",
                "library_prefix": prefix,
                "geo_sample_accession": library["geo_sample_accession"],
                "biological_sample_id": library["biological_sample_id"],
                "sample_name": library["sample_name"],
                "group": library["group"],
                "genotype": library["genotype"],
                "dev_stage_pcw": library["dev_stage_pcw"],
                "sex": library["sex"],
                "platform_refs": library["platform_refs"],
                "gene_symbol": symbol,
                "target_class": target.get("target_class", ""),
                "tier": target.get("tier", ""),
                "candidate_group": target.get("candidate_group", ""),
                "mechanism_class": target.get("mechanism_class", ""),
                "ensembl_id": window.get("ensembl_id", ""),
                "chrom": window.get("chrom", ""),
                "gene_start": window.get("gene_start", ""),
                "gene_end": window.get("gene_end", ""),
                "window_start": window.get("window_start", ""),
                "window_end": window.get("window_end", ""),
                "window_bp": args.window_bp,
                "n_window_peak_features": window.get("n_window_peak_features", 0),
                "n_cells": matrix_stats["n_cells"],
                "total_peak_counts": total_peak_counts,
                "target_peak_counts": count,
                "target_peak_cpm": cpm,
                "target_log2_peak_cpm": math.log2(cpm + 1.0) if math.isfinite(cpm) else math.nan,
                "target_peak_nonzero_entries": peak_nonzero.get(symbol, 0),
            }
            library_rows.append(row)
            if idx == 1:
                window_rows.append({key: row.get(key, "") for key in ["gene_symbol", "target_class", "tier", "candidate_group", "mechanism_class", "ensembl_id", "chrom", "gene_start", "gene_end", "window_start", "window_end", "window_bp", "n_window_peak_features"]})

    return library_rows, readiness_rows, window_rows


def aggregate_biological_samples(library_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in library_rows:
        key = (
            str(row.get("biological_sample_id", "")),
            str(row.get("gene_symbol", "")),
            str(row.get("group", "")),
        )
        grouped[key].append(row)

    sample_rows: list[dict[str, object]] = []
    for (sample_id, gene, group), rows in sorted(grouped.items()):
        total_peak_counts = sum(safe_float(row.get("total_peak_counts")) for row in rows)
        target_peak_counts = sum(safe_float(row.get("target_peak_counts")) for row in rows)
        n_cells = sum(int(safe_float(row.get("n_cells"))) for row in rows if math.isfinite(safe_float(row.get("n_cells"))))
        cpm = target_peak_counts / total_peak_counts * 1_000_000 if total_peak_counts else math.nan
        first = rows[0]
        sample_rows.append(
            {
                "dataset": "GSE305146",
                "biological_sample_id": sample_id,
                "group": group,
                "genotype": first.get("genotype", ""),
                "dev_stage_pcw": first.get("dev_stage_pcw", ""),
                "sex": first.get("sex", ""),
                "gene_symbol": gene,
                "target_class": first.get("target_class", ""),
                "tier": first.get("tier", ""),
                "candidate_group": first.get("candidate_group", ""),
                "mechanism_class": first.get("mechanism_class", ""),
                "library_count": len(rows),
                "n_cells": n_cells,
                "n_window_peak_features_total": sum(int(safe_float(row.get("n_window_peak_features"))) for row in rows if math.isfinite(safe_float(row.get("n_window_peak_features")))),
                "total_peak_counts": total_peak_counts,
                "target_peak_counts": target_peak_counts,
                "target_peak_cpm": cpm,
                "target_log2_peak_cpm": math.log2(cpm + 1.0) if math.isfinite(cpm) else math.nan,
                "target_peak_nonzero_entries": sum(int(safe_float(row.get("target_peak_nonzero_entries"))) for row in rows if math.isfinite(safe_float(row.get("target_peak_nonzero_entries")))),
            }
        )
    return sample_rows


def summarize_by_keys(rows: list[dict[str, object]], keys: list[str]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, ...], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[tuple(str(row.get(field, "")) for field in keys)].append(row)
    summaries: list[dict[str, object]] = []
    for key, group_rows in sorted(grouped.items()):
        log_values = [safe_float(row.get("target_log2_peak_cpm")) for row in group_rows]
        cpm_values = [safe_float(row.get("target_peak_cpm")) for row in group_rows]
        summary = dict(zip(keys, key, strict=False))
        summary.update(
            {
                "n": len(group_rows),
                "mean_peak_cpm": mean([v for v in cpm_values if math.isfinite(v)]),
                "median_peak_cpm": median([v for v in cpm_values if math.isfinite(v)]),
                "mean_log2_peak_cpm": mean([v for v in log_values if math.isfinite(v)]),
                "median_log2_peak_cpm": median([v for v in log_values if math.isfinite(v)]),
            }
        )
        summaries.append(summary)
    return summaries


def summarize_effects(rows: list[dict[str, object]], keys: list[str]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, ...], dict[str, list[dict[str, object]]]] = defaultdict(lambda: {"DS": [], "Control": []})
    for row in rows:
        key = tuple(str(row.get(field, "")) for field in keys)
        group = str(row.get("group", ""))
        if group in {"DS", "Control"}:
            grouped[key][group].append(row)

    summaries: list[dict[str, object]] = []
    pvalues: list[float] = []
    for key, values in sorted(grouped.items()):
        if not values["DS"] or not values["Control"]:
            continue
        ds_log = [safe_float(row.get("target_log2_peak_cpm")) for row in values["DS"]]
        control_log = [safe_float(row.get("target_log2_peak_cpm")) for row in values["Control"]]
        ds_cpm = [safe_float(row.get("target_peak_cpm")) for row in values["DS"]]
        control_cpm = [safe_float(row.get("target_peak_cpm")) for row in values["Control"]]
        first = (values["DS"] + values["Control"])[0]
        pvalue = safe_ttest(ds_log, control_log)
        pvalues.append(pvalue)
        summary = dict(zip(keys, key, strict=False))
        summary.update(
            {
                "dataset": "GSE305146",
                "target_class": first.get("target_class", ""),
                "tier": first.get("tier", ""),
                "candidate_group": first.get("candidate_group", ""),
                "mechanism_class": first.get("mechanism_class", ""),
                "n_ds": len(values["DS"]),
                "n_control": len(values["Control"]),
                "mean_peak_cpm_ds": mean([v for v in ds_cpm if math.isfinite(v)]),
                "mean_peak_cpm_control": mean([v for v in control_cpm if math.isfinite(v)]),
                "mean_log2_peak_cpm_ds": mean([v for v in ds_log if math.isfinite(v)]),
                "mean_log2_peak_cpm_control": mean([v for v in control_log if math.isfinite(v)]),
                "delta_log2_peak_cpm_ds_minus_control": mean([v for v in ds_log if math.isfinite(v)]) - mean([v for v in control_log if math.isfinite(v)]),
                "welch_pvalue": pvalue,
            }
        )
        summaries.append(summary)

    for summary, fdr in zip(summaries, bh_fdr(pvalues), strict=False):
        summary["bh_fdr"] = fdr
    return summaries


def write_note(path: Path, effect_rows: list[dict[str, object]], pcw_effect_rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    core = [row for row in effect_rows if str(row.get("target_class", "")).startswith("mechanism_core")]
    primary = [row for row in effect_rows if str(row.get("tier", "")).startswith("1_")]
    ranked = sorted(effect_rows, key=lambda row: abs(safe_float(row.get("delta_log2_peak_cpm_ds_minus_control"))), reverse=True)
    with path.open("w") as handle:
        handle.write("# GSE305146 scATAC Accessibility Proxy First Pass\n\n")
        handle.write("## Question\n\n")
        handle.write("Do fetal DS multiome samples show altered chromatin accessibility around DNMT3L/DNMT3A/DNMT3B, STAT1/STAT3, or final candidate genes?\n\n")
        handle.write("## Method\n\n")
        handle.write("This is a pseudo-bulk accessibility proxy, not chromVAR. For each final-analysis GSE305146 multiome library, peak counts were summed across peaks overlapping each target gene body plus 100 kb on either side. Counts were normalized by total peak counts per library, duplicate libraries were aggregated to biological samples, and DS-control differences were tested on log2 peak CPM values.\n\n")
        handle.write("## Core Genes\n\n")
        for row in sorted(core, key=lambda item: item.get("gene_symbol", "")):
            handle.write(
                f"- {row['gene_symbol']}: delta log2 peak CPM DS-Control {fmt(row['delta_log2_peak_cpm_ds_minus_control'])}; FDR {fmt(row['bh_fdr'])}; n DS {row['n_ds']}; n control {row['n_control']}.\n"
            )
        handle.write("\n## Primary Candidate Genes\n\n")
        for row in sorted(primary, key=lambda item: item.get("gene_symbol", "")):
            handle.write(
                f"- {row['gene_symbol']}: delta log2 peak CPM DS-Control {fmt(row['delta_log2_peak_cpm_ds_minus_control'])}; FDR {fmt(row['bh_fdr'])}.\n"
            )
        handle.write("\n## Largest Absolute Accessibility Shifts\n\n")
        for row in ranked[:10]:
            handle.write(
                f"- {row['gene_symbol']}: delta {fmt(row['delta_log2_peak_cpm_ds_minus_control'])}; FDR {fmt(row['bh_fdr'])}; class {row.get('target_class', '')}.\n"
            )
        handle.write("\n## Stage-Specific Caution\n\n")
        handle.write("Per-PCW effects are also written, but several stages have small DS/control sample counts. Treat those rows as exploratory direction checks rather than formal stage-resolved inference.\n")
        handle.write("\n## Interpretation\n\n")
        handle.write("This pass tests whether candidate loci are accessibility-shifted in DS fetal cortex, which is closer to a TF/regulatory mechanism than to direct methylation measurement. A true STAT/TF activity analysis still requires motif-level scoring such as chromVAR or a full differential-peak plus motif-enrichment workflow.\n")


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    library_rows, readiness_rows, window_rows = build_library_rows(args)
    sample_rows = aggregate_biological_samples(library_rows)
    group_summary = summarize_by_keys(sample_rows, ["dataset", "gene_symbol", "group"])
    effect_summary = summarize_effects(sample_rows, ["gene_symbol"])
    pcw_effect_summary = summarize_effects(sample_rows, ["gene_symbol", "dev_stage_pcw"])

    common_fields = [
        "dataset",
        "gene_symbol",
        "target_class",
        "tier",
        "candidate_group",
        "mechanism_class",
    ]
    write_tsv(
        outdir / "gse305146_accessibility_proxy_library_values.tsv",
        library_rows,
        common_fields
        + [
            "library_prefix",
            "geo_sample_accession",
            "biological_sample_id",
            "sample_name",
            "group",
            "genotype",
            "dev_stage_pcw",
            "sex",
            "platform_refs",
            "ensembl_id",
            "chrom",
            "gene_start",
            "gene_end",
            "window_start",
            "window_end",
            "window_bp",
            "n_window_peak_features",
            "n_cells",
            "total_peak_counts",
            "target_peak_counts",
            "target_peak_cpm",
            "target_log2_peak_cpm",
            "target_peak_nonzero_entries",
        ],
    )
    write_tsv(
        outdir / "gse305146_accessibility_proxy_biological_sample_values.tsv",
        sample_rows,
        common_fields
        + [
            "biological_sample_id",
            "group",
            "genotype",
            "dev_stage_pcw",
            "sex",
            "library_count",
            "n_cells",
            "n_window_peak_features_total",
            "total_peak_counts",
            "target_peak_counts",
            "target_peak_cpm",
            "target_log2_peak_cpm",
            "target_peak_nonzero_entries",
        ],
    )
    write_tsv(
        outdir / "gse305146_accessibility_proxy_group_summary.tsv",
        group_summary,
        ["dataset", "gene_symbol", "group", "n", "mean_peak_cpm", "median_peak_cpm", "mean_log2_peak_cpm", "median_log2_peak_cpm"],
    )
    write_tsv(
        outdir / "gse305146_accessibility_proxy_effect_summary.tsv",
        effect_summary,
        [
            "dataset",
            "gene_symbol",
            "target_class",
            "tier",
            "candidate_group",
            "mechanism_class",
            "n_ds",
            "n_control",
            "mean_peak_cpm_ds",
            "mean_peak_cpm_control",
            "mean_log2_peak_cpm_ds",
            "mean_log2_peak_cpm_control",
            "delta_log2_peak_cpm_ds_minus_control",
            "welch_pvalue",
            "bh_fdr",
        ],
    )
    write_tsv(
        outdir / "gse305146_accessibility_proxy_pcw_effect_summary.tsv",
        pcw_effect_summary,
        [
            "dataset",
            "gene_symbol",
            "dev_stage_pcw",
            "target_class",
            "tier",
            "candidate_group",
            "mechanism_class",
            "n_ds",
            "n_control",
            "mean_peak_cpm_ds",
            "mean_peak_cpm_control",
            "mean_log2_peak_cpm_ds",
            "mean_log2_peak_cpm_control",
            "delta_log2_peak_cpm_ds_minus_control",
            "welch_pvalue",
            "bh_fdr",
        ],
    )
    write_tsv(
        outdir / "gse305146_accessibility_proxy_readiness.tsv",
        readiness_rows,
        [
            "library_prefix",
            "biological_sample_id",
            "group",
            "dev_stage_pcw",
            "n_features",
            "n_cells",
            "n_entries",
            "last_gene_expression_row",
            "n_peak_features",
            "target_genes_with_coordinates",
            "target_peak_rows",
            "total_peak_counts",
        ],
    )
    write_tsv(
        outdir / "gse305146_accessibility_proxy_gene_windows.tsv",
        window_rows,
        [
            "gene_symbol",
            "target_class",
            "tier",
            "candidate_group",
            "mechanism_class",
            "ensembl_id",
            "chrom",
            "gene_start",
            "gene_end",
            "window_start",
            "window_end",
            "window_bp",
            "n_window_peak_features",
        ],
    )
    write_note(Path(args.note), effect_summary, pcw_effect_summary)
    print(f"Wrote GSE305146 accessibility proxy summaries to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
