#!/usr/bin/env python3
"""Add fetal-brain ATAC and histone-mark context to the final candidate matrix."""

from __future__ import annotations

import argparse
import csv
import gzip
from bisect import bisect_left
from collections import defaultdict
from glob import glob
from pathlib import Path


DEFAULT_MATRIX = "results/final_candidate_evidence/final_candidate_evidence_matrix.tsv"
DEFAULT_HUMAN_LINKS = "results/human_methylation_convergence/human_methylation_probe_gene_links.tsv"
DEFAULT_MOUSE_LINKS = "results/integrated_candidate_genes/mouse_triple_hit_dmr_gene_links.tsv"
DEFAULT_OUTDIR = "results/final_candidate_evidence"
DEFAULT_NOTE = "notes/final_candidate_regulatory_context.md"


HUMAN_ANNOTATIONS = {
    "human_fetal_ocr": ["data/geo/GSE149268/supplementary/GSE149268_annotation-ocr-hg38.bed.gz"],
    "human_fetal_pre": ["data/geo/GSE149268/supplementary/GSE149268_annotation-pre-hg38.bed.gz"],
    "human_fetal_atac_peak": ["data/geo/GSE149268/supplementary/GSM*_atac-*-optimal_idr-hg38.narrowPeak.gz"],
    "human_fetal_h3k27ac": ["data/geo/GSE149268/supplementary/GSM*_chip-*_H3K27ac-hg38.bed.gz"],
    "human_fetal_h3k4me1": ["data/geo/GSE149268/supplementary/GSM*_chip-*_H3K4me1-hg38.bed.gz"],
    "human_fetal_h3k27me3": ["data/geo/GSE149268/supplementary/GSM*_chip-*_H3K27me3-hg38.bed.gz"],
}


MOUSE_ANNOTATIONS = {
    "mouse_e14_5_forebrain_atac": ["data/geo/GSE172703/supplementary/*peaks_mm10.bed.gz"],
    "mouse_embryonic_forebrain_h3k27ac": ["data/geo/GSE201400/supplementary/GSE201400_*_H3K27AC.bed.gz"],
    "mouse_embryonic_forebrain_h3k27me3": ["data/geo/GSE201488/supplementary/GSE201488_*_H3K27ME3.bed.gz"],
    "mouse_embryonic_forebrain_h3k4me3": ["data/geo/GSE201487/supplementary/GSE201487_*_H3K4ME3.bed.gz"],
}


ADDED_FIELDS = [
    "human_linked_dnmt3l_sensitive_probes",
    "human_strict_fetal_positive_probe_ids",
    "human_strict_fetal_positive_probe_count",
    "human_any_regulatory_probe_count",
    "human_strict_regulatory_probe_count",
    "human_fetal_ocr_probe_count",
    "human_fetal_pre_probe_count",
    "human_fetal_atac_peak_probe_count",
    "human_fetal_h3k27ac_probe_count",
    "human_fetal_h3k4me1_probe_count",
    "human_fetal_h3k27me3_probe_count",
    "human_fetal_ocr_strict_probe_count",
    "human_fetal_atac_peak_strict_probe_count",
    "human_fetal_h3k27ac_strict_probe_count",
    "human_fetal_h3k27me3_strict_probe_count",
    "mouse_linked_triple_hit_dmr_ids",
    "mouse_linked_triple_hit_dmr_count",
    "mouse_any_regulatory_dmr_count",
    "mouse_e14_5_forebrain_atac_dmr_count",
    "mouse_embryonic_forebrain_h3k27ac_dmr_count",
    "mouse_embryonic_forebrain_h3k27me3_dmr_count",
    "mouse_embryonic_forebrain_h3k4me3_dmr_count",
    "regulatory_context_score",
    "regulatory_context_class",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", default=DEFAULT_MATRIX)
    parser.add_argument("--human-links", default=DEFAULT_HUMAN_LINKS)
    parser.add_argument("--mouse-links", default=DEFAULT_MOUSE_LINKS)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--note", default=DEFAULT_NOTE)
    return parser.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="", errors="replace")
    return path.open(newline="", errors="replace")


def clean_key(value: str | None) -> str:
    return (value or "").strip().replace("\r", "")


def is_yes(value: str | None) -> bool:
    return (value or "").strip().lower() == "yes"


def read_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = [clean_key(field) for field in reader.fieldnames or []]
        rows = [{clean_key(k): (v or "").strip() for k, v in row.items()} for row in reader]
    return rows, fields


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_bed_line(line: str) -> tuple[str, int, int] | None:
    if not line.strip() or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
        return None
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 3:
        return None
    try:
        start = int(float(parts[1]))
        end = int(float(parts[2]))
    except ValueError:
        return None
    if end <= start:
        return None
    return parts[0], start, end


def paths_from_globs(patterns: list[str]) -> list[Path]:
    paths: list[Path] = []
    for pattern in patterns:
        paths.extend(Path(path) for path in sorted(glob(pattern)))
    return paths


def load_annotation_group(patterns: list[str]) -> dict[str, tuple[list[int], list[int]]]:
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for path in paths_from_globs(patterns):
        with open_text(path) as handle:
            for line in handle:
                parsed = parse_bed_line(line)
                if parsed is None:
                    continue
                chrom, start, end = parsed
                intervals[chrom].append((start, end))

    indexed: dict[str, tuple[list[int], list[int]]] = {}
    for chrom, values in intervals.items():
        values.sort()
        starts = [start for start, _end in values]
        prefix_max_ends: list[int] = []
        current = 0
        for _start, end in values:
            current = max(current, end)
            prefix_max_ends.append(current)
        indexed[chrom] = (starts, prefix_max_ends)
    return indexed


def has_overlap(index: dict[str, tuple[list[int], list[int]]], chrom: str, start: int, end: int) -> bool:
    starts, prefix_max_ends = index.get(chrom, ([], []))
    if not starts:
        return False
    pos = bisect_left(starts, end) - 1
    return pos >= 0 and prefix_max_ends[pos] > start


def load_annotation_indexes(groups: dict[str, list[str]]) -> dict[str, dict[str, tuple[list[int], list[int]]]]:
    return {name: load_annotation_group(patterns) for name, patterns in groups.items()}


def read_human_probe_links(path: Path, candidate_symbols: set[str]) -> dict[str, list[dict[str, object]]]:
    rows: dict[str, dict[str, dict[str, object]]] = defaultdict(dict)
    with path.open(newline="", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        reader.fieldnames = [clean_key(field) for field in reader.fieldnames or []]
        for raw in reader:
            row = {clean_key(k): (v or "").strip() for k, v in raw.items()}
            symbol = row.get("gene_symbol", "").upper()
            if symbol not in candidate_symbols:
                continue
            probe_id = row.get("probe_id", "")
            chrom = row.get("hg38_chrom", "")
            if not probe_id or not chrom:
                continue
            try:
                start = int(float(row["probe_start"]))
                end = int(float(row["probe_end"]))
            except (KeyError, ValueError):
                continue
            if end <= start:
                end = start + 1
            existing = rows[symbol].get(probe_id)
            strict = is_yes(row.get("human_triple_phase3_fetal_mean_positive"))
            any_ds = is_yes(row.get("human_triple_phase3_any_ds"))
            if existing is None:
                rows[symbol][probe_id] = {
                    "id": probe_id,
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strict": strict,
                    "any_ds": any_ds,
                }
            else:
                existing["strict"] = bool(existing["strict"] or strict)
                existing["any_ds"] = bool(existing["any_ds"] or any_ds)
    return {symbol: list(values.values()) for symbol, values in rows.items()}


def read_mouse_dmr_links(path: Path, candidate_symbols: set[str]) -> dict[str, list[dict[str, object]]]:
    rows: dict[str, dict[str, dict[str, object]]] = defaultdict(dict)
    with path.open(newline="", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        reader.fieldnames = [clean_key(field) for field in reader.fieldnames or []]
        for raw in reader:
            row = {clean_key(k): (v or "").strip() for k, v in raw.items()}
            symbol = row.get("gene_symbol", "").upper()
            if symbol not in candidate_symbols:
                continue
            dmr_id = row.get("dmr_id", "")
            chrom = row.get("chrom", "")
            if not dmr_id or not chrom:
                continue
            try:
                start = int(float(row["dmr_start"]))
                end = int(float(row["dmr_end"]))
            except (KeyError, ValueError):
                continue
            rows[symbol][dmr_id] = {"id": dmr_id, "chrom": chrom, "start": start, "end": end}
    return {symbol: list(values.values()) for symbol, values in rows.items()}


def ids_with_overlap(intervals: list[dict[str, object]], index: dict[str, tuple[list[int], list[int]]]) -> set[str]:
    overlapping: set[str] = set()
    for interval in intervals:
        if has_overlap(index, str(interval["chrom"]), int(interval["start"]), int(interval["end"])):
            overlapping.add(str(interval["id"]))
    return overlapping


def summarize_human(symbol: str, probes_by_gene: dict[str, list[dict[str, object]]], indexes: dict[str, dict[str, tuple[list[int], list[int]]]]) -> dict[str, str]:
    probes = probes_by_gene.get(symbol, [])
    strict_probes = [probe for probe in probes if probe.get("strict")]
    overlap_sets: dict[str, set[str]] = {
        name: ids_with_overlap(probes, index)
        for name, index in indexes.items()
    }
    strict_overlap_sets: dict[str, set[str]] = {
        name: ids_with_overlap(strict_probes, index)
        for name, index in indexes.items()
    }
    any_regulatory = set().union(*overlap_sets.values()) if overlap_sets else set()
    strict_regulatory = set().union(*strict_overlap_sets.values()) if strict_overlap_sets else set()
    return {
        "human_linked_dnmt3l_sensitive_probes": str(len(probes)),
        "human_strict_fetal_positive_probe_ids": ",".join(sorted(str(probe["id"]) for probe in strict_probes)),
        "human_strict_fetal_positive_probe_count": str(len(strict_probes)),
        "human_any_regulatory_probe_count": str(len(any_regulatory)),
        "human_strict_regulatory_probe_count": str(len(strict_regulatory)),
        "human_fetal_ocr_probe_count": str(len(overlap_sets["human_fetal_ocr"])),
        "human_fetal_pre_probe_count": str(len(overlap_sets["human_fetal_pre"])),
        "human_fetal_atac_peak_probe_count": str(len(overlap_sets["human_fetal_atac_peak"])),
        "human_fetal_h3k27ac_probe_count": str(len(overlap_sets["human_fetal_h3k27ac"])),
        "human_fetal_h3k4me1_probe_count": str(len(overlap_sets["human_fetal_h3k4me1"])),
        "human_fetal_h3k27me3_probe_count": str(len(overlap_sets["human_fetal_h3k27me3"])),
        "human_fetal_ocr_strict_probe_count": str(len(strict_overlap_sets["human_fetal_ocr"])),
        "human_fetal_atac_peak_strict_probe_count": str(len(strict_overlap_sets["human_fetal_atac_peak"])),
        "human_fetal_h3k27ac_strict_probe_count": str(len(strict_overlap_sets["human_fetal_h3k27ac"])),
        "human_fetal_h3k27me3_strict_probe_count": str(len(strict_overlap_sets["human_fetal_h3k27me3"])),
    }


def summarize_mouse(symbol: str, dmrs_by_gene: dict[str, list[dict[str, object]]], indexes: dict[str, dict[str, tuple[list[int], list[int]]]]) -> dict[str, str]:
    dmrs = dmrs_by_gene.get(symbol, [])
    overlap_sets: dict[str, set[str]] = {
        name: ids_with_overlap(dmrs, index)
        for name, index in indexes.items()
    }
    any_regulatory = set().union(*overlap_sets.values()) if overlap_sets else set()
    return {
        "mouse_linked_triple_hit_dmr_ids": ",".join(sorted(str(dmr["id"]) for dmr in dmrs)),
        "mouse_linked_triple_hit_dmr_count": str(len(dmrs)),
        "mouse_any_regulatory_dmr_count": str(len(any_regulatory)),
        "mouse_e14_5_forebrain_atac_dmr_count": str(len(overlap_sets["mouse_e14_5_forebrain_atac"])),
        "mouse_embryonic_forebrain_h3k27ac_dmr_count": str(len(overlap_sets["mouse_embryonic_forebrain_h3k27ac"])),
        "mouse_embryonic_forebrain_h3k27me3_dmr_count": str(len(overlap_sets["mouse_embryonic_forebrain_h3k27me3"])),
        "mouse_embryonic_forebrain_h3k4me3_dmr_count": str(len(overlap_sets["mouse_embryonic_forebrain_h3k4me3"])),
    }


def regulatory_class(row: dict[str, str]) -> tuple[str, str]:
    human_strict_reg = int(row["human_strict_regulatory_probe_count"])
    human_any_reg = int(row["human_any_regulatory_probe_count"])
    mouse_reg = int(row["mouse_any_regulatory_dmr_count"])
    mouse_total = int(row["mouse_linked_triple_hit_dmr_count"])
    mouse_atac = int(row["mouse_e14_5_forebrain_atac_dmr_count"])
    mouse_h3k27ac = int(row["mouse_embryonic_forebrain_h3k27ac_dmr_count"])
    mouse_h3k27me3 = int(row["mouse_embryonic_forebrain_h3k27me3_dmr_count"])

    score = 0
    if human_strict_reg > 0:
        score += 3
    elif human_any_reg > 0:
        score += 1
    if mouse_reg > 0:
        score += 2
    if mouse_total and mouse_atac / mouse_total >= 0.5:
        score += 1
    if mouse_h3k27ac > 0:
        score += 1
    if mouse_h3k27me3 > 0:
        score += 1

    if human_strict_reg > 0 and mouse_reg > 0:
        label = "human_mouse_regulatory_supported"
    elif human_any_reg > 0 and mouse_reg > 0:
        label = "human_mouse_regulatory_possible"
    elif mouse_reg > 0:
        label = "mouse_regulatory_supported"
    elif human_any_reg > 0:
        label = "human_regulatory_supported"
    else:
        label = "no_regulatory_overlap_detected"
    return str(score), label


def write_note(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    class_counts: dict[str, int] = defaultdict(int)
    class_genes: dict[str, list[str]] = defaultdict(list)
    for row in rows:
        label = row["regulatory_context_class"]
        class_counts[label] += 1
        class_genes[label].append(row["gene_symbol"])

    top_rows = sorted(rows, key=lambda row: (-int(row["regulatory_context_score"]), row["tier"], row["gene_symbol"]))
    lines = [
        "# Final candidate regulatory context",
        "",
        "This pass adds candidate-level fetal brain ATAC and histone-mark context to the final 24-gene evidence matrix. It uses human GSE149268 hg38 fetal telencephalon OCR/pRE/ATAC/histone annotations and mouse GSE172703/GSE201400/GSE201488/GSE201487 mm10 embryonic forebrain annotations.",
        "",
        "Important interpretation: these overlaps support regulatory context, not causal methylation mediation by themselves.",
        "",
        "## Regulatory class summary",
        "",
        "| Regulatory class | n | Genes |",
        "|---|---:|---|",
    ]
    for label, count in sorted(class_counts.items()):
        lines.append(f"| {label} | {count} | {', '.join(class_genes[label])} |")
    lines.extend(
        [
            "",
            "## Highest regulatory-context scores",
            "",
            "| Gene | Tier | Class | Score | Human strict regulatory probes | Mouse regulatory DMRs | Mouse ATAC DMRs | Mouse H3K27ac DMRs | Mouse H3K27me3 DMRs |",
            "|---|---|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in top_rows[:12]:
        lines.append(
            "| {gene_symbol} | {tier} | {regulatory_context_class} | {regulatory_context_score} | "
            "{human_strict_regulatory_probe_count} | {mouse_any_regulatory_dmr_count} | "
            "{mouse_e14_5_forebrain_atac_dmr_count} | {mouse_embryonic_forebrain_h3k27ac_dmr_count} | "
            "{mouse_embryonic_forebrain_h3k27me3_dmr_count} |".format(**row)
        )
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            "- `results/final_candidate_evidence/final_candidate_evidence_matrix_with_regulatory_context.tsv`",
            "- `results/final_candidate_evidence/final_candidate_regulatory_context_summary.tsv`",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    rows, matrix_fields = read_tsv(Path(args.matrix))
    candidate_symbols = {row["gene_symbol"].upper() for row in rows}

    human_probe_links = read_human_probe_links(Path(args.human_links), candidate_symbols)
    mouse_dmr_links = read_mouse_dmr_links(Path(args.mouse_links), {row.get("mouse_symbol", "").upper() for row in rows})

    human_indexes = load_annotation_indexes(HUMAN_ANNOTATIONS)
    mouse_indexes = load_annotation_indexes(MOUSE_ANNOTATIONS)

    augmented_rows: list[dict[str, str]] = []
    for row in rows:
        symbol = row["gene_symbol"].upper()
        mouse_symbol = row.get("mouse_symbol", "").upper()
        augmented = {
            **row,
            **summarize_human(symbol, human_probe_links, human_indexes),
            **summarize_mouse(mouse_symbol, mouse_dmr_links, mouse_indexes),
        }
        score, label = regulatory_class(augmented)
        augmented["regulatory_context_score"] = score
        augmented["regulatory_context_class"] = label
        augmented_rows.append(augmented)

    outdir = Path(args.outdir)
    write_tsv(
        outdir / "final_candidate_evidence_matrix_with_regulatory_context.tsv",
        augmented_rows,
        matrix_fields + ADDED_FIELDS,
    )

    summary_rows = []
    for row in augmented_rows:
        summary_rows.append(
            {
                "tier": row["tier"],
                "gene_symbol": row["gene_symbol"],
                "evidence_class": row["evidence_class"],
                "regulatory_context_class": row["regulatory_context_class"],
                "regulatory_context_score": row["regulatory_context_score"],
                "human_strict_regulatory_probe_count": row["human_strict_regulatory_probe_count"],
                "human_any_regulatory_probe_count": row["human_any_regulatory_probe_count"],
                "mouse_any_regulatory_dmr_count": row["mouse_any_regulatory_dmr_count"],
                "mouse_e14_5_forebrain_atac_dmr_count": row["mouse_e14_5_forebrain_atac_dmr_count"],
                "mouse_embryonic_forebrain_h3k27ac_dmr_count": row["mouse_embryonic_forebrain_h3k27ac_dmr_count"],
                "mouse_embryonic_forebrain_h3k27me3_dmr_count": row["mouse_embryonic_forebrain_h3k27me3_dmr_count"],
            }
        )
    summary_rows.sort(key=lambda row: (-int(row["regulatory_context_score"]), row["tier"], row["gene_symbol"]))
    write_tsv(
        outdir / "final_candidate_regulatory_context_summary.tsv",
        summary_rows,
        [
            "tier",
            "gene_symbol",
            "evidence_class",
            "regulatory_context_class",
            "regulatory_context_score",
            "human_strict_regulatory_probe_count",
            "human_any_regulatory_probe_count",
            "mouse_any_regulatory_dmr_count",
            "mouse_e14_5_forebrain_atac_dmr_count",
            "mouse_embryonic_forebrain_h3k27ac_dmr_count",
            "mouse_embryonic_forebrain_h3k27me3_dmr_count",
        ],
    )
    write_note(Path(args.note), augmented_rows)

    print(f"Wrote {len(augmented_rows)} augmented candidate rows")
    print(f"Wrote {outdir / 'final_candidate_evidence_matrix_with_regulatory_context.tsv'}")
    print(f"Wrote {outdir / 'final_candidate_regulatory_context_summary.tsv'}")
    print(f"Wrote {args.note}")


if __name__ == "__main__":
    main()
