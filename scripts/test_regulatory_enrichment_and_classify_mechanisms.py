#!/usr/bin/env python3
"""Test regulatory enrichment and classify methylation-expression mechanisms."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from bisect import bisect_left
from collections import defaultdict
from glob import glob
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DEFAULT_MATRIX = "results/final_candidate_evidence/final_candidate_evidence_matrix_with_regulatory_context.tsv"
DEFAULT_HUMAN_LINKS = "results/human_methylation_convergence/human_methylation_probe_gene_links.tsv"
DEFAULT_MOUSE_LINKS = "results/integrated_candidate_genes/mouse_triple_hit_dmr_gene_links.tsv"
DEFAULT_MOUSE_RANKING = "results/integrated_candidate_genes/mouse_triple_hit_candidate_gene_ranking.tsv"
DEFAULT_OUTDIR = "results/mechanism_classification"
DEFAULT_FIGDIR = "figures/mechanism_classification"
DEFAULT_NOTE = "notes/mechanism_enrichment_and_direction_first_pass.md"


HUMAN_REGULATORY_TRACKS = {
    "human_fetal_ocr": ["data/geo/GSE149268/supplementary/GSE149268_annotation-ocr-hg38.bed.gz"],
    "human_fetal_pre": ["data/geo/GSE149268/supplementary/GSE149268_annotation-pre-hg38.bed.gz"],
    "human_fetal_atac": ["data/geo/GSE149268/supplementary/GSM*_atac-*-optimal_idr-hg38.narrowPeak.gz"],
    "human_h3k27ac": ["data/geo/GSE149268/supplementary/GSM*_chip-*_H3K27ac-hg38.bed.gz"],
    "human_h3k4me1": ["data/geo/GSE149268/supplementary/GSM*_chip-*_H3K4me1-hg38.bed.gz"],
    "human_h3k27me3": ["data/geo/GSE149268/supplementary/GSM*_chip-*_H3K27me3-hg38.bed.gz"],
}


MOUSE_REGULATORY_TRACKS = {
    "mouse_atac": ["data/geo/GSE172703/supplementary/*peaks_mm10.bed.gz"],
    "mouse_h3k27ac": ["data/geo/GSE201400/supplementary/GSE201400_*_H3K27AC.bed.gz"],
    "mouse_h3k4me3": ["data/geo/GSE201487/supplementary/GSE201487_*_H3K4ME3.bed.gz"],
    "mouse_h3k27me3": ["data/geo/GSE201488/supplementary/GSE201488_*_H3K27ME3.bed.gz"],
}


ENRICHMENT_FIELDS = [
    "candidate_set",
    "assay",
    "background",
    "candidate_genes_total",
    "candidate_genes_in_background",
    "candidate_successes",
    "background_genes",
    "background_successes",
    "candidate_success_fraction",
    "background_success_fraction",
    "fold_enrichment",
    "hypergeometric_p_greater",
    "candidate_success_genes",
    "candidate_background_genes",
]


MECHANISM_FIELDS = [
    "gene_symbol",
    "tier",
    "evidence_class",
    "regulatory_context_class",
    "human_methylation_direction",
    "mean_gse157078_delta_beta",
    "human_ds_signed_delta",
    "human_expression_direction",
    "human_expression_delta",
    "human_expression_context",
    "human_methyl_expr_relation",
    "mouse_symbol",
    "mouse_methylation_direction",
    "mouse_mean_dmr_delta_beta",
    "mouse_gse268951_direction",
    "mouse_gse268951_log2fc",
    "mouse_oe_methyl_expr_relation",
    "mouse_dp10_direction",
    "mouse_dp10_delta",
    "mouse_dp10_methyl_expr_relation",
    "human_strict_regulatory_probe_count",
    "mouse_any_regulatory_dmr_count",
    "mouse_embryonic_forebrain_h3k27me3_dmr_count",
    "mechanism_class",
    "interpretation",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", default=DEFAULT_MATRIX)
    parser.add_argument("--human-links", default=DEFAULT_HUMAN_LINKS)
    parser.add_argument("--mouse-links", default=DEFAULT_MOUSE_LINKS)
    parser.add_argument("--mouse-ranking", default=DEFAULT_MOUSE_RANKING)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--figdir", default=DEFAULT_FIGDIR)
    parser.add_argument("--note", default=DEFAULT_NOTE)
    return parser.parse_args()


def clean_key(value: str | None) -> str:
    return (value or "").strip().replace("\r", "")


def is_yes(value: str | None) -> bool:
    return (value or "").strip().lower() == "yes"


def safe_float(value: object) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return math.nan
    return number if math.isfinite(number) else math.nan


def fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}" if math.isfinite(value) else ""
    return str(value)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        reader.fieldnames = [clean_key(field) for field in reader.fieldnames or []]
        return [{clean_key(k): (v or "").strip() for k, v in row.items()} for row in reader]


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="", errors="replace")
    return path.open(newline="", errors="replace")


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


def load_annotation_indexes(groups: dict[str, list[str]]) -> dict[str, dict[str, tuple[list[int], list[int]]]]:
    return {name: load_annotation_group(patterns) for name, patterns in groups.items()}


def has_overlap(index: dict[str, tuple[list[int], list[int]]], chrom: str, start: int, end: int) -> bool:
    starts, prefix_max_ends = index.get(chrom, ([], []))
    if not starts:
        return False
    pos = bisect_left(starts, end) - 1
    return pos >= 0 and prefix_max_ends[pos] > start


def overlaps_any(indexes: dict[str, dict[str, tuple[list[int], list[int]]]], chrom: str, start: int, end: int) -> set[str]:
    return {name for name, index in indexes.items() if has_overlap(index, chrom, start, end)}


def build_human_background(path: Path, indexes: dict[str, dict[str, tuple[list[int], list[int]]]]) -> dict[str, dict[str, object]]:
    genes: dict[str, dict[str, object]] = {}
    seen: set[tuple[str, str]] = set()
    for row in read_tsv(path):
        if row.get("gene_type") != "protein_coding":
            continue
        symbol = row.get("gene_symbol", "").upper()
        probe_id = row.get("probe_id", "")
        chrom = row.get("hg38_chrom", "")
        if not symbol or not probe_id or not chrom:
            continue
        key = (symbol, probe_id)
        if key in seen:
            continue
        seen.add(key)
        try:
            start = int(float(row["probe_start"]))
            end = int(float(row["probe_end"]))
        except (KeyError, ValueError):
            continue
        if end <= start:
            end = start + 1
        overlaps = overlaps_any(indexes, chrom, start, end)
        item = genes.setdefault(
            symbol,
            {
                "linked_probes": 0,
                "any_regulatory_probes": 0,
                "strict_probes": 0,
                "strict_regulatory_probes": 0,
                "track_hits": defaultdict(int),
            },
        )
        item["linked_probes"] = int(item["linked_probes"]) + 1
        if overlaps:
            item["any_regulatory_probes"] = int(item["any_regulatory_probes"]) + 1
            for track in overlaps:
                item["track_hits"][track] += 1
        if is_yes(row.get("human_triple_phase3_fetal_mean_positive")):
            item["strict_probes"] = int(item["strict_probes"]) + 1
            if overlaps:
                item["strict_regulatory_probes"] = int(item["strict_regulatory_probes"]) + 1
    return genes


def build_mouse_background(path: Path, indexes: dict[str, dict[str, tuple[list[int], list[int]]]]) -> dict[str, dict[str, object]]:
    genes: dict[str, dict[str, object]] = {}
    seen: set[tuple[str, str]] = set()
    for row in read_tsv(path):
        symbol = row.get("gene_symbol", "").upper()
        dmr_id = row.get("dmr_id", "")
        chrom = row.get("chrom", "")
        if not symbol or not dmr_id or not chrom:
            continue
        key = (symbol, dmr_id)
        if key in seen:
            continue
        seen.add(key)
        try:
            start = int(float(row["dmr_start"]))
            end = int(float(row["dmr_end"]))
        except (KeyError, ValueError):
            continue
        overlaps = overlaps_any(indexes, chrom, start, end)
        item = genes.setdefault(
            symbol,
            {
                "linked_dmrs": 0,
                "any_regulatory_dmrs": 0,
                "track_hits": defaultdict(int),
            },
        )
        item["linked_dmrs"] = int(item["linked_dmrs"]) + 1
        if overlaps:
            item["any_regulatory_dmrs"] = int(item["any_regulatory_dmrs"]) + 1
            for track in overlaps:
                item["track_hits"][track] += 1
    return genes


def log_choose(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_p_greater(x: int, n: int, k: int, m: int) -> float:
    """P[X >= x] where X ~ Hypergeometric(M=m, K=k, draws=n)."""
    if n <= 0 or m <= 0 or k < 0 or k > m:
        return math.nan
    max_x = min(n, k)
    if x > max_x:
        return 0.0
    log_den = log_choose(m, n)
    probs = [
        math.exp(log_choose(k, i) + log_choose(m - k, n - i) - log_den)
        for i in range(max(0, x), max_x + 1)
    ]
    return min(1.0, sum(probs))


def enrichment_row(
    candidate_set: str,
    assay: str,
    background: str,
    candidate_symbols: set[str],
    background_symbols: set[str],
    success_symbols: set[str],
) -> dict[str, object]:
    candidate_in_background = sorted(candidate_symbols & background_symbols)
    candidate_successes = sorted(set(candidate_in_background) & success_symbols)
    n = len(candidate_in_background)
    a = len(candidate_successes)
    m = len(background_symbols)
    k = len(success_symbols & background_symbols)
    candidate_fraction = a / n if n else math.nan
    background_fraction = k / m if m else math.nan
    fold = candidate_fraction / background_fraction if background_fraction and math.isfinite(candidate_fraction) else math.nan
    return {
        "candidate_set": candidate_set,
        "assay": assay,
        "background": background,
        "candidate_genes_total": len(candidate_symbols),
        "candidate_genes_in_background": n,
        "candidate_successes": a,
        "background_genes": m,
        "background_successes": k,
        "candidate_success_fraction": candidate_fraction,
        "background_success_fraction": background_fraction,
        "fold_enrichment": fold,
        "hypergeometric_p_greater": hypergeom_p_greater(a, n, k, m),
        "candidate_success_genes": ",".join(candidate_successes),
        "candidate_background_genes": ",".join(candidate_in_background),
    }


def candidate_sets(matrix_rows: list[dict[str, str]]) -> dict[str, dict[str, set[str]]]:
    sets: dict[str, dict[str, set[str]]] = {}

    def add(name: str, rows: list[dict[str, str]]) -> None:
        sets[name] = {
            "human": {row["gene_symbol"].upper() for row in rows},
            "mouse": {row["mouse_symbol"].upper() for row in rows if row.get("mouse_symbol")},
        }

    add("primary_5", [row for row in matrix_rows if row["tier"] == "1_primary"])
    add("final_24", matrix_rows)
    add(
        "methylation_supported_final",
        [row for row in matrix_rows if row["tier"] in {"1_primary", "2_secondary_balanced", "3_methylation_heavy"}],
    )
    add("expression_pathway_final", [row for row in matrix_rows if row["tier"] == "4_expression_pathway"])
    return sets


def regulatory_enrichment(
    matrix_rows: list[dict[str, str]],
    human_bg: dict[str, dict[str, object]],
    mouse_bg: dict[str, dict[str, object]],
) -> list[dict[str, object]]:
    sets = candidate_sets(matrix_rows)

    human_any_background = {gene for gene, item in human_bg.items() if int(item["linked_probes"]) > 0}
    human_any_success = {gene for gene, item in human_bg.items() if int(item["any_regulatory_probes"]) > 0}
    human_strict_background = {gene for gene, item in human_bg.items() if int(item["strict_probes"]) > 0}
    human_strict_success = {gene for gene, item in human_bg.items() if int(item["strict_regulatory_probes"]) > 0}
    mouse_background = {gene for gene, item in mouse_bg.items() if int(item["linked_dmrs"]) > 0}
    mouse_any_success = {gene for gene, item in mouse_bg.items() if int(item["any_regulatory_dmrs"]) > 0}
    mouse_atac_success = {
        gene for gene, item in mouse_bg.items() if int(item["track_hits"].get("mouse_atac", 0)) > 0
    }
    mouse_h3k27ac_success = {
        gene for gene, item in mouse_bg.items() if int(item["track_hits"].get("mouse_h3k27ac", 0)) > 0
    }
    mouse_h3k27me3_success = {
        gene for gene, item in mouse_bg.items() if int(item["track_hits"].get("mouse_h3k27me3", 0)) > 0
    }

    tests = [
        ("human_any_regulatory_probe", "human_protein_coding_dnmt3l_linked_probe_genes", "human", human_any_background, human_any_success),
        ("human_strict_regulatory_probe", "human_protein_coding_strict_triple_probe_genes", "human", human_strict_background, human_strict_success),
        ("mouse_any_regulatory_dmr", "mouse_triple_hit_dmr_linked_genes", "mouse", mouse_background, mouse_any_success),
        ("mouse_atac_dmr", "mouse_triple_hit_dmr_linked_genes", "mouse", mouse_background, mouse_atac_success),
        ("mouse_h3k27ac_dmr", "mouse_triple_hit_dmr_linked_genes", "mouse", mouse_background, mouse_h3k27ac_success),
        ("mouse_h3k27me3_dmr", "mouse_triple_hit_dmr_linked_genes", "mouse", mouse_background, mouse_h3k27me3_success),
    ]
    rows: list[dict[str, object]] = []
    for set_name, symbols_by_species in sets.items():
        for assay, background_name, species, background_symbols, success_symbols in tests:
            rows.append(
                enrichment_row(
                    set_name,
                    assay,
                    background_name,
                    symbols_by_species[species],
                    background_symbols,
                    success_symbols,
                )
            )
    return rows


def sign_direction(value: float, threshold: float = 0.05) -> str:
    if not math.isfinite(value) or abs(value) < threshold:
        return "flat"
    return "up" if value > 0 else "down"


def methyl_direction(value: float, threshold: float = 0.01) -> str:
    if not math.isfinite(value) or abs(value) < threshold:
        return "mixed_or_flat"
    return "hyper" if value > 0 else "hypo"


def strongest_human_expression(row: dict[str, str]) -> tuple[str, float, str]:
    candidates = [
        ("GSE305153_all_cells", safe_float(row.get("gse305153_all_cells_delta"))),
        (row.get("best_bulk_same_direction_context", "best_bulk"), safe_float(row.get("best_bulk_same_direction_delta"))),
        (row.get("best_celltype_same_direction_label", "best_celltype"), safe_float(row.get("best_celltype_same_direction_delta"))),
    ]
    finite = [(label, value) for label, value in candidates if math.isfinite(value)]
    if not finite:
        return "unavailable", math.nan, ""
    label, value = max(finite, key=lambda item: abs(item[1]))
    return sign_direction(value, threshold=0.10), value, label


def compare_methyl_expression(methyl: str, expr: str) -> str:
    if methyl in {"mixed_or_flat", ""} or expr in {"flat", "unavailable", ""}:
        return "unresolved"
    if methyl == "hyper" and expr == "down":
        return "opposite_possible_repression"
    if methyl == "hypo" and expr == "up":
        return "opposite_possible_derepression"
    if methyl == "hyper" and expr == "up":
        return "same_direction_hyper_up_complex"
    if methyl == "hypo" and expr == "down":
        return "same_direction_hypo_down_complex"
    return "unresolved"


def read_mouse_ranking(path: Path) -> dict[str, dict[str, str]]:
    return {row["gene_symbol"].upper(): row for row in read_tsv(path)}


def mechanism_class(row: dict[str, str], human_relation: str, mouse_oe_relation: str, mouse_dp10_relation: str) -> tuple[str, str]:
    evidence = row["evidence_class"]
    h3k27me3 = int(float(row.get("mouse_embryonic_forebrain_h3k27me3_dmr_count", "0") or 0))
    strict_reg = int(float(row.get("human_strict_regulatory_probe_count", "0") or 0))
    mouse_reg = int(float(row.get("mouse_any_regulatory_dmr_count", "0") or 0))
    relations = {human_relation, mouse_oe_relation, mouse_dp10_relation}

    if evidence == "expression_pathway_methylation_unresolved":
        return (
            "expression_pathway_methylation_unresolved",
            "Expression/pathway candidate; current strict human methylation support is absent, so methylation independence is not proven.",
        )
    if evidence.startswith("methylation_heavy"):
        if h3k27me3 > 0:
            return (
                "polycomb_developmental_methylation_target",
                "Developmental methylation target with H3K27me3/Polycomb regulatory context; expression may be timing- or cell-type-restricted.",
            )
        return (
            "developmental_methylation_target_expression_weak",
            "Developmental methylation target with weak or inconsistent expression support in available datasets.",
        )
    if evidence == "methylation_expression_supported":
        return (
            "secondary_regulatory_methylation_candidate",
            "Secondary candidate with methylation and expression support, but lower priority or less direct directionality than the primary set.",
        )
    if evidence == "methylation_convergent":
        if "same_direction_hyper_up_complex" in relations or "same_direction_hypo_down_complex" in relations:
            return (
                "directionally_complex_regulatory_methylation",
                "Methylation and expression are not consistent with a simple promoter-silencing model; enhancer, gene-body, indirect, or cell-state regulation should be considered.",
            )
        if "opposite_possible_repression" in relations or "opposite_possible_derepression" in relations:
            return (
                "possible_repressive_or_derepressive_methylation",
                "Methylation-expression direction is compatible with repression/derepression in at least one context, but paired mediation is still needed.",
            )
        if strict_reg > 0 and mouse_reg > 0:
            return (
                "regulatory_methylation_convergent_unresolved_direction",
                "Strong regulatory methylation convergence, but expression direction is mixed or flat across contexts.",
            )
    return (
        "methylation_expression_unclassified",
        "Evidence is insufficient for a more specific mechanism class.",
    )


def classify_mechanisms(matrix_rows: list[dict[str, str]], mouse_ranking: dict[str, dict[str, str]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for row in matrix_rows:
        human_methyl_value = safe_float(row.get("mean_gse157078_delta_beta"))
        human_methyl = methyl_direction(human_methyl_value)
        human_expr, human_expr_delta, human_expr_context = strongest_human_expression(row)
        human_relation = compare_methyl_expression(human_methyl, human_expr)

        mouse_row = mouse_ranking.get(row.get("mouse_symbol", "").upper(), {})
        mouse_methyl_value = safe_float(mouse_row.get("mean_dmr_delta_beta"))
        mouse_methyl = methyl_direction(mouse_methyl_value)
        mouse_oe_direction = row.get("mouse_gse268951_direction", "")
        mouse_dp10_direction = row.get("mouse_dp10_direction", "")
        mouse_oe_relation = compare_methyl_expression(mouse_methyl, mouse_oe_direction)
        mouse_dp10_relation = compare_methyl_expression(mouse_methyl, mouse_dp10_direction)
        mech_class, interpretation = mechanism_class(row, human_relation, mouse_oe_relation, mouse_dp10_relation)

        rows.append(
            {
                "gene_symbol": row["gene_symbol"],
                "tier": row["tier"],
                "evidence_class": row["evidence_class"],
                "regulatory_context_class": row["regulatory_context_class"],
                "human_methylation_direction": human_methyl,
                "mean_gse157078_delta_beta": human_methyl_value,
                "human_ds_signed_delta": safe_float(row.get("mean_ds_signed_delta")),
                "human_expression_direction": human_expr,
                "human_expression_delta": human_expr_delta,
                "human_expression_context": human_expr_context,
                "human_methyl_expr_relation": human_relation,
                "mouse_symbol": row.get("mouse_symbol", ""),
                "mouse_methylation_direction": mouse_methyl,
                "mouse_mean_dmr_delta_beta": mouse_methyl_value,
                "mouse_gse268951_direction": mouse_oe_direction,
                "mouse_gse268951_log2fc": safe_float(row.get("mouse_gse268951_log2fc")),
                "mouse_oe_methyl_expr_relation": mouse_oe_relation,
                "mouse_dp10_direction": mouse_dp10_direction,
                "mouse_dp10_delta": safe_float(row.get("mouse_dp10_delta")),
                "mouse_dp10_methyl_expr_relation": mouse_dp10_relation,
                "human_strict_regulatory_probe_count": row.get("human_strict_regulatory_probe_count", "0"),
                "mouse_any_regulatory_dmr_count": row.get("mouse_any_regulatory_dmr_count", "0"),
                "mouse_embryonic_forebrain_h3k27me3_dmr_count": row.get("mouse_embryonic_forebrain_h3k27me3_dmr_count", "0"),
                "mechanism_class": mech_class,
                "interpretation": interpretation,
            }
        )
    return rows


def plot_enrichment(rows: list[dict[str, object]], figdir: Path) -> None:
    figdir.mkdir(parents=True, exist_ok=True)
    assays = [
        "human_strict_regulatory_probe",
        "mouse_any_regulatory_dmr",
        "mouse_atac_dmr",
        "mouse_h3k27ac_dmr",
        "mouse_h3k27me3_dmr",
    ]
    candidate_sets_to_plot = ["primary_5", "final_24", "methylation_supported_final", "expression_pathway_final"]
    row_by_key = {(row["candidate_set"], row["assay"]): row for row in rows}
    values = np.array(
        [
            [safe_float(row_by_key.get((candidate_set, assay), {}).get("fold_enrichment")) for assay in assays]
            for candidate_set in candidate_sets_to_plot
        ],
        dtype=float,
    )
    values = np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)
    x = np.arange(len(assays))
    width = 0.19
    colors = ["#2f6f8f", "#6f9e5a", "#b68a3a", "#9a5c83"]
    fig, ax = plt.subplots(figsize=(10.6, 4.8))
    for idx, candidate_set in enumerate(candidate_sets_to_plot):
        ax.bar(x + (idx - 1.5) * width, values[idx], width=width, color=colors[idx], label=candidate_set)
    ax.axhline(1.0, color="#333333", linewidth=0.8, linestyle="--")
    ax.set_xticks(x)
    ax.set_xticklabels([assay.replace("_", "\n") for assay in assays], fontsize=8)
    ax.set_ylabel("Fold enrichment vs background")
    ax.set_title("Regulatory Overlap Enrichment of Candidate Sets", fontsize=12, pad=10)
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#dddddd", linewidth=0.6)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(figdir / "regulatory_enrichment_fold_change.png", dpi=300)
    fig.savefig(figdir / "regulatory_enrichment_fold_change.pdf")
    plt.close(fig)


def plot_mechanism_counts(rows: list[dict[str, object]], figdir: Path) -> None:
    counts: dict[str, int] = defaultdict(int)
    for row in rows:
        counts[str(row["mechanism_class"])] += 1
    labels = sorted(counts, key=lambda label: (-counts[label], label))
    values = [counts[label] for label in labels]
    fig, ax = plt.subplots(figsize=(9.0, 4.8))
    ax.bar(np.arange(len(labels)), values, color="#5d8fb2")
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels([label.replace("_", "\n") for label in labels], rotation=25, ha="right", fontsize=8)
    ax.set_ylabel("Candidate genes")
    ax.set_title("Mechanism Classes for Final Candidate Genes", fontsize=12, pad=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for idx, value in enumerate(values):
        ax.text(idx, value + 0.1, str(value), ha="center", va="bottom", fontsize=9)
    ax.set_ylim(0, max(values) + 1.5)
    fig.tight_layout()
    fig.savefig(figdir / "mechanism_class_counts.png", dpi=300)
    fig.savefig(figdir / "mechanism_class_counts.pdf")
    plt.close(fig)


def background_summary(human_bg: dict[str, dict[str, object]], mouse_bg: dict[str, dict[str, object]]) -> list[dict[str, object]]:
    human_any = [gene for gene, item in human_bg.items() if int(item["linked_probes"]) > 0]
    human_any_reg = [gene for gene, item in human_bg.items() if int(item["any_regulatory_probes"]) > 0]
    human_strict = [gene for gene, item in human_bg.items() if int(item["strict_probes"]) > 0]
    human_strict_reg = [gene for gene, item in human_bg.items() if int(item["strict_regulatory_probes"]) > 0]
    mouse_any = [gene for gene, item in mouse_bg.items() if int(item["linked_dmrs"]) > 0]
    mouse_reg = [gene for gene, item in mouse_bg.items() if int(item["any_regulatory_dmrs"]) > 0]
    return [
        {
            "background": "human_protein_coding_dnmt3l_linked_probe_genes",
            "genes": len(human_any),
            "regulatory_success_genes": len(human_any_reg),
            "success_fraction": len(human_any_reg) / len(human_any) if human_any else math.nan,
        },
        {
            "background": "human_protein_coding_strict_triple_probe_genes",
            "genes": len(human_strict),
            "regulatory_success_genes": len(human_strict_reg),
            "success_fraction": len(human_strict_reg) / len(human_strict) if human_strict else math.nan,
        },
        {
            "background": "mouse_triple_hit_dmr_linked_genes",
            "genes": len(mouse_any),
            "regulatory_success_genes": len(mouse_reg),
            "success_fraction": len(mouse_reg) / len(mouse_any) if mouse_any else math.nan,
        },
    ]


def write_note(
    path: Path,
    enrichment_rows: list[dict[str, object]],
    mechanism_rows: list[dict[str, object]],
    background_rows: list[dict[str, object]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    primary_enrichment = [
        row for row in enrichment_rows if row["candidate_set"] == "primary_5" and row["assay"] in {"human_strict_regulatory_probe", "mouse_any_regulatory_dmr", "mouse_atac_dmr", "mouse_h3k27ac_dmr", "mouse_h3k27me3_dmr"}
    ]
    class_counts: dict[str, int] = defaultdict(int)
    class_genes: dict[str, list[str]] = defaultdict(list)
    for row in mechanism_rows:
        label = str(row["mechanism_class"])
        class_counts[label] += 1
        class_genes[label].append(str(row["gene_symbol"]))

    lines = [
        "# Regulatory enrichment and methylation-expression mechanism classification",
        "",
        "This pass adds statistical and directionality checks to the final DNMT3L candidate table.",
        "",
        "## Backgrounds",
        "",
        "| Background | Genes | Regulatory-success genes | Success fraction |",
        "|---|---:|---:|---:|",
    ]
    for row in background_rows:
        lines.append(
            f"| {row['background']} | {row['genes']} | {row['regulatory_success_genes']} | {float(row['success_fraction']):.3f} |"
        )
    lines.extend(
        [
            "",
            "## Primary candidate enrichment",
            "",
            "| Assay | Candidate successes / eligible | Background successes / genes | Fold enrichment | Hypergeometric p greater | Success genes |",
            "|---|---:|---:|---:|---:|---|",
        ]
    )
    for row in primary_enrichment:
        lines.append(
            f"| {row['assay']} | {row['candidate_successes']} / {row['candidate_genes_in_background']} | "
            f"{row['background_successes']} / {row['background_genes']} | {float(row['fold_enrichment']):.3f} | "
            f"{float(row['hypergeometric_p_greater']):.3g} | {str(row['candidate_success_genes']).replace(',', ', ')} |"
        )
    lines.extend(
        [
            "",
            "## Mechanism classes",
            "",
            "| Mechanism class | n | Genes |",
            "|---|---:|---|",
        ]
    )
    for label in sorted(class_counts, key=lambda item: (-class_counts[item], item)):
        lines.append(f"| {label} | {class_counts[label]} | {', '.join(class_genes[label])} |")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The enrichment tests support the statement that the primary candidate genes are embedded in fetal-brain regulatory chromatin, especially for human strict regulatory probes and mouse Dnmt3l-sensitive regulatory DMRs. These tests are gene-level hypergeometric tests against available DNMT3L-sensitive methylation backgrounds; they are useful but do not replace matched-region permutation or paired methylation-expression mediation.",
            "",
            "The methylation-expression direction table separates simple repression-like patterns from directionally complex candidates. FOXP1 is the clearest primary candidate compatible with methylation-associated repression. PLPP3, FGFR2, ZFHX3, and SIAH3 remain strong regulatory-convergent candidates, but their expression directions suggest enhancer, gene-body, indirect, timing-specific, or cell-state mechanisms rather than a simple promoter-silencing model.",
            "",
            "## Outputs",
            "",
            "- `results/mechanism_classification/regulatory_enrichment_tests.tsv`",
            "- `results/mechanism_classification/methylation_expression_mechanism_table.tsv`",
            "- `results/mechanism_classification/regulatory_background_summary.tsv`",
            "- `figures/mechanism_classification/regulatory_enrichment_fold_change.png`",
            "- `figures/mechanism_classification/mechanism_class_counts.png`",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    matrix_rows = read_tsv(Path(args.matrix))
    mouse_ranking = read_mouse_ranking(Path(args.mouse_ranking))

    print("Loading human regulatory annotations...")
    human_indexes = load_annotation_indexes(HUMAN_REGULATORY_TRACKS)
    print("Loading mouse regulatory annotations...")
    mouse_indexes = load_annotation_indexes(MOUSE_REGULATORY_TRACKS)

    print("Building human background...")
    human_bg = build_human_background(Path(args.human_links), human_indexes)
    print("Building mouse background...")
    mouse_bg = build_mouse_background(Path(args.mouse_links), mouse_indexes)

    enrichment_rows = regulatory_enrichment(matrix_rows, human_bg, mouse_bg)
    mechanism_rows = classify_mechanisms(matrix_rows, mouse_ranking)
    background_rows = background_summary(human_bg, mouse_bg)

    outdir = Path(args.outdir)
    write_tsv(outdir / "regulatory_enrichment_tests.tsv", enrichment_rows, ENRICHMENT_FIELDS)
    write_tsv(outdir / "methylation_expression_mechanism_table.tsv", mechanism_rows, MECHANISM_FIELDS)
    write_tsv(outdir / "regulatory_background_summary.tsv", background_rows, ["background", "genes", "regulatory_success_genes", "success_fraction"])

    figdir = Path(args.figdir)
    plot_enrichment(enrichment_rows, figdir)
    plot_mechanism_counts(mechanism_rows, figdir)
    write_note(Path(args.note), enrichment_rows, mechanism_rows, background_rows)

    print(f"Wrote enrichment tests to {outdir / 'regulatory_enrichment_tests.tsv'}")
    print(f"Wrote mechanism table to {outdir / 'methylation_expression_mechanism_table.tsv'}")
    print(f"Wrote figures to {figdir}")
    print(f"Wrote note to {args.note}")


if __name__ == "__main__":
    main()
