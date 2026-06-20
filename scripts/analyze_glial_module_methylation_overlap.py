#!/usr/bin/env python3
"""Intersect a Dnmt3l OE developmental module with methylation-linked genes."""

from __future__ import annotations

import argparse
import bisect
import csv
import math
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib.pyplot as plt
from scipy import stats

import rank_mouse_methylation_expression_candidate_genes as dmr_utils
import test_developmental_module_shifts as module_utils


DEFAULT_EXPR = "results/neurogenesis_mechanism/gse268951_rederived_gene_expression.tsv"
DEFAULT_OVERLAP = "results/neurogenesis_mechanism/dnmt3l_vs_dnmt3a_transcriptomic_overlap.tsv"
DEFAULT_FINAL = "results/neurogenesis_mechanism/final_candidate_neurogenesis_annotation.tsv"
DEFAULT_RAW_RNA = "data/geo/GSE268951/supplementary/GSE268951_all.fpkm_mRNA_anno.txt.gz"
DEFAULT_DMRS = "results/mouse_candidate_dmrs_sex_consistent/GSE271949_candidate_dmrs.tsv"
DEFAULT_TRIPLE_LINKS = "results/integrated_candidate_genes/mouse_triple_hit_dmr_gene_links.tsv"
DEFAULT_OUTDIR = "results/neurogenesis_mechanism"
DEFAULT_FIGDIR = "figures/neurogenesis_mechanism"
DEFAULT_MODULE = "glial_maturation_myelination"
DEFAULT_OUTPUT_PREFIX = "glial_module"

MODULE_LABELS = {
    "early_neurogenesis_neuron_differentiation": "early neurogenesis/neuron differentiation",
    "cortical_patterning_region_identity": "cortical patterning/region identity",
    "axon_dendrite_neuron_projection": "axon/dendrite/neuron projection",
    "synapse_formation_maturation": "synapse formation/maturation",
    "glial_maturation_myelination": "glial maturation/myelination",
    "neurodevelopment_signaling": "neurodevelopment signaling",
    "cell_cycle_progenitor_proliferation": "cell-cycle/progenitor proliferation",
    "final24_all_candidates": "final 24 candidates",
    "final16_neurodevelopment_candidates": "final 16 neurodevelopment candidates",
    "final_synapse_window_candidates": "final synapse-window candidates",
    "primary5_candidates": "primary five candidates",
}
BIN_SIZE = 100_000


@dataclass
class GeneMethylationSupport:
    any_dmr_ids: set[str] = field(default_factory=set)
    promoter_dmr_ids: set[str] = field(default_factory=set)
    gene_body_dmr_ids: set[str] = field(default_factory=set)
    nearest_dmr_ids: set[str] = field(default_factory=set)
    hyper_dmr_ids: set[str] = field(default_factory=set)
    hypo_dmr_ids: set[str] = field(default_factory=set)
    delta_by_dmr: dict[str, float] = field(default_factory=dict)
    min_abs_distance_to_tss: float = math.inf

    def add(
        self,
        dmr_id: str,
        relationship: str,
        direction: str,
        delta_beta: float,
        abs_distance_to_tss: float,
    ) -> None:
        self.any_dmr_ids.add(dmr_id)
        if relationship == "promoter_2kb":
            self.promoter_dmr_ids.add(dmr_id)
        elif relationship == "gene_body":
            self.gene_body_dmr_ids.add(dmr_id)
        elif relationship == "nearest_100kb":
            self.nearest_dmr_ids.add(dmr_id)
        if direction == "hyper":
            self.hyper_dmr_ids.add(dmr_id)
        elif direction == "hypo":
            self.hypo_dmr_ids.add(dmr_id)
        if math.isfinite(delta_beta):
            self.delta_by_dmr[dmr_id] = delta_beta
        self.min_abs_distance_to_tss = min(self.min_abs_distance_to_tss, abs_distance_to_tss)

    @property
    def direct_dmr_ids(self) -> set[str]:
        return self.promoter_dmr_ids | self.gene_body_dmr_ids

    @property
    def mean_delta_beta(self) -> float:
        values = list(self.delta_by_dmr.values())
        return sum(values) / len(values) if values else math.nan


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expr", default=DEFAULT_EXPR)
    parser.add_argument("--overlap", default=DEFAULT_OVERLAP)
    parser.add_argument("--final", default=DEFAULT_FINAL)
    parser.add_argument("--raw-rna", default=DEFAULT_RAW_RNA)
    parser.add_argument("--dmrs", default=DEFAULT_DMRS)
    parser.add_argument("--triple-links", default=DEFAULT_TRIPLE_LINKS)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--figdir", default=DEFAULT_FIGDIR)
    parser.add_argument("--module", default=DEFAULT_MODULE)
    parser.add_argument("--module-label", default="")
    parser.add_argument("--output-prefix", default=DEFAULT_OUTPUT_PREFIX)
    parser.add_argument("--note", default="")
    parser.add_argument("--promoter-window", type=int, default=2000)
    parser.add_argument("--nearest-window", type=int, default=100000)
    parser.add_argument("--min-fpkm-max", type=float, default=0.5)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [{k: (v or "").strip() for k, v in row.items()} for row in reader]


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})


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
    if isinstance(value, int):
        return str(value)
    return str(value)


def yes_no(value: bool) -> str:
    return "yes" if value else "no"


def overlap(start_a: int, end_a: int, start_b: int, end_b: int) -> bool:
    return min(end_a, end_b) > max(start_a, start_b)


def interval_bins(start: int, end: int) -> range:
    return range(max(0, start // BIN_SIZE), max(0, (max(start, end - 1)) // BIN_SIZE) + 1)


def expression_universe(
    expression_rows: dict[str, dict[str, object]],
    min_fpkm_max: float,
) -> set[str]:
    universe = set()
    for symbol, row in expression_rows.items():
        logfc = safe_float(row.get("gse268951_log2fc_oe_vs_control"))
        mean_fpkm_max = safe_float(row.get("mean_fpkm_max"))
        if math.isfinite(logfc) and mean_fpkm_max >= min_fpkm_max:
            universe.add(symbol)
    return universe


def expression_details(path: Path) -> dict[str, dict[str, str]]:
    return {row["gene_symbol"]: row for row in read_tsv(path) if row.get("gene_symbol")}


def build_direct_interval_bins(
    genes: dict[str, dmr_utils.Gene],
    universe: set[str],
    promoter_window: int,
) -> dict[str, dict[int, list[tuple[int, int, str, str, int]]]]:
    bins: dict[str, dict[int, list[tuple[int, int, str, str, int]]]] = defaultdict(lambda: defaultdict(list))
    for gene in genes.values():
        if gene.gene_symbol not in universe:
            continue
        promoter_start = max(0, gene.tss - promoter_window)
        promoter_end = gene.tss + promoter_window
        intervals = [
            (promoter_start, promoter_end, gene.gene_symbol, "promoter_2kb", gene.tss),
            (gene.start, gene.end, gene.gene_symbol, "gene_body", gene.tss),
        ]
        for start, end, symbol, relationship, tss in intervals:
            for bin_id in interval_bins(start, end):
                bins[gene.chrom][bin_id].append((start, end, symbol, relationship, tss))
    return bins


def build_tss_index(
    genes: dict[str, dmr_utils.Gene],
) -> tuple[dict[str, list[int]], dict[str, list[dmr_utils.Gene]]]:
    by_chrom: dict[str, list[dmr_utils.Gene]] = defaultdict(list)
    for gene in genes.values():
        by_chrom[gene.chrom].append(gene)
    tss_positions: dict[str, list[int]] = {}
    tss_genes: dict[str, list[dmr_utils.Gene]] = {}
    for chrom, chrom_genes in by_chrom.items():
        ordered = sorted(chrom_genes, key=lambda gene: (gene.tss, gene.gene_id))
        tss_positions[chrom] = [gene.tss for gene in ordered]
        tss_genes[chrom] = ordered
    return tss_positions, tss_genes


def nearest_gene_for_dmr(
    dmr: dmr_utils.DMR,
    tss_positions: dict[str, list[int]],
    tss_genes: dict[str, list[dmr_utils.Gene]],
) -> tuple[dmr_utils.Gene | None, float]:
    positions = tss_positions.get(dmr.chrom, [])
    genes = tss_genes.get(dmr.chrom, [])
    if not positions:
        return None, math.inf
    midpoint = (dmr.start + dmr.end) // 2
    idx = bisect.bisect_left(positions, midpoint)
    best_gene = None
    best_distance = math.inf
    for candidate_idx in {idx - 2, idx - 1, idx, idx + 1, idx + 2}:
        if candidate_idx < 0 or candidate_idx >= len(positions):
            continue
        distance = abs(midpoint - positions[candidate_idx])
        if distance < best_distance:
            best_distance = distance
            best_gene = genes[candidate_idx]
    return best_gene, best_distance


def summarize_dnmt3l_dmr_support(
    dmrs: list[dmr_utils.DMR],
    genes: dict[str, dmr_utils.Gene],
    universe: set[str],
    promoter_window: int,
    nearest_window: int,
) -> dict[str, GeneMethylationSupport]:
    support: dict[str, GeneMethylationSupport] = defaultdict(GeneMethylationSupport)
    direct_bins = build_direct_interval_bins(genes, universe, promoter_window)
    tss_positions, tss_genes = build_tss_index(genes)

    for dmr in dmrs:
        seen_direct: set[tuple[str, str]] = set()
        for bin_id in interval_bins(dmr.start, dmr.end):
            for start, end, symbol, relationship, tss in direct_bins.get(dmr.chrom, {}).get(bin_id, []):
                if (symbol, relationship) in seen_direct:
                    continue
                if not overlap(dmr.start, dmr.end, start, end):
                    continue
                seen_direct.add((symbol, relationship))
                distance = abs((dmr.start + dmr.end) // 2 - tss)
                support[symbol].add(dmr.dmr_id, relationship, dmr.direction, dmr.mean_delta_beta, distance)

        nearest_gene, nearest_distance = nearest_gene_for_dmr(dmr, tss_positions, tss_genes)
        if nearest_gene and nearest_gene.gene_symbol in universe and nearest_distance <= nearest_window:
            support[nearest_gene.gene_symbol].add(
                dmr.dmr_id,
                "nearest_100kb",
                dmr.direction,
                dmr.mean_delta_beta,
                nearest_distance,
            )
    return support


def summarize_triple_hit_support(path: Path) -> dict[str, GeneMethylationSupport]:
    support: dict[str, GeneMethylationSupport] = defaultdict(GeneMethylationSupport)
    for row in read_tsv(path):
        symbol = row.get("gene_symbol", "")
        if not symbol:
            continue
        relationship = row.get("relationship", "")
        support[symbol].add(
            row.get("dmr_id", ""),
            relationship,
            row.get("dmr_direction", ""),
            safe_float(row.get("dmr_mean_delta_beta")),
            safe_float(row.get("abs_distance_to_tss")),
        )
    return support


def dominant_relation(support: GeneMethylationSupport) -> str:
    if support.promoter_dmr_ids:
        return "promoter_2kb"
    if support.gene_body_dmr_ids:
        return "gene_body"
    if support.nearest_dmr_ids:
        return "nearest_100kb"
    return "none"


def methylation_expression_class(
    expression_direction: str,
    support: GeneMethylationSupport,
) -> str:
    if not support.any_dmr_ids:
        return "no_local_dmr_support"
    if expression_direction not in {"up", "down"}:
        return "dmr_linked_expression_subthreshold"
    if support.promoter_dmr_ids:
        mean_delta = support.mean_delta_beta
        if not math.isfinite(mean_delta):
            return "promoter_dmr_direction_unclear"
        if mean_delta > 0 and expression_direction == "down":
            return "promoter_hyper_expression_down"
        if mean_delta < 0 and expression_direction == "up":
            return "promoter_hypo_expression_up"
        if mean_delta > 0 and expression_direction == "up":
            return "promoter_hyper_expression_up"
        if mean_delta < 0 and expression_direction == "down":
            return "promoter_hypo_expression_down"
    if support.gene_body_dmr_ids:
        return "gene_body_dmr_expression_changed"
    return "nearest_dmr_expression_changed"


def set_for_class(
    class_name: str,
    dnmt3l_support: dict[str, GeneMethylationSupport],
    triple_support: dict[str, GeneMethylationSupport],
    expression: dict[str, dict[str, str]],
) -> set[str]:
    def changed(threshold: float) -> set[str]:
        return {
            symbol
            for symbol, row in expression.items()
            if abs(safe_float(row.get("gse268951_log2fc_oe_vs_control"))) >= threshold
        }

    if class_name == "any_dnmt3l_dmr_within_gene_or_100kb":
        return {symbol for symbol, item in dnmt3l_support.items() if item.any_dmr_ids}
    if class_name == "direct_dnmt3l_dmr_promoter_or_gene_body":
        return {symbol for symbol, item in dnmt3l_support.items() if item.direct_dmr_ids}
    if class_name == "promoter_2kb_dnmt3l_dmr":
        return {symbol for symbol, item in dnmt3l_support.items() if item.promoter_dmr_ids}
    if class_name == "gene_body_dnmt3l_dmr":
        return {symbol for symbol, item in dnmt3l_support.items() if item.gene_body_dmr_ids}
    if class_name == "triple_hit_dmr_any":
        return {symbol for symbol, item in triple_support.items() if item.any_dmr_ids}
    if class_name == "triple_hit_dmr_direct":
        return {symbol for symbol, item in triple_support.items() if item.direct_dmr_ids}
    if class_name == "dnmt3l_oe_de_p05_abslog2fc1":
        return {symbol for symbol, row in expression.items() if row.get("gse268951_de_p05_fc2") == "yes"}
    if class_name == "abs_log2fc_ge_0.25":
        return changed(0.25)
    if class_name == "abs_log2fc_ge_0.58":
        return changed(0.58)
    if class_name == "abs_log2fc_ge_0.25_and_any_dnmt3l_dmr":
        any_dmr = set_for_class("any_dnmt3l_dmr_within_gene_or_100kb", dnmt3l_support, triple_support, expression)
        return changed(0.25) & any_dmr
    if class_name == "abs_log2fc_ge_0.25_and_direct_dnmt3l_dmr":
        direct = set_for_class("direct_dnmt3l_dmr_promoter_or_gene_body", dnmt3l_support, triple_support, expression)
        return changed(0.25) & direct
    if class_name == "abs_log2fc_ge_0.58_and_any_dnmt3l_dmr":
        any_dmr = set_for_class("any_dnmt3l_dmr_within_gene_or_100kb", dnmt3l_support, triple_support, expression)
        return changed(0.58) & any_dmr
    if class_name == "abs_log2fc_ge_0.58_and_direct_dnmt3l_dmr":
        direct = set_for_class("direct_dnmt3l_dmr_promoter_or_gene_body", dnmt3l_support, triple_support, expression)
        return changed(0.58) & direct
    if class_name == "abs_log2fc_ge_0.58_and_triple_hit_dmr":
        triple = set_for_class("triple_hit_dmr_any", dnmt3l_support, triple_support, expression)
        return changed(0.58) & triple
    if class_name == "de_and_any_dnmt3l_dmr":
        de = set_for_class("dnmt3l_oe_de_p05_abslog2fc1", dnmt3l_support, triple_support, expression)
        any_dmr = set_for_class("any_dnmt3l_dmr_within_gene_or_100kb", dnmt3l_support, triple_support, expression)
        return de & any_dmr
    if class_name == "de_and_direct_dnmt3l_dmr":
        de = set_for_class("dnmt3l_oe_de_p05_abslog2fc1", dnmt3l_support, triple_support, expression)
        direct = set_for_class("direct_dnmt3l_dmr_promoter_or_gene_body", dnmt3l_support, triple_support, expression)
        return de & direct
    if class_name == "de_and_triple_hit_dmr":
        de = set_for_class("dnmt3l_oe_de_p05_abslog2fc1", dnmt3l_support, triple_support, expression)
        triple = set_for_class("triple_hit_dmr_any", dnmt3l_support, triple_support, expression)
        return de & triple
    raise ValueError(f"Unsupported class: {class_name}")


def fisher_row(class_name: str, hits: set[str], module_genes: set[str], universe: set[str]) -> dict[str, object]:
    background = universe - module_genes
    module_hits = hits & module_genes
    background_hits = hits & background
    table = [
        [len(module_hits), len(module_genes) - len(module_hits)],
        [len(background_hits), len(background) - len(background_hits)],
    ]
    odds_ratio, pvalue = stats.fisher_exact(table, alternative="greater")
    return {
        "feature_class": class_name,
        "module_genes": len(module_genes),
        "module_hits": len(module_hits),
        "module_fraction": len(module_hits) / len(module_genes) if module_genes else math.nan,
        "background_genes": len(background),
        "background_hits": len(background_hits),
        "background_fraction": len(background_hits) / len(background) if background else math.nan,
        "odds_ratio_module_enrichment": odds_ratio,
        "fisher_p_greater": pvalue,
        "module_hit_symbols": ",".join(sorted(module_hits)),
    }


def make_gene_rows(
    module_genes: set[str],
    expression: dict[str, dict[str, str]],
    dnmt3l_support: dict[str, GeneMethylationSupport],
    triple_support: dict[str, GeneMethylationSupport],
) -> list[dict[str, object]]:
    rows = []
    for symbol in sorted(module_genes):
        expr = expression.get(symbol, {})
        dmr = dnmt3l_support.get(symbol, GeneMethylationSupport())
        triple = triple_support.get(symbol, GeneMethylationSupport())
        logfc = safe_float(expr.get("gse268951_log2fc_oe_vs_control"))
        expression_direction = expr.get("gse268951_direction", "")
        if expression_direction not in {"up", "down", "flat"}:
            if logfc >= 0.25:
                expression_direction = "up"
            elif logfc <= -0.25:
                expression_direction = "down"
            else:
                expression_direction = "flat"
        rows.append(
            {
                "gene_symbol": symbol,
                "description": expr.get("description", ""),
                "gse268951_control_mean_fpkm": expr.get("gse268951_control_mean_fpkm", ""),
                "gse268951_dnmt3l_oe_mean_fpkm": expr.get("gse268951_dnmt3l_oe_mean_fpkm", ""),
                "gse268951_log2fc_oe_vs_control": logfc,
                "gse268951_effect_source": expr.get("gse268951_effect_source", ""),
                "gse268951_primary_pvalue": expr.get("gse268951_primary_pvalue", ""),
                "gse268951_primary_fdr": expr.get("gse268951_primary_fdr", ""),
                "gse268951_welch_pvalue": expr.get("gse268951_welch_pvalue", ""),
                "gse268951_bh_fdr": expr.get("gse268951_bh_fdr", ""),
                "gse268951_de_p05_fc2": expr.get("gse268951_de_p05_fc2", ""),
                "gse268951_de_fdr05_fc2": expr.get("gse268951_de_fdr05_fc2", ""),
                "gse268951_direction": expression_direction,
                "dnmt3l_any_dmr_count": len(dmr.any_dmr_ids),
                "dnmt3l_direct_dmr_count": len(dmr.direct_dmr_ids),
                "dnmt3l_promoter_2kb_dmr_count": len(dmr.promoter_dmr_ids),
                "dnmt3l_gene_body_dmr_count": len(dmr.gene_body_dmr_ids),
                "dnmt3l_nearest_100kb_dmr_count": len(dmr.nearest_dmr_ids),
                "dnmt3l_hyper_dmr_count": len(dmr.hyper_dmr_ids),
                "dnmt3l_hypo_dmr_count": len(dmr.hypo_dmr_ids),
                "dnmt3l_mean_dmr_delta_beta": dmr.mean_delta_beta,
                "dnmt3l_min_abs_distance_to_tss": "" if math.isinf(dmr.min_abs_distance_to_tss) else dmr.min_abs_distance_to_tss,
                "dnmt3l_dominant_relation": dominant_relation(dmr),
                "triple_hit_any_dmr_count": len(triple.any_dmr_ids),
                "triple_hit_direct_dmr_count": len(triple.direct_dmr_ids),
                "methylation_expression_class": methylation_expression_class(expression_direction, dmr),
                "dnmt3l_dmr_ids": ",".join(sorted(dmr.any_dmr_ids)[:30]) + ("..." if len(dmr.any_dmr_ids) > 30 else ""),
                "triple_hit_dmr_ids": ",".join(sorted(triple.any_dmr_ids)[:30]) + ("..." if len(triple.any_dmr_ids) > 30 else ""),
                "go_bp_terms": expr.get("go_bp_terms", ""),
            }
        )
    rows.sort(
        key=lambda row: (
            row["gse268951_de_p05_fc2"] != "yes",
            -int(row["triple_hit_any_dmr_count"]),
            -int(row["dnmt3l_direct_dmr_count"]),
            -int(row["dnmt3l_any_dmr_count"]),
            -abs(safe_float(row["gse268951_log2fc_oe_vs_control"])),
            str(row["gene_symbol"]),
        )
    )
    return rows


def module_label(module: str, override: str = "") -> str:
    return override or MODULE_LABELS.get(module, module.replace("_", " "))


def make_plot(summary_rows: list[dict[str, object]], figdir: Path, label: str, output_prefix: str) -> None:
    figdir.mkdir(parents=True, exist_ok=True)
    wanted = [
        "any_dnmt3l_dmr_within_gene_or_100kb",
        "direct_dnmt3l_dmr_promoter_or_gene_body",
        "triple_hit_dmr_any",
        "abs_log2fc_ge_0.25",
        "abs_log2fc_ge_0.25_and_direct_dnmt3l_dmr",
        "abs_log2fc_ge_0.58",
        "abs_log2fc_ge_0.58_and_direct_dnmt3l_dmr",
        "dnmt3l_oe_de_p05_abslog2fc1",
    ]
    labels = [
        "Any DMR\n<=100 kb",
        "Promoter/body\nDMR",
        "Triple-hit\nDMR",
        "Abs logFC\n>=0.25",
        "0.25 +\nprom/body",
        "Abs logFC\n>=0.58",
        "0.58 +\nprom/body",
        "DE\np<0.05 FC>=2",
    ]
    rows_by_name = {row["feature_class"]: row for row in summary_rows}
    module_fracs = [safe_float(rows_by_name[name]["module_fraction"]) for name in wanted]
    background_fracs = [safe_float(rows_by_name[name]["background_fraction"]) for name in wanted]
    x = range(len(wanted))

    fig, ax = plt.subplots(figsize=(9.2, 4.7))
    width = 0.38
    ax.bar([i - width / 2 for i in x], module_fracs, width=width, color="#1F77B4", label="Module")
    ax.bar([i + width / 2 for i in x], background_fracs, width=width, color="#9CA3AF", label="Background")
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("Fraction of genes")
    ax.set_title(f"Methylation-linked genes in Dnmt3l OE {label} module")
    ax.legend(frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#E5E7EB", linewidth=0.8)
    fig.tight_layout()
    fig.savefig(figdir / f"{output_prefix}_methylation_overlap_counts.png", dpi=260, bbox_inches="tight")
    fig.savefig(figdir / f"{output_prefix}_methylation_overlap_counts.pdf", bbox_inches="tight")
    plt.close(fig)


def write_note(
    path: Path,
    summary_rows: list[dict[str, object]],
    direction_rows: list[dict[str, object]],
    label: str,
    output_prefix: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    by_name = {row["feature_class"]: row for row in summary_rows}

    def count(name: str) -> str:
        row = by_name[name]
        return f"{row['module_hits']}/{row['module_genes']} ({safe_float(row['module_fraction']) * 100:.1f}%)"

    lines = [
        f"# {label.title()} Module Methylation Overlap",
        "",
        "Date: 2026-06-09",
        "",
        "## Question",
        "",
        f"The Dnmt3l overexpression rank test showed a shift in the {label} module. This analysis asks whether those module genes are linked to Dnmt3l-sensitive DNA methylation regions.",
        "",
        "## Definitions",
        "",
        "- `any DMR` means a sex-consistent GSE271949 Dnmt3l-sensitive DMR overlapping a gene body/promoter or assigned as the nearest gene within 100 kb.",
        "- `direct DMR` means promoter 2 kb or gene-body overlap.",
        "- `triple-hit DMR` means a stricter Dnmt3l-sensitive DMR that also overlaps the DNMT3A-dependent/postnatal neuronal methylation reference used in the mouse candidate-gene analysis.",
        "",
        "## Key Counts",
        "",
        f"- Any Dnmt3l DMR-linked module genes: {count('any_dnmt3l_dmr_within_gene_or_100kb')}.",
        f"- Direct promoter/body DMR-linked module genes: {count('direct_dnmt3l_dmr_promoter_or_gene_body')}.",
        f"- Promoter DMR-linked module genes: {count('promoter_2kb_dnmt3l_dmr')}.",
        f"- Gene-body DMR-linked module genes: {count('gene_body_dnmt3l_dmr')}.",
        f"- Triple-hit DMR-linked module genes: {count('triple_hit_dmr_any')}.",
        f"- Moderate shifted module genes with abs(log2FC) >= 0.25: {count('abs_log2fc_ge_0.25')}.",
        f"- Moderate shifted plus direct DMR-linked module genes: {count('abs_log2fc_ge_0.25_and_direct_dnmt3l_dmr')}.",
        f"- Larger shifted module genes with abs(log2FC) >= 0.58: {count('abs_log2fc_ge_0.58')}.",
        f"- Larger shifted plus direct DMR-linked module genes: {count('abs_log2fc_ge_0.58_and_direct_dnmt3l_dmr')}.",
        f"- Strong DE module genes: {count('dnmt3l_oe_de_p05_abslog2fc1')}.",
        f"- Strong DE plus direct DMR-linked module genes: {count('de_and_direct_dnmt3l_dmr')}.",
        "",
        "## Direction Classes",
        "",
    ]
    for row in direction_rows:
        lines.append(f"- {row['methylation_expression_class']}: {row['gene_count']} genes.")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            f"A high fraction of {label} module genes have nearby Dnmt3l-sensitive DMRs, and the direct promoter/body subset is the most useful methylation-expression subset. This supports a methylation-associated contribution to the module expression shift, but it does not prove direct methylation mediation gene by gene.",
            "",
            "For mechanistic interpretation, promoter hypermethylation with downregulation or promoter hypomethylation with upregulation is the cleanest simple repression model. Gene-body or nearest-DMR links, and promoter hypermethylation with upregulation, should be interpreted as enhancer/gene-body methylation, cell-state shifts, or indirect regulation until tested with matched methylation-expression mediation.",
            "",
            "## Outputs",
            "",
            f"- `results/neurogenesis_mechanism/{output_prefix}_methylation_overlap_summary.tsv`",
            f"- `results/neurogenesis_mechanism/{output_prefix}_methylation_gene_overlap.tsv`",
            f"- `results/neurogenesis_mechanism/{output_prefix}_methylation_direction_summary.tsv`",
            f"- `figures/neurogenesis_mechanism/{output_prefix}_methylation_overlap_counts.png`",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    figdir = Path(args.figdir)
    label = module_label(args.module, args.module_label)
    output_prefix = args.output_prefix
    note_path = Path(args.note) if args.note else Path("notes") / f"{output_prefix}_methylation_overlap.md"

    expression_rows = module_utils.load_expression_rows(Path(args.expr), Path(args.overlap))
    final_rows = module_utils.read_tsv(Path(args.final))
    modules = module_utils.build_modules(expression_rows, final_rows)
    if args.module not in modules:
        available = ", ".join(sorted(modules))
        raise SystemExit(f"Module {args.module!r} was not found. Available modules: {available}")
    universe = expression_universe(expression_rows, args.min_fpkm_max)
    module_genes = modules[args.module] & universe

    expression = expression_details(Path(args.expr))
    raw_genes = dmr_utils.read_gene_models(Path(args.raw_rna))
    dmrs = dmr_utils.read_dmrs(Path(args.dmrs))
    dnmt3l_support = summarize_dnmt3l_dmr_support(
        dmrs,
        raw_genes,
        universe,
        args.promoter_window,
        args.nearest_window,
    )
    triple_support = summarize_triple_hit_support(Path(args.triple_links))

    feature_classes = [
        "any_dnmt3l_dmr_within_gene_or_100kb",
        "direct_dnmt3l_dmr_promoter_or_gene_body",
        "promoter_2kb_dnmt3l_dmr",
        "gene_body_dnmt3l_dmr",
        "triple_hit_dmr_any",
        "triple_hit_dmr_direct",
        "abs_log2fc_ge_0.25",
        "abs_log2fc_ge_0.25_and_any_dnmt3l_dmr",
        "abs_log2fc_ge_0.25_and_direct_dnmt3l_dmr",
        "abs_log2fc_ge_0.58",
        "abs_log2fc_ge_0.58_and_any_dnmt3l_dmr",
        "abs_log2fc_ge_0.58_and_direct_dnmt3l_dmr",
        "abs_log2fc_ge_0.58_and_triple_hit_dmr",
        "dnmt3l_oe_de_p05_abslog2fc1",
        "de_and_any_dnmt3l_dmr",
        "de_and_direct_dnmt3l_dmr",
        "de_and_triple_hit_dmr",
    ]
    summary_rows = [
        fisher_row(
            class_name,
            set_for_class(class_name, dnmt3l_support, triple_support, expression),
            module_genes,
            universe,
        )
        for class_name in feature_classes
    ]

    gene_rows = make_gene_rows(module_genes, expression, dnmt3l_support, triple_support)
    direction_counts: dict[str, int] = defaultdict(int)
    for row in gene_rows:
        direction_counts[str(row["methylation_expression_class"])] += 1
    direction_rows = [
        {"methylation_expression_class": key, "gene_count": value}
        for key, value in sorted(direction_counts.items(), key=lambda item: (-item[1], item[0]))
    ]

    summary_fields = [
        "feature_class",
        "module_genes",
        "module_hits",
        "module_fraction",
        "background_genes",
        "background_hits",
        "background_fraction",
        "odds_ratio_module_enrichment",
        "fisher_p_greater",
        "module_hit_symbols",
    ]
    gene_fields = [
        "gene_symbol",
        "description",
        "gse268951_control_mean_fpkm",
        "gse268951_dnmt3l_oe_mean_fpkm",
        "gse268951_log2fc_oe_vs_control",
        "gse268951_effect_source",
        "gse268951_primary_pvalue",
        "gse268951_primary_fdr",
        "gse268951_welch_pvalue",
        "gse268951_bh_fdr",
        "gse268951_de_p05_fc2",
        "gse268951_de_fdr05_fc2",
        "gse268951_direction",
        "dnmt3l_any_dmr_count",
        "dnmt3l_direct_dmr_count",
        "dnmt3l_promoter_2kb_dmr_count",
        "dnmt3l_gene_body_dmr_count",
        "dnmt3l_nearest_100kb_dmr_count",
        "dnmt3l_hyper_dmr_count",
        "dnmt3l_hypo_dmr_count",
        "dnmt3l_mean_dmr_delta_beta",
        "dnmt3l_min_abs_distance_to_tss",
        "dnmt3l_dominant_relation",
        "triple_hit_any_dmr_count",
        "triple_hit_direct_dmr_count",
        "methylation_expression_class",
        "dnmt3l_dmr_ids",
        "triple_hit_dmr_ids",
        "go_bp_terms",
    ]
    summary_path = outdir / f"{output_prefix}_methylation_overlap_summary.tsv"
    gene_path = outdir / f"{output_prefix}_methylation_gene_overlap.tsv"
    direction_path = outdir / f"{output_prefix}_methylation_direction_summary.tsv"
    fig_path = figdir / f"{output_prefix}_methylation_overlap_counts.png"
    write_tsv(summary_path, summary_rows, summary_fields)
    write_tsv(gene_path, gene_rows, gene_fields)
    write_tsv(direction_path, direction_rows, ["methylation_expression_class", "gene_count"])
    make_plot(summary_rows, figdir, label, output_prefix)
    write_note(note_path, summary_rows, direction_rows, label, output_prefix)

    print(f"Wrote {summary_path}")
    print(f"Wrote {gene_path}")
    print(f"Wrote {direction_path}")
    print(f"Wrote {fig_path}")
    print(f"Wrote {note_path}")


if __name__ == "__main__":
    main()
