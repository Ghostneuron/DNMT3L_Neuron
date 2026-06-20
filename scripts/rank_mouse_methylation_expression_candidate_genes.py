#!/usr/bin/env python3
"""Rank mouse candidate genes by Dnmt3l-sensitive DMRs and expression support."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from urllib.parse import unquote


DEFAULT_DMRS = "results/mouse_dnmt3l_dnmt3a_neuron_overlap/gse271949_dmr_priority_triple_hits.tsv"
DEFAULT_RNA = "data/geo/GSE268951/supplementary/GSE268951_all.fpkm_mRNA_anno.txt.gz"
DEFAULT_COUNT_DE = "results/mouse_gse268951_count_de/gse268951_deseq2_edger_count_de_results.tsv"
DEFAULT_DP10 = "results/gse214273_dp10_scrna/gse214273_gene_effect_summary.tsv"
DEFAULT_OUTDIR = "results/integrated_candidate_genes"

CONTROL_SAMPLES = ["A-1", "A-2", "A-3", "B-1", "B-2", "B-3"]
OE_SAMPLES = ["C-1", "C-2", "C-3", "D-1", "D-2", "D-3"]
MALE_CONTROL_SAMPLES = ["A-1", "A-2", "A-3"]
MALE_OE_SAMPLES = ["C-1", "C-2", "C-3"]
FEMALE_CONTROL_SAMPLES = ["B-1", "B-2", "B-3"]
FEMALE_OE_SAMPLES = ["D-1", "D-2", "D-3"]
ALL_SAMPLES = CONTROL_SAMPLES + OE_SAMPLES
TARGET_GENES = ["Dnmt3l", "Dnmt3a", "Dnmt3b", "Stat1", "Stat3"]


@dataclass(frozen=True)
class DMR:
    dmr_id: str
    chrom: str
    start: int
    end: int
    direction: str
    n_cpgs: int
    mean_delta_beta: float
    row: dict[str, str]


@dataclass
class Gene:
    gene_id: str
    gene_symbol: str
    chrom: str
    start: int
    end: int
    strand: str
    description: str = ""
    transcript_count: int = 0
    fpkm: dict[str, float] = field(default_factory=lambda: defaultdict(float))

    @property
    def tss(self) -> int:
        return self.end if self.strand == "-" else self.start

    @property
    def length(self) -> int:
        return max(0, self.end - self.start)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dmrs", default=DEFAULT_DMRS)
    parser.add_argument("--rna", default=DEFAULT_RNA)
    parser.add_argument("--count-de", default=DEFAULT_COUNT_DE)
    parser.add_argument("--dp10", default=DEFAULT_DP10)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--promoter-window", type=int, default=2000)
    parser.add_argument("--nearest-window", type=int, default=100000)
    parser.add_argument("--log2fc-threshold", type=float, default=0.25)
    parser.add_argument("--strong-log2fc-threshold", type=float, default=0.58)
    parser.add_argument("--fpkm-pseudocount", type=float, default=0.1)
    return parser.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="", errors="replace")
    return path.open(newline="", errors="replace")


def safe_float(value: object) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return math.nan
    return number if math.isfinite(number) else math.nan


def safe_int(value: object) -> int:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return 0


def fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}" if math.isfinite(value) else ""
    if isinstance(value, int):
        return str(value)
    return str(value)


def mean(values: list[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    return sum(finite) / len(finite) if finite else math.nan


def log2_ratio(numerator: float, denominator: float, pseudocount: float) -> float:
    return math.log2((numerator + pseudocount) / (denominator + pseudocount))


def expression_direction(log2fc: float, threshold: float) -> str:
    if not math.isfinite(log2fc):
        return "not_available"
    if log2fc >= threshold:
        return "up"
    if log2fc <= -threshold:
        return "down"
    return "flat"


def normalize_chrom(chrom: str) -> str:
    value = str(chrom).strip()
    if not value:
        return ""
    if value.startswith("chr"):
        return value
    if value in {"M", "MT", "Mt", "mt"}:
        return "chrM"
    return f"chr{value}"


def chrom_key(chrom: str) -> tuple[int, int | str]:
    name = chrom.removeprefix("chr")
    if name.isdigit():
        return (0, int(name))
    special = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    return (0, special[name]) if name in special else (1, name)


def read_dmrs(path: Path) -> list[DMR]:
    dmrs: list[DMR] = []
    with open_text(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            dmrs.append(
                DMR(
                    dmr_id=row["dmr_id"],
                    chrom=row["chrom"],
                    start=int(row["start"]),
                    end=int(row["end"]),
                    direction=row.get("direction", ""),
                    n_cpgs=safe_int(row.get("n_cpgs")),
                    mean_delta_beta=safe_float(row.get("mean_delta_beta")),
                    row=row,
                )
            )
    return sorted(dmrs, key=lambda item: (chrom_key(item.chrom), item.start, item.end, item.dmr_id))


def read_gene_models(path: Path) -> dict[str, Gene]:
    genes: dict[str, Gene] = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene_id = row.get("GeneID", "").strip()
            if not gene_id:
                continue
            chrom = normalize_chrom(row.get("Chr", ""))
            start = safe_int(row.get("Start"))
            end = safe_int(row.get("End"))
            if not chrom or start <= 0 or end <= 0:
                continue
            if start > end:
                start, end = end, start
            symbol = row.get("GeneSymbol", "").strip() or gene_id
            description = unquote(row.get("Description", "").strip())
            strand = row.get("Strand", "").strip() or "."

            gene = genes.get(gene_id)
            if gene is None:
                gene = Gene(
                    gene_id=gene_id,
                    gene_symbol=symbol,
                    chrom=chrom,
                    start=start,
                    end=end,
                    strand=strand,
                    description=description,
                )
                genes[gene_id] = gene
            else:
                if gene.chrom == chrom:
                    gene.start = min(gene.start, start)
                    gene.end = max(gene.end, end)
                if gene.gene_symbol == gene.gene_id and symbol:
                    gene.gene_symbol = symbol
                if not gene.description and description:
                    gene.description = description
                if gene.strand == "." and strand:
                    gene.strand = strand

            gene.transcript_count += 1
            for sample in ALL_SAMPLES:
                gene.fpkm[sample] += safe_float(row.get(f"{sample}_FPKM"))
    return genes


def read_dp10_summary(path: Path) -> tuple[dict[str, dict[str, str]], dict[str, dict[str, str]]]:
    by_id: dict[str, dict[str, str]] = {}
    by_symbol: dict[str, dict[str, str]] = {}
    if not path.exists():
        return by_id, by_symbol
    with open_text(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            gene_id = row.get("gene_id", "")
            symbol = row.get("gene_symbol", "")
            if gene_id:
                by_id[gene_id] = row
            if symbol:
                by_symbol[symbol] = row
    return by_id, by_symbol


def read_gse268951_count_de(path: Path) -> tuple[dict[str, dict[str, str]], dict[str, dict[str, str]]]:
    by_id: dict[str, dict[str, str]] = {}
    by_symbol: dict[str, dict[str, str]] = {}
    if not path.exists():
        return by_id, by_symbol
    with open_text(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            gene_id = row.get("gene_id", "").strip()
            symbol = row.get("gene_symbol", "").strip() or row.get("GeneSymbol", "").strip()
            if "deseq2_log2fc" in row or "edger_log2fc" in row:
                log2fc = row.get("deseq2_log2fc") or row.get("edger_log2fc", "")
                fdr = row.get("deseq2_padj") or row.get("edger_fdr", "")
                pvalue = row.get("deseq2_pvalue") or row.get("edger_pvalue", "")
                row = {
                    **row,
                    "gene_symbol": symbol,
                    "gse268951_control_mean_count": row.get("control_mean_count", ""),
                    "gse268951_dnmt3l_oe_mean_count": row.get("dnmt3l_oe_mean_count", ""),
                    "gse268951_log2fc_oe_vs_control": log2fc,
                    "gse268951_primary_pvalue": pvalue,
                    "gse268951_primary_fdr": fdr,
                    "gse268951_edgeR_logCPM": row.get("edger_logCPM", ""),
                    "gse268951_effect_source": "public_GEO_processed_expected_count_DESeq2",
                    "gse268951_representative_transcript_id": "",
                    "gse268951_de_fdr05_fc2": (
                        "yes" if str(row.get("deseq2_padj05_abslog2fc1", "")).upper() == "TRUE" else "no"
                    ),
                }
            if gene_id:
                by_id[gene_id] = row
            if symbol:
                by_symbol[symbol] = row
    return by_id, by_symbol


def gene_expression_summary(
    gene: Gene,
    pseudocount: float,
    threshold: float,
    count_de_by_id: dict[str, dict[str, str]],
    count_de_by_symbol: dict[str, dict[str, str]],
) -> dict[str, object]:
    control_mean = mean([gene.fpkm[sample] for sample in CONTROL_SAMPLES])
    oe_mean = mean([gene.fpkm[sample] for sample in OE_SAMPLES])
    male_control_mean = mean([gene.fpkm[sample] for sample in MALE_CONTROL_SAMPLES])
    male_oe_mean = mean([gene.fpkm[sample] for sample in MALE_OE_SAMPLES])
    female_control_mean = mean([gene.fpkm[sample] for sample in FEMALE_CONTROL_SAMPLES])
    female_oe_mean = mean([gene.fpkm[sample] for sample in FEMALE_OE_SAMPLES])
    fpkm_log2fc = log2_ratio(oe_mean, control_mean, pseudocount)
    male_log2fc = log2_ratio(male_oe_mean, male_control_mean, pseudocount)
    female_log2fc = log2_ratio(female_oe_mean, female_control_mean, pseudocount)
    de = count_de_by_id.get(gene.gene_id) or count_de_by_symbol.get(gene.gene_symbol) or {}
    primary_log2fc = safe_float(de.get("gse268951_log2fc_oe_vs_control"))
    if not math.isfinite(primary_log2fc):
        primary_log2fc = fpkm_log2fc
    primary_fdr = safe_float(de.get("gse268951_primary_fdr"))
    primary_pvalue = safe_float(de.get("gse268951_primary_pvalue"))
    direction = de.get("gse268951_direction", "") or expression_direction(primary_log2fc, threshold)
    male_direction = expression_direction(male_log2fc, threshold)
    female_direction = expression_direction(female_log2fc, threshold)
    sex_consistent = (
        "yes"
        if male_direction == female_direction and male_direction in {"up", "down"}
        else "no"
    )
    return {
        "gse268951_control_mean_fpkm": control_mean,
        "gse268951_dnmt3l_oe_mean_fpkm": oe_mean,
        "gse268951_control_mean_count": safe_float(de.get("gse268951_control_mean_count")),
        "gse268951_dnmt3l_oe_mean_count": safe_float(de.get("gse268951_dnmt3l_oe_mean_count")),
        "gse268951_log2fc_oe_vs_control": primary_log2fc,
        "gse268951_fpkm_log2fc": fpkm_log2fc,
        "gse268951_primary_pvalue": primary_pvalue,
        "gse268951_primary_fdr": primary_fdr,
        "gse268951_edgeR_logCPM": safe_float(de.get("gse268951_edgeR_logCPM")),
        "gse268951_effect_source": de.get("gse268951_effect_source", "fpkm_fallback_log2_ratio"),
        "gse268951_representative_transcript_id": de.get("representative_transcript_id", ""),
        "gse268951_de_fdr05_fc2": de.get(
            "gse268951_de_fdr05_fc2",
            "yes" if math.isfinite(primary_fdr) and primary_fdr <= 0.05 and abs(primary_log2fc) >= 1.0 else "no",
        ),
        "gse268951_expression_direction": direction,
        "gse268951_male_log2fc_oe_vs_control": male_log2fc,
        "gse268951_female_log2fc_oe_vs_control": female_log2fc,
        "gse268951_male_expression_direction": male_direction,
        "gse268951_female_expression_direction": female_direction,
        "gse268951_sex_consistent_direction": sex_consistent,
        "gse268951_max_group_mean_fpkm": max(control_mean, oe_mean),
    }


def overlap(start_a: int, end_a: int, start_b: int, end_b: int) -> bool:
    return min(end_a, end_b) > max(start_a, start_b)


def signed_distance_to_tss(dmr: DMR, gene: Gene) -> int:
    midpoint = (dmr.start + dmr.end) // 2
    return midpoint - gene.tss


def absolute_interval_distance(start_a: int, end_a: int, start_b: int, end_b: int) -> int:
    if overlap(start_a, end_a, start_b, end_b):
        return 0
    if end_a <= start_b:
        return start_b - end_a
    return start_a - end_b


def relationship_for(dmr: DMR, gene: Gene, promoter_window: int) -> dict[str, object]:
    promoter_start = gene.tss - promoter_window
    promoter_end = gene.tss + promoter_window
    is_promoter = overlap(dmr.start, dmr.end, promoter_start, promoter_end)
    is_gene_body = overlap(dmr.start, dmr.end, gene.start, gene.end)
    signed_tss_distance = signed_distance_to_tss(dmr, gene)
    body_distance = absolute_interval_distance(dmr.start, dmr.end, gene.start, gene.end)
    if is_promoter:
        relationship = "promoter_2kb"
    elif is_gene_body:
        relationship = "gene_body"
    else:
        relationship = "nearest_100kb"
    return {
        "relationship": relationship,
        "promoter_2kb": "yes" if is_promoter else "no",
        "gene_body": "yes" if is_gene_body else "no",
        "signed_distance_to_tss": signed_tss_distance,
        "abs_distance_to_tss": abs(signed_tss_distance),
        "distance_to_gene_body": body_distance,
    }


def link_dmrs_to_genes(
    dmrs: list[DMR],
    genes: dict[str, Gene],
    promoter_window: int,
    nearest_window: int,
) -> list[dict[str, object]]:
    genes_by_chrom: dict[str, list[Gene]] = defaultdict(list)
    for gene in genes.values():
        genes_by_chrom[gene.chrom].append(gene)
    for chrom in genes_by_chrom:
        genes_by_chrom[chrom].sort(key=lambda gene: (gene.start, gene.end, gene.gene_id))

    links_by_key: dict[tuple[str, str], dict[str, object]] = {}
    for dmr in dmrs:
        chrom_genes = genes_by_chrom.get(dmr.chrom, [])
        nearest_gene = None
        nearest_distance = math.inf

        for gene in chrom_genes:
            distance = abs(signed_distance_to_tss(dmr, gene))
            if distance < nearest_distance:
                nearest_gene = gene
                nearest_distance = distance

            promoter_start = gene.tss - promoter_window
            promoter_end = gene.tss + promoter_window
            if not overlap(dmr.start, dmr.end, promoter_start, promoter_end) and not overlap(
                dmr.start, dmr.end, gene.start, gene.end
            ):
                continue
            links_by_key[(dmr.dmr_id, gene.gene_id)] = build_link_row(
                dmr,
                gene,
                relationship_for(dmr, gene, promoter_window),
            )

        if nearest_gene is not None and nearest_distance <= nearest_window:
            key = (dmr.dmr_id, nearest_gene.gene_id)
            nearest_relation = relationship_for(dmr, nearest_gene, promoter_window)
            existing = links_by_key.get(key)
            if existing is None:
                links_by_key[key] = build_link_row(dmr, nearest_gene, nearest_relation)
            else:
                existing["nearest_100kb"] = "yes"

    return sorted(
        links_by_key.values(),
        key=lambda row: (
            chrom_key(str(row["chrom"])),
            int(row["dmr_start"]),
            int(row["dmr_end"]),
            str(row["gene_symbol"]),
        ),
    )


def build_link_row(dmr: DMR, gene: Gene, relation: dict[str, object]) -> dict[str, object]:
    return {
        "dmr_id": dmr.dmr_id,
        "chrom": dmr.chrom,
        "dmr_start": dmr.start,
        "dmr_end": dmr.end,
        "dmr_direction": dmr.direction,
        "dmr_n_cpgs": dmr.n_cpgs,
        "dmr_mean_delta_beta": dmr.mean_delta_beta,
        "gene_id": gene.gene_id,
        "gene_symbol": gene.gene_symbol,
        "gene_chrom": gene.chrom,
        "gene_start": gene.start,
        "gene_end": gene.end,
        "gene_strand": gene.strand,
        "gene_tss": gene.tss,
        "gene_length": gene.length,
        "transcript_count": gene.transcript_count,
        "description": gene.description,
        "nearest_100kb": "yes" if relation["relationship"] == "nearest_100kb" else "no",
        **relation,
    }


def add_expression_to_link_rows(
    links: list[dict[str, object]],
    genes: dict[str, Gene],
    count_de_by_id: dict[str, dict[str, str]],
    count_de_by_symbol: dict[str, dict[str, str]],
    dp10_by_id: dict[str, dict[str, str]],
    dp10_by_symbol: dict[str, dict[str, str]],
    pseudocount: float,
    threshold: float,
) -> list[dict[str, object]]:
    expression_cache = {
        gene_id: gene_expression_summary(gene, pseudocount, threshold, count_de_by_id, count_de_by_symbol)
        for gene_id, gene in genes.items()
    }
    for row in links:
        gene_id = str(row["gene_id"])
        symbol = str(row["gene_symbol"])
        row.update(expression_cache.get(gene_id, {}))
        dp10 = dp10_by_id.get(gene_id) or dp10_by_symbol.get(symbol) or {}
        add_dp10_fields(row, dp10)
    return links


def add_dp10_fields(row: dict[str, object], dp10: dict[str, str]) -> None:
    row["dp10_wt_mean_cpm"] = safe_float(dp10.get("wt_mean_cpm"))
    row["dp10_dp10_mean_cpm"] = safe_float(dp10.get("dp10_mean_cpm"))
    row["dp10_delta_log2_cpm_plus_1"] = safe_float(dp10.get("delta_log2_cpm_plus_1"))
    row["dp10_delta_detected_fraction"] = safe_float(dp10.get("delta_detected_fraction"))
    row["dp10_permutation_pvalue_log2_cpm"] = safe_float(dp10.get("permutation_pvalue_log2_cpm"))
    row["dp10_expression_direction"] = dp10.get("direction", "not_available") if dp10 else "not_available"


def rank_genes(
    links: list[dict[str, object]],
    genes: dict[str, Gene],
    count_de_by_id: dict[str, dict[str, str]],
    count_de_by_symbol: dict[str, dict[str, str]],
    dp10_by_id: dict[str, dict[str, str]],
    dp10_by_symbol: dict[str, dict[str, str]],
    pseudocount: float,
    threshold: float,
    strong_threshold: float,
) -> list[dict[str, object]]:
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in links:
        grouped[str(row["gene_id"])].append(row)

    ranked = []
    for gene_id, rows in grouped.items():
        gene = genes[gene_id]
        expr = gene_expression_summary(gene, pseudocount, threshold, count_de_by_id, count_de_by_symbol)
        dp10 = dp10_by_id.get(gene_id) or dp10_by_symbol.get(gene.gene_symbol) or {}

        unique_dmrs = {str(row["dmr_id"]): row for row in rows}
        n_promoter = len({str(row["dmr_id"]) for row in rows if row.get("promoter_2kb") == "yes"})
        n_body = len({str(row["dmr_id"]) for row in rows if row.get("gene_body") == "yes"})
        n_nearest = len({str(row["dmr_id"]) for row in rows if row.get("nearest_100kb") == "yes"})
        deltas = [safe_float(row["dmr_mean_delta_beta"]) for row in unique_dmrs.values()]
        abs_tss_distances = [safe_int(row["abs_distance_to_tss"]) for row in rows]
        dmr_ids = sorted(unique_dmrs)
        directions = [str(row["dmr_direction"]) for row in unique_dmrs.values()]

        relationship_score = min(9, 3 * n_promoter) + min(4, 2 * n_body) + min(2, n_nearest)
        expression_score = 0
        if safe_float(expr["gse268951_max_group_mean_fpkm"]) >= 1:
            expression_score += 1
        if abs(safe_float(expr["gse268951_log2fc_oe_vs_control"])) >= threshold:
            expression_score += 1
        if abs(safe_float(expr["gse268951_log2fc_oe_vs_control"])) >= strong_threshold:
            expression_score += 1
        if expr.get("gse268951_de_fdr05_fc2") == "yes":
            expression_score += 1
        if expr["gse268951_sex_consistent_direction"] == "yes":
            expression_score += 1

        dp10_delta = safe_float(dp10.get("delta_log2_cpm_plus_1"))
        dp10_direction = dp10.get("direction", "not_available") if dp10 else "not_available"
        dp10_score = 0
        if abs(dp10_delta) >= threshold or dp10_direction in {"up", "down"}:
            dp10_score += 1
        if dp10_direction in {"up", "down"} and dp10_direction == expr["gse268951_expression_direction"]:
            dp10_score += 1

        row = {
            "gene_id": gene_id,
            "gene_symbol": gene.gene_symbol,
            "chrom": gene.chrom,
            "start": gene.start,
            "end": gene.end,
            "strand": gene.strand,
            "tss": gene.tss,
            "gene_length": gene.length,
            "transcript_count": gene.transcript_count,
            "description": gene.description,
            "n_triple_hit_dmrs": len(unique_dmrs),
            "n_promoter_2kb_dmrs": n_promoter,
            "n_gene_body_dmrs": n_body,
            "n_nearest_100kb_dmrs": n_nearest,
            "n_hyper_dmrs": sum(1 for direction in directions if direction == "hyper"),
            "n_hypo_dmrs": sum(1 for direction in directions if direction == "hypo"),
            "min_abs_distance_to_tss": min(abs_tss_distances) if abs_tss_distances else "",
            "mean_dmr_delta_beta": mean(deltas),
            "mean_abs_dmr_delta_beta": mean([abs(delta) for delta in deltas]),
            "max_abs_dmr_delta_beta": max([abs(delta) for delta in deltas], default=math.nan),
            "dmr_ids": ",".join(dmr_ids[:20]) + ("..." if len(dmr_ids) > 20 else ""),
            "relationship_score": relationship_score,
            "expression_score": expression_score,
            "dp10_support_score": dp10_score,
            "evidence_score": relationship_score + expression_score + dp10_score,
            **expr,
        }
        add_dp10_fields(row, dp10)
        ranked.append(row)

    ranked.sort(
        key=lambda row: (
            -safe_float(row["evidence_score"]),
            -safe_float(row["n_promoter_2kb_dmrs"]),
            -safe_float(row["n_gene_body_dmrs"]),
            -abs(safe_float(row["gse268951_log2fc_oe_vs_control"])),
            str(row["gene_symbol"]),
        )
    )
    return ranked


def target_gene_context(
    genes: dict[str, Gene],
    ranked_by_symbol: dict[str, dict[str, object]],
    count_de_by_id: dict[str, dict[str, str]],
    count_de_by_symbol: dict[str, dict[str, str]],
    dp10_by_id: dict[str, dict[str, str]],
    dp10_by_symbol: dict[str, dict[str, str]],
    pseudocount: float,
    threshold: float,
) -> list[dict[str, object]]:
    rows = []
    genes_by_symbol = {gene.gene_symbol: gene for gene in genes.values()}
    for symbol in TARGET_GENES:
        gene = genes_by_symbol.get(symbol)
        if gene is None:
            rows.append({"gene_symbol": symbol, "found_in_gse268951_annotation": "no"})
            continue
        expr = gene_expression_summary(gene, pseudocount, threshold, count_de_by_id, count_de_by_symbol)
        ranked = ranked_by_symbol.get(symbol, {})
        dp10 = dp10_by_id.get(gene.gene_id) or dp10_by_symbol.get(symbol) or {}
        row = {
            "gene_symbol": symbol,
            "gene_id": gene.gene_id,
            "found_in_gse268951_annotation": "yes",
            "chrom": gene.chrom,
            "start": gene.start,
            "end": gene.end,
            "strand": gene.strand,
            "n_triple_hit_dmrs": ranked.get("n_triple_hit_dmrs", 0),
            "n_promoter_2kb_dmrs": ranked.get("n_promoter_2kb_dmrs", 0),
            "n_gene_body_dmrs": ranked.get("n_gene_body_dmrs", 0),
            "n_nearest_100kb_dmrs": ranked.get("n_nearest_100kb_dmrs", 0),
            **expr,
        }
        add_dp10_fields(row, dp10)
        rows.append(row)
    return rows


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)

    dmrs = read_dmrs(Path(args.dmrs))
    genes = read_gene_models(Path(args.rna))
    count_de_by_id, count_de_by_symbol = read_gse268951_count_de(Path(args.count_de))
    dp10_by_id, dp10_by_symbol = read_dp10_summary(Path(args.dp10))

    links = link_dmrs_to_genes(dmrs, genes, args.promoter_window, args.nearest_window)
    add_expression_to_link_rows(
        links,
        genes,
        count_de_by_id,
        count_de_by_symbol,
        dp10_by_id,
        dp10_by_symbol,
        args.fpkm_pseudocount,
        args.log2fc_threshold,
    )
    ranked = rank_genes(
        links,
        genes,
        count_de_by_id,
        count_de_by_symbol,
        dp10_by_id,
        dp10_by_symbol,
        args.fpkm_pseudocount,
        args.log2fc_threshold,
        args.strong_log2fc_threshold,
    )
    ranked_by_symbol = {str(row["gene_symbol"]): row for row in ranked}
    targets = target_gene_context(
        genes,
        ranked_by_symbol,
        count_de_by_id,
        count_de_by_symbol,
        dp10_by_id,
        dp10_by_symbol,
        args.fpkm_pseudocount,
        args.log2fc_threshold,
    )

    link_fields = [
        "dmr_id",
        "chrom",
        "dmr_start",
        "dmr_end",
        "dmr_direction",
        "dmr_n_cpgs",
        "dmr_mean_delta_beta",
        "gene_id",
        "gene_symbol",
        "gene_chrom",
        "gene_start",
        "gene_end",
        "gene_strand",
        "gene_tss",
        "gene_length",
        "transcript_count",
        "relationship",
        "promoter_2kb",
        "gene_body",
        "nearest_100kb",
        "signed_distance_to_tss",
        "abs_distance_to_tss",
        "distance_to_gene_body",
        "gse268951_control_mean_fpkm",
        "gse268951_dnmt3l_oe_mean_fpkm",
        "gse268951_control_mean_count",
        "gse268951_dnmt3l_oe_mean_count",
        "gse268951_log2fc_oe_vs_control",
        "gse268951_fpkm_log2fc",
        "gse268951_primary_pvalue",
        "gse268951_primary_fdr",
        "gse268951_edgeR_logCPM",
        "gse268951_effect_source",
        "gse268951_representative_transcript_id",
        "gse268951_de_fdr05_fc2",
        "gse268951_expression_direction",
        "gse268951_male_log2fc_oe_vs_control",
        "gse268951_female_log2fc_oe_vs_control",
        "gse268951_sex_consistent_direction",
        "dp10_wt_mean_cpm",
        "dp10_dp10_mean_cpm",
        "dp10_delta_log2_cpm_plus_1",
        "dp10_delta_detected_fraction",
        "dp10_permutation_pvalue_log2_cpm",
        "dp10_expression_direction",
        "description",
    ]
    ranking_fields = [
        "gene_id",
        "gene_symbol",
        "chrom",
        "start",
        "end",
        "strand",
        "tss",
        "gene_length",
        "transcript_count",
        "n_triple_hit_dmrs",
        "n_promoter_2kb_dmrs",
        "n_gene_body_dmrs",
        "n_nearest_100kb_dmrs",
        "n_hyper_dmrs",
        "n_hypo_dmrs",
        "min_abs_distance_to_tss",
        "mean_dmr_delta_beta",
        "mean_abs_dmr_delta_beta",
        "max_abs_dmr_delta_beta",
        "relationship_score",
        "expression_score",
        "dp10_support_score",
        "evidence_score",
        "gse268951_control_mean_fpkm",
        "gse268951_dnmt3l_oe_mean_fpkm",
        "gse268951_control_mean_count",
        "gse268951_dnmt3l_oe_mean_count",
        "gse268951_log2fc_oe_vs_control",
        "gse268951_fpkm_log2fc",
        "gse268951_primary_pvalue",
        "gse268951_primary_fdr",
        "gse268951_edgeR_logCPM",
        "gse268951_effect_source",
        "gse268951_representative_transcript_id",
        "gse268951_de_fdr05_fc2",
        "gse268951_expression_direction",
        "gse268951_male_log2fc_oe_vs_control",
        "gse268951_female_log2fc_oe_vs_control",
        "gse268951_male_expression_direction",
        "gse268951_female_expression_direction",
        "gse268951_sex_consistent_direction",
        "dp10_wt_mean_cpm",
        "dp10_dp10_mean_cpm",
        "dp10_delta_log2_cpm_plus_1",
        "dp10_delta_detected_fraction",
        "dp10_permutation_pvalue_log2_cpm",
        "dp10_expression_direction",
        "dmr_ids",
        "description",
    ]
    target_fields = [
        "gene_symbol",
        "gene_id",
        "found_in_gse268951_annotation",
        "chrom",
        "start",
        "end",
        "strand",
        "n_triple_hit_dmrs",
        "n_promoter_2kb_dmrs",
        "n_gene_body_dmrs",
        "n_nearest_100kb_dmrs",
        "gse268951_control_mean_fpkm",
        "gse268951_dnmt3l_oe_mean_fpkm",
        "gse268951_control_mean_count",
        "gse268951_dnmt3l_oe_mean_count",
        "gse268951_log2fc_oe_vs_control",
        "gse268951_fpkm_log2fc",
        "gse268951_primary_pvalue",
        "gse268951_primary_fdr",
        "gse268951_edgeR_logCPM",
        "gse268951_effect_source",
        "gse268951_representative_transcript_id",
        "gse268951_de_fdr05_fc2",
        "gse268951_expression_direction",
        "gse268951_male_log2fc_oe_vs_control",
        "gse268951_female_log2fc_oe_vs_control",
        "gse268951_sex_consistent_direction",
        "dp10_wt_mean_cpm",
        "dp10_dp10_mean_cpm",
        "dp10_delta_log2_cpm_plus_1",
        "dp10_delta_detected_fraction",
        "dp10_permutation_pvalue_log2_cpm",
        "dp10_expression_direction",
    ]

    write_tsv(outdir / "mouse_triple_hit_dmr_gene_links.tsv", links, link_fields)
    write_tsv(outdir / "mouse_triple_hit_candidate_gene_ranking.tsv", ranked, ranking_fields)
    write_tsv(outdir / "mouse_target_gene_context.tsv", targets, target_fields)

    linked_dmr_ids = {str(row["dmr_id"]) for row in links}
    promoter_linked_dmr_ids = {str(row["dmr_id"]) for row in links if row.get("promoter_2kb") == "yes"}
    gene_body_linked_dmr_ids = {str(row["dmr_id"]) for row in links if row.get("gene_body") == "yes"}
    same_direction_rows = [
        row
        for row in ranked
        if row.get("gse268951_expression_direction") in {"up", "down"}
        and row.get("gse268951_expression_direction") == row.get("dp10_expression_direction")
    ]
    summary_rows = [
        {"metric": "input_triple_hit_dmrs", "value": len(dmrs)},
        {"metric": "gse268951_genes_with_annotation", "value": len(genes)},
        {"metric": "triple_hit_dmrs_linked_to_gene_within_100kb", "value": len(linked_dmr_ids)},
        {"metric": "triple_hit_dmrs_without_gene_within_100kb", "value": len(dmrs) - len(linked_dmr_ids)},
        {"metric": "triple_hit_dmrs_with_promoter_2kb_link", "value": len(promoter_linked_dmr_ids)},
        {"metric": "triple_hit_dmrs_with_gene_body_link", "value": len(gene_body_linked_dmr_ids)},
        {"metric": "dmr_gene_links", "value": len(links)},
        {"metric": "ranked_candidate_genes", "value": len(ranked)},
        {
            "metric": "ranked_genes_with_promoter_2kb_dmr",
            "value": sum(1 for row in ranked if safe_int(row["n_promoter_2kb_dmrs"]) > 0),
        },
        {
            "metric": "ranked_genes_with_gene_body_dmr",
            "value": sum(1 for row in ranked if safe_int(row["n_gene_body_dmrs"]) > 0),
        },
        {
            "metric": "ranked_genes_with_gse268951_abs_log2fc_ge_0.25",
            "value": sum(
                1 for row in ranked if abs(safe_float(row["gse268951_log2fc_oe_vs_control"])) >= 0.25
            ),
        },
        {
            "metric": "ranked_genes_with_dp10_abs_delta_ge_0.25",
            "value": sum(1 for row in ranked if abs(safe_float(row["dp10_delta_log2_cpm_plus_1"])) >= 0.25),
        },
        {
            "metric": "ranked_genes_with_same_nonflat_gse268951_and_dp10_direction",
            "value": len(same_direction_rows),
        },
    ]
    write_tsv(outdir / "mouse_triple_hit_gene_set_summary.tsv", summary_rows, ["metric", "value"])


if __name__ == "__main__":
    main()
