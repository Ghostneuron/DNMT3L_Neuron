#!/usr/bin/env python3
"""Analyze neurogenesis function and Dnmt3l-vs-Dnmt3a transcriptomic overlap."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from bisect import bisect_left, bisect_right
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


DEFAULT_GSE268951 = "data/geo/GSE268951/supplementary/GSE268951_all.fpkm_mRNA_anno.txt.gz"
DEFAULT_GSE268951_COUNT_DE = "results/mouse_gse268951_count_de/gse268951_deseq2_edger_count_de_results.tsv"
DEFAULT_GSE141587_DE = "data/geo/GSE141587/supplementary/GSE141587_RNAseq_P39_DEgenes_edgeR_qlfTest_Dnmt3aCKO_vs_Control.txt.gz"
DEFAULT_DNMT3L_DMRS = "results/mouse_candidate_dmrs_sex_consistent/GSE271949_candidate_dmrs.tsv"
DEFAULT_TRIPLE_LINKS = "results/integrated_candidate_genes/mouse_triple_hit_dmr_gene_links.tsv"
DEFAULT_SUMMARY = "results/manuscript_summary/final_candidate_manuscript_summary.tsv"
DEFAULT_OUTDIR = "results/neurogenesis_mechanism"
DEFAULT_FIGDIR = "figures/neurogenesis_mechanism"
DEFAULT_NOTE = "notes/neurogenesis_and_dnmt3a_independence_first_pass.md"


CONTROL_COLS = ["A-1_FPKM", "A-2_FPKM", "A-3_FPKM", "B-1_FPKM", "B-2_FPKM", "B-3_FPKM"]
OE_COLS = ["C-1_FPKM", "C-2_FPKM", "C-3_FPKM", "D-1_FPKM", "D-2_FPKM", "D-3_FPKM"]
PSEUDOCOUNT = 0.1


NEUROGENESIS_KEYWORDS = {
    "core_neurogenesis": [
        "neurogenesis",
        "generation of neurons",
        "neuron differentiation",
        "neuronal differentiation",
        "neuron development",
        "neural precursor",
        "neural progenitor",
        "neural stem",
    ],
    "brain_region_patterning": [
        "brain development",
        "forebrain development",
        "telencephalon",
        "cerebral cortex",
        "cortical",
        "diencephalon",
        "midbrain",
        "hindbrain",
        "neural tube",
        "pattern specification",
    ],
    "neurite_synapse_migration": [
        "axon",
        "dendrite",
        "synapse",
        "synaptic",
        "neuron projection",
        "cell migration",
        "neuronal migration",
        "axon guidance",
    ],
    "gliogenesis_or_glial_context": [
        "gliogenesis",
        "glial",
        "astrocyte",
        "oligodendrocyte",
    ],
    "neurogenesis_signaling": [
        "wnt",
        "notch",
        "bmp",
        "fibroblast growth factor",
        "fgf",
        "retinoic acid",
        "tgf",
        "mapk",
        "pi3k",
    ],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gse268951", default=DEFAULT_GSE268951)
    parser.add_argument("--gse268951-count-de", default=DEFAULT_GSE268951_COUNT_DE)
    parser.add_argument("--gse141587-de", default=DEFAULT_GSE141587_DE)
    parser.add_argument("--dnmt3l-dmrs", default=DEFAULT_DNMT3L_DMRS)
    parser.add_argument("--triple-links", default=DEFAULT_TRIPLE_LINKS)
    parser.add_argument("--summary", default=DEFAULT_SUMMARY)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--figdir", default=DEFAULT_FIGDIR)
    parser.add_argument("--note", default=DEFAULT_NOTE)
    parser.add_argument("--pvalue-cutoff", type=float, default=0.05)
    parser.add_argument("--abs-log2fc-cutoff", type=float, default=1.0)
    return parser.parse_args()


def clean_key(value: str | None) -> str:
    return (value or "").strip().replace("\r", "")


def safe_float(value: object) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return math.nan
    return number if math.isfinite(number) else math.nan


def fmt(value: object, digits: int = 6) -> str:
    number = safe_float(value)
    if not math.isfinite(number):
        return ""
    return f"{number:.{digits}g}"


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        reader.fieldnames = [clean_key(field) for field in reader.fieldnames or []]
        return [{clean_key(k): (v or "").strip() for k, v in row.items()} for row in reader]


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def bh_fdr(pvalues: list[float]) -> list[float]:
    n = len(pvalues)
    order = sorted(range(n), key=lambda i: (math.inf if not math.isfinite(pvalues[i]) else pvalues[i]))
    qvalues = [math.nan] * n
    running = 1.0
    for rank, idx in reversed(list(enumerate(order, start=1))):
        p = pvalues[idx]
        if not math.isfinite(p):
            qvalues[idx] = math.nan
            continue
        running = min(running, p * n / rank)
        qvalues[idx] = min(running, 1.0)
    return qvalues


def split_go_terms(go_text: str) -> list[str]:
    terms: list[str] = []
    for item in (go_text or "").split(";"):
        item = item.strip()
        if not item:
            continue
        if ":" in item:
            # Keep the readable GO name after the GO ID.
            parts = item.split(":", 2)
            if len(parts) == 3 and parts[0] == "GO":
                item = parts[2]
        terms.append(item)
    return terms


def normalize_chrom(chrom: object) -> str:
    value = str(chrom or "").strip()
    if not value or value == "-":
        return ""
    return value if value.startswith("chr") else f"chr{value}"


def neurogenesis_flags(go_terms: list[str], description: str = "") -> dict[str, object]:
    haystack = " ; ".join(go_terms + [description]).lower()
    flags: dict[str, object] = {}
    matched_terms: list[str] = []
    for category, keywords in NEUROGENESIS_KEYWORDS.items():
        category_terms = [term for term in go_terms if any(keyword in term.lower() for keyword in keywords)]
        category_hit = bool(category_terms) or any(keyword in haystack for keyword in keywords)
        flags[category] = "yes" if category_hit else "no"
        matched_terms.extend(category_terms)
    unique_terms = sorted(set(matched_terms))
    flags["neurogenesis_category_count"] = sum(1 for category in NEUROGENESIS_KEYWORDS if flags[category] == "yes")
    flags["neurogenesis_related"] = "yes" if flags["neurogenesis_category_count"] else "no"
    flags["matched_neurogenesis_go_terms"] = "; ".join(unique_terms[:20])
    return flags


def load_gse268951_gene_expression(path: Path) -> tuple[list[dict[str, object]], dict[str, dict[str, object]]]:
    grouped: dict[str, dict[str, object]] = {}
    for row in read_tsv(path):
        gene_id = row.get("GeneID", "")
        symbol = row.get("GeneSymbol", "")
        if not gene_id or not symbol or symbol == "-":
            continue
        rec = grouped.setdefault(
            gene_id,
            {
                "gene_id": gene_id,
                "gene_symbol": symbol,
                "chrom": normalize_chrom(row.get("Chr", "")),
                "start": row.get("Start", ""),
                "end": row.get("End", ""),
                "strand": row.get("Strand", ""),
                "description": row.get("Description", ""),
                "go_bp_terms": set(),
                "values": {col: 0.0 for col in CONTROL_COLS + OE_COLS},
            },
        )
        for col in CONTROL_COLS + OE_COLS:
            rec["values"][col] += safe_float(row.get(col, 0)) if math.isfinite(safe_float(row.get(col, 0))) else 0.0
        for term in split_go_terms(row.get("GO_BP", "")):
            rec["go_bp_terms"].add(term)

    rows: list[dict[str, object]] = []
    pvalues: list[float] = []
    for rec in grouped.values():
        values = rec["values"]
        control = np.array([values[col] for col in CONTROL_COLS], dtype=float)
        oe = np.array([values[col] for col in OE_COLS], dtype=float)
        control_log = np.log2(control + PSEUDOCOUNT)
        oe_log = np.log2(oe + PSEUDOCOUNT)
        control_mean = float(np.mean(control))
        oe_mean = float(np.mean(oe))
        log2fc = float(np.log2((oe_mean + PSEUDOCOUNT) / (control_mean + PSEUDOCOUNT)))
        ttest = stats.ttest_ind(oe_log, control_log, equal_var=False, nan_policy="omit")
        pvalue = float(ttest.pvalue) if math.isfinite(float(ttest.pvalue)) else math.nan
        pvalues.append(pvalue)
        go_terms = sorted(rec["go_bp_terms"])
        flags = neurogenesis_flags(go_terms, str(rec["description"]))
        row = {
            "gene_id": rec["gene_id"],
            "gene_symbol": rec["gene_symbol"],
            "chrom": rec["chrom"],
            "start": rec["start"],
            "end": rec["end"],
            "strand": rec["strand"],
            "description": rec["description"],
            "gse268951_control_mean_fpkm": fmt(control_mean),
            "gse268951_dnmt3l_oe_mean_fpkm": fmt(oe_mean),
            "gse268951_log2fc_oe_vs_control": fmt(log2fc),
            "gse268951_effect_source": "processed_FPKM_Welch_fallback",
            "gse268951_primary_pvalue": fmt(pvalue),
            "gse268951_primary_fdr": "",
            "gse268951_welch_pvalue": fmt(pvalue),
            "gse268951_fpkm_log2fc": fmt(log2fc),
            "gse268951_direction": "up" if log2fc >= 0.25 else ("down" if log2fc <= -0.25 else "flat"),
            "go_bp_terms": "; ".join(go_terms[:80]),
            **flags,
        }
        rows.append(row)

    qvalues = bh_fdr(pvalues)
    for row, qvalue in zip(rows, qvalues):
        log2fc = safe_float(row["gse268951_log2fc_oe_vs_control"])
        pvalue = safe_float(row["gse268951_welch_pvalue"])
        row["gse268951_bh_fdr"] = fmt(qvalue)
        row["gse268951_primary_fdr"] = fmt(qvalue)
        row["gse268951_de_p05_fc2"] = "yes" if pvalue <= 0.05 and abs(log2fc) >= 1.0 else "no"
        row["gse268951_de_fdr05_fc2"] = "yes" if qvalue <= 0.05 and abs(log2fc) >= 1.0 else "no"

    by_symbol: dict[str, dict[str, object]] = {}
    for row in sorted(rows, key=lambda r: abs(safe_float(r["gse268951_log2fc_oe_vs_control"])), reverse=True):
        by_symbol.setdefault(str(row["gene_symbol"]), row)
    return rows, by_symbol


def load_gse268951_count_de(path: Path) -> tuple[list[dict[str, object]], dict[str, dict[str, object]]]:
    rows: list[dict[str, object]] = []
    for raw in read_tsv(path):
        gene_id = raw.get("gene_id", "")
        symbol = raw.get("gene_symbol") or raw.get("GeneSymbol", "")
        if not gene_id or not symbol or symbol == "-":
            continue
        if raw.get("gse268951_log2fc_oe_vs_control", "") != "":
            log2fc = safe_float(raw.get("gse268951_log2fc_oe_vs_control"))
            pvalue = safe_float(raw.get("gse268951_primary_pvalue"))
            fdr = safe_float(raw.get("gse268951_primary_fdr"))
            fpkm_log2fc = safe_float(raw.get("gse268951_fpkm_log2fc"))
            description = raw.get("description", "")
            go_text = raw.get("go_bp_terms", "")
            row = {
                "gene_id": gene_id,
                "representative_transcript_id": raw.get("representative_transcript_id", ""),
                "transcript_count_for_gene": raw.get("transcript_count_for_gene", ""),
                "gene_symbol": symbol,
                "chrom": normalize_chrom(raw.get("chrom", "")),
                "start": raw.get("start", ""),
                "end": raw.get("end", ""),
                "strand": raw.get("strand", ""),
                "description": description,
                "gse268951_control_mean_count": raw.get("gse268951_control_mean_count", ""),
                "gse268951_dnmt3l_oe_mean_count": raw.get("gse268951_dnmt3l_oe_mean_count", ""),
                "gse268951_control_mean_fpkm": raw.get("gse268951_control_mean_fpkm", ""),
                "gse268951_dnmt3l_oe_mean_fpkm": raw.get("gse268951_dnmt3l_oe_mean_fpkm", ""),
                "gse268951_log2fc_oe_vs_control": fmt(log2fc),
                "gse268951_effect_source": raw.get("gse268951_effect_source", "validated_edgeR_best_transcript"),
                "gse268951_primary_pvalue": fmt(pvalue),
                "gse268951_primary_fdr": fmt(fdr),
                "gse268951_welch_pvalue": "",
                "gse268951_bh_fdr": fmt(fdr),
                "gse268951_fpkm_log2fc": fmt(fpkm_log2fc),
                "gse268951_edgeR_logFC": raw.get("gse268951_edgeR_logFC", fmt(log2fc)),
                "gse268951_edgeR_logCPM": raw.get("gse268951_edgeR_logCPM", ""),
                "gse268951_edgeR_pvalue": raw.get("gse268951_edgeR_pvalue", fmt(pvalue)),
                "gse268951_edgeR_fdr": raw.get("gse268951_edgeR_fdr", fmt(fdr)),
                "gse268951_direction": raw.get("gse268951_direction") or ("up" if log2fc >= 0.25 else ("down" if log2fc <= -0.25 else "flat")),
                "gse268951_de_p05_fc2": raw.get("gse268951_de_p05_fc2") or ("yes" if pvalue <= 0.05 and abs(log2fc) >= 1.0 else "no"),
                "gse268951_de_fdr05_fc2": raw.get("gse268951_de_fdr05_fc2") or ("yes" if fdr <= 0.05 and abs(log2fc) >= 1.0 else "no"),
            }
        else:
            log2fc = safe_float(raw.get("deseq2_log2fc"))
            pvalue = safe_float(raw.get("deseq2_pvalue"))
            fdr = safe_float(raw.get("deseq2_padj"))
            fpkm_log2fc = safe_float(raw.get("fpkm_log2fc"))
            description = raw.get("Description", "")
            go_text = raw.get("GO_BP", "")
            row = {
                "gene_id": gene_id,
                "representative_transcript_id": "",
                "transcript_count_for_gene": "",
                "gene_symbol": symbol,
                "chrom": normalize_chrom(raw.get("Chr", "")),
                "start": raw.get("Start", ""),
                "end": raw.get("End", ""),
                "strand": raw.get("Strand", ""),
                "description": description,
                "gse268951_control_mean_count": fmt(raw.get("control_mean_count")),
                "gse268951_dnmt3l_oe_mean_count": fmt(raw.get("dnmt3l_oe_mean_count")),
                "gse268951_control_mean_fpkm": fmt(raw.get("control_mean_fpkm")),
                "gse268951_dnmt3l_oe_mean_fpkm": fmt(raw.get("dnmt3l_oe_mean_fpkm")),
                "gse268951_log2fc_oe_vs_control": fmt(log2fc),
                "gse268951_effect_source": "processed_expected_count_DESeq2_sex_adjusted",
                "gse268951_primary_pvalue": fmt(pvalue),
                "gse268951_primary_fdr": fmt(fdr),
                "gse268951_welch_pvalue": fmt(pvalue),
                "gse268951_bh_fdr": fmt(fdr),
                "gse268951_fpkm_log2fc": fmt(fpkm_log2fc),
                "gse268951_count_deseq2_log2fc": fmt(log2fc),
                "gse268951_count_deseq2_pvalue": fmt(pvalue),
                "gse268951_count_deseq2_fdr": fmt(fdr),
                "gse268951_count_edger_log2fc": fmt(raw.get("edger_log2fc")),
                "gse268951_count_edger_fdr": fmt(raw.get("edger_fdr")),
                "gse268951_direction": "up" if log2fc >= 0.25 else ("down" if log2fc <= -0.25 else "flat"),
                "gse268951_de_p05_fc2": "yes" if pvalue <= 0.05 and abs(log2fc) >= 1.0 else "no",
                "gse268951_de_fdr05_fc2": "yes" if fdr <= 0.05 and abs(log2fc) >= 1.0 else "no",
            }
        go_terms = split_go_terms(go_text)
        flags = neurogenesis_flags(go_terms, description)
        row["go_bp_terms"] = "; ".join(go_terms[:80])
        row.update(flags)
        rows.append(row)

    by_symbol: dict[str, dict[str, object]] = {}
    for row in sorted(rows, key=lambda r: abs(safe_float(r["gse268951_log2fc_oe_vs_control"])), reverse=True):
        by_symbol.setdefault(str(row["gene_symbol"]), row)
    return rows, by_symbol


def load_dnmt3a_de(path: Path) -> dict[str, dict[str, object]]:
    rows: dict[str, dict[str, object]] = {}
    for row in read_tsv(path):
        symbol = row.get("geneName", "")
        if not symbol:
            continue
        logfc = safe_float(row.get("logFC"))
        fdr = safe_float(row.get("FDR"))
        rows[symbol] = {
            "dnmt3a_gene_id": row.get("geneID", ""),
            "dnmt3a_logfc_cko_vs_control": fmt(logfc),
            "dnmt3a_fdr": fmt(fdr),
            "dnmt3a_direction": "up" if logfc >= 0.25 else ("down" if logfc <= -0.25 else "flat"),
            "dnmt3a_de_fdr05": "yes" if fdr <= 0.05 else "no",
            "dnmt3a_de_fdr05_fc2": "yes" if fdr <= 0.05 and abs(logfc) >= 1 else "no",
        }
    return rows


def load_dmrs(path: Path) -> dict[str, list[dict[str, object]]]:
    by_chrom: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in read_tsv(path):
        chrom = normalize_chrom(row.get("chrom", ""))
        start = int(float(row.get("start", 0)))
        end = int(float(row.get("end", 0)))
        by_chrom[chrom].append(
            {
                "dmr_id": row.get("dmr_id", ""),
                "start": start,
                "end": end,
                "direction": row.get("direction", ""),
                "delta_beta": safe_float(row.get("mean_delta_beta")),
            }
        )
    for chrom in by_chrom:
        by_chrom[chrom].sort(key=lambda item: item["start"])
    return by_chrom


def dmr_support_for_gene(gene: dict[str, object], dmrs_by_chrom: dict[str, list[dict[str, object]]]) -> dict[str, object]:
    chrom = normalize_chrom(gene.get("chrom", ""))
    start = int(float(gene.get("start") or 0))
    end = int(float(gene.get("end") or 0))
    strand = str(gene.get("strand", ""))
    if not chrom or not start or not end:
        return {
            "nearby_dnmt3l_dmr_count_100kb": 0,
            "gene_body_dnmt3l_dmr_count": 0,
            "promoter_2kb_dnmt3l_dmr_count": 0,
            "nearby_dnmt3l_dmr_directions": "",
        }
    tss = end if strand == "-" else start
    promoter_start = max(0, tss - 2000)
    promoter_end = tss + 2000
    query_start = max(0, start - 100000)
    query_end = end + 100000
    dmrs = dmrs_by_chrom.get(chrom, [])
    starts = [int(dmr["start"]) for dmr in dmrs]
    left = bisect_left(starts, query_start)
    right = bisect_right(starts, query_end)
    nearby = []
    gene_body = []
    promoter = []
    for dmr in dmrs[left:right]:
        d_start = int(dmr["start"])
        d_end = int(dmr["end"])
        if d_end < query_start or d_start > query_end:
            continue
        nearby.append(dmr)
        if d_start < end and d_end > start:
            gene_body.append(dmr)
        if d_start < promoter_end and d_end > promoter_start:
            promoter.append(dmr)
    directions = sorted(set(str(dmr["direction"]) for dmr in nearby if dmr.get("direction")))
    return {
        "nearby_dnmt3l_dmr_count_100kb": len(nearby),
        "gene_body_dnmt3l_dmr_count": len(gene_body),
        "promoter_2kb_dnmt3l_dmr_count": len(promoter),
        "nearby_dnmt3l_dmr_directions": ",".join(directions),
    }


def load_triple_hit_support(path: Path) -> dict[str, dict[str, object]]:
    support: dict[str, dict[str, object]] = {}
    for row in read_tsv(path):
        symbol = row.get("gene_symbol", "")
        if not symbol:
            continue
        rec = support.setdefault(
            symbol,
            {
                "triple_hit_dmr_count": 0,
                "triple_hit_gene_body_count": 0,
                "triple_hit_promoter_2kb_count": 0,
                "triple_hit_dmr_ids": [],
            },
        )
        rec["triple_hit_dmr_count"] += 1
        rec["triple_hit_gene_body_count"] += 1 if row.get("gene_body") == "yes" else 0
        rec["triple_hit_promoter_2kb_count"] += 1 if row.get("promoter_2kb") == "yes" else 0
        rec["triple_hit_dmr_ids"].append(row.get("dmr_id", ""))
    for rec in support.values():
        rec["triple_hit_dmr_ids"] = ",".join(sorted(set(rec["triple_hit_dmr_ids"])))
    return support


def transcriptomic_class(dnmt3l_row: dict[str, object], dnmt3a_row: dict[str, object] | None) -> str:
    dnmt3l_de = dnmt3l_row.get("gse268951_de_fdr05_fc2") == "yes"
    dnmt3a_de = bool(dnmt3a_row and dnmt3a_row.get("dnmt3a_de_fdr05_fc2") == "yes")
    if dnmt3l_de and dnmt3a_de:
        d1 = dnmt3l_row.get("gse268951_direction")
        d2 = dnmt3a_row.get("dnmt3a_direction") if dnmt3a_row else ""
        return "dnmt3l_de_and_dnmt3a_de_same_direction" if d1 == d2 else "dnmt3l_de_and_dnmt3a_de_opposite_direction"
    if dnmt3l_de and not dnmt3a_de:
        return "dnmt3l_de_not_dnmt3a_de"
    if not dnmt3l_de and dnmt3a_de:
        return "dnmt3a_de_not_dnmt3l_de"
    return "not_de_in_either_first_pass"


def methylation_independence_label(row: dict[str, object]) -> str:
    if row.get("transcriptomic_overlap_class") != "dnmt3l_de_not_dnmt3a_de":
        return "not_applicable"
    nearby = int(row.get("nearby_dnmt3l_dmr_count_100kb", 0) or 0)
    triple = int(row.get("triple_hit_dmr_count", 0) or 0)
    if triple > 0:
        return "dnmt3l_specific_expression_with_dnmt3a_related_methylation_context"
    if nearby > 0:
        return "dnmt3l_specific_expression_with_nearby_dnmt3l_methylation"
    return "dnmt3l_specific_expression_no_nearby_dnmt3l_dmr_methylation_independent_or_indirect_candidate"


def build_overlap_rows(
    gse268951_rows: list[dict[str, object]],
    dnmt3a_de: dict[str, dict[str, object]],
    dmrs_by_chrom: dict[str, list[dict[str, object]]],
    triple_support: dict[str, dict[str, object]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for expr in gse268951_rows:
        symbol = str(expr["gene_symbol"])
        d3a = dnmt3a_de.get(symbol, {})
        dmr = dmr_support_for_gene(expr, dmrs_by_chrom)
        triple = triple_support.get(
            symbol,
            {
                "triple_hit_dmr_count": 0,
                "triple_hit_gene_body_count": 0,
                "triple_hit_promoter_2kb_count": 0,
                "triple_hit_dmr_ids": "",
            },
        )
        row = {
            "gene_symbol": symbol,
            "gene_id": expr["gene_id"],
            "representative_transcript_id": expr.get("representative_transcript_id", ""),
            "transcript_count_for_gene": expr.get("transcript_count_for_gene", ""),
            "description": expr["description"],
            "gse268951_control_mean_count": expr.get("gse268951_control_mean_count", ""),
            "gse268951_dnmt3l_oe_mean_count": expr.get("gse268951_dnmt3l_oe_mean_count", ""),
            "gse268951_control_mean_fpkm": expr["gse268951_control_mean_fpkm"],
            "gse268951_dnmt3l_oe_mean_fpkm": expr["gse268951_dnmt3l_oe_mean_fpkm"],
            "gse268951_log2fc_oe_vs_control": expr["gse268951_log2fc_oe_vs_control"],
            "gse268951_effect_source": expr.get("gse268951_effect_source", ""),
            "gse268951_primary_pvalue": expr.get("gse268951_primary_pvalue", expr.get("gse268951_welch_pvalue", "")),
            "gse268951_primary_fdr": expr.get("gse268951_primary_fdr", expr.get("gse268951_bh_fdr", "")),
            "gse268951_welch_pvalue": expr["gse268951_welch_pvalue"],
            "gse268951_bh_fdr": expr["gse268951_bh_fdr"],
            "gse268951_fpkm_log2fc": expr.get("gse268951_fpkm_log2fc", ""),
            "gse268951_count_deseq2_log2fc": expr.get("gse268951_count_deseq2_log2fc", ""),
            "gse268951_count_deseq2_pvalue": expr.get("gse268951_count_deseq2_pvalue", ""),
            "gse268951_count_deseq2_fdr": expr.get("gse268951_count_deseq2_fdr", ""),
            "gse268951_count_edger_log2fc": expr.get("gse268951_count_edger_log2fc", ""),
            "gse268951_count_edger_fdr": expr.get("gse268951_count_edger_fdr", ""),
            "gse268951_direction": expr["gse268951_direction"],
            "gse268951_de_p05_fc2": expr["gse268951_de_p05_fc2"],
            "gse268951_de_fdr05_fc2": expr.get("gse268951_de_fdr05_fc2", ""),
            "dnmt3a_logfc_cko_vs_control": d3a.get("dnmt3a_logfc_cko_vs_control", ""),
            "dnmt3a_fdr": d3a.get("dnmt3a_fdr", ""),
            "dnmt3a_direction": d3a.get("dnmt3a_direction", "not_in_table"),
            "dnmt3a_de_fdr05_fc2": d3a.get("dnmt3a_de_fdr05_fc2", "no"),
            **dmr,
            **triple,
            "neurogenesis_related": expr["neurogenesis_related"],
            "neurogenesis_category_count": expr["neurogenesis_category_count"],
            "matched_neurogenesis_go_terms": expr["matched_neurogenesis_go_terms"],
        }
        row["transcriptomic_overlap_class"] = transcriptomic_class(expr, d3a)
        row["methylation_independence_screen"] = methylation_independence_label(row)
        rows.append(row)
    rows.sort(
        key=lambda r: (
            r["transcriptomic_overlap_class"] != "dnmt3l_de_not_dnmt3a_de",
            r["neurogenesis_related"] != "yes",
            -abs(safe_float(r["gse268951_log2fc_oe_vs_control"])),
        )
    )
    return rows


def final_candidate_neurogenesis_rows(
    summary_path: Path,
    by_symbol: dict[str, dict[str, object]],
    overlap_by_symbol: dict[str, dict[str, object]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for row in read_tsv(summary_path):
        gene = row["gene_symbol"]
        mouse = gene[0] + gene[1:].lower() if gene else ""
        expr = by_symbol.get(mouse, {})
        overlap = overlap_by_symbol.get(mouse, {})
        out = {
            "gene_symbol": gene,
            "mouse_symbol": mouse,
            "tier": row.get("tier", ""),
            "candidate_group": row.get("candidate_group", ""),
            "mechanism_class": row.get("mechanism_class", ""),
            "proposed_neurogenesis_role": neurogenesis_role(gene, expr, row),
            "neurogenesis_related": expr.get("neurogenesis_related", "no"),
            "core_neurogenesis": expr.get("core_neurogenesis", "no"),
            "brain_region_patterning": expr.get("brain_region_patterning", "no"),
            "neurite_synapse_migration": expr.get("neurite_synapse_migration", "no"),
            "gliogenesis_or_glial_context": expr.get("gliogenesis_or_glial_context", "no"),
            "neurogenesis_signaling": expr.get("neurogenesis_signaling", "no"),
            "neurogenesis_category_count": expr.get("neurogenesis_category_count", 0),
            "matched_neurogenesis_go_terms": expr.get("matched_neurogenesis_go_terms", ""),
            "gse268951_representative_transcript_id": expr.get("representative_transcript_id", ""),
            "gse268951_log2fc_oe_vs_control": expr.get("gse268951_log2fc_oe_vs_control", ""),
            "gse268951_primary_fdr": expr.get("gse268951_primary_fdr", ""),
            "gse268951_effect_source": expr.get("gse268951_effect_source", ""),
            "gse268951_de_p05_fc2": expr.get("gse268951_de_p05_fc2", ""),
            "gse268951_de_fdr05_fc2": expr.get("gse268951_de_fdr05_fc2", ""),
            "dnmt3a_logfc_cko_vs_control": overlap.get("dnmt3a_logfc_cko_vs_control", ""),
            "dnmt3a_de_fdr05_fc2": overlap.get("dnmt3a_de_fdr05_fc2", ""),
            "transcriptomic_overlap_class": overlap.get("transcriptomic_overlap_class", ""),
            "methylation_independence_screen": overlap.get("methylation_independence_screen", ""),
            "nearby_dnmt3l_dmr_count_100kb": overlap.get("nearby_dnmt3l_dmr_count_100kb", ""),
            "triple_hit_dmr_count": overlap.get("triple_hit_dmr_count", ""),
            "short_interpretation": row.get("short_interpretation", ""),
        }
        rows.append(out)
    return rows


def neurogenesis_role(gene: str, expr: dict[str, object], summary_row: dict[str, str]) -> str:
    if expr.get("core_neurogenesis") == "yes":
        return "direct_GO_neurogenesis_or_neuron_development"
    if expr.get("brain_region_patterning") == "yes":
        return "brain_region_patterning_or_cortical_development"
    if expr.get("neurite_synapse_migration") == "yes":
        return "neurite_synapse_migration_or_neuronal_maturation"
    if expr.get("neurogenesis_signaling") == "yes":
        return "neurogenesis_related_signaling"
    if summary_row.get("mechanism_class") == "polycomb_developmental_methylation_target":
        return "developmental_chromatin_candidate_without_direct_GO_neurogenesis_term"
    return "not_directly_neurogenesis_annotated_in_GSE268951_GO"


def make_final_candidate_heatmap(rows: list[dict[str, object]], figdir: Path) -> None:
    figdir.mkdir(parents=True, exist_ok=True)
    columns = [
        "core_neurogenesis",
        "brain_region_patterning",
        "neurite_synapse_migration",
        "gliogenesis_or_glial_context",
        "neurogenesis_signaling",
        "gse268951_de_fdr05_fc2",
        "dnmt3a_de_fdr05_fc2",
    ]
    col_labels = [
        "Core\nneurogenesis",
        "Brain-region\npatterning",
        "Neurite/synapse/\nmigration",
        "Glial\ncontext",
        "Neurogenesis\nsignaling",
        "Dnmt3l OE\nDE",
        "Dnmt3a cKO\nDE",
    ]
    matrix = np.array([[1 if str(row.get(col, "")) == "yes" else 0 for col in columns] for row in rows], dtype=float)
    labels = [str(row["gene_symbol"]) for row in rows]
    fig, ax = plt.subplots(figsize=(9.8, 8.0))
    im = ax.imshow(matrix, cmap="YlGnBu", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(np.arange(len(columns)))
    ax.set_xticklabels(col_labels, fontsize=8.2)
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels, fontsize=8.4)
    ax.tick_params(axis="x", top=True, bottom=False, labeltop=True, labelbottom=False, length=0)
    ax.tick_params(axis="y", length=0)
    for i, row in enumerate(rows):
        tier = str(row.get("tier", ""))
        color = "#D85C4A" if tier == "1_primary" else ("#E3A33B" if tier.startswith("2_") else ("#6E5A9B" if tier.startswith("3_") else "#4A7C59"))
        ax.scatter(-0.62, i, marker="s", s=24, color=color, clip_on=False)
        for j in range(len(columns)):
            if matrix[i, j] == 1:
                ax.text(j, i, "yes", ha="center", va="center", fontsize=7.2, color="white")
    ax.set_title("Neurogenesis annotations for final DNMT3L-linked candidates", fontsize=12.5, fontweight="bold", pad=54)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks(np.arange(-0.5, len(columns), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(labels), 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=1.1)
    ax.tick_params(which="minor", bottom=False, left=False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.025)
    cbar.set_ticks([0, 1])
    cbar.set_ticklabels(["no", "yes"])
    cbar.ax.tick_params(labelsize=8)
    fig.subplots_adjust(left=0.16, right=0.91, top=0.80, bottom=0.08)
    fig.text(
        0.16,
        0.035,
        "Row markers: red = primary; yellow = secondary; purple = developmental/Polycomb; green = expression/pathway methylation-unresolved.",
        ha="left",
        va="center",
        fontsize=8.1,
        color="#4B5563",
    )
    fig.savefig(figdir / "final_candidate_neurogenesis_annotation_heatmap.png", dpi=260, bbox_inches="tight")
    fig.savefig(figdir / "final_candidate_neurogenesis_annotation_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)


def make_overlap_barplot(rows: list[dict[str, object]], figdir: Path) -> None:
    figdir.mkdir(parents=True, exist_ok=True)
    classes = [
        "dnmt3l_de_not_dnmt3a_de",
        "dnmt3l_de_and_dnmt3a_de_same_direction",
        "dnmt3l_de_and_dnmt3a_de_opposite_direction",
        "dnmt3a_de_not_dnmt3l_de",
    ]
    labels = [
        "Dnmt3l DE\nnot Dnmt3a DE",
        "Both DE\nsame direction",
        "Both DE\nopposite",
        "Dnmt3a DE\nnot Dnmt3l DE",
    ]
    all_counts = [sum(1 for row in rows if row["transcriptomic_overlap_class"] == cls) for cls in classes]
    neuro_counts = [
        sum(1 for row in rows if row["transcriptomic_overlap_class"] == cls and row["neurogenesis_related"] == "yes")
        for cls in classes
    ]
    x = np.arange(len(classes))
    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    ax.bar(x - 0.18, all_counts, width=0.36, color="#9CA3AF", label="All genes")
    ax.bar(x + 0.18, neuro_counts, width=0.36, color="#2F6F8F", label="Neurogenesis-related")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8.4)
    ax.set_ylabel("Gene count")
    ax.set_title("Transcriptomic overlap between Dnmt3l OE and Dnmt3a cKO", fontsize=12.5, fontweight="bold")
    ax.legend(frameon=False, fontsize=8.5)
    ax.spines[["top", "right"]].set_visible(False)
    for i, (a, n) in enumerate(zip(all_counts, neuro_counts)):
        ax.text(i - 0.18, a + max(all_counts) * 0.015, str(a), ha="center", va="bottom", fontsize=8)
        ax.text(i + 0.18, n + max(all_counts) * 0.015, str(n), ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    fig.savefig(figdir / "dnmt3l_dnmt3a_transcriptomic_overlap_counts.png", dpi=260, bbox_inches="tight")
    fig.savefig(figdir / "dnmt3l_dnmt3a_transcriptomic_overlap_counts.pdf", bbox_inches="tight")
    plt.close(fig)


def write_note(
    path: Path,
    final_rows: list[dict[str, object]],
    overlap_rows: list[dict[str, object]],
    dnmt3l_sig_neuro_rows: list[dict[str, object]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    n_dnmt3l_de = sum(1 for row in overlap_rows if row["gse268951_de_fdr05_fc2"] == "yes")
    n_dnmt3a_de = sum(1 for row in overlap_rows if row["dnmt3a_de_fdr05_fc2"] == "yes")
    n_dnmt3l_only = sum(1 for row in overlap_rows if row["transcriptomic_overlap_class"] == "dnmt3l_de_not_dnmt3a_de")
    n_dnmt3l_only_neuro = sum(
        1
        for row in overlap_rows
        if row["transcriptomic_overlap_class"] == "dnmt3l_de_not_dnmt3a_de" and row["neurogenesis_related"] == "yes"
    )
    n_mi = sum(
        1
        for row in overlap_rows
        if row["methylation_independence_screen"]
        == "dnmt3l_specific_expression_no_nearby_dnmt3l_dmr_methylation_independent_or_indirect_candidate"
    )
    final_neuro = [row for row in final_rows if row["neurogenesis_related"] == "yes"]
    top_final = ", ".join(f"`{row['gene_symbol']}`" for row in final_neuro[:20])
    top_dnmt3l_only = [
        row
        for row in overlap_rows
        if row["transcriptomic_overlap_class"] == "dnmt3l_de_not_dnmt3a_de" and row["neurogenesis_related"] == "yes"
    ][:20]
    top_dnmt3l_only_text = ", ".join(f"`{row['gene_symbol']}`" for row in top_dnmt3l_only)
    lines = [
        "# Neurogenesis And Dnmt3a-Independence First Pass",
        "",
        "Date: 2026-06-09",
        "",
        "This analysis addresses two questions: whether final DNMT3L-linked candidates are functionally connected to neurogenesis, and whether Dnmt3l-responsive transcripts that are not altered in a Dnmt3a cKO reference can be identified for methylation-independent or DNMT3A-independent follow-up.",
        "",
        "## Data used",
        "",
        "- `GSE268951`: P28 medial frontal cortex Dnmt3l overexpression RNA-seq. Primary expression effects use a public GEO processed expected-count reanalysis from `GSE268951_all.fpkm_mRNA_anno.txt.gz` with DESeq2 and edgeR; FPKM is retained only as abundance context. The local `CON-VS-OE_DE_anno.xls` table is excluded from manuscript evidence because it is not a GEO public supplementary file.",
        "- `GSE141587`: P39 frontal cortex excitatory-neuron Dnmt3a cKO RNA-seq edgeR table.",
        "- `GSE271949`: Dnmt3l overexpression methylation DMRs, used to mark nearby DMR support.",
        "- Final 24-gene manuscript candidate table.",
        "",
        "## Main counts",
        "",
        f"- Public expected-count Dnmt3l-OE DE genes using FDR <= 0.05 and absolute logFC >= 1: {n_dnmt3l_de}.",
        f"- Dnmt3a-cKO DE genes using FDR <= 0.05 and absolute logFC >= 1: {n_dnmt3a_de}.",
        f"- Dnmt3l-OE DE but not Dnmt3a-cKO DE genes: {n_dnmt3l_only}.",
        f"- Dnmt3l-only DE genes with neurogenesis-related GO/context terms: {n_dnmt3l_only_neuro}.",
        f"- Dnmt3l-only DE genes with no nearby Dnmt3l DMR within the first-pass 100 kb screen: {n_mi}.",
        "",
        "## Final candidate neurogenesis interpretation",
        "",
        f"Final candidates with neurogenesis-related GO/context annotation in the GSE268951 annotation include: {top_final or 'none by this first-pass GO keyword screen'}.",
        "",
        "The most neurogenesis-relevant final-candidate mechanisms are not limited to direct `neurogenesis` GO terms. Patterning genes such as `OTX1`, `OTX2`, `LHX2`, `PRDM16`, and `NR2F1`, neuron identity genes such as `BCL11B`, and signaling genes such as `FZD7`, `FGFR2`, `PLPP3`, `RXRA`, `CRABP2`, and `CYP26A1` should be treated as neurodevelopmental/neurogenesis-context candidates even when the exact GO term is broader.",
        "",
        "## Dnmt3l-not-Dnmt3a interpretation",
        "",
        f"Top Dnmt3l-OE-only neurogenesis-related genes from the public expected-count expression screen include: {top_dnmt3l_only_text or 'none'}.",
        "",
        "Genes altered by Dnmt3l OE but not by Dnmt3a cKO are useful as DNMT3A-independent candidates. They should not automatically be called methylation-independent. A stronger in-silico methylation-independent/indirect class requires both Dnmt3l-specific expression and absence of nearby Dnmt3l methylation support in the available methylome screen. Even then, this remains a negative local-DMR screen, not proof that methylation is irrelevant.",
        "",
        "## Generated outputs",
        "",
        "- `results/neurogenesis_mechanism/final_candidate_neurogenesis_annotation.tsv`",
        "- `results/neurogenesis_mechanism/gse268951_rederived_gene_expression.tsv`",
        "- `results/neurogenesis_mechanism/dnmt3l_oe_significant_neurogenesis_genes.tsv`",
        "- `results/neurogenesis_mechanism/dnmt3l_vs_dnmt3a_transcriptomic_overlap.tsv`",
        "- `results/neurogenesis_mechanism/dnmt3l_not_dnmt3a_neurogenesis_candidates.tsv`",
        "- `figures/neurogenesis_mechanism/final_candidate_neurogenesis_annotation_heatmap.png`",
        "- `figures/neurogenesis_mechanism/dnmt3l_dnmt3a_transcriptomic_overlap_counts.png`",
        "",
        "## Caveats",
        "",
        "The Dnmt3l dataset is an overexpression model in medial frontal cortex, whereas the Dnmt3a reference is a conditional knockout in P39 frontal cortex excitatory pyramidal neurons. OE and KO are not reciprocal perturbations, and the cell populations are not identical. Therefore, `Dnmt3l altered but not Dnmt3a altered` should be interpreted as DNMT3A-independent candidate evidence, not definitive DNMT3A independence. Methylation-independent claims require additional negative evidence from matched methylation and expression in the same cell type.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    figdir = Path(args.figdir)
    if Path(args.gse268951_count_de).exists():
        gse268951_rows, by_symbol = load_gse268951_count_de(Path(args.gse268951_count_de))
    else:
        gse268951_rows, by_symbol = load_gse268951_gene_expression(Path(args.gse268951))
    dnmt3a_de = load_dnmt3a_de(Path(args.gse141587_de))
    dmrs_by_chrom = load_dmrs(Path(args.dnmt3l_dmrs))
    triple_support = load_triple_hit_support(Path(args.triple_links))
    overlap_rows = build_overlap_rows(gse268951_rows, dnmt3a_de, dmrs_by_chrom, triple_support)
    overlap_by_symbol = {str(row["gene_symbol"]): row for row in overlap_rows}
    final_rows = final_candidate_neurogenesis_rows(Path(args.summary), by_symbol, overlap_by_symbol)

    gse268951_fields = [
        "gene_id",
        "representative_transcript_id",
        "transcript_count_for_gene",
        "gene_symbol",
        "chrom",
        "start",
        "end",
        "strand",
        "description",
        "gse268951_control_mean_count",
        "gse268951_dnmt3l_oe_mean_count",
        "gse268951_control_mean_fpkm",
        "gse268951_dnmt3l_oe_mean_fpkm",
        "gse268951_log2fc_oe_vs_control",
        "gse268951_effect_source",
        "gse268951_primary_pvalue",
        "gse268951_primary_fdr",
        "gse268951_edgeR_logFC",
        "gse268951_edgeR_logCPM",
        "gse268951_edgeR_pvalue",
        "gse268951_edgeR_fdr",
        "gse268951_welch_pvalue",
        "gse268951_bh_fdr",
        "gse268951_fpkm_log2fc",
        "gse268951_direction",
        "gse268951_de_p05_fc2",
        "gse268951_de_fdr05_fc2",
        "neurogenesis_related",
        "neurogenesis_category_count",
        *NEUROGENESIS_KEYWORDS.keys(),
        "matched_neurogenesis_go_terms",
        "go_bp_terms",
    ]
    overlap_fields = [
        "gene_symbol",
        "gene_id",
        "representative_transcript_id",
        "transcript_count_for_gene",
        "description",
        "gse268951_control_mean_count",
        "gse268951_dnmt3l_oe_mean_count",
        "gse268951_control_mean_fpkm",
        "gse268951_dnmt3l_oe_mean_fpkm",
        "gse268951_log2fc_oe_vs_control",
        "gse268951_effect_source",
        "gse268951_primary_pvalue",
        "gse268951_primary_fdr",
        "gse268951_edgeR_logFC",
        "gse268951_edgeR_logCPM",
        "gse268951_edgeR_pvalue",
        "gse268951_edgeR_fdr",
        "gse268951_welch_pvalue",
        "gse268951_bh_fdr",
        "gse268951_fpkm_log2fc",
        "gse268951_direction",
        "gse268951_de_p05_fc2",
        "gse268951_de_fdr05_fc2",
        "dnmt3a_logfc_cko_vs_control",
        "dnmt3a_fdr",
        "dnmt3a_direction",
        "dnmt3a_de_fdr05_fc2",
        "transcriptomic_overlap_class",
        "methylation_independence_screen",
        "nearby_dnmt3l_dmr_count_100kb",
        "gene_body_dnmt3l_dmr_count",
        "promoter_2kb_dnmt3l_dmr_count",
        "nearby_dnmt3l_dmr_directions",
        "triple_hit_dmr_count",
        "triple_hit_gene_body_count",
        "triple_hit_promoter_2kb_count",
        "triple_hit_dmr_ids",
        "neurogenesis_related",
        "neurogenesis_category_count",
        "matched_neurogenesis_go_terms",
    ]
    final_fields = [
        "gene_symbol",
        "mouse_symbol",
        "tier",
        "candidate_group",
        "mechanism_class",
        "proposed_neurogenesis_role",
        "neurogenesis_related",
        "core_neurogenesis",
        "brain_region_patterning",
        "neurite_synapse_migration",
        "gliogenesis_or_glial_context",
        "neurogenesis_signaling",
        "neurogenesis_category_count",
        "matched_neurogenesis_go_terms",
        "gse268951_representative_transcript_id",
        "gse268951_log2fc_oe_vs_control",
        "gse268951_primary_fdr",
        "gse268951_effect_source",
        "gse268951_de_p05_fc2",
        "gse268951_de_fdr05_fc2",
        "dnmt3a_logfc_cko_vs_control",
        "dnmt3a_de_fdr05_fc2",
        "transcriptomic_overlap_class",
        "methylation_independence_screen",
        "nearby_dnmt3l_dmr_count_100kb",
        "triple_hit_dmr_count",
        "short_interpretation",
    ]

    sig_neuro = [
        row
        for row in gse268951_rows
        if row["gse268951_de_fdr05_fc2"] == "yes" and row["neurogenesis_related"] == "yes"
    ]
    sig_neuro.sort(key=lambda r: abs(safe_float(r["gse268951_log2fc_oe_vs_control"])), reverse=True)
    dnmt3l_only_neuro = [
        row
        for row in overlap_rows
        if row["transcriptomic_overlap_class"] == "dnmt3l_de_not_dnmt3a_de" and row["neurogenesis_related"] == "yes"
    ]

    write_tsv(outdir / "gse268951_rederived_gene_expression.tsv", gse268951_rows, gse268951_fields)
    write_tsv(outdir / "dnmt3l_oe_significant_neurogenesis_genes.tsv", sig_neuro, gse268951_fields)
    write_tsv(outdir / "dnmt3l_vs_dnmt3a_transcriptomic_overlap.tsv", overlap_rows, overlap_fields)
    write_tsv(outdir / "dnmt3l_not_dnmt3a_neurogenesis_candidates.tsv", dnmt3l_only_neuro, overlap_fields)
    write_tsv(outdir / "final_candidate_neurogenesis_annotation.tsv", final_rows, final_fields)
    make_final_candidate_heatmap(final_rows, figdir)
    make_overlap_barplot(overlap_rows, figdir)
    write_note(Path(args.note), final_rows, overlap_rows, sig_neuro)

    print(f"Wrote {outdir / 'gse268951_rederived_gene_expression.tsv'}")
    print(f"Wrote {outdir / 'dnmt3l_oe_significant_neurogenesis_genes.tsv'}")
    print(f"Wrote {outdir / 'dnmt3l_vs_dnmt3a_transcriptomic_overlap.tsv'}")
    print(f"Wrote {outdir / 'dnmt3l_not_dnmt3a_neurogenesis_candidates.tsv'}")
    print(f"Wrote {outdir / 'final_candidate_neurogenesis_annotation.tsv'}")
    print(f"Wrote {figdir / 'final_candidate_neurogenesis_annotation_heatmap.png'}")
    print(f"Wrote {figdir / 'dnmt3l_dnmt3a_transcriptomic_overlap_counts.png'}")
    print(f"Wrote {args.note}")


if __name__ == "__main__":
    main()
