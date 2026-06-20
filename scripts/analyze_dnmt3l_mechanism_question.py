#!/usr/bin/env python3
"""Build a focused evidence table for the DNMT3L mechanism question."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
from pathlib import Path

from scipy import stats


GSE157078_SAMPLES = [
    {"sample": "Z1", "condition": "control", "meth": "SAMPLE 1 Methylated signal", "unmeth": "SAMPLE 1 Unmethylated Signal", "pval": "SAMPLE 1 Detection Pval"},
    {"sample": "Z2", "condition": "control", "meth": "SAMPLE 2 Methylated signal", "unmeth": "SAMPLE 2 Unmethylated Signal", "pval": "SAMPLE 2 Detection Pval"},
    {"sample": "Z3", "condition": "control", "meth": "SAMPLE 3 Methylated signal", "unmeth": "SAMPLE 3 Unmethylated Signal", "pval": "SAMPLE 3 Detection Pval"},
    {"sample": "D1", "condition": "DNMT3L_overexpression", "meth": "SAMPLE 4 Methylated signal", "unmeth": "SAMPLE 4 Unmethylated Signal", "pval": "SAMPLE 4 Detection Pval"},
    {"sample": "D2", "condition": "DNMT3L_overexpression", "meth": "SAMPLE 5 Methylated signal", "unmeth": "SAMPLE 5 Unmethylated Signal", "pval": "SAMPLE 5 Detection Pval"},
    {"sample": "D3", "condition": "DNMT3L_overexpression", "meth": "SAMPLE 6 Methylated signal", "unmeth": "SAMPLE 6 Unmethylated Signal", "pval": "SAMPLE 6 Detection Pval"},
]
GSE157078_CONTROL = ["Z1", "Z2", "Z3"]
GSE157078_CASE = ["D1", "D2", "D3"]

GSE157080_SAMPLES = [
    {"sample": "Z1", "condition": "control", "column": "Z1_FPKM"},
    {"sample": "Z2", "condition": "control", "column": "Z2_FPKM"},
    {"sample": "Z3", "condition": "control", "column": "Z3_FPKM"},
    {"sample": "D1", "condition": "DNMT3L_overexpression", "column": "D1_FPKM"},
    {"sample": "D2", "condition": "DNMT3L_overexpression", "column": "D2_FPKM"},
    {"sample": "D3", "condition": "DNMT3L_overexpression", "column": "D3_FPKM"},
]

TARGET_AXIS_GENES = ["DNMT3L", "DNMT3A", "DNMT3B", "STAT1", "STAT3"]

GENE_ALIASES = {
    "PLPP3": ["PLPP3", "PPAP2B"],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gse157078-signals", default="data/geo/GSE157078/supplementary/GSE157078_signals.txt.gz")
    parser.add_argument("--epic-manifest", default="data/geo/GPL21145/GPL21145_MethylationEPIC_15073387_v-1-0.csv.gz")
    parser.add_argument("--gse157080-expression", default="data/geo/GSE157080/supplementary/GSE157080_all.fpkm_anno.txt.gz")
    parser.add_argument(
        "--gse157080-author-de",
        default="data/geo/GSE157080/supplementary/GSE157080_Table_S3_author_count_de.tsv",
    )
    parser.add_argument("--final-summary", default="results/manuscript_summary/final_candidate_manuscript_summary.tsv")
    parser.add_argument("--mechanism-table", default="results/mechanism_classification/methylation_expression_mechanism_table.tsv")
    parser.add_argument("--dnmt3a-overlap", default="results/neurogenesis_mechanism/dnmt3l_vs_dnmt3a_transcriptomic_overlap.tsv")
    parser.add_argument("--module-tests", default="results/neurogenesis_mechanism/developmental_module_rank_tests.tsv")
    parser.add_argument("--ds-target-summary", default="results/ds_expression_comparison/ds_target_gene_effect_summary.tsv")
    parser.add_argument("--outdir", default="results/mechanism_question")
    parser.add_argument("--note", default="notes/dnmt3l_mechanism_question_first_pass.md")
    parser.add_argument("--detection-p", type=float, default=0.01)
    parser.add_argument("--paper-delta", type=float, default=0.10)
    parser.add_argument("--paper-p", type=float, default=0.05)
    parser.add_argument("--beta-offset", type=float, default=100.0)
    parser.add_argument("--expression-pseudocount", type=float, default=0.1)
    parser.add_argument(
        "--exclude-snp-probes",
        action="store_true",
        help="Exclude all EPIC probes with any SNP annotation. Off by default because this is more stringent than the paper-level QC.",
    )
    return parser.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="", errors="replace")
    return path.open(newline="", errors="replace")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", errors="replace") as handle:
        return [{k: (v or "").strip() for k, v in row.items()} for row in csv.DictReader(handle, delimiter="\t")]


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


def mean(values: list[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    return sum(finite) / len(finite) if finite else math.nan


def fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}" if math.isfinite(value) else ""
    return str(value)


def direction(value: float, threshold: float = 0.1) -> str:
    if not math.isfinite(value) or abs(value) < threshold:
        return "flat"
    return "up" if value > 0 else "down"


def yes_no(value: bool) -> str:
    return "yes" if value else "no"


def chrom_name(value: str) -> str:
    value = (value or "").strip()
    if not value:
        return ""
    return value if value.startswith("chr") else f"chr{value}"


def read_epic_annotation(path: Path) -> dict[str, dict[str, str]]:
    annotations: dict[str, dict[str, str]] = {}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("IlmnID,"):
                header = next(csv.reader([line]))
                break
        else:
            raise RuntimeError(f"Could not find EPIC header in {path}")
        reader = csv.DictReader(handle, fieldnames=header)
        for row in reader:
            probe = row.get("Name") or row.get("IlmnID")
            chrom = chrom_name(row.get("CHR", ""))
            pos = (row.get("MAPINFO") or "").strip()
            if not probe or not chrom or not pos.isdigit():
                continue
            annotations[probe] = {
                "chrom": chrom,
                "position": pos,
                "gene_symbol": row.get("UCSC_RefGene_Name", ""),
                "gene_group": row.get("UCSC_RefGene_Group", ""),
                "cpg_relation": row.get("Relation_to_UCSC_CpG_Island", ""),
                "snp_id": row.get("SNP_ID", ""),
                "snp_maf": row.get("SNP_MinorAlleleFrequency", ""),
            }
    return annotations


def beta_values(row: dict[str, str], detection_p: float, beta_offset: float) -> dict[str, float] | None:
    values: dict[str, float] = {}
    for sample in GSE157078_SAMPLES:
        pvalue = safe_float(row.get(sample["pval"], ""))
        if not math.isfinite(pvalue) or pvalue > detection_p:
            return None
        methylated = safe_float(row.get(sample["meth"], ""))
        unmethylated = safe_float(row.get(sample["unmeth"], ""))
        denominator = methylated + unmethylated + beta_offset
        if not math.isfinite(denominator) or denominator <= 0:
            return None
        values[sample["sample"]] = methylated / denominator
    return values


def call_gse157078_paper_like(args: argparse.Namespace) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    annotations = read_epic_annotation(Path(args.epic_manifest))
    counts = {
        "signals_seen": 0,
        "passed_detection": 0,
        "annotated": 0,
        "autosomal": 0,
        "snp_annotated": 0,
        "analysis_background": 0,
        "paper_like_candidates": 0,
        "hyper": 0,
        "hypo": 0,
    }
    candidates: list[dict[str, object]] = []
    with open_text(Path(args.gse157078_signals)) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            counts["signals_seen"] += 1
            probe = row.get("ID_REF", "")
            values = beta_values(row, args.detection_p, args.beta_offset)
            if values is None:
                continue
            counts["passed_detection"] += 1
            anno = annotations.get(probe)
            if not anno:
                continue
            counts["annotated"] += 1
            if anno["chrom"] in {"chrX", "chrY"}:
                continue
            counts["autosomal"] += 1
            has_snp = bool(anno.get("snp_id") or anno.get("snp_maf"))
            if has_snp:
                counts["snp_annotated"] += 1
            if args.exclude_snp_probes and has_snp:
                continue
            counts["analysis_background"] += 1

            control = [values[sample] for sample in GSE157078_CONTROL]
            case = [values[sample] for sample in GSE157078_CASE]
            delta = mean(case) - mean(control)
            _stat, pvalue = stats.ttest_ind(case, control, equal_var=True)
            if not math.isfinite(pvalue):
                continue
            if pvalue <= args.paper_p and abs(delta) >= args.paper_delta:
                call = "hyper" if delta > 0 else "hypo"
                counts["paper_like_candidates"] += 1
                counts[call] += 1
                candidates.append(
                    {
                        "probe_id": probe,
                        "chrom": anno["chrom"],
                        "start": int(anno["position"]) - 1,
                        "end": int(anno["position"]),
                        "direction": call,
                        "delta_beta": delta,
                        "pvalue_student_t": pvalue,
                        "control_mean_beta": mean(control),
                        "dnmt3l_oe_mean_beta": mean(case),
                        "gene_symbol": anno.get("gene_symbol", ""),
                        "gene_group": anno.get("gene_group", ""),
                        "cpg_relation": anno.get("cpg_relation", ""),
                    }
                )

    summary = [{"metric": key, "value": value} for key, value in counts.items()]
    return candidates, summary


def expression_values(row: dict[str, str], pseudocount: float) -> tuple[list[float], list[float], dict[str, float]]:
    sample_values = {sample["sample"]: safe_float(row.get(sample["column"], "")) for sample in GSE157080_SAMPLES}
    control = [math.log2(max(sample_values[sample], 0.0) + pseudocount) for sample in GSE157078_CONTROL]
    case = [math.log2(max(sample_values[sample], 0.0) + pseudocount) for sample in GSE157078_CASE]
    return control, case, sample_values


def read_gse157080_author_de(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    rows = read_tsv(path)
    out: dict[str, dict[str, str]] = {}
    for row in rows:
        symbol = (row.get("GeneSymbol") or "").strip().upper()
        gene_id = (row.get("gene_id") or "").strip().upper()
        if symbol and symbol not in out:
            out[symbol] = row
        if gene_id and gene_id not in out:
            out[gene_id] = row
    return out


def summarize_gse157080_expression(args: argparse.Namespace, genes: list[str]) -> list[dict[str, object]]:
    wanted = {gene.upper() for gene in genes}
    alias_to_gene = {gene: gene for gene in wanted}
    for gene, aliases in GENE_ALIASES.items():
        if gene in wanted:
            for alias in aliases:
                alias_to_gene[alias.upper()] = gene
    selected: dict[str, dict[str, str]] = {}
    sample_cols = [sample["column"] for sample in GSE157080_SAMPLES]
    with open_text(Path(args.gse157080_expression)) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            symbol = (row.get("GeneSymbol") or "").strip().upper()
            gene = alias_to_gene.get(symbol)
            if gene is None:
                continue
            current = selected.get(gene)
            current_mean = mean([safe_float(current.get(col, "")) for col in sample_cols]) if current else -1
            row_mean = mean([safe_float(row.get(col, "")) for col in sample_cols])
            if current is None or row_mean > current_mean:
                selected[gene] = row

    author_de = read_gse157080_author_de(Path(args.gse157080_author_de))
    rows: list[dict[str, object]] = []
    for gene in sorted(wanted):
        row = selected.get(gene)
        if row is None:
            rows.append({"gene_symbol": gene, "found": "no"})
            continue
        control, case, sample_values = expression_values(row, args.expression_pseudocount)
        fpkm_delta = mean(case) - mean(control)
        _stat, pvalue = stats.ttest_ind(case, control, equal_var=False)
        author = author_de.get(str(row.get("gene_id", "")).upper()) or author_de.get(gene) or {}
        author_delta = safe_float(author.get("log2FoldChange"))
        author_pvalue = safe_float(author.get("pvalue"))
        author_qvalue = safe_float(author.get("qvalue"))
        if math.isfinite(author_delta):
            primary_delta = author_delta
            primary_pvalue = author_pvalue
            primary_fdr = author_qvalue
            source = "author_Table_S3_count_based_DE"
        else:
            primary_delta = fpkm_delta
            primary_pvalue = pvalue
            primary_fdr = math.nan
            source = "processed_FPKM_Welch_fallback"
        rows.append(
            {
                "gene_symbol": gene,
                "found": "yes",
                "gene_id": row.get("gene_id", ""),
                "control_mean_fpkm": mean([sample_values[sample] for sample in GSE157078_CONTROL]),
                "dnmt3l_oe_mean_fpkm": mean([sample_values[sample] for sample in GSE157078_CASE]),
                "control_mean_log2_fpkm": mean(control),
                "dnmt3l_oe_mean_log2_fpkm": mean(case),
                "log2_delta_oe_vs_control": primary_delta,
                "expression_pvalue": primary_pvalue,
                "expression_fdr": primary_fdr,
                "expression_effect_source": source,
                "author_log2FoldChange": author_delta,
                "author_pvalue": author_pvalue,
                "author_qvalue": author_qvalue,
                "fpkm_log2_delta_oe_vs_control": fpkm_delta,
                "fpkm_welch_pvalue_log2": pvalue,
                "direction": direction(primary_delta),
                "description": row.get("Description", ""),
            }
        )
    return rows


def first_number_after(pattern: str, value: str) -> float:
    match = re.search(pattern, value or "")
    return safe_float(match.group(1)) if match else math.nan


def parse_summary_value(rows: list[dict[str, str]], metric: str) -> str:
    for row in rows:
        if row.get("metric") == metric:
            return row.get("value", "")
    return ""


def load_by_gene(path: Path, field: str = "gene_symbol") -> dict[str, dict[str, str]]:
    rows = read_tsv(path)
    return {(row.get(field) or "").upper(): row for row in rows if row.get(field)}


def relation_label(methyl_direction: str, expr_direction: str) -> str:
    if methyl_direction == "hyper" and expr_direction == "down":
        return "repression_like"
    if methyl_direction == "hypo" and expr_direction == "up":
        return "derepression_like"
    if methyl_direction in {"hyper", "hypo"} and expr_direction in {"up", "down"}:
        return "directionally_complex"
    return "unresolved"


def build_mechanism_table(
    args: argparse.Namespace,
    expression_rows: list[dict[str, object]],
    paper_summary: list[dict[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    final_rows = read_tsv(Path(args.final_summary))
    mechanism_by_gene = load_by_gene(Path(args.mechanism_table))
    dnmt3a_by_mouse_gene = load_by_gene(Path(args.dnmt3a_overlap))
    expr_by_gene = {str(row.get("gene_symbol", "")).upper(): row for row in expression_rows}
    module_rows = read_tsv(Path(args.module_tests))
    ds_target_rows = read_tsv(Path(args.ds_target_summary))

    axis_expr = {gene: expr_by_gene.get(gene, {}) for gene in TARGET_AXIS_GENES}
    stat3_ds_contexts = [row for row in ds_target_rows if row.get("gene_symbol") == "STAT3"]
    stat3_ds_positive = sum(1 for row in stat3_ds_contexts if safe_float(row.get("delta", row.get("effect", ""))) > 0)
    if stat3_ds_positive == 0:
        stat3_ds_positive = sum(1 for row in stat3_ds_contexts if (row.get("direction") or "").lower() == "up")

    table: list[dict[str, object]] = []
    for row in final_rows:
        gene = row.get("gene_symbol", "")
        gene_key = gene.upper()
        mech = mechanism_by_gene.get(gene_key, {})
        expr = expr_by_gene.get(gene_key, {})
        mouse_symbol = mech.get("mouse_symbol") or row.get("gene_symbol", "")
        dnmt3a = dnmt3a_by_mouse_gene.get(mouse_symbol.upper(), {})

        human_strict = safe_float(row.get("human_strict_regulatory_probe_count") or mech.get("human_strict_regulatory_probe_count"))
        human_probes = safe_float(row.get("human_strict_fetal_positive_probe_count") or row.get("human_strict_fetal_positive_probes"))
        if not math.isfinite(human_probes):
            human_probes = first_number_after(r"strict probes\s+([0-9.]+)", row.get("human_methylation", ""))
        mouse_triple = safe_float(row.get("mouse_methylation", "").split("DMRs ")[-1].split(";")[0]) if "DMRs " in row.get("mouse_methylation", "") else safe_float(row.get("mouse_linked_triple_hit_dmr_count"))
        mouse_reg = safe_float(row.get("mouse regulatory DMRs", "") or mech.get("mouse_any_regulatory_dmr_count"))
        if not math.isfinite(mouse_reg):
            mouse_reg = safe_float(row.get("regulatory_context", "").split("mouse regulatory DMRs ")[-1]) if "mouse regulatory DMRs " in row.get("regulatory_context", "") else math.nan

        methyl_score = 0
        methyl_score += 1 if human_probes > 0 else 0
        methyl_score += 1 if human_strict > 0 else 0
        methyl_score += 1 if mouse_triple > 0 else 0
        methyl_score += 1 if mouse_reg > 0 else 0
        methyl_score += 1 if "repression" in (mech.get("human_methyl_expr_relation", "") + mech.get("mouse_oe_methyl_expr_relation", "")) else 0

        expr_direction = str(expr.get("direction", ""))
        direct_delta = safe_float(expr.get("log2_delta_oe_vs_control"))
        direct_p = safe_float(expr.get("expression_pvalue"))
        if not math.isfinite(direct_p):
            direct_p = safe_float(expr.get("fpkm_welch_pvalue_log2"))
        if not math.isfinite(direct_p):
            direct_p = safe_float(expr.get("welch_pvalue_log2"))
        direct_fdr = safe_float(expr.get("expression_fdr"))
        direct_source = str(expr.get("expression_effect_source", ""))
        direct_expr_signal = expr_direction in {"up", "down"} and (
            (math.isfinite(direct_fdr) and direct_fdr <= 0.05)
            or (not math.isfinite(direct_fdr) and math.isfinite(direct_p) and direct_p <= 0.1)
        )
        dnmt3l_not_dnmt3a = dnmt3a.get("transcriptomic_overlap_class") == "dnmt3l_de_not_dnmt3a_de"
        no_nearby_dmr = "no_nearby" in dnmt3a.get("methylation_independence_screen", "")
        complex_direction = "directionally_complex" in row.get("mechanism_class", "") or "directionally_complex" in mech.get("mechanism_class", "")

        tf_indirect_score = 0
        tf_indirect_score += 1 if direct_expr_signal else 0
        tf_indirect_score += 1 if dnmt3l_not_dnmt3a else 0
        tf_indirect_score += 1 if no_nearby_dmr else 0
        tf_indirect_score += 1 if complex_direction else 0
        tf_indirect_score += 1 if "TF" in row.get("human_expression", "") or "transcription" in row.get("short_interpretation", "").lower() else 0

        if methyl_score >= 4 and tf_indirect_score >= 2:
            mechanism_call = "both_methylation_and_tf_indirect_plausible"
        elif methyl_score >= 4:
            mechanism_call = "methylation_linked_priority"
        elif tf_indirect_score >= 3:
            mechanism_call = "tf_indirect_priority_methylation_unresolved"
        else:
            mechanism_call = "mixed_or_unresolved"

        methyl_direction = mech.get("human_methylation_direction", "")
        direct_relation = relation_label(methyl_direction, expr_direction)
        table.append(
            {
                "gene_symbol": gene,
                "tier": row.get("tier", ""),
                "candidate_group": row.get("candidate_group", ""),
                "previous_mechanism_class": row.get("mechanism_class", ""),
                "focused_mechanism_call": mechanism_call,
                "methylation_axis_score": methyl_score,
                "tf_indirect_axis_score": tf_indirect_score,
                "human_strict_fetal_positive_probes": row.get("human_methylation", ""),
                "human_regulatory_probe_count": human_strict,
                "mouse_triple_hit_dmr_count": mouse_triple,
                "mouse_regulatory_dmr_count": mouse_reg,
                "human_direct_dnmt3l_oe_log2fc": direct_delta,
                "human_direct_dnmt3l_oe_pvalue": direct_p,
                "human_direct_dnmt3l_oe_fdr": direct_fdr,
                "human_direct_expression_source": direct_source,
                "human_direct_dnmt3l_oe_direction": expr_direction,
                "human_direct_methyl_expr_relation": direct_relation,
                "mouse_gse268951_direction": mech.get("mouse_gse268951_direction", ""),
                "mouse_gse268951_log2fc": mech.get("mouse_gse268951_log2fc", ""),
                "dnmt3a_transcriptomic_overlap_class": dnmt3a.get("transcriptomic_overlap_class", ""),
                "methylation_independence_screen": dnmt3a.get("methylation_independence_screen", ""),
                "nearby_dnmt3l_dmr_count_100kb": dnmt3a.get("nearby_dnmt3l_dmr_count_100kb", ""),
                "stage_celltype_interpretation": stage_interpretation(row, mechanism_call),
            }
        )

    registry = [
        {
            "test": "GSE157078 paper-like DMC QC",
            "status": "completed_first_pass",
            "main_output": "results/mechanism_question/gse157078_paper_like_dmcs.tsv",
            "readout": f"{parse_summary_value(paper_summary, 'paper_like_candidates')} paper-like CpGs; {parse_summary_value(paper_summary, 'hyper')} hyper; {parse_summary_value(paper_summary, 'hypo')} hypo",
            "answers_question_by": "validates the direct human DNMT3L-OE methylation foreground before mechanism claims",
        },
        {
            "test": "GSE157080 direct expression of candidates and STAT axis",
            "status": "completed_first_pass",
            "main_output": "results/mechanism_question/gse157080_direct_expression_final_candidates.tsv",
            "readout": f"STAT1 log2FC {fmt(axis_expr.get('STAT1', {}).get('log2_delta_oe_vs_control', ''))}; STAT3 log2FC {fmt(axis_expr.get('STAT3', {}).get('log2_delta_oe_vs_control', ''))}",
            "answers_question_by": "tests direct transcriptional response in the same human DNMT3L-OE system",
        },
        {
            "test": "Mouse Dnmt3l DMR overlap with neuronal DNMT3A-dependent methylation",
            "status": "completed_existing",
            "main_output": "results/mouse_dnmt3l_dnmt3a_neuron_overlap/gse271949_dmr_priority_triple_hits.tsv",
            "readout": "used as mouse triple-hit methylation evidence in the focused table",
            "answers_question_by": "supports or weakens DNMT3A-mediated methylation interpretation",
        },
        {
            "test": "Dnmt3l-OE transcriptome versus Dnmt3a-cKO transcriptome",
            "status": "completed_existing",
            "main_output": "results/neurogenesis_mechanism/dnmt3l_vs_dnmt3a_transcriptomic_overlap.tsv",
            "readout": "flags Dnmt3l-responsive genes that do not pass Dnmt3a-cKO transcriptomic threshold",
            "answers_question_by": "prioritizes DNMT3A-independent or indirect transcriptional candidates",
        },
        {
            "test": "SH-SY5Y phase-specific DNMT3L WGBS validation",
            "status": "completed_existing",
            "main_output": "results/shsy5y_gse168276/gse168276_gse157078_target_concordance_summary.tsv",
            "readout": "Phase 3 shows strongest positive targeted concordance but modest per-probe agreement",
            "answers_question_by": "tests whether DNMT3L methylation effects strengthen with neuronal differentiation",
        },
        {
            "test": "DS brain methylation concordance",
            "status": "completed_existing",
            "main_output": "results/ds_methylation_dnmt3l_concordance/ds_dnmt3l_candidate_probe_concordance_summary.tsv",
            "readout": "fetal and adult DS methylation concordance is context-dependent",
            "answers_question_by": "tests whether DNMT3L-sensitive loci resemble DS brain methylation",
        },
        {
            "test": "DS cell-type expression and candidate pseudobulk",
            "status": "completed_existing",
            "main_output": "results/cross_species_candidate_expression/gse305153_candidate_celltype_effects.tsv",
            "readout": "cell-type context is available for final candidates and target genes",
            "answers_question_by": "places DNMT3L/STAT/candidate effects into fetal cortical cell contexts",
        },
        {
            "test": "GSE74486 OXBS 5mC versus 5hmC separation",
            "status": "not_started_high_value",
            "main_output": "planned: results/mechanism_question/gse74486_oxbs_5mc_5hmc_candidate_summary.tsv",
            "readout": "not yet run",
            "answers_question_by": "distinguishes true methylation from hydroxymethylation in DS adult brain",
        },
        {
            "test": "scATAC/chromVAR TF activity in DS fetal cortex",
            "status": "not_started_high_value",
            "main_output": "planned: results/mechanism_question/gse305146_tf_activity_target_motifs.tsv",
            "readout": "not yet run",
            "answers_question_by": "directly tests STAT/TF pathway activity independent of DNA methylation",
        },
    ]

    return table, registry


def stage_interpretation(row: dict[str, str], mechanism_call: str) -> str:
    group = row.get("candidate_group", "")
    gene = row.get("gene_symbol", "")
    if "methylation_linked" in mechanism_call and gene == "FOXP1":
        return "human fetal/DS and mouse postnatal evidence support a methylation-repression candidate, with cortical excitatory-neuron relevance"
    if "methylation_linked" in mechanism_call:
        return "methylation-linked evidence is substantial, but paired methylation-expression mediation is still needed"
    if "both" in mechanism_call:
        return "methylation-linked regulatory evidence is present, but expression direction or DNMT3A-independence suggests additional TF, enhancer, timing, or cell-state mechanisms"
    if "TF" in mechanism_call or "tf_indirect" in mechanism_call:
        return "best treated as a transcriptional/pathway candidate until stronger local methylation mediation is shown"
    if "Methylation-heavy" in group:
        return "stronger as a developmental methylation target than as an expression biomarker in current windows"
    return "mixed evidence; needs paired methylation-expression or perturbation validation"


def write_note(
    path: Path,
    paper_summary: list[dict[str, object]],
    table: list[dict[str, object]],
    registry: list[dict[str, object]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    calls: dict[str, int] = {}
    for row in table:
        calls[str(row["focused_mechanism_call"])] = calls.get(str(row["focused_mechanism_call"]), 0) + 1
    primary = [row for row in table if str(row.get("tier", "")).startswith("1_")]
    with path.open("w") as handle:
        handle.write("# DNMT3L Mechanism Question First Pass\n\n")
        handle.write("Question: does DNMT3L affect neurodevelopment through DNMT3A/DNMT3B-mediated methylation, through STAT/TF pathways independent of methylation, or both depending on stage/cell type?\n\n")
        handle.write("## What Was Added\n\n")
        handle.write("1. A paper-like `GSE157078` DMC QC call using detection p <= 0.01, autosomal EPIC probes, Student t-test p <= 0.05, and absolute delta beta >= 0.10. SNP annotations are counted, but not excluded by default because excluding every SNP-annotated probe is more stringent than the paper-level QC.\n")
        handle.write("2. A direct `GSE157080` expression summary for final candidates plus DNMT3L/DNMT3A/DNMT3B/STAT1/STAT3, using author supplementary Table S3 count-based DE statistics for primary log2FC/p/q values and GEO FPKM only as abundance context.\n")
        handle.write("3. A focused mechanism table combining existing methylation, DNMT3A-reference, expression, DS, and regulatory evidence.\n\n")
        handle.write("## Paper-Like GSE157078 QC\n\n")
        for row in paper_summary:
            handle.write(f"- {row['metric']}: {row['value']}\n")
        handle.write("\nThis is a QC foreground, not a claim that we exactly reproduced every original preprocessing decision.\n\n")
        handle.write("## Focused Mechanism Calls\n\n")
        for call, count in sorted(calls.items()):
            handle.write(f"- {call}: {count}\n")
        handle.write("\nPrimary candidates:\n\n")
        for row in primary:
            handle.write(f"- {row['gene_symbol']}: {row['focused_mechanism_call']}; methylation score {row['methylation_axis_score']}; TF/indirect score {row['tf_indirect_axis_score']}.\n")
        handle.write("\n## Test Registry\n\n")
        for row in registry:
            handle.write(f"- {row['test']}: {row['status']}. {row['readout']}\n")
        handle.write("\n## Working Interpretation\n\n")
        handle.write("The strongest current answer is a dual-axis model. Several candidates carry substantial methylation-linked evidence, especially human fetal regulatory methylation and mouse DNMT3A-linked DMR support. However, many expression directions are not compatible with simple promoter silencing, and the direct STAT/TF pathway remains an important parallel mechanism, especially for early human neural differentiation. The next highest-value missing tests are OXBS 5mC/5hmC separation in `GSE74486` and TF activity/accessibility analysis in `GSE305146`.\n")


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)

    paper_dmcs, paper_summary = call_gse157078_paper_like(args)
    write_tsv(
        outdir / "gse157078_paper_like_dmcs.tsv",
        paper_dmcs,
        [
            "probe_id",
            "chrom",
            "start",
            "end",
            "direction",
            "delta_beta",
            "pvalue_student_t",
            "control_mean_beta",
            "dnmt3l_oe_mean_beta",
            "gene_symbol",
            "gene_group",
            "cpg_relation",
        ],
    )
    write_tsv(outdir / "gse157078_paper_like_dmc_summary.tsv", paper_summary, ["metric", "value"])

    final_rows = read_tsv(Path(args.final_summary))
    final_genes = [row["gene_symbol"] for row in final_rows if row.get("gene_symbol")]
    expression_rows = summarize_gse157080_expression(args, sorted(set(final_genes + TARGET_AXIS_GENES)))
    write_tsv(
        outdir / "gse157080_direct_expression_final_candidates.tsv",
        expression_rows,
        [
            "gene_symbol",
            "found",
            "gene_id",
            "control_mean_fpkm",
            "dnmt3l_oe_mean_fpkm",
            "control_mean_log2_fpkm",
            "dnmt3l_oe_mean_log2_fpkm",
            "log2_delta_oe_vs_control",
            "expression_pvalue",
            "expression_fdr",
            "expression_effect_source",
            "author_log2FoldChange",
            "author_pvalue",
            "author_qvalue",
            "fpkm_log2_delta_oe_vs_control",
            "fpkm_welch_pvalue_log2",
            "direction",
            "description",
        ],
    )

    mechanism_table, registry = build_mechanism_table(args, expression_rows, paper_summary)
    write_tsv(
        outdir / "dnmt3l_mechanism_evidence_table.tsv",
        mechanism_table,
        [
            "gene_symbol",
            "tier",
            "candidate_group",
            "previous_mechanism_class",
            "focused_mechanism_call",
            "methylation_axis_score",
            "tf_indirect_axis_score",
            "human_strict_fetal_positive_probes",
            "human_regulatory_probe_count",
            "mouse_triple_hit_dmr_count",
            "mouse_regulatory_dmr_count",
            "human_direct_dnmt3l_oe_log2fc",
            "human_direct_dnmt3l_oe_pvalue",
            "human_direct_dnmt3l_oe_fdr",
            "human_direct_expression_source",
            "human_direct_dnmt3l_oe_direction",
            "human_direct_methyl_expr_relation",
            "mouse_gse268951_direction",
            "mouse_gse268951_log2fc",
            "dnmt3a_transcriptomic_overlap_class",
            "methylation_independence_screen",
            "nearby_dnmt3l_dmr_count_100kb",
            "stage_celltype_interpretation",
        ],
    )
    write_tsv(
        outdir / "dnmt3l_mechanism_test_registry.tsv",
        registry,
        ["test", "status", "main_output", "readout", "answers_question_by"],
    )
    write_note(Path(args.note), paper_summary, mechanism_table, registry)


if __name__ == "__main__":
    main()
