#!/usr/bin/env python3
"""Build manuscript-ready candidate summary tables and draft results text."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


DEFAULT_MATRIX = "results/final_candidate_evidence/final_candidate_evidence_matrix_with_regulatory_context.tsv"
DEFAULT_MECHANISM = "results/mechanism_classification/methylation_expression_mechanism_table.tsv"
DEFAULT_PERMUTATION = "results/mechanism_classification/matched_regulatory_permutation_tests.tsv"
DEFAULT_OUTDIR = "results/manuscript_summary"
DEFAULT_NOTE = "notes/manuscript_results_candidate_mechanism.md"


SUMMARY_FIELDS = [
    "gene_symbol",
    "tier",
    "candidate_group",
    "mechanism_class",
    "human_methylation",
    "human_expression",
    "mouse_methylation",
    "mouse_expression",
    "regulatory_context",
    "directionality",
    "paper_role",
    "short_interpretation",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", default=DEFAULT_MATRIX)
    parser.add_argument("--mechanism", default=DEFAULT_MECHANISM)
    parser.add_argument("--permutation", default=DEFAULT_PERMUTATION)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--note", default=DEFAULT_NOTE)
    return parser.parse_args()


def clean_key(value: str | None) -> str:
    return (value or "").strip().replace("\r", "")


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
        writer.writerows(rows)


def safe_float(value: object) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return math.nan
    return number if math.isfinite(number) else math.nan


def fmt_num(value: object, digits: int = 2) -> str:
    number = safe_float(value)
    if not math.isfinite(number):
        return ""
    return f"{number:.{digits}f}"


def tier_label(tier: str) -> str:
    labels = {
        "1_primary": "Primary methylation-convergent",
        "2_secondary_balanced": "Secondary methylation + expression",
        "3_methylation_heavy": "Methylation-heavy developmental",
        "4_expression_pathway": "Expression/pathway methylation-unresolved",
    }
    return labels.get(tier, tier)


def paper_role(tier: str, mechanism_class: str) -> str:
    if tier == "1_primary":
        return "Main candidate panel"
    if mechanism_class == "polycomb_developmental_methylation_target":
        return "Developmental/Polycomb mechanism panel"
    if tier == "2_secondary_balanced":
        return "Secondary validation table"
    return "Pathway/context table"


def short_interpretation(row: dict[str, str], mechanism: dict[str, str]) -> str:
    gene = row["gene_symbol"]
    cls = mechanism.get("mechanism_class", "")
    if gene == "FOXP1":
        return "Best primary candidate for repression-like methylation-expression direction."
    if cls == "directionally_complex_regulatory_methylation":
        return "Strong regulatory-convergent candidate; direction suggests enhancer, gene-body, indirect, timing-specific, or cell-state mechanism."
    if cls == "polycomb_developmental_methylation_target":
        return "Developmental methylation target with Polycomb/H3K27me3 context and weak or timing-restricted expression."
    if cls == "secondary_regulatory_methylation_candidate":
        return "Secondary methylation-expression candidate for validation."
    if cls == "expression_pathway_methylation_unresolved":
        return "Expression/pathway candidate; methylation independence is not proven."
    return mechanism.get("interpretation", "")


def build_summary_rows(matrix_rows: list[dict[str, str]], mechanism_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    mechanism_by_gene = {row["gene_symbol"]: row for row in mechanism_rows}
    rows: list[dict[str, str]] = []
    for row in matrix_rows:
        gene = row["gene_symbol"]
        mech = mechanism_by_gene.get(gene, {})
        human_methyl = (
            f"{mech.get('human_methylation_direction', '')}; "
            f"strict probes {row.get('human_strict_fetal_positive_probe_count', row.get('human_strict_fetal_positive_probes', ''))}; "
            f"DS signed delta {fmt_num(mech.get('human_ds_signed_delta'), 3)}"
        )
        human_expression = (
            f"{mech.get('human_expression_direction', '')} "
            f"{fmt_num(mech.get('human_expression_delta'), 2)} "
            f"({mech.get('human_expression_context', '')})"
        )
        mouse_methyl = (
            f"{mech.get('mouse_methylation_direction', '')}; "
            f"DMRs {row.get('mouse_linked_triple_hit_dmr_count', '')}; "
            f"mean delta {fmt_num(mech.get('mouse_mean_dmr_delta_beta'), 3)}"
        )
        mouse_expression = (
            f"OE {mech.get('mouse_gse268951_direction', '')} {fmt_num(mech.get('mouse_gse268951_log2fc'), 2)}; "
            f"Dp10 {mech.get('mouse_dp10_direction', '')} {fmt_num(mech.get('mouse_dp10_delta'), 2)}"
        )
        regulatory = (
            f"{row.get('regulatory_context_class', '')}; "
            f"human strict-reg probes {row.get('human_strict_regulatory_probe_count', '')}; "
            f"mouse regulatory DMRs {row.get('mouse_any_regulatory_dmr_count', '')}"
        )
        directionality = (
            f"human {mech.get('human_methyl_expr_relation', '')}; "
            f"mouse OE {mech.get('mouse_oe_methyl_expr_relation', '')}; "
            f"mouse Dp10 {mech.get('mouse_dp10_methyl_expr_relation', '')}"
        )
        rows.append(
            {
                "gene_symbol": gene,
                "tier": row["tier"],
                "candidate_group": tier_label(row["tier"]),
                "mechanism_class": mech.get("mechanism_class", ""),
                "human_methylation": human_methyl,
                "human_expression": human_expression,
                "mouse_methylation": mouse_methyl,
                "mouse_expression": mouse_expression,
                "regulatory_context": regulatory,
                "directionality": directionality,
                "paper_role": paper_role(row["tier"], mech.get("mechanism_class", "")),
                "short_interpretation": short_interpretation(row, mech),
            }
        )
    return rows


def select_permutation_rows(rows: list[dict[str, str]]) -> dict[tuple[str, str], dict[str, str]]:
    assays = {
        "human_atac_or_ocr_probe",
        "human_h3k27ac_probe",
        "human_h3k27me3_probe",
        "human_strict_any_regulatory_probe",
        "mouse_any_regulatory_dmr",
        "mouse_atac_dmr",
        "mouse_h3k27ac_dmr",
        "mouse_h3k27me3_dmr",
    }
    sets = {"primary_5", "final_24", "methylation_supported_final", "expression_pathway_final"}
    return {
        (row["candidate_set"], row["assay"]): row
        for row in rows
        if row.get("candidate_set") in sets and row.get("assay") in assays
    }


def mechanism_counts(mechanism_rows: list[dict[str, str]]) -> dict[str, list[str]]:
    out: dict[str, list[str]] = defaultdict(list)
    for row in mechanism_rows:
        out[row["mechanism_class"]].append(row["gene_symbol"])
    return dict(out)


def write_results_note(
    path: Path,
    summary_rows: list[dict[str, str]],
    mechanism_rows: list[dict[str, str]],
    permutation_rows: list[dict[str, str]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    perm = select_permutation_rows(permutation_rows)
    classes = mechanism_counts(mechanism_rows)
    primary = [row["gene_symbol"] for row in summary_rows if row["tier"] == "1_primary"]
    direction_complex = classes.get("directionally_complex_regulatory_methylation", [])
    repression_like = classes.get("possible_repressive_or_derepressive_methylation", [])
    polycomb = classes.get("polycomb_developmental_methylation_target", [])
    expression_unresolved = classes.get("expression_pathway_methylation_unresolved", [])

    def gene_list(genes: list[str]) -> str:
        if not genes:
            return ""
        if len(genes) == 1:
            return genes[0]
        return ", ".join(genes[:-1]) + ", and " + genes[-1]

    def ptext(set_name: str, assay: str) -> str:
        row = perm[(set_name, assay)]
        return (
            f"{row['candidate_successes']}/{row['candidate_records']}, "
            f"fold {float(row['fold_over_null']):.2f}, empirical p={float(row['empirical_p_greater']):.3g}"
        )

    lines = [
        "# Manuscript draft: candidate mechanism results",
        "",
        "## Suggested Results Subsection Title",
        "",
        "Integrated methylome, transcriptome, and chromatin analyses identify DNMT3L-linked regulatory candidates in neurodevelopment",
        "",
        "## Draft Results Text",
        "",
        "To prioritize putative downstream targets of altered DNMT3L dosage, we integrated mouse Dnmt3l-sensitive methylation regions, human DNMT3L-sensitive CpGs, DS-associated methylation concordance, fetal brain expression changes, and fetal brain chromatin annotations. This analysis identified a primary methylation-convergent candidate set consisting of "
        + gene_list(primary)
        + ". These genes were not selected because they are themselves methyltransferases or STAT-family genes; rather, they represent downstream developmental regulators or signaling genes near DNMT3L-linked methylation signals.",
        "",
        "Candidate methylation sites were strongly enriched for human fetal-brain regulatory annotations after relationship- and direction-matched permutation testing. In the primary-five set, human DNMT3L-linked CpG-gene records overlapped ATAC/OCR annotations above matched expectation ("
        + ptext("primary_5", "human_atac_or_ocr_probe")
        + "), H3K27ac ("
        + ptext("primary_5", "human_h3k27ac_probe")
        + "), and H3K27me3 ("
        + ptext("primary_5", "human_h3k27me3_probe")
        + "). The broader final 24-gene set showed similar human regulatory enrichment, including strict human triple-support probes ("
        + ptext("final_24", "human_strict_any_regulatory_probe")
        + "). In contrast, mouse Dnmt3l-sensitive DMRs supported regulatory context at individual loci but were not enriched above matched mouse DMR background for ATAC or histone-mark overlap.",
        "",
        "Methylation-expression directionality suggested more than one regulatory mechanism. "
        + gene_list(repression_like)
        + " showed the clearest direction compatible with methylation-associated repression. In contrast, "
        + gene_list(direction_complex)
        + " showed strong regulatory convergence but methylation and expression changed in the same direction in at least one major context, arguing against a simple promoter-silencing model. These genes may instead reflect enhancer methylation, gene-body methylation, indirect regulation, developmental timing, or cell-state shifts.",
        "",
        "A separate group of methylation-heavy developmental genes, including "
        + gene_list(polycomb)
        + ", showed strong regulatory or H3K27me3/Polycomb context but weaker or more context-restricted expression changes. These genes are best interpreted as developmental methylation targets rather than simple expression biomarkers. Finally, "
        + gene_list(expression_unresolved)
        + " remained expression/pathway-associated candidates whose methylation dependence is unresolved rather than proven independent.",
        "",
        "Together, these results support a model in which altered DNMT3L dosage is linked to fetal-brain regulatory methylation programs near neurodevelopmental genes. The strongest current evidence supports regulatory convergence, while causal methylation mediation remains to be tested by paired methylation-expression modeling or perturbation/rescue data.",
        "",
        "## Suggested Main Figures",
        "",
        "- Figure 1: normal and DS expression landscape of DNMT3L/DNMT3A/DNMT3B/STAT1/STAT3.",
        "- Figure 2: final candidate evidence heatmap with methylation, expression, and regulatory context.",
        "- Figure 3: locus panels for PLPP3, FOXP1, FGFR2, ZFHX3, and SIAH3.",
        "- Figure 4: matched regulatory permutation test and mechanism-class summary.",
        "- Figure 5: working model of DNMT3L-linked regulatory convergence in fetal neurodevelopment.",
        "",
        "## Suggested Tables",
        "",
        "- Table 1: datasets and analysis modules.",
        "- Table 2: manuscript candidate summary table.",
        "- Table 3: methylation-expression mechanism classification.",
        "- Supplementary Table: full final 24-gene evidence matrix with regulatory context.",
        "",
        "## Generated Output",
        "",
        "- `results/manuscript_summary/final_candidate_manuscript_summary.tsv`",
        "- `notes/manuscript_results_candidate_mechanism.md`",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    matrix_rows = read_tsv(Path(args.matrix))
    mechanism_rows = read_tsv(Path(args.mechanism))
    permutation_rows = read_tsv(Path(args.permutation))
    summary_rows = build_summary_rows(matrix_rows, mechanism_rows)

    outdir = Path(args.outdir)
    write_tsv(outdir / "final_candidate_manuscript_summary.tsv", summary_rows, SUMMARY_FIELDS)
    write_results_note(Path(args.note), summary_rows, mechanism_rows, permutation_rows)
    print(f"Wrote {len(summary_rows)} candidate manuscript rows")
    print(f"Wrote {outdir / 'final_candidate_manuscript_summary.tsv'}")
    print(f"Wrote {args.note}")


if __name__ == "__main__":
    main()
