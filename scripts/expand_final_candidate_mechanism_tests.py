#!/usr/bin/env python3
"""Expand candidate summaries with explicit in-silico mechanism tests."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DEFAULT_SUMMARY = "results/manuscript_summary/final_candidate_manuscript_summary.tsv"
DEFAULT_MECHANISM = "results/mechanism_classification/methylation_expression_mechanism_table.tsv"
DEFAULT_MATRIX = "results/final_candidate_evidence/final_candidate_evidence_matrix_with_regulatory_context.tsv"
DEFAULT_OUTDIR = "results/manuscript_summary"
DEFAULT_FIGDIR = "figures/manuscript_summary"
DEFAULT_NOTE = "notes/expanded_candidate_in_silico_mechanism_tests.md"


EXTRA_FIELDS = [
    "proposed_in_silico_mechanism",
    "test_priority",
    "core_testable_hypothesis",
    "methylation_tests",
    "chromatin_tests",
    "expression_celltype_tests",
    "cross_species_tests",
    "pathway_network_tests",
    "expected_supporting_pattern",
    "key_counter_evidence",
    "main_caveat",
    "recommended_manuscript_use",
]


MECHANISM_COLUMNS = [
    "repression_or_derepression",
    "enhancer_gene_body_complex",
    "polycomb_h3k27me3",
    "developmental_timing",
    "pathway_network",
    "cell_type_specificity",
    "cross_species_support",
    "regulatory_chromatin",
]


PATHWAY_GROUPS = {
    "FZD7": "Wnt/progenitor signaling",
    "TGIF2": "TGF-beta/SMAD and progenitor signaling",
    "PRDM8": "neuronal maturation and transcriptional repression",
    "CRABP2": "retinoic acid transport/signaling",
    "BCL11B": "corticofugal/deep-layer neuronal identity",
    "RXRA": "retinoid receptor signaling",
    "CYP26A1": "retinoic acid catabolism",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", default=DEFAULT_SUMMARY)
    parser.add_argument("--mechanism", default=DEFAULT_MECHANISM)
    parser.add_argument("--matrix", default=DEFAULT_MATRIX)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--figdir", default=DEFAULT_FIGDIR)
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


def has_repression_relation(row: dict[str, str]) -> bool:
    relations = " ".join(
        [
            row.get("human_methyl_expr_relation", ""),
            row.get("mouse_oe_methyl_expr_relation", ""),
            row.get("mouse_dp10_methyl_expr_relation", ""),
        ]
    )
    return "opposite_possible_repression" in relations or "opposite_possible_derepression" in relations


def has_complex_relation(row: dict[str, str]) -> bool:
    relations = " ".join(
        [
            row.get("human_methyl_expr_relation", ""),
            row.get("mouse_oe_methyl_expr_relation", ""),
            row.get("mouse_dp10_methyl_expr_relation", ""),
        ]
    )
    return "same_direction" in relations or "complex" in relations


def has_celltype_context(row: dict[str, str]) -> bool:
    context = row.get("human_expression_context", "")
    if not context:
        return False
    return context not in {"GSE305153_all_cells", "fetal_brain_final"}


def regulatory_score(row: dict[str, str], matrix: dict[str, str]) -> float:
    context = row.get("regulatory_context_class", "") or matrix.get("regulatory_context_class", "")
    if "human_mouse_regulatory_supported" in context:
        return 1.0
    if "human_regulatory_supported" in context or "mouse_regulatory_supported" in context:
        return 0.75
    if "possible" in context:
        return 0.5
    return 0.0


def cross_species_score(row: dict[str, str], matrix: dict[str, str]) -> float:
    context = row.get("regulatory_context_class", "") or matrix.get("regulatory_context_class", "")
    if "human_mouse_regulatory_supported" in context:
        return 1.0
    if row.get("human_methylation_direction") and row.get("mouse_methylation_direction"):
        return 0.5
    return 0.0


def priority_for(row: dict[str, str], matrix: dict[str, str]) -> str:
    tier = row["tier"]
    if tier == "1_primary":
        return "primary"
    if tier == "2_secondary_balanced":
        return "high_secondary"
    if row["mechanism_class"] == "polycomb_developmental_methylation_target":
        return "high_developmental"
    if tier == "4_expression_pathway" and has_repression_relation(row):
        return "medium_pathway_repression_like"
    return "medium_context"


def mechanism_for(gene: str, row: dict[str, str], matrix: dict[str, str]) -> str:
    cls = row["mechanism_class"]
    if cls == "possible_repressive_or_derepressive_methylation":
        return "repression_like_methylation"
    if cls == "directionally_complex_regulatory_methylation":
        return "enhancer_gene_body_or_cell_state_methylation"
    if cls == "secondary_regulatory_methylation_candidate":
        if has_repression_relation(row) and not has_complex_relation(row):
            return "secondary_repression_like_methylation"
        if "opposite_possible_derepression" in row.get("human_methyl_expr_relation", ""):
            return "secondary_context_dependent_derepression"
        if has_repression_relation(row) and has_complex_relation(row):
            return "secondary_mixed_repression_and_complex_regulation"
        return "secondary_enhancer_gene_body_or_cell_state_methylation"
    if cls == "polycomb_developmental_methylation_target":
        return "polycomb_developmental_timing_methylation"
    if cls == "developmental_methylation_target_expression_weak":
        return "developmental_methylation_priming_with_weak_expression"
    if cls == "expression_pathway_methylation_unresolved":
        return f"pathway_expression_candidate_{PATHWAY_GROUPS.get(gene, 'methylation_unresolved')}"
    return cls


def hypothesis_for(gene: str, row: dict[str, str], mechanism: str) -> str:
    if "repression_like" in mechanism:
        return (
            f"{gene} may be downregulated or derepressed when nearby DNMT3L-linked methylation changes occur "
            "at promoter-proximal or regulatory CpG/DMR features."
        )
    if "polycomb" in mechanism:
        return (
            f"{gene} may be part of a Polycomb-associated developmental regulatory program in which DNMT3L-linked "
            "methylation marks poised or timing-sensitive chromatin rather than producing a large steady-state expression shift."
        )
    if "priming" in mechanism:
        return (
            f"{gene} may carry developmental methylation changes that prime later expression or cell-state potential, "
            "with weak expression effects in currently available datasets."
        )
    if "pathway_expression" in mechanism:
        pathway = PATHWAY_GROUPS.get(gene, "neurodevelopmental pathway")
        return (
            f"{gene} may report a {pathway} response downstream of DS/DNMT3L perturbation; current data should test "
            "whether nearby methylation support is present, absent, or only mouse-specific."
        )
    return (
        f"{gene} may be regulated through enhancer, gene-body, developmental timing, or cell-composition effects rather "
        "than a simple promoter methylation repression model."
    )


def methylation_tests_for(row: dict[str, str], matrix: dict[str, str]) -> str:
    tests = [
        "map linked CpGs/DMRs to promoter_2kb, gene_body, and nearest_100kb classes",
        "compare DNMT3L/Dnmt3l perturbation direction with DS signed methylation direction",
    ]
    if row.get("human_methylation_direction") in {"hyper", "hypo"}:
        tests.append("stratify human DNMT3L-sensitive probes by delta-beta magnitude")
    if row.get("mouse_methylation_direction") in {"hyper", "hypo"}:
        tests.append("summarize mouse DMR delta-beta and triple-hit overlap strength")
    if matrix.get("human_strict_fetal_positive_probe_count", "0") != "0":
        tests.append("test strict fetal-positive probe concordance")
    return "; ".join(tests)


def chromatin_tests_for(row: dict[str, str], matrix: dict[str, str]) -> str:
    tests = [
        "intersect linked features with fetal brain OCR/ATAC and H3K27ac",
        "test H3K27me3/Polycomb overlap",
    ]
    if row["mechanism_class"] in {"polycomb_developmental_methylation_target", "developmental_methylation_target_expression_weak"}:
        tests.append("prioritize bivalent or poised promoter/enhancer interpretation")
    if has_complex_relation(row):
        tests.append("separate promoter features from enhancer/gene-body features")
    return "; ".join(tests)


def expression_tests_for(row: dict[str, str]) -> str:
    context = row.get("human_expression_context", "")
    tests = [
        "compare DS/control direction within each dataset instead of pooling absolute expression units",
        "test donor/sample-level pseudobulk effects in fetal cortex single-cell data",
    ]
    if context:
        tests.append(f"focus on strongest current human context: {context}")
    if row.get("mouse_gse268951_direction") and row.get("mouse_gse268951_direction") != "flat":
        tests.append("compare mouse Dnmt3l OE expression direction")
    if row.get("mouse_dp10_direction") and row.get("mouse_dp10_direction") != "flat":
        tests.append("compare Dp10 expression direction as dosage-model context")
    return "; ".join(tests)


def cross_species_tests_for(row: dict[str, str], matrix: dict[str, str]) -> str:
    tests = [
        "compare human gene symbol with mouse ortholog symbol",
        "check whether methylation direction and expression direction agree across human DS, mouse Dnmt3l OE, and Dp10",
    ]
    if "human_mouse_regulatory_supported" in (row.get("regulatory_context_class", "") or matrix.get("regulatory_context_class", "")):
        tests.append("require both human and mouse regulatory overlap for strongest support")
    return "; ".join(tests)


def pathway_tests_for(gene: str, row: dict[str, str], mechanism: str) -> str:
    tests = [
        "test co-expression with candidate set and target-panel genes",
        "run pathway or gene-set enrichment in DS and DNMT3L/Dnmt3l perturbation contrasts",
    ]
    if gene in PATHWAY_GROUPS:
        tests.append(f"test {PATHWAY_GROUPS[gene]} gene-set connectivity")
    if "polycomb" in mechanism:
        tests.append("compare with cortical patterning and neuronal differentiation gene sets")
    return "; ".join(tests)


def expected_pattern_for(gene: str, row: dict[str, str], mechanism: str) -> str:
    if "repression_like" in mechanism:
        return "opposite methylation-expression direction, regulatory overlap, and consistent effect direction in at least two relevant contexts"
    if "polycomb" in mechanism:
        return "strong methylation plus H3K27me3/poised chromatin overlap, with expression restricted to developmental stage or cell type"
    if "priming" in mechanism:
        return "methylation and regulatory support without a large bulk expression shift, but with stage or cell-type restricted expression"
    if "pathway_expression" in mechanism:
        return "coherent expression/pathway shift, with methylation support treated as present, mouse-specific, or unresolved"
    return "methylation features in enhancer/gene-body/regulatory chromatin plus same-direction or context-dependent expression changes"


def counter_evidence_for(row: dict[str, str], mechanism: str) -> str:
    if "repression_like" in mechanism:
        return "same-direction methylation-expression change, lack of promoter/regulatory overlap, or inconsistent donor-level effect"
    if "polycomb" in mechanism or "priming" in mechanism:
        return "absence of H3K27me3/poised chromatin support or no developmental/cell-type expression pattern"
    if "pathway_expression" in mechanism:
        return "pathway expression shift without any linked methylation or regulatory support in human brain data"
    return "only bulk expression change without linked methylation/regulatory support, or effects explained by cell composition alone"


def caveat_for(row: dict[str, str], mechanism: str) -> str:
    if row["tier"] == "4_expression_pathway":
        return "methylation dependence is unresolved; do not describe as methylation-independent without negative methylation evidence in matched cells"
    if "polycomb" in mechanism or "priming" in mechanism:
        return "weak expression does not exclude developmental function, but it limits direct expression-mediated claims"
    if row["tier"] == "2_secondary_balanced":
        return "secondary priority because cross-species directionality, regulatory overlap, or expression context is less complete than primary genes"
    return "requires paired methylation-expression or perturbation validation for causal interpretation"


def manuscript_use_for(row: dict[str, str], mechanism: str) -> str:
    if row["tier"] == "1_primary":
        return "main candidate mechanism panel"
    if row["tier"] == "2_secondary_balanced":
        return "secondary validation or supplementary candidate table"
    if "polycomb" in mechanism:
        return "developmental/Polycomb mechanism subsection"
    if row["tier"] == "4_expression_pathway":
        return "pathway/methylation-unresolved subsection"
    return "supplementary mechanism table"


def mechanism_scores(row: dict[str, str], matrix: dict[str, str], mechanism: str) -> dict[str, float]:
    h3k27me3 = safe_float(row.get("mouse_embryonic_forebrain_h3k27me3_dmr_count"))
    if not math.isfinite(h3k27me3):
        h3k27me3 = safe_float(matrix.get("mouse_embryonic_forebrain_h3k27me3_dmr_count"))
    return {
        "repression_or_derepression": 1.0 if has_repression_relation(row) else 0.0,
        "enhancer_gene_body_complex": 1.0 if has_complex_relation(row) else 0.0,
        "polycomb_h3k27me3": 1.0
        if row["mechanism_class"] == "polycomb_developmental_methylation_target"
        else (0.5 if h3k27me3 > 0 else 0.0),
        "developmental_timing": 1.0
        if row["tier"] == "3_methylation_heavy" or "developmental" in mechanism
        else (0.5 if row["tier"] in {"1_primary", "2_secondary_balanced"} else 0.0),
        "pathway_network": 1.0 if row["tier"] == "4_expression_pathway" else (0.5 if row["tier"] == "3_methylation_heavy" else 0.0),
        "cell_type_specificity": 1.0 if has_celltype_context(row) else 0.5,
        "cross_species_support": cross_species_score(row, matrix),
        "regulatory_chromatin": regulatory_score(row, matrix),
    }


def build_rows(
    summary_rows: list[dict[str, str]],
    mechanism_rows: list[dict[str, str]],
    matrix_rows: list[dict[str, str]],
) -> tuple[list[dict[str, object]], list[dict[str, float]]]:
    mechanism_by_gene = {row["gene_symbol"]: row for row in mechanism_rows}
    matrix_by_gene = {row["gene_symbol"]: row for row in matrix_rows}
    expanded: list[dict[str, object]] = []
    score_rows: list[dict[str, float]] = []

    for summary in summary_rows:
        gene = summary["gene_symbol"]
        mechanism_row = mechanism_by_gene.get(gene, {})
        matrix = matrix_by_gene.get(gene, {})
        merged = {**summary, **mechanism_row}
        mechanism = mechanism_for(gene, merged, matrix)
        extra = {
            "proposed_in_silico_mechanism": mechanism,
            "test_priority": priority_for(merged, matrix),
            "core_testable_hypothesis": hypothesis_for(gene, merged, mechanism),
            "methylation_tests": methylation_tests_for(merged, matrix),
            "chromatin_tests": chromatin_tests_for(merged, matrix),
            "expression_celltype_tests": expression_tests_for(merged),
            "cross_species_tests": cross_species_tests_for(merged, matrix),
            "pathway_network_tests": pathway_tests_for(gene, merged, mechanism),
            "expected_supporting_pattern": expected_pattern_for(gene, merged, mechanism),
            "key_counter_evidence": counter_evidence_for(merged, mechanism),
            "main_caveat": caveat_for(merged, mechanism),
            "recommended_manuscript_use": manuscript_use_for(merged, mechanism),
        }
        out = {**summary, **extra}
        expanded.append(out)
        scores = mechanism_scores(merged, matrix, mechanism)
        score_rows.append({"gene_symbol": gene, **scores})

    return expanded, score_rows


def make_heatmap(score_rows: list[dict[str, float]], expanded_rows: list[dict[str, object]], figdir: Path) -> None:
    figdir.mkdir(parents=True, exist_ok=True)
    expanded_by_gene = {str(row["gene_symbol"]): row for row in expanded_rows}
    rows = [row for row in score_rows if expanded_by_gene[row["gene_symbol"]]["tier"] != "1_primary"]
    tier_order = {"2_secondary_balanced": 0, "3_methylation_heavy": 1, "4_expression_pathway": 2}
    rows.sort(key=lambda row: (tier_order.get(str(expanded_by_gene[row["gene_symbol"]]["tier"]), 9), row["gene_symbol"]))

    matrix = np.array([[float(row[col]) for col in MECHANISM_COLUMNS] for row in rows], dtype=float)
    labels = [row["gene_symbol"] for row in rows]
    col_labels = [
        "Repression/\nderepression",
        "Enhancer/\ngene body",
        "Polycomb/\nH3K27me3",
        "Developmental\ntiming",
        "Pathway/\nnetwork",
        "Cell type\nspecificity",
        "Cross-species\nsupport",
        "Regulatory\nchromatin",
    ]

    fig, ax = plt.subplots(figsize=(11.4, 7.4))
    im = ax.imshow(matrix, vmin=0, vmax=1, cmap="YlGnBu", aspect="auto")
    ax.set_xticks(np.arange(len(col_labels)))
    ax.set_xticklabels(col_labels, fontsize=8.5)
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels, fontsize=9)
    ax.tick_params(axis="x", top=True, bottom=False, labeltop=True, labelbottom=False, length=0)
    ax.tick_params(axis="y", length=0)

    for i, gene in enumerate(labels):
        tier = str(expanded_by_gene[gene]["tier"])
        if tier == "2_secondary_balanced":
            color = "#D85C4A"
        elif tier == "3_methylation_heavy":
            color = "#6E5A9B"
        else:
            color = "#4A7C59"
        ax.scatter(-0.62, i, marker="s", s=24, color=color, clip_on=False)

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            if value >= 0.99:
                text = "yes"
            elif value > 0:
                text = "partial"
            else:
                text = ""
            color = "white" if value >= 0.75 else "#111827"
            ax.text(j, i, text, ha="center", va="center", fontsize=7.2, color=color)

    ax.set_title("In-silico mechanism tests for non-primary candidates", fontsize=13, fontweight="bold", pad=64)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks(np.arange(-0.5, len(col_labels), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(labels), 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=1.2)
    ax.tick_params(which="minor", bottom=False, left=False)

    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.025)
    cbar.set_ticks([0, 0.5, 1])
    cbar.set_ticklabels(["no", "partial", "yes"])
    cbar.ax.tick_params(labelsize=8)

    fig.subplots_adjust(left=0.13, right=0.91, top=0.78, bottom=0.11)
    fig.text(
        0.13,
        0.055,
        "Row markers: red = secondary methylation-expression; purple = developmental/Polycomb; green = expression/pathway methylation-unresolved.",
        ha="left",
        va="center",
        fontsize=8.3,
        color="#4B5563",
    )
    fig.savefig(figdir / "nonprimary_candidate_in_silico_mechanism_test_map.png", dpi=260, bbox_inches="tight")
    fig.savefig(figdir / "nonprimary_candidate_in_silico_mechanism_test_map.pdf", bbox_inches="tight")
    plt.close(fig)


def grouped_gene_list(rows: list[dict[str, object]], tier: str) -> str:
    genes = [str(row["gene_symbol"]) for row in rows if row["tier"] == tier]
    return ", ".join(f"`{gene}`" for gene in genes)


def write_note(path: Path, expanded_rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    secondary = [row for row in expanded_rows if row["tier"] == "2_secondary_balanced"]
    developmental = [row for row in expanded_rows if row["tier"] == "3_methylation_heavy"]
    pathway = [row for row in expanded_rows if row["tier"] == "4_expression_pathway"]

    lines = [
        "# Expanded In-Silico Mechanism Tests For Final Candidates",
        "",
        "Date: 2026-06-09",
        "",
        "This note expands the existing `final_candidate_manuscript_summary.tsv` into explicit computational mechanism tests. It does not create a new candidate set; it adds testable hypotheses, readouts, and caveats for the same final 24 genes.",
        "",
        "Generated outputs:",
        "",
        "- `results/manuscript_summary/final_candidate_in_silico_mechanism_tests.tsv`",
        "- `results/manuscript_summary/nonprimary_candidate_in_silico_mechanism_tests.tsv`",
        "- `figures/manuscript_summary/nonprimary_candidate_in_silico_mechanism_test_map.png`",
        "- `figures/manuscript_summary/nonprimary_candidate_in_silico_mechanism_test_map.pdf`",
        "",
        "## Secondary methylation-expression candidates",
        "",
        f"Genes: {grouped_gene_list(expanded_rows, '2_secondary_balanced')}.",
        "",
        "These genes have methylation and expression support, but are lower priority than the five primary genes because cross-species directionality, regulatory support, or expression context is less complete. The main in-silico tests are methylation-expression direction checks, promoter versus enhancer/gene-body localization, fetal regulatory chromatin overlap, and donor/sample-level expression consistency.",
        "",
        "| Gene | Proposed mechanism | Highest-value in-silico test | Main caveat |",
        "| --- | --- | --- | --- |",
    ]
    for row in secondary:
        lines.append(
            f"| `{row['gene_symbol']}` | {row['proposed_in_silico_mechanism']} | {row['expected_supporting_pattern']} | {row['main_caveat']} |"
        )

    lines.extend(
        [
            "",
            "## Developmental and Polycomb methylation candidates",
            "",
            f"Genes: {grouped_gene_list(expanded_rows, '3_methylation_heavy')}.",
            "",
            "These genes are best handled as developmental regulatory targets rather than simple expression biomarkers. The main in-silico tests are H3K27me3/Polycomb overlap, bivalent or poised chromatin context, fetal stage or pseudotime restriction, and cell-type-specific expression.",
            "",
            "| Gene | Proposed mechanism | Highest-value in-silico test | Main caveat |",
            "| --- | --- | --- | --- |",
        ]
    )
    for row in developmental:
        lines.append(
            f"| `{row['gene_symbol']}` | {row['proposed_in_silico_mechanism']} | {row['expected_supporting_pattern']} | {row['main_caveat']} |"
        )

    lines.extend(
        [
            "",
            "## Expression/pathway methylation-unresolved candidates",
            "",
            f"Genes: {grouped_gene_list(expanded_rows, '4_expression_pathway')}.",
            "",
            "These genes should not be called methylation-independent. The stronger statement is that their methylation dependence is unresolved in the current data. The main in-silico tests are pathway connectivity, co-expression with DS/DNMT3L signatures, fetal cell-type specificity, and whether human methylation support emerges when less stringent or cell-type-matched methylation data are used.",
            "",
            "| Gene | Proposed mechanism | Highest-value in-silico test | Main caveat |",
            "| --- | --- | --- | --- |",
        ]
    )
    for row in pathway:
        lines.append(
            f"| `{row['gene_symbol']}` | {row['proposed_in_silico_mechanism']} | {row['expected_supporting_pattern']} | {row['main_caveat']} |"
        )

    lines.extend(
        [
            "",
            "## Practical interpretation",
            "",
            "The expanded table is most useful as a supplementary mechanism roadmap. It tells the reader what each non-primary gene is being used for: secondary methylation-expression validation, developmental/Polycomb timing, or expression/pathway context. For the paper, the safest language is that these candidates define testable computational mechanisms and validation priorities, not proven downstream targets.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    summary_rows = read_tsv(Path(args.summary))
    mechanism_rows = read_tsv(Path(args.mechanism))
    matrix_rows = read_tsv(Path(args.matrix))
    expanded_rows, score_rows = build_rows(summary_rows, mechanism_rows, matrix_rows)

    summary_fields = list(summary_rows[0].keys())
    full_fields = summary_fields + EXTRA_FIELDS
    outdir = Path(args.outdir)
    full_path = outdir / "final_candidate_in_silico_mechanism_tests.tsv"
    nonprimary_path = outdir / "nonprimary_candidate_in_silico_mechanism_tests.tsv"
    score_path = outdir / "candidate_mechanism_test_scores.tsv"

    write_tsv(full_path, expanded_rows, full_fields)
    write_tsv(nonprimary_path, [row for row in expanded_rows if row["tier"] != "1_primary"], full_fields)
    write_tsv(score_path, score_rows, ["gene_symbol", *MECHANISM_COLUMNS])
    make_heatmap(score_rows, expanded_rows, Path(args.figdir))
    write_note(Path(args.note), expanded_rows)

    print(f"Wrote {full_path}")
    print(f"Wrote {nonprimary_path}")
    print(f"Wrote {score_path}")
    print(f"Wrote {Path(args.figdir) / 'nonprimary_candidate_in_silico_mechanism_test_map.png'}")
    print(f"Wrote {Path(args.figdir) / 'nonprimary_candidate_in_silico_mechanism_test_map.pdf'}")
    print(f"Wrote {args.note}")


if __name__ == "__main__":
    main()
