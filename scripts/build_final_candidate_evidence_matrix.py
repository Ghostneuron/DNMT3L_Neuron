#!/usr/bin/env python3
"""Build a final evidence matrix for prioritized DNMT3L-neurodevelopment genes."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


DEFAULT_HUMAN_METHYL = "results/human_methylation_convergence/human_methylation_gene_support.tsv"
DEFAULT_CROSS_SPECIES = "results/cross_species_candidate_expression/cross_species_candidate_priority.tsv"
DEFAULT_MOUSE = "results/integrated_candidate_genes/mouse_triple_hit_candidate_gene_ranking.tsv"
DEFAULT_OUTDIR = "results/final_candidate_evidence"
DEFAULT_NOTES = "notes/final_candidate_gene_evidence_matrix.md"


CANDIDATES: list[dict[str, str]] = [
    {
        "gene_symbol": "PLPP3",
        "tier": "1_primary",
        "evidence_class": "methylation_convergent",
        "manuscript_use": "Primary candidate",
        "common_theme": "developmental regulatory/signaling target",
        "special_feature": "Balanced mouse methylation, strict human methylation, and DS expression support; lipid-phosphate/vascular-glial signaling candidate.",
    },
    {
        "gene_symbol": "FOXP1",
        "tier": "1_primary",
        "evidence_class": "methylation_convergent",
        "manuscript_use": "Primary candidate",
        "common_theme": "developmental regulatory/signaling target",
        "special_feature": "Neurodevelopmental transcription factor with convergent mouse-human methylation support and expression decrease in several contexts.",
    },
    {
        "gene_symbol": "FGFR2",
        "tier": "1_primary",
        "evidence_class": "methylation_convergent",
        "manuscript_use": "Primary candidate",
        "common_theme": "developmental regulatory/signaling target",
        "special_feature": "FGF receptor connected to progenitor patterning and differentiation; supported by strict human methylation and expression evidence.",
    },
    {
        "gene_symbol": "ZFHX3",
        "tier": "1_primary",
        "evidence_class": "methylation_convergent",
        "manuscript_use": "Primary candidate",
        "common_theme": "developmental transcriptional regulation",
        "special_feature": "High integrated score with multiple strict human methylation probes; less obvious prior candidate, so useful for novelty.",
    },
    {
        "gene_symbol": "SIAH3",
        "tier": "1_primary",
        "evidence_class": "methylation_convergent",
        "manuscript_use": "Primary candidate",
        "common_theme": "developmental regulatory/signaling target",
        "special_feature": "High integrated score with promoter-linked human methylation support; potentially novel E3-ligase-family mechanism.",
    },
    {
        "gene_symbol": "FAM222A",
        "tier": "2_secondary_balanced",
        "evidence_class": "methylation_expression_supported",
        "manuscript_use": "Secondary validation candidate",
        "common_theme": "methylation-linked expression candidate",
        "special_feature": "Strong expression support plus strict human methylation support, but mouse/human directions are less simple than top-tier genes.",
    },
    {
        "gene_symbol": "COL26A1",
        "tier": "2_secondary_balanced",
        "evidence_class": "methylation_expression_supported",
        "manuscript_use": "Secondary validation candidate",
        "common_theme": "methylation-linked expression candidate",
        "special_feature": "Extracellular-matrix/cell-environment candidate with strict human methylation and strong expression support.",
    },
    {
        "gene_symbol": "KIRREL3",
        "tier": "2_secondary_balanced",
        "evidence_class": "methylation_expression_supported",
        "manuscript_use": "Secondary validation candidate",
        "common_theme": "methylation-linked expression candidate",
        "special_feature": "Synaptic/cell-adhesion gene with human methylation support and moderate expression support.",
    },
    {
        "gene_symbol": "CAMTA1",
        "tier": "2_secondary_balanced",
        "evidence_class": "methylation_expression_supported",
        "manuscript_use": "Secondary validation candidate",
        "common_theme": "methylation-linked expression candidate",
        "special_feature": "Calmodulin-responsive neurodevelopmental transcription factor with human-mouse methylation support.",
    },
    {
        "gene_symbol": "LMX1B",
        "tier": "2_secondary_balanced",
        "evidence_class": "methylation_expression_supported",
        "manuscript_use": "Secondary validation candidate",
        "common_theme": "methylation-linked expression candidate",
        "special_feature": "LIM-homeobox developmental patterning gene with strict human methylation and moderate expression support.",
    },
    {
        "gene_symbol": "OTX1",
        "tier": "3_methylation_heavy",
        "evidence_class": "methylation_heavy_expression_weak",
        "manuscript_use": "Developmental methylation target",
        "common_theme": "patterning/fate gene near methylation-sensitive regions",
        "special_feature": "Very strong strict human methylation support, but expression change is weak or not consistent in available DS data.",
    },
    {
        "gene_symbol": "OTX2",
        "tier": "3_methylation_heavy",
        "evidence_class": "methylation_heavy_expression_conflicting",
        "manuscript_use": "Developmental methylation target",
        "common_theme": "patterning/fate gene near methylation-sensitive regions",
        "special_feature": "Strong human methylation support and classic developmental patterning role; expression direction is conflicting.",
    },
    {
        "gene_symbol": "LHX2",
        "tier": "3_methylation_heavy",
        "evidence_class": "methylation_heavy_expression_weak",
        "manuscript_use": "Developmental methylation target",
        "common_theme": "patterning/fate gene near methylation-sensitive regions",
        "special_feature": "Forebrain/cortical selector gene with strict human methylation support but weak expression evidence.",
    },
    {
        "gene_symbol": "DRAXIN",
        "tier": "3_methylation_heavy",
        "evidence_class": "methylation_heavy_expression_weak",
        "manuscript_use": "Developmental methylation target",
        "common_theme": "axon guidance/developmental target near methylation-sensitive regions",
        "special_feature": "Axon guidance candidate with multiple strict human methylation probes; expression change is weak.",
    },
    {
        "gene_symbol": "NFIA",
        "tier": "3_methylation_heavy",
        "evidence_class": "methylation_heavy_expression_weak",
        "manuscript_use": "Developmental methylation target",
        "common_theme": "gliogenesis/developmental target near methylation-sensitive regions",
        "special_feature": "Glial-development transcription factor with human methylation support but limited expression evidence.",
    },
    {
        "gene_symbol": "PRDM16",
        "tier": "3_methylation_heavy",
        "evidence_class": "methylation_heavy_expression_weak",
        "manuscript_use": "Developmental methylation target",
        "common_theme": "epigenetic/developmental transcriptional regulation",
        "special_feature": "PR-domain developmental regulator with methylation support but weak expression convergence.",
    },
    {
        "gene_symbol": "NR2F1",
        "tier": "3_methylation_heavy",
        "evidence_class": "methylation_heavy_expression_weak",
        "manuscript_use": "Developmental methylation target",
        "common_theme": "arealization/developmental target near methylation-sensitive regions",
        "special_feature": "Cortical arealization nuclear receptor with strict human methylation support but weak expression convergence.",
    },
    {
        "gene_symbol": "FZD7",
        "tier": "4_expression_pathway",
        "evidence_class": "expression_pathway_methylation_unresolved",
        "manuscript_use": "Pathway/context candidate",
        "common_theme": "expression/pathway-associated developmental signal",
        "special_feature": "Strong Wnt receptor expression signal, including regional/cell-type effects; strict human methylation support was not detected by current filters.",
    },
    {
        "gene_symbol": "TGIF2",
        "tier": "4_expression_pathway",
        "evidence_class": "expression_pathway_methylation_unresolved",
        "manuscript_use": "Pathway/context candidate",
        "common_theme": "expression/pathway-associated developmental signal",
        "special_feature": "TGF-beta-related homeobox candidate with strong expression support; methylation mechanism remains unresolved.",
    },
    {
        "gene_symbol": "PRDM8",
        "tier": "4_expression_pathway",
        "evidence_class": "expression_pathway_methylation_unresolved",
        "manuscript_use": "Pathway/context candidate",
        "common_theme": "expression/pathway-associated developmental signal",
        "special_feature": "Neuronal subtype/cortical identity candidate with strong mouse and expression support but no strict human methylation triple support.",
    },
    {
        "gene_symbol": "CRABP2",
        "tier": "4_expression_pathway",
        "evidence_class": "expression_pathway_methylation_unresolved",
        "manuscript_use": "Pathway/context candidate",
        "common_theme": "retinoic-acid pathway signal",
        "special_feature": "Retinoic-acid binding candidate with mixed regional expression; methylation support is not strict in current human convergence filters.",
    },
    {
        "gene_symbol": "BCL11B",
        "tier": "4_expression_pathway",
        "evidence_class": "expression_pathway_methylation_unresolved",
        "manuscript_use": "Pathway/context candidate",
        "common_theme": "neuronal subtype identity signal",
        "special_feature": "CTIP2/layer 5-6 neuron identity gene with moderate expression support; current strict human methylation support is absent.",
    },
    {
        "gene_symbol": "RXRA",
        "tier": "4_expression_pathway",
        "evidence_class": "expression_pathway_methylation_unresolved",
        "manuscript_use": "Pathway/context candidate",
        "common_theme": "retinoic-acid pathway signal",
        "special_feature": "Retinoid receptor with mouse methylation/expression support; strict human methylation convergence was not detected.",
    },
    {
        "gene_symbol": "CYP26A1",
        "tier": "4_expression_pathway",
        "evidence_class": "expression_pathway_methylation_unresolved",
        "manuscript_use": "Pathway/context candidate",
        "common_theme": "retinoic-acid pathway signal",
        "special_feature": "Retinoic-acid catabolism gene with strong mouse evidence and moderate expression support; human methylation convergence is unresolved.",
    },
]


FIELDS = [
    "tier",
    "gene_symbol",
    "evidence_class",
    "manuscript_use",
    "mechanistic_inference",
    "mouse_symbol",
    "mouse_rank",
    "mouse_evidence_score",
    "mouse_triple_hit_dmrs",
    "mouse_promoter_2kb_dmrs",
    "mouse_gene_body_dmrs",
    "mouse_gse268951_direction",
    "mouse_gse268951_log2fc",
    "mouse_dp10_direction",
    "mouse_dp10_delta",
    "human_strict_fetal_positive_probes",
    "human_triple_any_ds_probes",
    "human_triple_fetal_ds_probes",
    "human_methylation_support_score",
    "mean_gse157078_delta_beta",
    "mean_shsy_phase3_signed_delta",
    "mean_ds_signed_delta",
    "human_expression_support_level",
    "human_expression_support_score",
    "cross_species_expression_priority_score",
    "bulk_contexts_same_direction",
    "bulk_contexts_opposite_direction",
    "celltypes_same_direction",
    "celltypes_opposite_direction",
    "gse305153_all_cells_direction",
    "gse305153_all_cells_delta",
    "best_bulk_same_direction_context",
    "best_bulk_same_direction_delta",
    "best_celltype_same_direction_label",
    "best_celltype_same_direction_delta",
    "common_theme",
    "special_feature",
    "description",
]


SUMMARY_FIELDS = [
    "evidence_class",
    "n_genes",
    "genes",
    "interpretation",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--human-methylation", default=DEFAULT_HUMAN_METHYL)
    parser.add_argument("--cross-species", default=DEFAULT_CROSS_SPECIES)
    parser.add_argument("--mouse-ranking", default=DEFAULT_MOUSE)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--notes", default=DEFAULT_NOTES)
    return parser.parse_args()


def clean_key(value: str | None) -> str:
    return (value or "").strip().replace("\r", "")


def read_tsv_by_key(path: Path, key: str) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    with path.open(newline="", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        reader.fieldnames = [clean_key(field) for field in reader.fieldnames or []]
        for raw in reader:
            row = {clean_key(k): (v or "").strip() for k, v in raw.items()}
            value = row.get(key, "")
            if value:
                rows[value.upper()] = row
    return rows


def to_float(value: str) -> float:
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


def first_value(*values: str) -> str:
    for value in values:
        if value not in ("", None):
            return value
    return ""


def infer_mechanism(row: dict[str, str]) -> str:
    strict = to_float(row.get("human_strict_fetal_positive_probes", ""))
    human_any = to_float(row.get("human_triple_any_ds_probes", ""))
    expression_level = row.get("human_expression_support_level", "")
    evidence_class = row.get("evidence_class", "")
    if evidence_class == "methylation_convergent":
        return "Strong methylation-linked/convergent candidate; paired methylation-expression mediation or perturbation validation is still needed."
    if strict > 0 and expression_level in {"strong", "moderate"}:
        return "Methylation-supported candidate with useful expression evidence, but lower priority than Tier 1."
    if strict > 0 or human_any > 0:
        return "Methylation-supported developmental target; expression effect may be subtle, cell-type-specific, or developmentally transient."
    return "Methylation-unresolved expression/pathway candidate; do not call methylation-independent without additional tests."


def build_rows(
    human_methyl: dict[str, dict[str, str]],
    cross_species: dict[str, dict[str, str]],
    mouse: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for candidate in CANDIDATES:
        gene = candidate["gene_symbol"].upper()
        hm = human_methyl.get(gene, {})
        cs = cross_species.get(gene, {})
        mouse_symbol = first_value(hm.get("mouse_gene_symbol", ""), cs.get("mouse_gene_symbol", ""))
        m = mouse.get(mouse_symbol.upper(), {}) if mouse_symbol else {}

        row = {
            **candidate,
            "mouse_symbol": mouse_symbol,
            "mouse_rank": first_value(hm.get("mouse_rank", ""), cs.get("mouse_rank", "")),
            "mouse_evidence_score": first_value(hm.get("mouse_evidence_score", ""), cs.get("mouse_evidence_score", ""), m.get("evidence_score", "")),
            "mouse_triple_hit_dmrs": m.get("n_triple_hit_dmrs", ""),
            "mouse_promoter_2kb_dmrs": m.get("n_promoter_2kb_dmrs", ""),
            "mouse_gene_body_dmrs": m.get("n_gene_body_dmrs", ""),
            "mouse_gse268951_direction": first_value(hm.get("mouse_gse268951_direction", ""), cs.get("mouse_gse268951_direction", ""), m.get("gse268951_expression_direction", "")),
            "mouse_gse268951_log2fc": first_value(cs.get("mouse_gse268951_log2fc", ""), m.get("gse268951_log2fc_oe_vs_control", "")),
            "mouse_dp10_direction": first_value(hm.get("mouse_dp10_direction", ""), cs.get("mouse_dp10_direction", ""), m.get("dp10_expression_direction", "")),
            "mouse_dp10_delta": first_value(cs.get("mouse_dp10_delta", ""), m.get("dp10_delta_log2_cpm_plus_1", "")),
            "human_strict_fetal_positive_probes": hm.get("n_human_triple_phase3_fetal_mean_positive_probes", "0"),
            "human_triple_any_ds_probes": hm.get("n_human_triple_phase3_any_ds_probes", "0"),
            "human_triple_fetal_ds_probes": hm.get("n_human_triple_phase3_fetal_ds_probes", "0"),
            "human_methylation_support_score": hm.get("human_methylation_support_score", "0"),
            "mean_gse157078_delta_beta": hm.get("mean_gse157078_delta_beta", ""),
            "mean_shsy_phase3_signed_delta": hm.get("mean_shsy_phase3_signed_delta", ""),
            "mean_ds_signed_delta": hm.get("mean_ds_signed_delta", ""),
            "human_expression_support_level": first_value(hm.get("human_expression_support_level", ""), cs.get("human_support_level", "")),
            "human_expression_support_score": first_value(hm.get("human_expression_support_score", ""), cs.get("human_support_score", "")),
            "cross_species_expression_priority_score": first_value(hm.get("cross_species_expression_priority_score", ""), cs.get("cross_species_priority_score", "")),
            "bulk_contexts_same_direction": cs.get("bulk_contexts_same_direction", ""),
            "bulk_contexts_opposite_direction": cs.get("bulk_contexts_opposite_direction", ""),
            "celltypes_same_direction": cs.get("celltypes_same_direction", ""),
            "celltypes_opposite_direction": cs.get("celltypes_opposite_direction", ""),
            "gse305153_all_cells_direction": cs.get("gse305153_all_cells_direction", ""),
            "gse305153_all_cells_delta": cs.get("gse305153_all_cells_delta", ""),
            "best_bulk_same_direction_context": cs.get("best_bulk_same_direction_context", ""),
            "best_bulk_same_direction_delta": cs.get("best_bulk_same_direction_delta", ""),
            "best_celltype_same_direction_label": cs.get("best_celltype_same_direction_label", ""),
            "best_celltype_same_direction_delta": cs.get("best_celltype_same_direction_delta", ""),
            "description": first_value(hm.get("description", ""), cs.get("description", ""), m.get("description", "")),
        }
        row["mechanistic_inference"] = infer_mechanism(row)
        rows.append(row)
    return rows


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})


def summarize(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    interpretations = {
        "methylation_convergent": "Best current evidence for DNMT3L/DNMT3A/DNMT3B methylation-linked regulation.",
        "methylation_expression_supported": "Good secondary genes with both methylation and expression support, but weaker or less direct than Tier 1.",
        "methylation_heavy_expression_weak": "Likely developmental methylation targets; expression may require cell-type, region, or timing-specific tests.",
        "methylation_heavy_expression_conflicting": "Methylation-supported developmental target with conflicting expression direction.",
        "expression_pathway_methylation_unresolved": "Expression/pathway candidates; current filters do not prove methylation independence.",
    }
    by_class: dict[str, list[str]] = {}
    for row in rows:
        by_class.setdefault(row["evidence_class"], []).append(row["gene_symbol"])
    return [
        {
            "evidence_class": evidence_class,
            "n_genes": str(len(genes)),
            "genes": ",".join(genes),
            "interpretation": interpretations.get(evidence_class, ""),
        }
        for evidence_class, genes in by_class.items()
    ]


def write_notes(path: Path, rows: list[dict[str, str]], summary_rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    primary = [row["gene_symbol"] for row in rows if row["tier"] == "1_primary"]
    secondary = [row["gene_symbol"] for row in rows if row["tier"] == "2_secondary_balanced"]
    methylation_heavy = [row["gene_symbol"] for row in rows if row["tier"] == "3_methylation_heavy"]
    expression = [row["gene_symbol"] for row in rows if row["tier"] == "4_expression_pathway"]

    lines = [
        "# Final candidate gene evidence matrix",
        "",
        "This note summarizes the current computational prioritization of DNMT3L-related neurodevelopmental candidate genes. The matrix separates methylation-convergent genes from expression/pathway genes because absence of support in the current methylation filters does not prove methylation independence.",
        "",
        "## Final tiers",
        "",
        f"- Primary methylation-convergent candidates: {', '.join(primary)}.",
        f"- Secondary methylation + expression candidates: {', '.join(secondary)}.",
        f"- Methylation-heavy developmental targets: {', '.join(methylation_heavy)}.",
        f"- Expression/pathway-associated, methylation-unresolved candidates: {', '.join(expression)}.",
        "",
        "## Common features",
        "",
        "- Most candidates are downstream neurodevelopmental regulators rather than DNMT3L, DNMT3A, DNMT3B, STAT1, or STAT3 themselves.",
        "- The recurring biology is developmental patterning, progenitor/neuron fate, axon guidance, glial/vascular or extracellular-matrix context, and Wnt/FGF/TGF/retinoic-acid signaling.",
        "- The methylation-convergent candidates have mouse Dnmt3l-overexpression DMR support plus human DNMT3L-sensitive/DS methylation support.",
        "- The expression/pathway candidates should be called methylation-unresolved, not methylation-independent, until broader enhancer/TAD methylation scans and methylation-expression mediation analyses are done.",
        "",
        "## Evidence class summary",
        "",
        "| Evidence class | n | Genes | Interpretation |",
        "|---|---:|---|---|",
    ]
    for row in summary_rows:
        genes = row["genes"].replace(",", ", ")
        lines.append(f"| {row['evidence_class']} | {row['n_genes']} | {genes} | {row['interpretation']} |")
    lines.extend(
        [
            "",
            "## Practical manuscript wording",
            "",
            "A conservative conclusion is that DNMT3L dosage converges on developmentally methylated regulatory regions near neurodevelopmental regulators, with strongest current support for PLPP3, FOXP1, FGFR2, ZFHX3, and SIAH3. Additional expression/pathway candidates such as FZD7, TGIF2, PRDM8, and CRABP2 may reflect downstream developmental state changes, but methylation independence has not been proven.",
            "",
            "## Outputs",
            "",
            "- `results/final_candidate_evidence/final_candidate_evidence_matrix.tsv`",
            "- `results/final_candidate_evidence/final_candidate_evidence_summary.tsv`",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    human_methyl = read_tsv_by_key(Path(args.human_methylation), "gene_symbol")
    cross_species = read_tsv_by_key(Path(args.cross_species), "human_symbol")
    mouse = read_tsv_by_key(Path(args.mouse_ranking), "gene_symbol")

    rows = build_rows(human_methyl, cross_species, mouse)
    summary_rows = summarize(rows)

    outdir = Path(args.outdir)
    write_tsv(outdir / "final_candidate_evidence_matrix.tsv", rows, FIELDS)
    write_tsv(outdir / "final_candidate_evidence_summary.tsv", summary_rows, SUMMARY_FIELDS)
    write_notes(Path(args.notes), rows, summary_rows)

    print(f"Wrote {len(rows)} candidate rows to {outdir / 'final_candidate_evidence_matrix.tsv'}")
    print(f"Wrote {len(summary_rows)} summary rows to {outdir / 'final_candidate_evidence_summary.tsv'}")
    print(f"Wrote notes to {args.notes}")


if __name__ == "__main__":
    main()
