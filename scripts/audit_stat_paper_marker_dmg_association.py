#!/usr/bin/env python3
"""Audit DNMT3L-STAT paper markers against project methylation outputs."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "mechanism_question" / "stat_paper_marker_dmg_audit.tsv"

MARKERS = [
    {
        "human": "HES3",
        "mouse": "Hes3",
        "class": "bHLH TF / Notch-related differentiation factor",
        "paper_claim": "DNMT3L-STAT1/3 activated marker; reported methylation-independent in HNP/SH-SY5Y assays",
        "motif_evidence": "No specific HES3 HOMER/chromVAR motif support in current summaries; belongs to bHLH/Notch differentiation biology.",
    },
    {
        "human": "ASCL1",
        "mouse": "Ascl1",
        "class": "proneural bHLH TF",
        "paper_claim": "DNMT3L-STAT1/3 activated marker; reported methylation-independent in HNP/SH-SY5Y assays",
        "motif_evidence": "ASCL1-like de novo motif appears in human DNMT3L-OE hyper/all DMR HOMER output; Ascl1 motif also appears in adult DS hyper DMR and mouse Dnmt3l-OE hyper DMR raw HOMER outputs.",
    },
    {
        "human": "NEUROD2",
        "mouse": "Neurod2",
        "class": "neurogenic bHLH TF",
        "paper_claim": "DNMT3L-STAT1/3 activated marker; mostly methylation-independent, with partial decitabine sensitivity noted by authors",
        "motif_evidence": "NEUROD2 motif was included in fetal chromVAR-style target motif activity summaries; not a direct top known HOMER hit in DNMT3L-OE DMR summaries.",
    },
    {
        "human": "NEUROG2",
        "mouse": "Neurog2",
        "class": "neurogenic bHLH TF",
        "paper_claim": "DNMT3L-STAT1/3 activated marker; mostly methylation-independent, with partial decitabine sensitivity noted by authors",
        "motif_evidence": "NEUROG2 motif was included in fetal chromVAR-style target motif activity summaries; neurogenic bHLH motif family appears in fetal DS scATAC/HOMER summaries.",
    },
    {
        "human": "CDKN1A",
        "mouse": "Cdkn1a",
        "class": "CDK inhibitor / p21 cell-cycle exit marker",
        "paper_claim": "DNMT3L-STAT1/3 activated marker; reported methylation-independent in HNP/SH-SY5Y assays",
        "motif_evidence": "Not a TF motif target for HOMER/chromVAR; interpreted as downstream cell-cycle effector.",
    },
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def split_gene_field(value: str) -> set[str]:
    genes: set[str] = set()
    for chunk in (value or "").replace(",", ";").split(";"):
        token = chunk.strip()
        if token:
            genes.add(token)
    return genes


def count_human_hits(rows: list[dict[str, str]], gene: str) -> int:
    count = 0
    for row in rows:
        names = set()
        for field in ("ucsc_refgene_name", "gencode_name"):
            names.update(split_gene_field(row.get(field, "")))
        if gene in names:
            count += 1
    return count


def main() -> None:
    human_dmcs = read_tsv(ROOT / "results" / "human_candidate_dmrs" / "GSE157078_candidate_dmcs.tsv")
    human_dmrs = read_tsv(ROOT / "results" / "human_candidate_dmrs" / "GSE157078_candidate_dmrs.tsv")

    mouse_overlap_rows = {
        row["gene_symbol"]: row
        for row in read_tsv(ROOT / "results" / "neurogenesis_mechanism" / "dnmt3l_vs_dnmt3a_transcriptomic_overlap.tsv")
    }

    modules_by_gene: dict[str, list[str]] = defaultdict(list)
    for row in read_tsv(ROOT / "results" / "neurogenesis_mechanism" / "developmental_module_gene_membership.tsv"):
        gene = row["gene_symbol"]
        module = row["module"]
        if module not in modules_by_gene[gene]:
            modules_by_gene[gene].append(module)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "human_gene",
        "mouse_gene",
        "marker_class",
        "paper_claim",
        "mouse_model_context",
        "human_dnmt3l_oe_direct_candidate_dmc_count",
        "human_dnmt3l_oe_direct_candidate_dmr_count",
        "mouse_dnmt3l_oe_log2fc",
        "mouse_dnmt3l_oe_fdr",
        "mouse_dnmt3l_oe_direction",
        "mouse_dnmt3a_cko_logfc",
        "mouse_dnmt3a_cko_fdr",
        "mouse_nearby_dnmt3l_dmr_count_100kb",
        "mouse_gene_body_dnmt3l_dmr_count",
        "mouse_promoter_2kb_dnmt3l_dmr_count",
        "mouse_nearby_dnmt3l_dmr_directions",
        "mouse_triple_hit_dmr_count",
        "developmental_modules",
        "motif_or_chromatin_evidence",
        "interpretation",
    ]
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for marker in MARKERS:
            mouse = marker["mouse"]
            overlap = mouse_overlap_rows.get(mouse, {})
            direct_human_dmc = count_human_hits(human_dmcs, marker["human"])
            direct_human_dmr = count_human_hits(human_dmrs, marker["human"])
            nearby = overlap.get("nearby_dnmt3l_dmr_count_100kb", "")
            gene_body = overlap.get("gene_body_dnmt3l_dmr_count", "")
            promoter = overlap.get("promoter_2kb_dnmt3l_dmr_count", "")
            human_direct = direct_human_dmc or direct_human_dmr
            mouse_local = any(int(x or 0) > 0 for x in (nearby, gene_body, promoter))
            mouse_context = (
                "GSE271949/GSE268951 P28 medial frontal cortex, Emx1 lineage; postnatal cortical maturation/synapse-glia window, "
                "not the same stage as fetal HNP neurogenesis."
            )
            if not human_direct and mouse_local:
                interpretation = (
                    "No direct human DNMT3L-OE HNP candidate DMC/DMR at this gene. Mouse P28 Dnmt3l-OE has nearby hyper-DMR support, "
                    "but that is a later maturation-stage association; expression is not FDR-significant, so it is not proof that this gene is "
                    "methylation-controlled during neurogenesis."
                )
            elif not human_direct:
                interpretation = (
                    "No direct human DNMT3L-OE HNP candidate DMC/DMR and no strong local mouse DMR-expression coupling in current tables."
                )
            else:
                interpretation = (
                    "Direct human DNMT3L-OE candidate methylation hit detected; would require locus-specific follow-up."
                )

            writer.writerow(
                {
                    "human_gene": marker["human"],
                    "mouse_gene": mouse,
                    "marker_class": marker["class"],
                    "paper_claim": marker["paper_claim"],
                    "human_dnmt3l_oe_direct_candidate_dmc_count": direct_human_dmc,
                    "human_dnmt3l_oe_direct_candidate_dmr_count": direct_human_dmr,
                    "mouse_model_context": mouse_context,
                    "mouse_dnmt3l_oe_log2fc": overlap.get("gse268951_log2fc_oe_vs_control", ""),
                    "mouse_dnmt3l_oe_fdr": overlap.get("gse268951_primary_fdr", overlap.get("gse268951_bh_fdr", "")),
                    "mouse_dnmt3l_oe_direction": overlap.get("gse268951_direction", ""),
                    "mouse_dnmt3a_cko_logfc": overlap.get("dnmt3a_logfc_cko_vs_control", ""),
                    "mouse_dnmt3a_cko_fdr": overlap.get("dnmt3a_fdr", ""),
                    "mouse_nearby_dnmt3l_dmr_count_100kb": nearby,
                    "mouse_gene_body_dnmt3l_dmr_count": gene_body,
                    "mouse_promoter_2kb_dnmt3l_dmr_count": promoter,
                    "mouse_nearby_dnmt3l_dmr_directions": overlap.get("nearby_dnmt3l_dmr_directions", ""),
                    "mouse_triple_hit_dmr_count": overlap.get("triple_hit_dmr_count", ""),
                    "developmental_modules": ";".join(modules_by_gene.get(mouse, [])),
                    "motif_or_chromatin_evidence": marker["motif_evidence"],
                    "interpretation": interpretation,
                }
            )

    print(f"Wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
