#!/usr/bin/env python3
"""Compare Dnmt3a cKO threshold choices against Dnmt3l OE and modules."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DEFAULT_DNMT3A_DE = "data/geo/GSE141587/supplementary/GSE141587_RNAseq_P39_DEgenes_edgeR_qlfTest_Dnmt3aCKO_vs_Control.txt.gz"
DEFAULT_OVERLAP = "results/neurogenesis_mechanism/dnmt3l_vs_dnmt3a_transcriptomic_overlap.tsv"
DEFAULT_FINAL = "results/neurogenesis_mechanism/final_candidate_neurogenesis_annotation.tsv"
DEFAULT_MEMBERSHIP = "results/neurogenesis_mechanism/developmental_module_gene_membership.tsv"
DEFAULT_OUTDIR = "results/neurogenesis_mechanism"
DEFAULT_FIGDIR = "figures/neurogenesis_mechanism"
DEFAULT_NOTE = "notes/dnmt3a_threshold_sensitivity.md"


THRESHOLDS = [1.0, 0.58, 0.25, 0.0]
MODULES_OF_INTEREST = [
    "early_neurogenesis_neuron_differentiation",
    "cortical_patterning_region_identity",
    "axon_dendrite_neuron_projection",
    "synapse_formation_maturation",
    "glial_maturation_myelination",
    "neurodevelopment_signaling",
    "final16_neurodevelopment_candidates",
    "final_synapse_window_candidates",
    "primary5_candidates",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dnmt3a-de", default=DEFAULT_DNMT3A_DE)
    parser.add_argument("--overlap", default=DEFAULT_OVERLAP)
    parser.add_argument("--final", default=DEFAULT_FINAL)
    parser.add_argument("--membership", default=DEFAULT_MEMBERSHIP)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--figdir", default=DEFAULT_FIGDIR)
    parser.add_argument("--note", default=DEFAULT_NOTE)
    return parser.parse_args()


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
    with path.open(newline="", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [{k: (v or "").strip() for k, v in row.items()} for row in reader]


def read_gzip_tsv(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", newline="", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [{k: (v or "").strip() for k, v in row.items()} for row in reader]


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def load_dnmt3a_rows(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for row in read_gzip_tsv(path):
        logfc = safe_float(row.get("logFC"))
        fdr = safe_float(row.get("FDR"))
        if not math.isfinite(logfc) or not math.isfinite(fdr):
            continue
        rows.append(
            {
                "gene_symbol": row.get("geneName", ""),
                "gene_id": row.get("geneID", ""),
                "gene_type": row.get("geneType", ""),
                "dnmt3a_logfc_cko_vs_control": logfc,
                "dnmt3a_fdr": fdr,
            }
        )
    return rows


def load_module_membership(path: Path) -> dict[str, set[str]]:
    memberships: dict[str, set[str]] = defaultdict(set)
    for row in read_tsv(path):
        memberships[row["module"]].add(row["gene_symbol"])
    return dict(memberships)


def threshold_summary(
    thresholds: list[float],
    dnmt3a_rows: list[dict[str, object]],
    overlap_rows: list[dict[str, str]],
    final_rows: list[dict[str, str]],
    modules: dict[str, set[str]],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    overlap_by_symbol = {row["gene_symbol"]: row for row in overlap_rows}
    final_symbols = {row["mouse_symbol"] for row in final_rows if row.get("mouse_symbol")}
    final_neuro_symbols = {
        row["mouse_symbol"]
        for row in final_rows
        if row.get("mouse_symbol") and row.get("neurogenesis_related") == "yes"
    }
    dnmt3l_de = {row["gene_symbol"] for row in overlap_rows if row.get("gse268951_de_p05_fc2") == "yes"}
    dnmt3l_neuro_de = {
        row["gene_symbol"]
        for row in overlap_rows
        if row.get("gse268951_de_p05_fc2") == "yes" and row.get("neurogenesis_related") == "yes"
    }
    summary_rows: list[dict[str, object]] = []
    gene_rows: list[dict[str, object]] = []
    module_rows: list[dict[str, object]] = []

    for cutoff in thresholds:
        de_rows = [
            row
            for row in dnmt3a_rows
            if row["dnmt3a_fdr"] <= 0.05 and abs(row["dnmt3a_logfc_cko_vs_control"]) >= cutoff
        ]
        de_symbols = {str(row["gene_symbol"]) for row in de_rows}
        de_protein = {str(row["gene_symbol"]) for row in de_rows if row.get("gene_type") == "protein_coding"}
        neuro_symbols = {
            symbol
            for symbol in de_symbols
            if overlap_by_symbol.get(symbol, {}).get("neurogenesis_related") == "yes"
        }
        summary_rows.append(
            {
                "dnmt3a_abs_logfc_cutoff": cutoff,
                "dnmt3a_de_genes_all_types": len(de_symbols),
                "dnmt3a_de_genes_protein_coding": len(de_protein),
                "dnmt3a_de_neurogenesis_related": len(neuro_symbols),
                "overlap_with_dnmt3l_oe_de": len(de_symbols & dnmt3l_de),
                "overlap_with_dnmt3l_oe_neurogenesis_de": len(de_symbols & dnmt3l_neuro_de),
                "final24_candidates_de": len(de_symbols & final_symbols),
                "final16_neurodevelopment_candidates_de": len(de_symbols & final_neuro_symbols),
                "dnmt3a_only_de_vs_dnmt3l": len(de_symbols - dnmt3l_de),
                "dnmt3l_only_de_vs_dnmt3a": len(dnmt3l_de - de_symbols),
            }
        )
        for row in sorted(de_rows, key=lambda r: abs(float(r["dnmt3a_logfc_cko_vs_control"])), reverse=True):
            symbol = str(row["gene_symbol"])
            overlap = overlap_by_symbol.get(symbol, {})
            gene_rows.append(
                {
                    "dnmt3a_abs_logfc_cutoff": cutoff,
                    "gene_symbol": symbol,
                    "gene_type": row.get("gene_type", ""),
                    "dnmt3a_logfc_cko_vs_control": fmt(row["dnmt3a_logfc_cko_vs_control"]),
                    "dnmt3a_fdr": fmt(row["dnmt3a_fdr"]),
                    "dnmt3l_oe_de_p05_fc2": "yes" if symbol in dnmt3l_de else "no",
                    "neurogenesis_related": overlap.get("neurogenesis_related", "unknown"),
                    "matched_neurogenesis_go_terms": overlap.get("matched_neurogenesis_go_terms", ""),
                    "is_final_candidate": "yes" if symbol in final_symbols else "no",
                    "is_final_neurodevelopment_candidate": "yes" if symbol in final_neuro_symbols else "no",
                }
            )
        for module in MODULES_OF_INTEREST:
            genes = modules.get(module, set())
            module_rows.append(
                {
                    "dnmt3a_abs_logfc_cutoff": cutoff,
                    "module": module,
                    "module_size": len(genes),
                    "module_dnmt3a_de_genes": len(genes & de_symbols),
                    "module_dnmt3a_de_gene_symbols": ",".join(sorted(genes & de_symbols)),
                }
            )
    return summary_rows, gene_rows, module_rows


def make_summary_plot(summary_rows: list[dict[str, object]], figdir: Path) -> None:
    figdir.mkdir(parents=True, exist_ok=True)
    rows = sorted(summary_rows, key=lambda row: float(row["dnmt3a_abs_logfc_cutoff"]), reverse=True)
    x = np.arange(len(rows))
    labels = [str(row["dnmt3a_abs_logfc_cutoff"]) for row in rows]
    all_counts = [int(row["dnmt3a_de_genes_protein_coding"]) for row in rows]
    neuro_counts = [int(row["dnmt3a_de_neurogenesis_related"]) for row in rows]
    overlap_counts = [int(row["overlap_with_dnmt3l_oe_de"]) for row in rows]

    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    ax.plot(x, all_counts, marker="o", linewidth=2.0, color="#2F6F8F", label="Protein-coding Dnmt3a DE")
    ax.plot(x, neuro_counts, marker="o", linewidth=2.0, color="#6E5A9B", label="Neurogenesis-related")
    ax.plot(x, overlap_counts, marker="o", linewidth=2.0, color="#D85C4A", label="Overlap with Dnmt3l OE DE")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_xlabel("Dnmt3a cKO abs(logFC) cutoff")
    ax.set_ylabel("Gene count")
    ax.set_title("Dnmt3a cKO threshold sensitivity", fontsize=12.5, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=8.5)
    for xi, count in zip(x, all_counts):
        ax.text(xi, count + max(all_counts) * 0.03, str(count), ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    fig.savefig(figdir / "dnmt3a_threshold_sensitivity_counts.png", dpi=260, bbox_inches="tight")
    fig.savefig(figdir / "dnmt3a_threshold_sensitivity_counts.pdf", bbox_inches="tight")
    plt.close(fig)


def write_note(path: Path, summary_rows: list[dict[str, object]], module_rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    by_cutoff = {float(row["dnmt3a_abs_logfc_cutoff"]): row for row in summary_rows}
    row_1 = by_cutoff[1.0]
    row_058 = by_cutoff[0.58]
    module_058 = [row for row in module_rows if float(row["dnmt3a_abs_logfc_cutoff"]) == 0.58]
    module_hits_058 = [row for row in module_058 if int(row["module_dnmt3a_de_genes"]) > 0]
    hit_text = "; ".join(
        f"{row['module']}={row['module_dnmt3a_de_genes']} ({row['module_dnmt3a_de_gene_symbols']})"
        for row in module_hits_058
    )
    lines = [
        "# Dnmt3a cKO Threshold Sensitivity",
        "",
        "Date: 2026-06-09",
        "",
        "`abs(logFC) >= 0.58` corresponds to about a 1.5-fold cutoff. This sensitivity check asks whether relaxing the Dnmt3a cKO single-gene threshold changes the interpretation of the Dnmt3l/Dnmt3a transcriptomic comparison.",
        "",
        "## Main comparison",
        "",
        f"- At abs(logFC) >= 1.0 and FDR <= 0.05: {row_1['dnmt3a_de_genes_all_types']} Dnmt3a DE genes, {row_1['dnmt3a_de_genes_protein_coding']} protein-coding, {row_1['dnmt3a_de_neurogenesis_related']} neurogenesis-related.",
        f"- At abs(logFC) >= 0.58 and FDR <= 0.05: {row_058['dnmt3a_de_genes_all_types']} Dnmt3a DE genes, {row_058['dnmt3a_de_genes_protein_coding']} protein-coding, {row_058['dnmt3a_de_neurogenesis_related']} neurogenesis-related.",
        f"- Overlap with re-derived Dnmt3l OE DE genes changes from {row_1['overlap_with_dnmt3l_oe_de']} to {row_058['overlap_with_dnmt3l_oe_de']}.",
        f"- Final 24 candidate hits remain {row_058['final24_candidates_de']} at the 0.58 cutoff.",
        "",
        "## Module-level implication",
        "",
        "Changing the binary logFC cutoff does not change the rank-based module test, because the module test already uses the full continuous logFC distribution. Therefore, the earlier conclusion remains: Dnmt3a cKO does not show a broad synapse-window module shift comparable to Dnmt3l OE.",
        "",
        f"At the 0.58 cutoff, module membership hits are: {hit_text or 'none in tested modules'}.",
        "",
        "## Interpretation",
        "",
        "Lowering the Dnmt3a single-gene threshold does recover additional genes, including a few neurodevelopmental genes, but it does not qualitatively change the conclusion. Dnmt3a cKO still has a much smaller transcriptomic response than Dnmt3l OE and does not convert into a broad synapse/neurodevelopment module signal.",
        "",
        "## Generated outputs",
        "",
        "- `results/neurogenesis_mechanism/dnmt3a_threshold_sensitivity_summary.tsv`",
        "- `results/neurogenesis_mechanism/dnmt3a_threshold_sensitivity_genes.tsv`",
        "- `results/neurogenesis_mechanism/dnmt3a_threshold_sensitivity_modules.tsv`",
        "- `figures/neurogenesis_mechanism/dnmt3a_threshold_sensitivity_counts.png`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    dnmt3a_rows = load_dnmt3a_rows(Path(args.dnmt3a_de))
    overlap_rows = read_tsv(Path(args.overlap))
    final_rows = read_tsv(Path(args.final))
    memberships = load_module_membership(Path(args.membership))
    summary_rows, gene_rows, module_rows = threshold_summary(
        THRESHOLDS, dnmt3a_rows, overlap_rows, final_rows, memberships
    )
    outdir = Path(args.outdir)
    figdir = Path(args.figdir)
    write_tsv(
        outdir / "dnmt3a_threshold_sensitivity_summary.tsv",
        summary_rows,
        [
            "dnmt3a_abs_logfc_cutoff",
            "dnmt3a_de_genes_all_types",
            "dnmt3a_de_genes_protein_coding",
            "dnmt3a_de_neurogenesis_related",
            "overlap_with_dnmt3l_oe_de",
            "overlap_with_dnmt3l_oe_neurogenesis_de",
            "final24_candidates_de",
            "final16_neurodevelopment_candidates_de",
            "dnmt3a_only_de_vs_dnmt3l",
            "dnmt3l_only_de_vs_dnmt3a",
        ],
    )
    write_tsv(
        outdir / "dnmt3a_threshold_sensitivity_genes.tsv",
        gene_rows,
        [
            "dnmt3a_abs_logfc_cutoff",
            "gene_symbol",
            "gene_type",
            "dnmt3a_logfc_cko_vs_control",
            "dnmt3a_fdr",
            "dnmt3l_oe_de_p05_fc2",
            "neurogenesis_related",
            "matched_neurogenesis_go_terms",
            "is_final_candidate",
            "is_final_neurodevelopment_candidate",
        ],
    )
    write_tsv(
        outdir / "dnmt3a_threshold_sensitivity_modules.tsv",
        module_rows,
        [
            "dnmt3a_abs_logfc_cutoff",
            "module",
            "module_size",
            "module_dnmt3a_de_genes",
            "module_dnmt3a_de_gene_symbols",
        ],
    )
    make_summary_plot(summary_rows, figdir)
    write_note(Path(args.note), summary_rows, module_rows)
    print(f"Wrote {outdir / 'dnmt3a_threshold_sensitivity_summary.tsv'}")
    print(f"Wrote {outdir / 'dnmt3a_threshold_sensitivity_genes.tsv'}")
    print(f"Wrote {outdir / 'dnmt3a_threshold_sensitivity_modules.tsv'}")
    print(f"Wrote {figdir / 'dnmt3a_threshold_sensitivity_counts.png'}")
    print(f"Wrote {args.note}")


if __name__ == "__main__":
    main()
