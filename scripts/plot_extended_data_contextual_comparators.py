#!/usr/bin/env python3
"""Build contextual Drosophila and postnatal DS PFC extended-data figure."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np
import pandas as pd


FIGDIR = Path("figures/contextual_comparators")
OUTBASE = FIGDIR / "extended_data_drosophila_postnatal_context"
SUMMARY = Path("results/manuscript_summary/drosophila_postnatal_context_summary.tsv")

DROSOPHILA = [
    {
        "accession": "GSE263069",
        "label": "GSE263069\nlarval brain RNA-seq",
        "neural_specificity": 3,
        "assay_directness": 3,
        "mammalian_methylation_relevance": 0,
        "pathway_comparator_value": 3,
        "note": "human DNMT3L-OE in elavC155-Gal4 larval brains",
    },
    {
        "accession": "GSE68861",
        "label": "GSE68861\nwhole-larva array",
        "neural_specificity": 1,
        "assay_directness": 1,
        "mammalian_methylation_relevance": 0,
        "pathway_comparator_value": 1,
        "note": "whole-larva transgenic/generational DNMT3L expression series",
    },
]


def read_optional(path: Path) -> pd.DataFrame:
    if path.exists():
        return pd.read_csv(path, sep="\t")
    return pd.DataFrame()


def write_summary(target_df: pd.DataFrame, deg_df: pd.DataFrame, audit_df: pd.DataFrame) -> None:
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)

    def audit_value(metric: str, fallback: str = "NA") -> str:
        if audit_df.empty:
            return fallback
        hits = audit_df.loc[audit_df["metric"].eq(metric), "value"]
        return str(hits.iloc[0]) if not hits.empty else fallback

    total_rows = "NA"
    total_matches = "NA"
    if not deg_df.empty:
        total_rows = str(int(pd.to_numeric(deg_df["rows_loaded"], errors="coerce").fillna(0).sum()))
        total_matches = str(
            int(pd.to_numeric(deg_df["target_or_candidate_matches"], errors="coerce").fillna(0).sum())
        )

    dnmt3l_summary = "NA"
    stat3_summary = "NA"
    if not target_df.empty:
        dnmt3l = target_df.loc[target_df["gene_symbol"].eq("DNMT3L"), "interpretation"]
        stat3 = target_df.loc[target_df["gene_symbol"].eq("STAT3"), "interpretation"]
        if not dnmt3l.empty:
            dnmt3l_summary = str(dnmt3l.iloc[0])
        if not stat3.empty:
            stat3_summary = str(stat3.iloc[0])

    rows = [
        {
            "context": "Drosophila DNMT3L-OE",
            "resource": "GSE263069",
            "assay_or_data": "third instar larval brain RNA-seq",
            "public_or_local_status": "candidate metadata in config; raw expression not yet processed locally",
            "use_in_manuscript": "Tier 4 methylation-independent neural developmental pathway comparator",
            "current_conclusion": "cleaner Drosophila comparator because it uses larval brains and neuronal-driver DNMT3L-OE",
            "caveat": "cannot test mammalian DNMT3A/DNMT3B-mediated DNA methylation; ortholog/pathway analysis remains a future step",
        },
        {
            "context": "Drosophila DNMT3L-OE",
            "resource": "GSE68861",
            "assay_or_data": "whole-larva Affymetrix expression array",
            "public_or_local_status": "candidate metadata in config; raw expression not yet processed locally",
            "use_in_manuscript": "Lower-priority methylation-independent pathway comparator",
            "current_conclusion": "less direct for neurodevelopment because it is whole-larva and generational",
            "caveat": "tissue mixture and design make it unsuitable for primary brain/neuron claims",
        },
        {
            "context": "Risgaard postnatal DS PFC",
            "resource": audit_value("current_dbgap_accession", "phs004098.v2.p1"),
            "assay_or_data": "single-cell/single-nucleus genomic data",
            "public_or_local_status": audit_value("current_dbgap_status", "controlled_access_metadata_ready"),
            "use_in_manuscript": "Contextual postnatal PFC resource; full reanalysis requires dbGaP authorized access",
            "current_conclusion": "metadata support dataset relevance, but public files do not provide full local matrices",
            "caveat": "do not claim donor-level or cell-type statistical reanalysis from controlled-access data",
        },
        {
            "context": "Risgaard postnatal DS PFC",
            "resource": f"Zenodo {audit_value('zenodo_record_id', '17527842')}",
            "assay_or_data": "public code/readme/figure-generation archive",
            "public_or_local_status": (
                f"{audit_value('downloaded_file_count', '4')} files downloaded; "
                f"{audit_value('available_not_downloaded_file_count', '1')} listed but not downloaded"
            ),
            "use_in_manuscript": "Methods and availability audit; not a replacement for the controlled genomic matrices",
            "current_conclusion": (
                f"code sample audit detected {audit_value('code_sample_ids_detected', '51')} sample IDs, "
                f"including {audit_value('ts21_code_sample_ids_detected', '24')} Ts21 and "
                f"{audit_value('control_code_sample_ids_detected', '27')} control IDs"
            ),
            "caveat": "public code does not by itself expose all expression or ATAC count matrices",
        },
        {
            "context": "Risgaard postnatal DS PFC",
            "resource": "DS_PFC Shiny target-gene plot",
            "assay_or_data": "public visual expression summary for target genes",
            "public_or_local_status": "visual audit completed from saved Shiny plot",
            "use_in_manuscript": "Contextual target-gene expression reference",
            "current_conclusion": f"DNMT3L: {dnmt3l_summary}; STAT3: {stat3_summary}",
            "caveat": "visual inspection only; raw per-donor values and model output require controlled data or downloadable summaries",
        },
        {
            "context": "Risgaard postnatal DS PFC",
            "resource": "DS_PFC Shiny strict DEG tables",
            "assay_or_data": "public Shiny DEG output under strict default filters",
            "public_or_local_status": f"{total_rows} strict DEG rows loaded across audited queries",
            "use_in_manuscript": "Negative/limited public-audit context for target genes and final candidates",
            "current_conclusion": f"{total_matches} target-gene or 24-candidate matches under strict default filters",
            "caveat": "strict default filters do not rule out weaker, subtype-specific, or donor-level effects",
        },
    ]
    pd.DataFrame(rows).to_csv(SUMMARY, sep="\t", index=False)


def score_target_expression(summary: str) -> tuple[float, str]:
    text = str(summary).lower()
    if "very low" in text or "near-baseline" in text:
        return 0.12, "very low"
    if "rare" in text:
        return 0.20, "rare"
    if "broad" in text and "stronger" in text:
        return 0.72, "broad, selected high"
    if "broad" in text or "clearly detectable" in text or "detected broadly" in text:
        return 0.65, "broad"
    return 0.35, "detected"


def main() -> None:
    FIGDIR.mkdir(parents=True, exist_ok=True)
    target_df = read_optional(Path("results/postnatal_pfc_science_aea1549/shiny_target_expression_visual_audit.tsv"))
    deg_df = read_optional(Path("results/postnatal_pfc_science_aea1549/shiny_deg_strict_query_summary.tsv"))
    audit_df = read_optional(Path("results/postnatal_pfc_science_aea1549/dataset_audit_summary.tsv"))
    write_summary(target_df, deg_df, audit_df)

    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.facecolor": "white",
        }
    )

    fig = plt.figure(figsize=(13.5, 8.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.05, 1.0], height_ratios=[1.0, 1.0])

    ax = fig.add_subplot(gs[0, 0])
    score_cols = [
        "neural_specificity",
        "assay_directness",
        "mammalian_methylation_relevance",
        "pathway_comparator_value",
    ]
    labels = [
        "Neural\nspecificity",
        "Assay\nmatch",
        "Mammalian\nmethylation",
        "Pathway\ncomparator",
    ]
    dros = pd.DataFrame(DROSOPHILA)
    mat = dros[score_cols].to_numpy(dtype=float)
    cmap = ListedColormap(["#F2F2F2", "#C7D9DD", "#73A9AD", "#2F6F8F"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)
    ax.imshow(mat, aspect="auto", cmap=cmap, norm=norm)
    ax.set_xticks(np.arange(len(labels)), labels)
    ax.set_yticks(np.arange(len(dros)), dros["label"])
    ax.set_title("A  Drosophila DNMT3L-OE comparator value", loc="left", fontweight="bold")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            ax.text(j, i, f"{int(mat[i, j])}", ha="center", va="center", fontsize=9, fontweight="bold")
    ax.set_xlabel(
        "Scores: 0 absent/not useful, 3 strongest; no mammalian DNMT3A/DNMT3B methylation system",
        fontsize=8,
        color="0.25",
    )

    ax = fig.add_subplot(gs[0, 1])
    resources = [
        ("dbGaP metadata", 1, "ready"),
        ("Full count/ATAC matrices", 0, "controlled"),
        ("Zenodo code archive", 1, "ready"),
        ("Shiny target-gene plot", 1, "ready"),
        ("Strict Shiny DEG audit", 1, "ready"),
    ]
    status_mat = np.array([[r[1] for r in resources]])
    status_cmap = ListedColormap(["#D9B26F", "#4C78A8"])
    ax.imshow(status_mat, aspect="auto", cmap=status_cmap, vmin=0, vmax=1)
    ax.set_yticks([0], ["Risgaard\npostnatal PFC"])
    ax.set_xticks(np.arange(len(resources)), [r[0] for r in resources], rotation=30, ha="right")
    ax.set_title("B  Public availability boundary", loc="left", fontweight="bold")
    for j, (_, value, label) in enumerate(resources):
        ax.text(j, 0, label, ha="center", va="center", fontsize=8, fontweight="bold", color="white" if value else "black")
    ax.set_ylim(0.5, -0.8)
    ax.set_xlabel(
        "Public Shiny/code are contextual; full matrices require dbGaP access",
        fontsize=8,
        color="0.25",
    )

    ax = fig.add_subplot(gs[1, 0])
    genes = ["DNMT3L", "DNMT3A", "DNMT3B", "STAT1", "STAT3"]
    if target_df.empty:
        scores = [np.nan] * len(genes)
        labels = ["not audited"] * len(genes)
    else:
        lookup = target_df.set_index("gene_symbol")["visual_expression_summary"].to_dict()
        scored = [score_target_expression(lookup.get(g, "")) for g in genes]
        scores = [s for s, _ in scored]
        labels = [lab for _, lab in scored]
    y = np.arange(len(genes))
    colors = ["#C44E52" if g == "STAT3" else "#4C78A8" if g.startswith("DNMT") else "#59A14F" for g in genes]
    ax.barh(y, scores, color=colors)
    ax.set_yticks(y, genes)
    ax.set_xlim(0, 1.0)
    ax.set_xlabel("Qualitative detection score from public Shiny plot")
    ax.set_title("C  Postnatal PFC target-gene visual audit", loc="left", fontweight="bold")
    ax.invert_yaxis()
    for i, (score, label) in enumerate(zip(scores, labels)):
        ax.text(min(score + 0.03, 0.92), i, label, va="center", fontsize=8)

    ax = fig.add_subplot(gs[1, 1])
    if deg_df.empty:
        query_ids = ["no DEG audit"]
        rows = [0]
        matches = [0]
    else:
        query_ids = [q.replace("_default_strict", "").replace("_", "\n") for q in deg_df["query_id"]]
        rows = pd.to_numeric(deg_df["rows_loaded"], errors="coerce").fillna(0).to_numpy()
        matches = pd.to_numeric(deg_df["target_or_candidate_matches"], errors="coerce").fillna(0).to_numpy()
    x = np.arange(len(query_ids))
    ax.bar(x, rows, color="#8E6C8A", label="strict DEG rows loaded")
    ax.scatter(x, matches, color="#E45756", s=42, zorder=3, label="target/candidate matches")
    ax.set_xticks(x, query_ids)
    ax.set_ylabel("Rows or matches")
    ax.set_title("D  Public Shiny strict DEG audit", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=8)
    for i, value in enumerate(matches):
        ax.text(i, max(rows) * 0.06 if len(rows) else 0.5, f"matches {int(value)}", ha="center", fontsize=8)

    fig.suptitle(
        "Extended Data Figure 9. Drosophila and early postnatal DS PFC contextual resources",
        fontweight="bold",
    )
    fig.savefig(OUTBASE.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(OUTBASE.with_suffix(".pdf"), bbox_inches="tight")
    print(OUTBASE.with_suffix(".png"))
    print(OUTBASE.with_suffix(".pdf"))
    print(SUMMARY)


if __name__ == "__main__":
    main()
