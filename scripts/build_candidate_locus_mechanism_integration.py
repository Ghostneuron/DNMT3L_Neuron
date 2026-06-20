#!/usr/bin/env python3
"""Integrate candidate-locus evidence for the DNMT3L mechanism question.

This script intentionally keeps three evidence layers separate:
1. candidate-linked methylation/expression evidence from existing gene/probe tables,
2. local fetal DS scATAC support from hg38 peak/window overlap,
3. context-level HOMER motif-family enrichment from altered-region sets.

The HOMER layer is not treated as exact motif-instance evidence at a candidate
locus unless a separate motif-scanning step is performed.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap


FINAL_SUMMARY = Path("results/manuscript_summary/final_candidate_manuscript_summary.tsv")
FINAL_EVIDENCE = Path("results/final_candidate_evidence/final_candidate_evidence_matrix_with_regulatory_context.tsv")
FINAL_TESTS = Path("results/manuscript_summary/final_candidate_in_silico_mechanism_tests.tsv")
HUMAN_DMR = Path("results/human_candidate_dmrs/GSE157078_candidate_dmrs.tsv")
HUMAN_DMC = Path("results/human_candidate_dmrs/GSE157078_candidate_dmcs.tsv")
DS_PROBE_VALUES = Path("results/ds_methylation_dnmt3l_concordance/ds_dnmt3l_candidate_probe_values.tsv")
MOUSE_DMR_LINKS = Path("results/integrated_candidate_genes/mouse_triple_hit_dmr_gene_links.tsv")
MOUSE_DNMT3A_PRIORITY = Path("results/mouse_dnmt3l_dnmt3a_neuron_overlap/gse271949_dmr_priority_triple_hits.tsv")
ATAC_WINDOWS = Path("results/mechanism_question/gse305146_accessibility_proxy_gene_windows.tsv")
ATAC_EFFECTS = Path("results/mechanism_question/gse305146_accessibility_proxy_effect_summary.tsv")
DS_OPEN_PEAKS = Path("results/gse305146_scATAC_homer_motif_enrichment/tables/GSE305146_DS_open_recurrent_peak_calls.tsv")
CONTROL_OPEN_PEAKS = Path(
    "results/gse305146_scATAC_homer_motif_enrichment/tables/GSE305146_Control_open_recurrent_peak_calls.tsv"
)
MOTIF_GROUPS = Path("results/mechanism_question/homer_cross_context_motif_family_group_summary.tsv")

DEVELOPMENTAL_FAMILIES = {
    "HOMEBOX_ANY",
    "LHX_HOMEBOX",
    "DLX_HOMEBOX",
    "NKX_HOMEBOX",
    "GSX_HOMEBOX",
    "EMX_HOMEBOX",
    "ISL_HOMEBOX",
    "POU_HOMEBOX",
    "NEURO_BHLH",
    "BHLH_ANY",
    "OLIG_BHLH",
    "SOX_HMG",
    "RFX_XBOX",
    "MEF2",
}
CPG_CHROMATIN_FAMILIES = {"KLF_SP_GC", "CTCF_BORIS", "ZINC_FINGER_OTHER"}
STAT_FAMILIES = {"STAT"}

ACTIVE_SUPPORT_CALLS = {"suggestive_only", "significant_single_context", "strong_single_context", "recurrent_significant"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="results/mechanism_question")
    parser.add_argument("--figdir", default="figures/mechanism_question")
    return parser.parse_args()


def read_tsv(path: Path, dtype: object | None = None) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t", dtype=dtype)


def fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        return f"{value:.6g}"
    if pd.isna(value):
        return ""
    return str(value)


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fieldnames})


def to_float(value: object) -> float:
    if value is None or pd.isna(value):
        return float("nan")
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def to_int(value: object) -> int:
    if value is None or pd.isna(value) or value == "":
        return 0
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def split_list(value: object) -> list[str]:
    if value is None or pd.isna(value):
        return []
    text = str(value).strip()
    if not text or text.lower() == "nan" or text in {"0", "0.0"}:
        return []
    return [part.strip() for part in re.split(r"[;,]", text) if part.strip()]


def annotation_tokens(row: pd.Series, columns: list[str]) -> set[str]:
    tokens: set[str] = set()
    for column in columns:
        if column not in row or pd.isna(row[column]):
            continue
        text = str(row[column])
        text = text.replace("///", ";").replace("|", ";").replace(",", ";")
        for token in text.split(";"):
            token = token.strip()
            if token and token.lower() not in {"nan", "na", "none"}:
                tokens.add(token.upper())
    return tokens


def add_gene_tokens(df: pd.DataFrame) -> pd.DataFrame:
    columns = [column for column in ["ucsc_refgene_name", "gencode_name"] if column in df.columns]
    out = df.copy()
    out["_gene_tokens"] = [annotation_tokens(row, columns) for _, row in out.iterrows()]
    return out


def filter_gene_annotated(df: pd.DataFrame, gene: str) -> pd.DataFrame:
    gene_key = gene.upper()
    return df[df["_gene_tokens"].map(lambda tokens: gene_key in tokens)].copy()


def direction_counts(df: pd.DataFrame) -> tuple[int, int, str]:
    if "direction" not in df.columns or df.empty:
        return 0, 0, ""
    counts = df["direction"].fillna("unknown").astype(str).value_counts().to_dict()
    hyper = int(counts.get("hyper", 0))
    hypo = int(counts.get("hypo", 0))
    summary = ";".join(f"{key}:{counts[key]}" for key in sorted(counts))
    return hyper, hypo, summary


def summarize_human_dnmt3l(df: pd.DataFrame, gene: str, id_col: str, delta_col: str) -> dict[str, object]:
    matched = filter_gene_annotated(df, gene)
    hyper, hypo, dir_summary = direction_counts(matched)
    ids = sorted({str(value) for value in matched[id_col].dropna().tolist()}) if id_col in matched.columns else []
    return {
        "count": len(ids) if ids else int(len(matched)),
        "hyper_count": hyper,
        "hypo_count": hypo,
        "direction_summary": dir_summary,
        "mean_delta_beta": matched[delta_col].astype(float).mean() if delta_col in matched.columns and not matched.empty else np.nan,
        "max_abs_delta_beta": matched[delta_col].astype(float).abs().max()
        if delta_col in matched.columns and not matched.empty
        else np.nan,
        "ids": ";".join(ids[:25]),
    }


def summarize_ds_probes(df: pd.DataFrame, gene: str, stage_regex: str) -> dict[str, object]:
    matched = filter_gene_annotated(df, gene)
    stage = matched[matched["context"].astype(str).str.contains(stage_regex, case=False, regex=True, na=False)].copy()
    if stage.empty:
        return {
            "probe_count": 0,
            "context_count": 0,
            "mean_ds_delta_beta": np.nan,
            "mean_signed_ds_delta_by_dnmt3l_direction": np.nan,
            "concordant_fraction": np.nan,
            "positive_context_count": 0,
            "negative_context_count": 0,
            "contexts": "",
        }
    context_means = stage.groupby("context")["ds_delta_beta"].mean()
    concordant = stage["concordant"].astype(str).str.lower().eq("yes")
    contexts = sorted(stage["context"].dropna().astype(str).unique().tolist())
    return {
        "probe_count": int(stage["probe_id"].nunique()),
        "context_count": int(stage["context"].nunique()),
        "mean_ds_delta_beta": stage["ds_delta_beta"].astype(float).mean(),
        "mean_signed_ds_delta_by_dnmt3l_direction": stage["signed_ds_delta_by_dnmt3l_direction"].astype(float).mean(),
        "concordant_fraction": float(concordant.mean()),
        "positive_context_count": int((context_means > 0).sum()),
        "negative_context_count": int((context_means < 0).sum()),
        "contexts": ";".join(contexts[:20]),
    }


def summarize_mouse_gene(
    links: pd.DataFrame, priority: pd.DataFrame, mouse_symbol: str
) -> dict[str, object]:
    if not mouse_symbol:
        return {}
    matched = links[links["gene_symbol"].astype(str).str.upper() == mouse_symbol.upper()].copy()
    priority_ids = set(priority["dmr_id"].astype(str).tolist()) if not priority.empty else set()
    dmr_ids = sorted({str(value) for value in matched["dmr_id"].dropna().tolist()})
    in_priority = [dmr_id for dmr_id in dmr_ids if dmr_id in priority_ids]
    priority_subset = priority[priority["dmr_id"].astype(str).isin(dmr_ids)].copy() if not priority.empty else pd.DataFrame()
    dnmt3a_ref_col = "dnmt3a_p39_cko_hypo_cg_dmrs_overlap"
    h3k27me3_col = "h3k27me3_p39_cko_up_overlap"
    return {
        "mouse_dnmt3l_dmr_count": len(dmr_ids),
        "mouse_dnmt3l_promoter_2kb_dmr_count": int(matched.get("promoter_2kb", pd.Series(dtype=str)).astype(str).eq("yes").sum()),
        "mouse_dnmt3l_gene_body_dmr_count": int(matched.get("gene_body", pd.Series(dtype=str)).astype(str).eq("yes").sum()),
        "mouse_dnmt3l_nearest_100kb_dmr_count": int(matched.get("nearest_100kb", pd.Series(dtype=str)).astype(str).eq("yes").sum()),
        "mouse_dnmt3l_mean_delta_beta": matched["dmr_mean_delta_beta"].astype(float).mean()
        if "dmr_mean_delta_beta" in matched.columns and not matched.empty
        else np.nan,
        "mouse_dnmt3a_reference_overlap_dmr_count": len(in_priority),
        "mouse_dnmt3a_p39_cko_hypo_overlap_dmr_count": int(
            priority_subset.get(dnmt3a_ref_col, pd.Series(dtype=str)).astype(str).eq("yes").sum()
        ),
        "mouse_h3k27me3_p39_cko_up_overlap_dmr_count": int(
            priority_subset.get(h3k27me3_col, pd.Series(dtype=str)).astype(str).eq("yes").sum()
        ),
        "mouse_dnmt3l_dmr_ids": ";".join(dmr_ids[:25]),
    }


def overlap_peak_summary(peaks: pd.DataFrame, window: pd.Series) -> dict[str, object]:
    if window is None or window.empty:
        return {"count": 0, "mean_delta_frequency": np.nan, "max_abs_delta_frequency": np.nan, "intervals": ""}
    chrom = str(window["chrom"])
    start = int(window["window_start"])
    end = int(window["window_end"])
    subset = peaks[(peaks["chrom"].astype(str) == chrom) & (peaks["end"].astype(int) >= start) & (peaks["start"].astype(int) <= end)]
    intervals = [
        f"{row.chrom}:{int(row.start)}-{int(row.end)}"
        for row in subset.sort_values("delta_frequency_ds_minus_control", ascending=False).head(20).itertuples()
    ]
    return {
        "count": int(len(subset)),
        "mean_delta_frequency": subset["delta_frequency_ds_minus_control"].astype(float).mean() if not subset.empty else np.nan,
        "max_abs_delta_frequency": subset["delta_frequency_ds_minus_control"].astype(float).abs().max() if not subset.empty else np.nan,
        "intervals": ";".join(intervals),
    }


def supported_motif_rows(motifs: pd.DataFrame, analysis_group: str) -> pd.DataFrame:
    subset = motifs[
        (motifs["analysis_group"].astype(str) == analysis_group)
        & (motifs["support_call"].astype(str).isin(ACTIVE_SUPPORT_CALLS))
        & (motifs["motif_family"].astype(str) != "OTHER")
    ].copy()
    if subset.empty:
        return subset
    return subset.sort_values(["best_q_value", "best_rank"], ascending=[True, True])


def active_contexts(row: dict[str, object]) -> dict[str, str]:
    contexts: dict[str, str] = {}
    direct_hyper = to_int(row["human_dnmt3l_oe_hyper_dmr_count"]) + to_int(row["human_dnmt3l_oe_hyper_dmc_count"])
    direct_hypo = to_int(row["human_dnmt3l_oe_hypo_dmr_count"]) + to_int(row["human_dnmt3l_oe_hypo_dmc_count"])
    curated_dnmt3l = to_int(row["human_dnmt3l_oe_curated_probe_count"])
    curated_dnmt3l_delta = to_float(row["human_dnmt3l_oe_curated_mean_delta_beta"])

    if direct_hyper > 0:
        contexts["human_DNMT3L_OE_methylation|hyper"] = "candidate annotated to DNMT3L-OE hyper-DMR/DMC"
    if direct_hypo > 0:
        contexts["human_DNMT3L_OE_methylation|hypo"] = "candidate annotated to DNMT3L-OE hypo-DMR/DMC"
    if direct_hyper + direct_hypo == 0 and curated_dnmt3l > 0:
        contexts[
            "human_DNMT3L_OE_methylation|hyper" if curated_dnmt3l_delta >= 0 else "human_DNMT3L_OE_methylation|hypo"
        ] = "candidate linked by curated DNMT3L-sensitive probe evidence matrix"

    fetal_delta = to_float(row["human_fetal_ds_mean_delta_beta"])
    if to_int(row["human_fetal_ds_effective_probe_count"]) > 0:
        if math.isnan(fetal_delta):
            fetal_delta = to_float(row["human_ds_curated_mean_signed_delta_beta"])
        contexts[
            "human_fetal_DS_methylation|hyper" if fetal_delta >= 0 else "human_fetal_DS_methylation|hypo"
        ] = "candidate-linked fetal DS methylation probes or curated evidence matrix"

    adult_delta = to_float(row["human_adult_ds_mean_delta_beta"])
    if to_int(row["human_adult_ds_probe_count"]) > 0:
        contexts[
            "human_adult_DS_methylation|hyper" if adult_delta >= 0 else "human_adult_DS_methylation|hypo"
        ] = "candidate-linked adult DS methylation probes"

    if to_int(row["mouse_dnmt3l_dmr_count"]) > 0:
        contexts["mouse_Dnmt3l_OE_methylation|hyper"] = "mouse ortholog linked to Dnmt3l-OE hyper-DMRs"
    if to_int(row["gse305146_ds_open_peak_count_100kb"]) > 0:
        contexts["human_fetal_DS_scATAC|DS_open"] = "candidate hg38 window overlaps fetal DS-open recurrent peaks"
    if to_int(row["gse305146_control_open_peak_count_100kb"]) > 0:
        contexts["human_fetal_DS_scATAC|Control_open"] = "candidate hg38 window overlaps fetal control-open recurrent peaks"
    return contexts


def motif_family_flags(motif_rows: list[dict[str, object]]) -> tuple[str, str, str, str]:
    families = {str(row["motif_family"]) for row in motif_rows}
    developmental = bool(families & DEVELOPMENTAL_FAMILIES)
    cpg = bool(families & CPG_CHROMATIN_FAMILIES)
    stat = bool(families & STAT_FAMILIES)
    family_text = ";".join(sorted(families))
    return ("yes" if developmental else "no", "yes" if cpg else "no", "yes" if stat else "no", family_text)


def mechanism_call(row: dict[str, object]) -> str:
    human_methyl = (
        to_int(row["human_dnmt3l_oe_effective_probe_or_region_count"])
        + to_int(row["human_fetal_ds_effective_probe_count"])
        + to_int(row["human_adult_ds_probe_count"])
        > 0
    )
    mouse_methyl = to_int(row["mouse_dnmt3l_dmr_count"]) > 0
    ds_atac = to_int(row["gse305146_ds_open_peak_count_100kb"]) > 0 or to_int(row["gse305146_control_open_peak_count_100kb"]) > 0
    dev_tf = row.get("developmental_tf_context_support") == "yes"
    cpg_tf = row.get("cpg_ctcf_gc_context_support") == "yes"
    directionality = str(row.get("directionality", "")).lower()
    base_class = str(row.get("mechanism_class", "")).lower()

    if "opposite_possible_repression" in directionality and human_methyl and mouse_methyl:
        return "cross_species_repression_like_methylation_candidate"
    if human_methyl and mouse_methyl and dev_tf and cpg_tf:
        return "cross_species_methylation_with_developmental_and_cpg_tf_context"
    if human_methyl and ds_atac and dev_tf:
        return "human_stage_specific_methylation_accessibility_tf_context"
    if mouse_methyl and not human_methyl:
        return "mouse_methylation_supported_human_methylation_unresolved"
    if "expression_pathway" in base_class:
        return "expression_pathway_candidate_methylation_not_excluded"
    if human_methyl or mouse_methyl:
        return "methylation_associated_without_local_atac_support"
    return "mechanism_unresolved_by_current_tables"


def methylation_inference(row: dict[str, object]) -> str:
    if to_int(row["mouse_dnmt3a_p39_cko_hypo_overlap_dmr_count"]) > 0:
        return "mouse Dnmt3l-OE DMRs overlap postnatal Dnmt3a-cKO hypomethylated CG regions; supports a DNMT3A-linked methylation component"
    if to_int(row["mouse_dnmt3l_dmr_count"]) > 0 or to_int(row["human_dnmt3l_oe_effective_probe_or_region_count"]) > 0:
        return "methylation-associated; current candidate-linked data do not prove DNMT3A dependence"
    return "methylation dependence unresolved; expression/pathway evidence should not be called methylation-independent"


def stat_inference(row: dict[str, object]) -> str:
    if row.get("stat_context_support") == "yes":
        return "STAT motif family appears in an active altered-region context; this is context-level, not locus-specific binding proof"
    return "no candidate-linked STAT motif-enrichment support detected in the active HOMER context layer"


def build_heatmap(rows: list[dict[str, object]], figdir: Path) -> None:
    figdir.mkdir(parents=True, exist_ok=True)
    columns = [
        ("Human DNMT3L\nmethyl.", lambda row: to_int(row["human_dnmt3l_oe_effective_probe_or_region_count"]) > 0),
        ("Fetal DS\nmethyl.", lambda row: to_int(row["human_fetal_ds_effective_probe_count"]) > 0),
        ("Adult DS\nmethyl.", lambda row: to_int(row["human_adult_ds_probe_count"]) > 0),
        ("Mouse Dnmt3l\nDMR", lambda row: to_int(row["mouse_dnmt3l_dmr_count"]) > 0),
        ("Dnmt3a ref.\noverlap", lambda row: to_int(row["mouse_dnmt3a_p39_cko_hypo_overlap_dmr_count"]) > 0),
        ("Fetal DS-open\nATAC", lambda row: to_int(row["gse305146_ds_open_peak_count_100kb"]) > 0),
        ("Control-open\nATAC", lambda row: to_int(row["gse305146_control_open_peak_count_100kb"]) > 0),
        ("Dev. TF\ncontext", lambda row: row.get("developmental_tf_context_support") == "yes"),
        ("CpG/CTCF\ncontext", lambda row: row.get("cpg_ctcf_gc_context_support") == "yes"),
        ("STAT\ncontext", lambda row: row.get("stat_context_support") == "yes"),
    ]
    matrix = np.array([[1 if fn(row) else 0 for _, fn in columns] for row in rows], dtype=float)
    fig_height = max(6.0, 0.34 * len(rows) + 1.8)
    fig, ax = plt.subplots(figsize=(10.2, fig_height))
    cmap = ListedColormap(["#f4f4f2", "#2b6cb0"])
    im = ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(range(len(columns)))
    ax.set_xticklabels([label for label, _ in columns], rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([str(row["gene_symbol"]) for row in rows], fontsize=9)
    ax.set_title("Candidate-locus mechanism evidence", fontsize=12, pad=12)
    ax.set_xticks(np.arange(-0.5, len(columns), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(rows), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.0)
    ax.tick_params(which="minor", bottom=False, left=False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02, ticks=[0.25, 0.75])
    cbar.ax.set_yticklabels(["absent", "present"])
    fig.tight_layout()
    fig.savefig(figdir / "candidate_locus_mechanism_evidence_heatmap.png", dpi=300)
    fig.savefig(figdir / "candidate_locus_mechanism_evidence_heatmap.pdf")
    plt.close(fig)


def write_note(path: Path, candidate_rows: list[dict[str, object]], motif_rows: list[dict[str, object]]) -> None:
    calls = pd.Series([row["candidate_locus_mechanism_call"] for row in candidate_rows]).value_counts()
    dev_count = sum(row["developmental_tf_context_support"] == "yes" for row in candidate_rows)
    cpg_count = sum(row["cpg_ctcf_gc_context_support"] == "yes" for row in candidate_rows)
    stat_count = sum(row["stat_context_support"] == "yes" for row in candidate_rows)
    primary = [row["gene_symbol"] for row in candidate_rows if str(row["tier"]) == "1_primary"]
    lines = [
        "# Candidate-locus mechanism integration",
        "",
        "This analysis integrates the 24 final candidates with local methylation, fetal DS scATAC, DNMT3A-overlap reference, and HOMER motif-family context.",
        "",
        "Important scope rule: the motif layer is context-level HOMER enrichment for altered-region sets. It is not exact motif-instance evidence at a candidate locus unless a separate motif-scanning step is run.",
        "",
        "## Main output",
        "",
        "- Candidate table: `results/mechanism_question/candidate_locus_mechanism_integration.tsv`",
        "- Motif-support long table: `results/mechanism_question/candidate_locus_context_motif_support.tsv`",
        "- Summary table: `results/mechanism_question/candidate_locus_mechanism_summary.tsv`",
        "- Heatmap: `figures/mechanism_question/candidate_locus_mechanism_evidence_heatmap.png`",
        "",
        "## High-level result",
        "",
        f"- Primary candidates retained in this integration: {', '.join(primary)}.",
        f"- Candidates with developmental TF context support: {dev_count}/{len(candidate_rows)}.",
        f"- Candidates with CpG/CTCF/GC-rich motif context support: {cpg_count}/{len(candidate_rows)}.",
        f"- Candidates with STAT context support: {stat_count}/{len(candidate_rows)}.",
        f"- Context-level motif-support rows written: {len(motif_rows)}.",
        "",
        "## Mechanism-call counts",
        "",
    ]
    for call, count in calls.items():
        lines.append(f"- {call}: {count}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The current data support a combined model rather than a clean either/or model. Many candidates have methylation support in human DNMT3L-OE or DS methylation data and mouse Dnmt3l-OE DMR data, while the HOMER layer repeatedly points to developmental TF families and CpG/CTCF/GC-rich context. Local fetal DS scATAC support is present for a subset of candidate windows, but these ATAC proxy tests were not FDR-significant at the gene-window level, so they are best treated as locus context rather than expression proof.",
            "",
            "The STAT branch remains weaker in these candidate-locus outputs: STAT1/STAT3 may still matter upstream or through expression/network state, but these HOMER region-set summaries do not provide strong candidate-local STAT binding evidence.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines))


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    figdir = Path(args.figdir)
    outdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    final_summary = read_tsv(FINAL_SUMMARY)
    final_evidence = read_tsv(FINAL_EVIDENCE)
    final_tests = read_tsv(FINAL_TESTS)
    human_dmr = add_gene_tokens(read_tsv(HUMAN_DMR))
    human_dmc = add_gene_tokens(read_tsv(HUMAN_DMC))
    ds_values = add_gene_tokens(read_tsv(DS_PROBE_VALUES))
    mouse_links = read_tsv(MOUSE_DMR_LINKS)
    mouse_priority = read_tsv(MOUSE_DNMT3A_PRIORITY)
    windows = read_tsv(ATAC_WINDOWS)
    atac_effects = read_tsv(ATAC_EFFECTS)
    ds_open = read_tsv(DS_OPEN_PEAKS)
    control_open = read_tsv(CONTROL_OPEN_PEAKS)
    motifs = read_tsv(MOTIF_GROUPS)

    evidence_by_gene = final_evidence.set_index("gene_symbol").to_dict(orient="index")
    tests_by_gene = final_tests.set_index("gene_symbol").to_dict(orient="index")
    windows_by_gene = windows.set_index("gene_symbol")
    atac_effects_by_gene = atac_effects.set_index("gene_symbol").to_dict(orient="index")

    candidate_rows: list[dict[str, object]] = []
    motif_support_rows: list[dict[str, object]] = []

    for summary_row in final_summary.to_dict(orient="records"):
        gene = str(summary_row["gene_symbol"])
        evidence = evidence_by_gene.get(gene, {})
        tests = tests_by_gene.get(gene, {})
        mouse_symbol = str(evidence.get("mouse_symbol", "") or "")

        dmr_summary = summarize_human_dnmt3l(human_dmr, gene, "dmr_id", "mean_delta_beta")
        dmc_summary = summarize_human_dnmt3l(human_dmc, gene, "probe_id", "delta_beta")
        fetal_ds = summarize_ds_probes(ds_values, gene, "fetal")
        adult_ds = summarize_ds_probes(ds_values, gene, "adult")
        mouse = summarize_mouse_gene(mouse_links, mouse_priority, mouse_symbol)

        if gene in windows_by_gene.index:
            window_row = windows_by_gene.loc[gene]
            if isinstance(window_row, pd.DataFrame):
                window_row = window_row.iloc[0]
            ds_peak = overlap_peak_summary(ds_open, window_row)
            control_peak = overlap_peak_summary(control_open, window_row)
        else:
            ds_peak = {"count": 0, "mean_delta_frequency": np.nan, "max_abs_delta_frequency": np.nan, "intervals": ""}
            control_peak = {"count": 0, "mean_delta_frequency": np.nan, "max_abs_delta_frequency": np.nan, "intervals": ""}

        atac_effect = atac_effects_by_gene.get(gene, {})

        row: dict[str, object] = {
            "gene_symbol": gene,
            "mouse_symbol": mouse_symbol,
            "tier": summary_row.get("tier", ""),
            "candidate_group": summary_row.get("candidate_group", ""),
            "mechanism_class": summary_row.get("mechanism_class", ""),
            "paper_role": summary_row.get("paper_role", ""),
            "directionality": summary_row.get("directionality", ""),
            "base_short_interpretation": summary_row.get("short_interpretation", ""),
            "human_dnmt3l_oe_dmr_count": dmr_summary["count"],
            "human_dnmt3l_oe_hyper_dmr_count": dmr_summary["hyper_count"],
            "human_dnmt3l_oe_hypo_dmr_count": dmr_summary["hypo_count"],
            "human_dnmt3l_oe_dmr_direction_summary": dmr_summary["direction_summary"],
            "human_dnmt3l_oe_dmr_mean_delta_beta": dmr_summary["mean_delta_beta"],
            "human_dnmt3l_oe_dmr_max_abs_delta_beta": dmr_summary["max_abs_delta_beta"],
            "human_dnmt3l_oe_dmr_ids": dmr_summary["ids"],
            "human_dnmt3l_oe_dmc_count": dmc_summary["count"],
            "human_dnmt3l_oe_hyper_dmc_count": dmc_summary["hyper_count"],
            "human_dnmt3l_oe_hypo_dmc_count": dmc_summary["hypo_count"],
            "human_dnmt3l_oe_dmc_mean_delta_beta": dmc_summary["mean_delta_beta"],
            "human_dnmt3l_oe_dmc_ids": dmc_summary["ids"],
            "human_dnmt3l_oe_curated_probe_count": max(
                to_int(evidence.get("human_strict_fetal_positive_probes", 0)),
                to_int(evidence.get("human_strict_fetal_positive_probe_count", 0)),
                len(split_list(evidence.get("human_linked_dnmt3l_sensitive_probes", ""))),
            ),
            "human_dnmt3l_oe_curated_mean_delta_beta": evidence.get("mean_gse157078_delta_beta", np.nan),
            "human_dnmt3l_oe_effective_probe_or_region_count": max(
                to_int(dmr_summary["count"]) + to_int(dmc_summary["count"]),
                to_int(evidence.get("human_strict_fetal_positive_probes", 0)),
                to_int(evidence.get("human_strict_fetal_positive_probe_count", 0)),
                len(split_list(evidence.get("human_linked_dnmt3l_sensitive_probes", ""))),
            ),
            "human_fetal_ds_probe_count": fetal_ds["probe_count"],
            "human_fetal_ds_curated_probe_count": max(
                to_int(evidence.get("human_triple_fetal_ds_probes", 0)),
                to_int(evidence.get("human_strict_fetal_positive_probes", 0)),
                to_int(evidence.get("human_strict_fetal_positive_probe_count", 0)),
            ),
            "human_fetal_ds_effective_probe_count": max(
                to_int(fetal_ds["probe_count"]),
                to_int(evidence.get("human_triple_fetal_ds_probes", 0)),
                to_int(evidence.get("human_strict_fetal_positive_probes", 0)),
                to_int(evidence.get("human_strict_fetal_positive_probe_count", 0)),
            ),
            "human_any_ds_curated_probe_count": to_int(evidence.get("human_triple_any_ds_probes", 0)),
            "human_ds_curated_mean_signed_delta_beta": evidence.get("mean_ds_signed_delta", np.nan),
            "human_fetal_ds_context_count": fetal_ds["context_count"],
            "human_fetal_ds_mean_delta_beta": fetal_ds["mean_ds_delta_beta"],
            "human_fetal_ds_mean_signed_delta_by_dnmt3l_direction": fetal_ds[
                "mean_signed_ds_delta_by_dnmt3l_direction"
            ],
            "human_fetal_ds_concordant_fraction": fetal_ds["concordant_fraction"],
            "human_fetal_ds_positive_context_count": fetal_ds["positive_context_count"],
            "human_fetal_ds_negative_context_count": fetal_ds["negative_context_count"],
            "human_fetal_ds_contexts": fetal_ds["contexts"],
            "human_adult_ds_probe_count": adult_ds["probe_count"],
            "human_adult_ds_context_count": adult_ds["context_count"],
            "human_adult_ds_mean_delta_beta": adult_ds["mean_ds_delta_beta"],
            "human_adult_ds_mean_signed_delta_by_dnmt3l_direction": adult_ds[
                "mean_signed_ds_delta_by_dnmt3l_direction"
            ],
            "human_adult_ds_concordant_fraction": adult_ds["concordant_fraction"],
            "human_adult_ds_positive_context_count": adult_ds["positive_context_count"],
            "human_adult_ds_negative_context_count": adult_ds["negative_context_count"],
            "human_adult_ds_contexts": adult_ds["contexts"],
            "mouse_dnmt3l_dmr_count": mouse.get("mouse_dnmt3l_dmr_count", 0),
            "mouse_dnmt3l_promoter_2kb_dmr_count": mouse.get("mouse_dnmt3l_promoter_2kb_dmr_count", 0),
            "mouse_dnmt3l_gene_body_dmr_count": mouse.get("mouse_dnmt3l_gene_body_dmr_count", 0),
            "mouse_dnmt3l_nearest_100kb_dmr_count": mouse.get("mouse_dnmt3l_nearest_100kb_dmr_count", 0),
            "mouse_dnmt3l_mean_delta_beta": mouse.get("mouse_dnmt3l_mean_delta_beta", np.nan),
            "mouse_dnmt3a_reference_overlap_dmr_count": mouse.get("mouse_dnmt3a_reference_overlap_dmr_count", 0),
            "mouse_dnmt3a_p39_cko_hypo_overlap_dmr_count": mouse.get(
                "mouse_dnmt3a_p39_cko_hypo_overlap_dmr_count", 0
            ),
            "mouse_h3k27me3_p39_cko_up_overlap_dmr_count": mouse.get(
                "mouse_h3k27me3_p39_cko_up_overlap_dmr_count", 0
            ),
            "mouse_dnmt3l_dmr_ids": mouse.get("mouse_dnmt3l_dmr_ids", ""),
            "gse305146_ds_open_peak_count_100kb": ds_peak["count"],
            "gse305146_ds_open_mean_delta_frequency": ds_peak["mean_delta_frequency"],
            "gse305146_ds_open_max_abs_delta_frequency": ds_peak["max_abs_delta_frequency"],
            "gse305146_ds_open_peak_intervals": ds_peak["intervals"],
            "gse305146_control_open_peak_count_100kb": control_peak["count"],
            "gse305146_control_open_mean_delta_frequency": control_peak["mean_delta_frequency"],
            "gse305146_control_open_max_abs_delta_frequency": control_peak["max_abs_delta_frequency"],
            "gse305146_control_open_peak_intervals": control_peak["intervals"],
            "gse305146_window_delta_log2_peak_cpm_ds_minus_control": atac_effect.get(
                "delta_log2_peak_cpm_ds_minus_control", np.nan
            ),
            "gse305146_window_bh_fdr": atac_effect.get("bh_fdr", np.nan),
            "human_any_regulatory_probe_count": evidence.get("human_any_regulatory_probe_count", 0),
            "human_strict_regulatory_probe_count": evidence.get("human_strict_regulatory_probe_count", 0),
            "human_fetal_h3k27me3_probe_count": evidence.get("human_fetal_h3k27me3_probe_count", 0),
            "mouse_any_regulatory_dmr_count": evidence.get("mouse_any_regulatory_dmr_count", 0),
            "mouse_embryonic_forebrain_h3k27me3_dmr_count": evidence.get(
                "mouse_embryonic_forebrain_h3k27me3_dmr_count", 0
            ),
            "regulatory_context_class": evidence.get("regulatory_context_class", ""),
            "proposed_in_silico_mechanism": tests.get("proposed_in_silico_mechanism", ""),
            "core_testable_hypothesis": tests.get("core_testable_hypothesis", ""),
            "motif_support_scope": "context_level_region_set_enrichment_not_candidate_motif_instance",
        }

        context_map = active_contexts(row)
        candidate_motif_rows: list[dict[str, object]] = []
        for analysis_group, local_anchor in context_map.items():
            for motif_row in supported_motif_rows(motifs, analysis_group).to_dict(orient="records"):
                support_row = {
                    "gene_symbol": gene,
                    "tier": row["tier"],
                    "analysis_group": analysis_group,
                    "candidate_local_anchor": local_anchor,
                    "motif_family": motif_row.get("motif_family", ""),
                    "support_call": motif_row.get("support_call", ""),
                    "best_q_value": motif_row.get("best_q_value", np.nan),
                    "best_rank": motif_row.get("best_rank", np.nan),
                    "best_motif_name": motif_row.get("best_motif_name", ""),
                    "support_scope": "context_level_region_set_enrichment_not_candidate_motif_instance",
                }
                motif_support_rows.append(support_row)
                candidate_motif_rows.append(support_row)

        dev_flag, cpg_flag, stat_flag, families = motif_family_flags(candidate_motif_rows)
        row["developmental_tf_context_support"] = dev_flag
        row["cpg_ctcf_gc_context_support"] = cpg_flag
        row["stat_context_support"] = stat_flag
        row["supported_context_motif_families"] = families
        row["candidate_locus_mechanism_call"] = mechanism_call(row)
        row["methylation_dependence_inference"] = methylation_inference(row)
        row["stat_tf_independence_inference"] = stat_inference(row)
        candidate_rows.append(row)

    candidate_fields = [
        "gene_symbol",
        "mouse_symbol",
        "tier",
        "candidate_group",
        "mechanism_class",
        "paper_role",
        "candidate_locus_mechanism_call",
        "developmental_tf_context_support",
        "cpg_ctcf_gc_context_support",
        "stat_context_support",
        "human_dnmt3l_oe_dmr_count",
        "human_dnmt3l_oe_hyper_dmr_count",
        "human_dnmt3l_oe_hypo_dmr_count",
        "human_dnmt3l_oe_dmr_mean_delta_beta",
        "human_dnmt3l_oe_dmc_count",
        "human_dnmt3l_oe_curated_probe_count",
        "human_dnmt3l_oe_effective_probe_or_region_count",
        "human_dnmt3l_oe_curated_mean_delta_beta",
        "human_fetal_ds_probe_count",
        "human_fetal_ds_curated_probe_count",
        "human_fetal_ds_effective_probe_count",
        "human_any_ds_curated_probe_count",
        "human_fetal_ds_context_count",
        "human_fetal_ds_mean_delta_beta",
        "human_ds_curated_mean_signed_delta_beta",
        "human_fetal_ds_concordant_fraction",
        "human_adult_ds_probe_count",
        "human_adult_ds_context_count",
        "human_adult_ds_mean_delta_beta",
        "human_adult_ds_concordant_fraction",
        "mouse_dnmt3l_dmr_count",
        "mouse_dnmt3l_promoter_2kb_dmr_count",
        "mouse_dnmt3l_gene_body_dmr_count",
        "mouse_dnmt3l_nearest_100kb_dmr_count",
        "mouse_dnmt3l_mean_delta_beta",
        "mouse_dnmt3a_p39_cko_hypo_overlap_dmr_count",
        "mouse_h3k27me3_p39_cko_up_overlap_dmr_count",
        "gse305146_ds_open_peak_count_100kb",
        "gse305146_control_open_peak_count_100kb",
        "gse305146_window_delta_log2_peak_cpm_ds_minus_control",
        "gse305146_window_bh_fdr",
        "human_any_regulatory_probe_count",
        "human_strict_regulatory_probe_count",
        "human_fetal_h3k27me3_probe_count",
        "mouse_any_regulatory_dmr_count",
        "mouse_embryonic_forebrain_h3k27me3_dmr_count",
        "regulatory_context_class",
        "supported_context_motif_families",
        "motif_support_scope",
        "methylation_dependence_inference",
        "stat_tf_independence_inference",
        "proposed_in_silico_mechanism",
        "core_testable_hypothesis",
        "directionality",
        "base_short_interpretation",
        "human_dnmt3l_oe_dmr_direction_summary",
        "human_dnmt3l_oe_dmr_ids",
        "human_dnmt3l_oe_dmc_ids",
        "human_fetal_ds_contexts",
        "human_adult_ds_contexts",
        "mouse_dnmt3l_dmr_ids",
        "gse305146_ds_open_peak_intervals",
        "gse305146_control_open_peak_intervals",
    ]
    write_tsv(outdir / "candidate_locus_mechanism_integration.tsv", candidate_rows, candidate_fields)

    motif_fields = [
        "gene_symbol",
        "tier",
        "analysis_group",
        "candidate_local_anchor",
        "motif_family",
        "support_call",
        "best_q_value",
        "best_rank",
        "best_motif_name",
        "support_scope",
    ]
    write_tsv(outdir / "candidate_locus_context_motif_support.tsv", motif_support_rows, motif_fields)

    summary_df = pd.DataFrame(candidate_rows)
    summary_rows = []
    for call, group in summary_df.groupby("candidate_locus_mechanism_call", sort=False):
        summary_rows.append(
            {
                "candidate_locus_mechanism_call": call,
                "n_candidates": int(len(group)),
                "genes": ";".join(group["gene_symbol"].astype(str).tolist()),
                "n_human_dnmt3l_oe_dmr_supported": int((group["human_dnmt3l_oe_dmr_count"].astype(int) > 0).sum()),
                "n_human_dnmt3l_oe_effective_supported": int(
                    (group["human_dnmt3l_oe_effective_probe_or_region_count"].astype(int) > 0).sum()
                ),
                "n_fetal_ds_methylation_supported": int((group["human_fetal_ds_effective_probe_count"].astype(int) > 0).sum()),
                "n_adult_ds_methylation_supported": int((group["human_adult_ds_probe_count"].astype(int) > 0).sum()),
                "n_mouse_dnmt3l_dmr_supported": int((group["mouse_dnmt3l_dmr_count"].astype(int) > 0).sum()),
                "n_mouse_dnmt3a_overlap_supported": int(
                    (group["mouse_dnmt3a_p39_cko_hypo_overlap_dmr_count"].astype(int) > 0).sum()
                ),
                "n_ds_open_atac_window_supported": int((group["gse305146_ds_open_peak_count_100kb"].astype(int) > 0).sum()),
                "n_developmental_tf_context_supported": int((group["developmental_tf_context_support"] == "yes").sum()),
                "n_cpg_ctcf_gc_context_supported": int((group["cpg_ctcf_gc_context_support"] == "yes").sum()),
                "n_stat_context_supported": int((group["stat_context_support"] == "yes").sum()),
            }
        )
    write_tsv(
        outdir / "candidate_locus_mechanism_summary.tsv",
        summary_rows,
        [
            "candidate_locus_mechanism_call",
            "n_candidates",
            "genes",
            "n_human_dnmt3l_oe_dmr_supported",
            "n_human_dnmt3l_oe_effective_supported",
            "n_fetal_ds_methylation_supported",
            "n_adult_ds_methylation_supported",
            "n_mouse_dnmt3l_dmr_supported",
            "n_mouse_dnmt3a_overlap_supported",
            "n_ds_open_atac_window_supported",
            "n_developmental_tf_context_supported",
            "n_cpg_ctcf_gc_context_supported",
            "n_stat_context_supported",
        ],
    )

    build_heatmap(candidate_rows, figdir)
    write_note(outdir / "candidate_locus_mechanism_integration.md", candidate_rows, motif_support_rows)

    print(f"Wrote {len(candidate_rows)} candidate rows")
    print(f"Wrote {len(motif_support_rows)} candidate x context motif-support rows")
    print(outdir / "candidate_locus_mechanism_integration.tsv")


if __name__ == "__main__":
    main()
