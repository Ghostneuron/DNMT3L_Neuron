#!/usr/bin/env python3
"""Scan candidate-linked altered regions for HOMER motif instances.

This is a tighter follow-up to the context-level HOMER enrichment analysis. It
uses candidate-linked DMRs/probes/ATAC peaks as foreground regions and scans
them with motif files selected from existing HOMER enrichment outputs.

Important interpretation rule: an instance hit is sequence-level motif
compatibility. It is not TF binding, and unsupported STAT motif scans are kept
as negative-control/exploratory rows rather than enrichment evidence.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import subprocess
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap


FINAL_CANDIDATES = Path("results/manuscript_summary/final_candidate_manuscript_summary.tsv")
FINAL_EVIDENCE = Path("results/final_candidate_evidence/final_candidate_evidence_matrix_with_regulatory_context.tsv")
HUMAN_DMR = Path("results/human_candidate_dmrs/GSE157078_candidate_dmrs.tsv")
HUMAN_DMC = Path("results/human_candidate_dmrs/GSE157078_candidate_dmcs.tsv")
DS_PROBES = Path("results/ds_methylation_dnmt3l_concordance/ds_dnmt3l_candidate_probe_values.tsv")
MOUSE_DMR_LINKS = Path("results/integrated_candidate_genes/mouse_triple_hit_dmr_gene_links.tsv")
ATAC_WINDOWS = Path("results/mechanism_question/gse305146_accessibility_proxy_gene_windows.tsv")
DS_OPEN_PEAKS = Path("results/gse305146_scATAC_homer_motif_enrichment/tables/GSE305146_DS_open_recurrent_peak_calls.tsv")
CONTROL_OPEN_PEAKS = Path("results/gse305146_scATAC_homer_motif_enrichment/tables/GSE305146_Control_open_recurrent_peak_calls.tsv")
MOTIF_GROUPS = Path("results/mechanism_question/homer_cross_context_motif_family_group_summary.tsv")

HOMER_BIN = Path("/Users/jili/tools/homer/bin")
ANNOTATE_PEAKS = HOMER_BIN / "annotatePeaks.pl"
STAT_CONTROL_MOTIFS = [
    {
        "analysis_group": "known_library_STAT_control",
        "motif_path": HOMER_BIN.parent / "motifs" / "stat1.motif",
        "motif_name": "STAT1(Stat)/HelaS3-STAT1-ChIP-Seq(GSE12782)/Homer",
        "motif_families": "STAT",
        "best_q_value": np.nan,
        "support_call": "known_library_control_not_enriched",
        "context_enriched": "no",
        "stat_negative_control": "yes",
        "motif_category": "stat",
    },
    {
        "analysis_group": "known_library_STAT_control",
        "motif_path": HOMER_BIN.parent / "motifs" / "stat3.motif",
        "motif_name": "Stat3(Stat)/mES-Stat3-ChIP-Seq(GSE11431)/Homer",
        "motif_families": "STAT",
        "best_q_value": np.nan,
        "support_call": "known_library_control_not_enriched",
        "context_enriched": "no",
        "stat_negative_control": "yes",
        "motif_category": "stat",
    },
]

HOMER_BASES = {
    "human_DNMT3L_OE_methylation": Path("results/homer_motif_enrichment"),
    "human_fetal_DS_methylation": Path("results/ds_homer_motif_enrichment"),
    "human_adult_DS_methylation": Path("results/adult_human_ds_homer_motif_enrichment_sample25k_s8"),
    "mouse_Dnmt3l_OE_methylation": Path("results/mouse_homer_motif_enrichment"),
    "mouse_DS_model_methylation": Path("results/mouse_ds_homer_motif_enrichment"),
    "human_fetal_DS_scATAC": Path("results/gse305146_scATAC_homer_motif_enrichment"),
}

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
CPG_CTCF_FAMILIES = {"KLF_SP_GC", "CTCF_BORIS", "ZINC_FINGER_OTHER"}
SIGNIFICANT_SUPPORT = {"strong_single_context", "significant_single_context", "recurrent_significant"}

REGION_GROUPS = [
    {
        "name": "human_dnmt3l_oe_candidate_dmr_hyper_hg19",
        "genome": "hg19",
        "region_kind": "human_dnmt3l_dmr",
        "direction": "hyper",
        "analysis_groups": ["human_DNMT3L_OE_methylation|hyper"],
    },
    {
        "name": "human_dnmt3l_oe_candidate_dmr_hypo_hg19",
        "genome": "hg19",
        "region_kind": "human_dnmt3l_dmr",
        "direction": "hypo",
        "analysis_groups": ["human_DNMT3L_OE_methylation|hypo"],
    },
    {
        "name": "human_dnmt3l_oe_candidate_dmc_hyper_hg19",
        "genome": "hg19",
        "region_kind": "human_dnmt3l_dmc",
        "direction": "hyper",
        "analysis_groups": ["human_DNMT3L_OE_methylation|hyper"],
    },
    {
        "name": "human_dnmt3l_oe_candidate_dmc_hypo_hg19",
        "genome": "hg19",
        "region_kind": "human_dnmt3l_dmc",
        "direction": "hypo",
        "analysis_groups": ["human_DNMT3L_OE_methylation|hypo"],
    },
    {
        "name": "human_fetal_ds_candidate_probe_hyper_hg19",
        "genome": "hg19",
        "region_kind": "human_ds_probe_fetal",
        "direction": "hyper",
        "analysis_groups": ["human_fetal_DS_methylation|hyper"],
    },
    {
        "name": "human_fetal_ds_candidate_probe_hypo_hg19",
        "genome": "hg19",
        "region_kind": "human_ds_probe_fetal",
        "direction": "hypo",
        "analysis_groups": ["human_fetal_DS_methylation|hypo"],
    },
    {
        "name": "human_adult_ds_candidate_probe_hyper_hg19",
        "genome": "hg19",
        "region_kind": "human_ds_probe_adult",
        "direction": "hyper",
        "analysis_groups": ["human_adult_DS_methylation|hyper"],
    },
    {
        "name": "human_adult_ds_candidate_probe_hypo_hg19",
        "genome": "hg19",
        "region_kind": "human_ds_probe_adult",
        "direction": "hypo",
        "analysis_groups": ["human_adult_DS_methylation|hypo"],
    },
    {
        "name": "mouse_dnmt3l_oe_candidate_dmr_hyper_mm10",
        "genome": "mm10",
        "region_kind": "mouse_dnmt3l_dmr",
        "direction": "hyper",
        "analysis_groups": ["mouse_Dnmt3l_OE_methylation|hyper"],
    },
    {
        "name": "gse305146_ds_open_candidate_peaks_hg38",
        "genome": "hg38",
        "region_kind": "human_fetal_ds_atac_ds_open",
        "direction": "DS_open",
        "analysis_groups": ["human_fetal_DS_scATAC|DS_open"],
    },
    {
        "name": "gse305146_control_open_candidate_peaks_hg38",
        "genome": "hg38",
        "region_kind": "human_fetal_ds_atac_control_open",
        "direction": "Control_open",
        "analysis_groups": ["human_fetal_DS_scATAC|Control_open"],
    },
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="results/mechanism_question/candidate_locus_homer_instance_scan")
    parser.add_argument("--figdir", default="figures/mechanism_question")
    parser.add_argument("--q-threshold", type=float, default=0.05)
    parser.add_argument("--probe-window-bp", type=int, default=200)
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
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fieldnames})


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
        text = str(row[column]).replace("///", ";").replace("|", ";").replace(",", ";")
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


def clean_id(text: object) -> str:
    text = str(text)
    text = re.sub(r"[^A-Za-z0-9_.:+-]+", "_", text)
    return text.strip("_")[:180]


def centered_window(start: object, end: object, width: int) -> tuple[int, int]:
    start_i = int(float(start))
    end_i = int(float(end))
    center = int(round((start_i + end_i) / 2))
    half = max(1, width // 2)
    return max(0, center - half), center + half


def append_region(
    rows: list[dict[str, object]],
    *,
    source_group: str,
    region_kind: str,
    genome: str,
    gene: str,
    region_id: str,
    chrom: str,
    start: int,
    end: int,
    direction: str,
    score: float = 0.0,
    extra: str = "",
) -> None:
    if end <= start:
        end = start + 1
    peak_id = clean_id(f"{source_group}|{gene}|{region_id}")
    rows.append(
        {
            "peak_id": peak_id,
            "source_group": source_group,
            "region_kind": region_kind,
            "genome": genome,
            "gene_symbol": gene,
            "region_id": region_id,
            "chrom": str(chrom),
            "start": int(start),
            "end": int(end),
            "direction": direction,
            "score": score,
            "strand": "+",
            "extra": extra,
        }
    )


def build_candidate_regions(
    final_genes: list[str],
    evidence: pd.DataFrame,
    probe_window_bp: int,
) -> dict[str, list[dict[str, object]]]:
    groups: dict[str, list[dict[str, object]]] = {group["name"]: [] for group in REGION_GROUPS}
    evidence_by_gene = evidence.set_index("gene_symbol").to_dict(orient="index")

    human_dmr = add_gene_tokens(read_tsv(HUMAN_DMR))
    human_dmc = add_gene_tokens(read_tsv(HUMAN_DMC))
    ds_probes = add_gene_tokens(read_tsv(DS_PROBES))
    mouse_links = read_tsv(MOUSE_DMR_LINKS)
    windows = read_tsv(ATAC_WINDOWS)
    ds_open = read_tsv(DS_OPEN_PEAKS)
    control_open = read_tsv(CONTROL_OPEN_PEAKS)

    for gene in final_genes:
        gene_key = gene.upper()
        ev = evidence_by_gene.get(gene, {})
        curated_dnmt3l_probes = set(split_list(ev.get("human_linked_dnmt3l_sensitive_probes", "")))
        curated_ds_probes = set(split_list(ev.get("human_strict_fetal_positive_probe_ids", "")))

        matched_dmr = human_dmr[human_dmr["_gene_tokens"].map(lambda tokens: gene_key in tokens)].copy()
        for row in matched_dmr.to_dict(orient="records"):
            direction = str(row.get("direction", ""))
            group_name = f"human_dnmt3l_oe_candidate_dmr_{direction}_hg19"
            if group_name not in groups:
                continue
            append_region(
                groups[group_name],
                source_group=group_name,
                region_kind="human_dnmt3l_dmr",
                genome="hg19",
                gene=gene,
                region_id=str(row["dmr_id"]),
                chrom=str(row["chrom"]),
                start=int(row["start"]),
                end=int(row["end"]),
                direction=direction,
                score=float(row.get("mean_delta_beta", 0.0)),
            )

        matched_dmc = human_dmc[
            human_dmc["_gene_tokens"].map(lambda tokens: gene_key in tokens)
            | human_dmc["probe_id"].astype(str).isin(curated_dnmt3l_probes)
        ].copy()
        for row in matched_dmc.to_dict(orient="records"):
            direction = str(row.get("direction", ""))
            group_name = f"human_dnmt3l_oe_candidate_dmc_{direction}_hg19"
            if group_name not in groups:
                continue
            start, end = centered_window(row["start"], row["end"], probe_window_bp)
            append_region(
                groups[group_name],
                source_group=group_name,
                region_kind="human_dnmt3l_dmc",
                genome="hg19",
                gene=gene,
                region_id=str(row["probe_id"]),
                chrom=str(row["chrom"]),
                start=start,
                end=end,
                direction=direction,
                score=float(row.get("delta_beta", 0.0)),
            )

        matched_ds = ds_probes[
            ds_probes["_gene_tokens"].map(lambda tokens: gene_key in tokens)
            | ds_probes["probe_id"].astype(str).isin(curated_ds_probes)
        ].copy()
        if not matched_ds.empty:
            matched_ds["stage"] = np.where(
                matched_ds["context"].astype(str).str.contains("adult", case=False, na=False), "adult", "fetal"
            )
            grouped = (
                matched_ds.groupby(["stage", "probe_id", "chrom", "start", "end"], dropna=False)
                .agg(mean_ds_delta_beta=("ds_delta_beta", "mean"), contexts=("context", lambda values: ";".join(sorted(set(map(str, values))))))
                .reset_index()
            )
            for row in grouped.to_dict(orient="records"):
                stage = str(row["stage"])
                direction = "hyper" if float(row["mean_ds_delta_beta"]) >= 0 else "hypo"
                group_name = f"human_{stage}_ds_candidate_probe_{direction}_hg19"
                if group_name not in groups:
                    continue
                start, end = centered_window(row["start"], row["end"], probe_window_bp)
                append_region(
                    groups[group_name],
                    source_group=group_name,
                    region_kind=f"human_ds_probe_{stage}",
                    genome="hg19",
                    gene=gene,
                    region_id=str(row["probe_id"]),
                    chrom=str(row["chrom"]),
                    start=start,
                    end=end,
                    direction=direction,
                    score=float(row["mean_ds_delta_beta"]),
                    extra=str(row["contexts"]),
                )

    mouse_symbol_by_gene = evidence.set_index("gene_symbol")["mouse_symbol"].dropna().astype(str).to_dict()
    for gene, mouse_symbol in mouse_symbol_by_gene.items():
        matched_mouse = mouse_links[mouse_links["gene_symbol"].astype(str).str.upper() == mouse_symbol.upper()].copy()
        for row in matched_mouse.to_dict(orient="records"):
            direction = str(row.get("dmr_direction", "hyper"))
            if direction != "hyper":
                continue
            append_region(
                groups["mouse_dnmt3l_oe_candidate_dmr_hyper_mm10"],
                source_group="mouse_dnmt3l_oe_candidate_dmr_hyper_mm10",
                region_kind="mouse_dnmt3l_dmr",
                genome="mm10",
                gene=gene,
                region_id=str(row["dmr_id"]),
                chrom=str(row["chrom"]),
                start=int(row["dmr_start"]),
                end=int(row["dmr_end"]),
                direction=direction,
                score=float(row.get("dmr_mean_delta_beta", 0.0)),
                extra=str(row.get("relationship", "")),
            )

    windows = windows[windows["gene_symbol"].isin(final_genes)].copy()
    for peak_label, peak_df, group_name in [
        ("DS_open", ds_open, "gse305146_ds_open_candidate_peaks_hg38"),
        ("Control_open", control_open, "gse305146_control_open_candidate_peaks_hg38"),
    ]:
        for window in windows.to_dict(orient="records"):
            gene = str(window["gene_symbol"])
            subset = peak_df[
                (peak_df["chrom"].astype(str) == str(window["chrom"]))
                & (peak_df["end"].astype(int) >= int(window["window_start"]))
                & (peak_df["start"].astype(int) <= int(window["window_end"]))
            ].copy()
            for row in subset.to_dict(orient="records"):
                region_id = f"{peak_label}:{row['chrom']}:{int(row['start'])}-{int(row['end'])}"
                append_region(
                    groups[group_name],
                    source_group=group_name,
                    region_kind=f"human_fetal_ds_atac_{peak_label}",
                    genome="hg38",
                    gene=gene,
                    region_id=region_id,
                    chrom=str(row["chrom"]),
                    start=int(row["start"]),
                    end=int(row["end"]),
                    direction=peak_label,
                    score=float(row.get("delta_frequency_ds_minus_control", 0.0)),
                )

    for group_name, rows in groups.items():
        deduped = {row["peak_id"]: row for row in rows}
        groups[group_name] = list(deduped.values())
    return groups


def motif_category(families: set[str]) -> str:
    if "STAT" in families:
        return "stat"
    if families & CPG_CTCF_FAMILIES:
        return "cpg_ctcf_gc"
    if families & DEVELOPMENTAL_FAMILIES:
        return "developmental_tf"
    return "other_supported"


def motif_path_for(row: pd.Series) -> Path | None:
    base = HOMER_BASES.get(str(row["broad_context"]))
    if base is None:
        return None
    rank = int(float(row["best_rank"]))
    path = base / "homer_out" / str(row["best_set_id"]) / "knownResults" / f"known{rank}.motif"
    return path if path.exists() else None


def select_motifs(q_threshold: float) -> dict[str, list[dict[str, object]]]:
    motif_groups = read_tsv(MOTIF_GROUPS)
    selected: dict[str, list[dict[str, object]]] = {}
    for analysis_group, group_df in motif_groups.groupby("analysis_group"):
        rows: dict[Path, dict[str, object]] = {}
        for _, row in group_df.iterrows():
            family = str(row["motif_family"])
            support = str(row["support_call"])
            q_value = float(row["best_q_value"]) if not pd.isna(row["best_q_value"]) else float("nan")
            include = family != "OTHER" and support in SIGNIFICANT_SUPPORT and not math.isnan(q_value) and q_value <= q_threshold
            include_stat_control = family == "STAT" and not math.isnan(q_value)
            if not include and not include_stat_control:
                continue
            motif_path = motif_path_for(row)
            if motif_path is None:
                continue
            record = rows.setdefault(
                motif_path,
                {
                    "analysis_group": analysis_group,
                    "motif_path": motif_path,
                    "motif_name": row["best_motif_name"],
                    "motif_families": set(),
                    "best_q_value": q_value,
                    "support_call": support,
                    "context_enriched": "yes" if include else "no",
                    "stat_negative_control": "yes" if include_stat_control and not include else "no",
                },
            )
            record["motif_families"].add(family)
            record["best_q_value"] = min(float(record["best_q_value"]), q_value)
            if include:
                record["context_enriched"] = "yes"
                record["support_call"] = support
        selected[str(analysis_group)] = [
            {
                **record,
                "motif_families": ";".join(sorted(record["motif_families"])),
                "motif_category": motif_category(set(record["motif_families"])),
            }
            for record in rows.values()
        ]
    return selected


def write_bed(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for row in rows:
            handle.write(
                "\t".join(
                    [
                        str(row["chrom"]),
                        str(int(row["start"])),
                        str(int(row["end"])),
                        str(row["peak_id"]),
                        str(int(round(abs(float(row.get("score", 0.0))) * 1000))),
                        str(row.get("strand", "+")),
                    ]
                )
                + "\n"
            )


def run_homer_scan(
    bed_path: Path,
    genome: str,
    motif_records: list[dict[str, object]],
    output_path: Path,
    log_path: Path,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["PATH"] = f"{HOMER_BIN}:{env.get('PATH', '')}"
    cmd = [
        str(ANNOTATE_PEAKS),
        str(bed_path),
        genome,
        "-size",
        "given",
        "-m",
        *[str(record["motif_path"]) for record in motif_records],
        "-nmotifs",
    ]
    with output_path.open("w") as stdout, log_path.open("w") as stderr:
        subprocess.run(cmd, check=True, cwd=Path.cwd(), env=env, stdout=stdout, stderr=stderr)


def parse_homer_output(
    output_path: Path,
    region_by_peak_id: dict[str, dict[str, object]],
    motif_records: list[dict[str, object]],
) -> list[dict[str, object]]:
    if not output_path.exists() or output_path.stat().st_size == 0:
        return []
    motif_by_name = {str(record["motif_name"]): record for record in motif_records}
    table = pd.read_csv(output_path, sep="\t")
    motif_columns: list[tuple[str, str]] = []
    for column in table.columns:
        if "/Homer" not in column:
            continue
        motif_name = column.split(" Distance From Peak", 1)[0]
        motif_columns.append((column, motif_name))

    hit_rows: list[dict[str, object]] = []
    for row in table.to_dict(orient="records"):
        peak_id = str(row[table.columns[0]])
        region = region_by_peak_id.get(peak_id)
        if region is None:
            continue
        for column, motif_name in motif_columns:
            try:
                count = int(float(row[column]))
            except (TypeError, ValueError):
                count = 0
            if count <= 0:
                continue
            motif = motif_by_name.get(motif_name)
            if motif is None:
                continue
            hit_rows.append(
                {
                    **{key: region.get(key, "") for key in [
                        "peak_id",
                        "source_group",
                        "region_kind",
                        "genome",
                        "gene_symbol",
                        "region_id",
                        "chrom",
                        "start",
                        "end",
                        "direction",
                    ]},
                    "motif_name": motif_name,
                    "motif_families": motif["motif_families"],
                    "motif_category": motif["motif_category"],
                    "motif_count": count,
                    "context_enriched": motif["context_enriched"],
                    "stat_negative_control": motif["stat_negative_control"],
                    "support_call": motif["support_call"],
                    "best_q_value": motif["best_q_value"],
                    "analysis_group": motif["analysis_group"],
                }
            )
    return hit_rows


def summarize_hits(
    final_genes: list[str],
    region_rows: list[dict[str, object]],
    hit_rows: list[dict[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    region_df = pd.DataFrame(region_rows)
    hit_df = pd.DataFrame(hit_rows)
    summary_rows: list[dict[str, object]] = []
    context_rows: list[dict[str, object]] = []

    for gene in final_genes:
        gene_regions = region_df[region_df["gene_symbol"] == gene] if not region_df.empty else pd.DataFrame()
        gene_hits = hit_df[hit_df["gene_symbol"] == gene] if not hit_df.empty else pd.DataFrame()
        row: dict[str, object] = {
            "gene_symbol": gene,
            "n_regions_scanned": int(len(gene_regions)),
            "n_human_dnmt3l_regions": int(gene_regions["region_kind"].astype(str).str.startswith("human_dnmt3l").sum()) if not gene_regions.empty else 0,
            "n_human_ds_probe_regions": int(gene_regions["region_kind"].astype(str).str.startswith("human_ds_probe").sum()) if not gene_regions.empty else 0,
            "n_mouse_dnmt3l_regions": int((gene_regions["region_kind"].astype(str) == "mouse_dnmt3l_dmr").sum()) if not gene_regions.empty else 0,
            "n_fetal_ds_atac_peak_regions": int(gene_regions["region_kind"].astype(str).str.startswith("human_fetal_ds_atac").sum()) if not gene_regions.empty else 0,
        }
        for category in ["developmental_tf", "cpg_ctcf_gc", "stat", "other_supported"]:
            subset = gene_hits[gene_hits["motif_category"] == category] if not gene_hits.empty else pd.DataFrame()
            row[f"n_regions_with_{category}_motif"] = int(subset["peak_id"].nunique()) if not subset.empty else 0
            row[f"{category}_source_groups"] = ";".join(sorted(subset["source_group"].dropna().astype(str).unique())) if not subset.empty else ""
            row[f"{category}_motif_families"] = ";".join(
                sorted({family for families in subset.get("motif_families", pd.Series(dtype=str)).dropna().astype(str) for family in families.split(";") if family})
            ) if not subset.empty else ""
        stat_hits = row["n_regions_with_stat_motif"]
        if stat_hits:
            row["stat_instance_interpretation"] = "STAT motif instances present, but selected STAT motifs were not enriched in their source HOMER contexts"
        else:
            row["stat_instance_interpretation"] = "no STAT motif instances detected in scanned candidate-linked regions"
        if row["n_regions_with_developmental_tf_motif"] and row["n_regions_with_cpg_ctcf_gc_motif"]:
            row["candidate_motif_instance_interpretation"] = "candidate-linked regions carry both developmental TF and CpG/CTCF/GC motif instances"
        elif row["n_regions_with_developmental_tf_motif"]:
            row["candidate_motif_instance_interpretation"] = "candidate-linked regions carry developmental TF motif instances"
        elif row["n_regions_with_cpg_ctcf_gc_motif"]:
            row["candidate_motif_instance_interpretation"] = "candidate-linked regions carry CpG/CTCF/GC motif instances"
        else:
            row["candidate_motif_instance_interpretation"] = "no selected enriched motif instances detected in scanned candidate-linked regions"
        summary_rows.append(row)

    if not region_df.empty:
        for (source_group, region_kind), group in region_df.groupby(["source_group", "region_kind"]):
            group_hits = hit_df[hit_df["source_group"] == source_group] if not hit_df.empty else pd.DataFrame()
            for category in ["developmental_tf", "cpg_ctcf_gc", "stat", "other_supported"]:
                subset = group_hits[group_hits["motif_category"] == category] if not group_hits.empty else pd.DataFrame()
                n_regions_with_motif = int(subset["peak_id"].nunique()) if not subset.empty and "peak_id" in subset else 0
                n_genes_with_motif = int(subset["gene_symbol"].nunique()) if not subset.empty and "gene_symbol" in subset else 0
                genes_with_motif = (
                    ";".join(sorted(subset["gene_symbol"].dropna().astype(str).unique()))
                    if not subset.empty and "gene_symbol" in subset
                    else ""
                )
                motif_families = (
                    ";".join(
                        sorted(
                            {
                                family
                                for families in subset.get("motif_families", pd.Series(dtype=str)).dropna().astype(str)
                                for family in families.split(";")
                                if family
                            }
                        )
                    )
                    if not subset.empty
                    else ""
                )
                context_rows.append(
                    {
                        "source_group": source_group,
                        "region_kind": region_kind,
                        "motif_category": category,
                        "n_regions_scanned": int(len(group)),
                        "n_regions_with_motif": n_regions_with_motif,
                        "fraction_regions_with_motif": (float(n_regions_with_motif) / len(group)) if len(group) else np.nan,
                        "n_candidate_genes_with_motif": n_genes_with_motif,
                        "genes_with_motif": genes_with_motif,
                        "motif_families": motif_families,
                    }
                )
    return summary_rows, context_rows


def make_heatmap(summary_rows: list[dict[str, object]], figdir: Path) -> None:
    figdir.mkdir(parents=True, exist_ok=True)
    columns = [
        ("Dev.\ninstances", "n_regions_with_developmental_tf_motif"),
        ("CpG/CTCF\ninstances", "n_regions_with_cpg_ctcf_gc_motif"),
        ("STAT ctrl\ninstances", "n_regions_with_stat_motif"),
        ("Human DNMT3L\nregions", "n_human_dnmt3l_regions"),
        ("Human DS\nprobes", "n_human_ds_probe_regions"),
        ("Mouse Dnmt3l\nDMRs", "n_mouse_dnmt3l_regions"),
        ("Fetal DS\nATAC peaks", "n_fetal_ds_atac_peak_regions"),
    ]
    matrix = np.array([[1 if int(row[column]) > 0 else 0 for _, column in columns] for row in summary_rows], dtype=float)
    fig_height = max(6.0, 0.34 * len(summary_rows) + 1.8)
    fig, ax = plt.subplots(figsize=(8.8, fig_height))
    cmap = ListedColormap(["#f4f4f2", "#3a7d44"])
    im = ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(range(len(columns)))
    ax.set_xticklabels([label for label, _ in columns], rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(len(summary_rows)))
    ax.set_yticklabels([str(row["gene_symbol"]) for row in summary_rows], fontsize=9)
    ax.set_title("Candidate-linked HOMER motif instances", fontsize=12, pad=12)
    ax.set_xticks(np.arange(-0.5, len(columns), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(summary_rows), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.0)
    ax.tick_params(which="minor", bottom=False, left=False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02, ticks=[0.25, 0.75])
    cbar.ax.set_yticklabels(["absent", "present"])
    fig.tight_layout()
    fig.savefig(figdir / "candidate_locus_homer_motif_instance_heatmap.png", dpi=300)
    fig.savefig(figdir / "candidate_locus_homer_motif_instance_heatmap.pdf")
    plt.close(fig)


def write_note(path: Path, summary_rows: list[dict[str, object]], context_rows: list[dict[str, object]]) -> None:
    n = len(summary_rows)
    dev = sum(int(row["n_regions_with_developmental_tf_motif"]) > 0 for row in summary_rows)
    cpg = sum(int(row["n_regions_with_cpg_ctcf_gc_motif"]) > 0 for row in summary_rows)
    stat = sum(int(row["n_regions_with_stat_motif"]) > 0 for row in summary_rows)
    lines = [
        "# Candidate-locus HOMER motif-instance scan",
        "",
        "This follow-up scans candidate-linked altered regions with HOMER motif files selected from the existing enrichment outputs.",
        "",
        "Important scope rule: motif instances indicate sequence compatibility only. They do not prove TF binding, and STAT motif instances are reported separately from STAT motif enrichment.",
        "",
        "## Main outputs",
        "",
        "- Gene summary: `results/mechanism_question/candidate_locus_homer_instance_scan/candidate_locus_homer_motif_instance_gene_summary.tsv`",
        "- Region hit table: `results/mechanism_question/candidate_locus_homer_instance_scan/candidate_locus_homer_motif_instance_region_hits.tsv`",
        "- Context summary: `results/mechanism_question/candidate_locus_homer_instance_scan/candidate_locus_homer_motif_instance_context_summary.tsv`",
        "- Heatmap: `figures/mechanism_question/candidate_locus_homer_motif_instance_heatmap.png`",
        "",
        "## High-level result",
        "",
        f"- Candidates with developmental TF motif instances: {dev}/{n}.",
        f"- Candidates with CpG/CTCF/GC motif instances: {cpg}/{n}.",
        f"- Candidates with STAT motif instances: {stat}/{n}.",
        "",
        "## Interpretation",
        "",
        "This strengthens the candidate-locus interpretation for developmental and CpG/CTCF motif context when motif instances occur inside the same candidate-linked DMR/probe/ATAC regions. STAT remains weaker: any STAT rows should be read as motif-presence checks from non-enriched STAT motif families, not as evidence for STAT-driven chromatin regulation.",
        "",
    ]
    context_df = pd.DataFrame(context_rows)
    if not context_df.empty:
        top = context_df.sort_values(["motif_category", "fraction_regions_with_motif"], ascending=[True, False]).head(12)
        lines.extend(["## Top context hit rates", ""])
        for row in top.to_dict(orient="records"):
            lines.append(
                f"- {row['source_group']} / {row['motif_category']}: {row['n_regions_with_motif']}/{row['n_regions_scanned']} regions"
            )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines))


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    figdir = Path(args.figdir)
    input_dir = outdir / "inputs"
    scan_dir = outdir / "annotatePeaks"
    outdir.mkdir(parents=True, exist_ok=True)
    scan_dir.mkdir(parents=True, exist_ok=True)

    final = read_tsv(FINAL_CANDIDATES)
    evidence = read_tsv(FINAL_EVIDENCE)
    final_genes = final["gene_symbol"].astype(str).tolist()

    region_groups = build_candidate_regions(final_genes, evidence, args.probe_window_bp)
    motifs_by_group = select_motifs(args.q_threshold)

    all_region_rows: list[dict[str, object]] = []
    all_hit_rows: list[dict[str, object]] = []
    scan_manifest: list[dict[str, object]] = []

    for group in REGION_GROUPS:
        group_name = group["name"]
        rows = region_groups.get(group_name, [])
        if not rows:
            continue
        motif_records: list[dict[str, object]] = []
        for analysis_group in group["analysis_groups"]:
            motif_records.extend(motifs_by_group.get(analysis_group, []))
        motif_records.extend([record for record in STAT_CONTROL_MOTIFS if Path(record["motif_path"]).exists()])
        dedup: dict[Path, dict[str, object]] = {Path(record["motif_path"]): record for record in motif_records}
        motif_records = list(dedup.values())
        if not motif_records:
            continue

        bed_path = input_dir / f"{group_name}.bed"
        output_path = scan_dir / f"{group_name}.annotatePeaks.nmotifs.tsv"
        log_path = scan_dir / f"{group_name}.annotatePeaks.log"
        write_bed(bed_path, rows)
        run_homer_scan(bed_path, group["genome"], motif_records, output_path, log_path)

        region_by_peak_id = {row["peak_id"]: row for row in rows}
        hits = parse_homer_output(output_path, region_by_peak_id, motif_records)
        all_region_rows.extend(rows)
        all_hit_rows.extend(hits)
        scan_manifest.append(
            {
                "source_group": group_name,
                "genome": group["genome"],
                "bed": bed_path,
                "annotatePeaks_output": output_path,
                "annotatePeaks_log": log_path,
                "n_regions": len(rows),
                "n_motif_files": len(motif_records),
                "n_region_motif_hits": len(hits),
                "analysis_groups": ";".join(group["analysis_groups"]),
            }
        )

    write_tsv(
        outdir / "candidate_locus_homer_motif_instance_scan_manifest.tsv",
        scan_manifest,
        [
            "source_group",
            "genome",
            "bed",
            "annotatePeaks_output",
            "annotatePeaks_log",
            "n_regions",
            "n_motif_files",
            "n_region_motif_hits",
            "analysis_groups",
        ],
    )
    write_tsv(
        outdir / "candidate_locus_homer_motif_instance_regions.tsv",
        all_region_rows,
        ["peak_id", "source_group", "region_kind", "genome", "gene_symbol", "region_id", "chrom", "start", "end", "direction", "score", "extra"],
    )
    write_tsv(
        outdir / "candidate_locus_homer_motif_instance_region_hits.tsv",
        all_hit_rows,
        [
            "peak_id",
            "source_group",
            "region_kind",
            "genome",
            "gene_symbol",
            "region_id",
            "chrom",
            "start",
            "end",
            "direction",
            "motif_name",
            "motif_families",
            "motif_category",
            "motif_count",
            "context_enriched",
            "stat_negative_control",
            "support_call",
            "best_q_value",
            "analysis_group",
        ],
    )

    gene_summary, context_summary = summarize_hits(final_genes, all_region_rows, all_hit_rows)
    write_tsv(
        outdir / "candidate_locus_homer_motif_instance_gene_summary.tsv",
        gene_summary,
        [
            "gene_symbol",
            "n_regions_scanned",
            "n_human_dnmt3l_regions",
            "n_human_ds_probe_regions",
            "n_mouse_dnmt3l_regions",
            "n_fetal_ds_atac_peak_regions",
            "n_regions_with_developmental_tf_motif",
            "developmental_tf_source_groups",
            "developmental_tf_motif_families",
            "n_regions_with_cpg_ctcf_gc_motif",
            "cpg_ctcf_gc_source_groups",
            "cpg_ctcf_gc_motif_families",
            "n_regions_with_stat_motif",
            "stat_source_groups",
            "stat_motif_families",
            "n_regions_with_other_supported_motif",
            "other_supported_source_groups",
            "other_supported_motif_families",
            "stat_instance_interpretation",
            "candidate_motif_instance_interpretation",
        ],
    )
    write_tsv(
        outdir / "candidate_locus_homer_motif_instance_context_summary.tsv",
        context_summary,
        [
            "source_group",
            "region_kind",
            "motif_category",
            "n_regions_scanned",
            "n_regions_with_motif",
            "fraction_regions_with_motif",
            "n_candidate_genes_with_motif",
            "genes_with_motif",
            "motif_families",
        ],
    )
    make_heatmap(gene_summary, figdir)
    write_note(outdir / "candidate_locus_homer_motif_instance_scan.md", gene_summary, context_summary)

    print(f"Scanned {len(all_region_rows)} candidate-linked regions")
    print(f"Wrote {len(all_hit_rows)} region x motif hit rows")
    print(outdir / "candidate_locus_homer_motif_instance_gene_summary.tsv")


if __name__ == "__main__":
    main()
