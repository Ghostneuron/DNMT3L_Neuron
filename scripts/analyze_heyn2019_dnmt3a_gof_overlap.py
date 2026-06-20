#!/usr/bin/env python3
"""Compare DNMT3L-linked candidates with Heyn 2019 DNMT3A-GOF DMR gene sets."""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


DEFAULT_OUTDIR = "results/heyn2019_dnmt3a_gof_overlap"
DEFAULT_FIGDIR = "figures/heyn2019_dnmt3a_gof_overlap"
DEFAULT_SUPPDIR = "manuscript/Supplementary Tables"


@dataclass(frozen=True)
class GeneTss:
    chrom: str
    tss: int
    start: int
    end: int
    strand: str
    gene_id: str
    gene_symbol: str
    gene_biotype: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--heyn-human-dmrs",
        default="data/heyn2019_dnmt3a_gof/raw/41588_2018_274_MOESM3_table2_human_hyper_dmr.txt",
    )
    parser.add_argument(
        "--heyn-mouse-dmrs",
        default="data/heyn2019_dnmt3a_gof/raw/41588_2018_274_MOESM5_table5_mouse_npc_hyper_dmr.txt",
    )
    parser.add_argument(
        "--human-gtf",
        default="data/reference/gencode/Homo_sapiens.GRCh37.75.gtf.gz",
    )
    parser.add_argument(
        "--mouse-gtf",
        default="data/reference/gencode/Mus_musculus.GRCm38.91.gtf.gz",
    )
    parser.add_argument(
        "--final-candidates",
        default="results/manuscript_summary/final_candidate_manuscript_summary.tsv",
    )
    parser.add_argument(
        "--candidate-integration",
        default="results/mechanism_question/candidate_locus_mechanism_integration.tsv",
    )
    parser.add_argument(
        "--human-methylation-genes",
        default="results/human_methylation_convergence/human_methylation_gene_support.tsv",
    )
    parser.add_argument(
        "--human-dnmt3l-dmrs",
        default="results/human_candidate_dmrs/GSE157078_candidate_dmrs.tsv",
    )
    parser.add_argument(
        "--mouse-dnmt3l-dmrs",
        default="results/mouse_candidate_dmrs/GSE271949_candidate_dmrs.tsv",
    )
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--figdir", default=DEFAULT_FIGDIR)
    parser.add_argument("--suppdir", default=DEFAULT_SUPPDIR)
    return parser.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", errors="replace")
    return path.open(errors="replace")


ATTR_RE = re.compile(r'([A-Za-z0-9_]+) "([^"]*)"')


def parse_attrs(attr_text: str) -> dict[str, str]:
    return {match.group(1): match.group(2) for match in ATTR_RE.finditer(attr_text)}


def normalize_chrom(chrom: str) -> str:
    chrom = str(chrom).strip()
    if not chrom:
        return chrom
    if chrom.startswith("chr"):
        return chrom
    return f"chr{chrom}"


def read_gene_tss(gtf_path: Path) -> tuple[dict[str, list[GeneTss]], set[str]]:
    genes: dict[str, list[GeneTss]] = defaultdict(list)
    symbols: set[str] = set()
    with open_text(gtf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = parse_attrs(parts[8])
            biotype = attrs.get("gene_biotype") or attrs.get("gene_type") or ""
            if biotype != "protein_coding":
                continue
            chrom = normalize_chrom(parts[0])
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            tss = start if strand == "+" else end
            symbol = attrs.get("gene_name", attrs.get("gene_id", ""))
            gene = GeneTss(
                chrom=chrom,
                tss=tss,
                start=start,
                end=end,
                strand=strand,
                gene_id=attrs.get("gene_id", ""),
                gene_symbol=symbol,
                gene_biotype=biotype,
            )
            genes[chrom].append(gene)
            symbols.add(symbol)
    for chrom in genes:
        genes[chrom].sort(key=lambda item: item.tss)
    return genes, symbols


def read_bed4(path: Path, source: str) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    with open_text(path) as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            rows.append(
                {
                    "source": source,
                    "dmr_id": parts[3] if len(parts) > 3 and parts[3] else f"{source}_{line_no:06d}",
                    "chrom": normalize_chrom(parts[0]),
                    "start": int(parts[1]),
                    "end": int(parts[2]),
                }
            )
    return pd.DataFrame(rows)


def nearest_genes_for_interval(
    chrom_genes: list[GeneTss],
    start: int,
    end: int,
) -> tuple[list[GeneTss], int]:
    if not chrom_genes:
        return [], -1
    tsses = [gene.tss for gene in chrom_genes]
    positions = {
        max(0, min(len(chrom_genes) - 1, bisect.bisect_left(tsses, start) + offset))
        for offset in (-2, -1, 0, 1, 2)
    }
    positions.update(
        max(0, min(len(chrom_genes) - 1, bisect.bisect_left(tsses, end) + offset))
        for offset in (-2, -1, 0, 1, 2)
    )
    best_distance = math.inf
    best: list[GeneTss] = []
    for pos in sorted(positions):
        gene = chrom_genes[pos]
        if start <= gene.tss <= end:
            distance = 0
        else:
            distance = min(abs(gene.tss - start), abs(gene.tss - end))
        if distance < best_distance:
            best_distance = distance
            best = [gene]
        elif distance == best_distance:
            best.append(gene)
    return best, int(best_distance)


def annotate_nearest_genes(dmrs: pd.DataFrame, genes_by_chrom: dict[str, list[GeneTss]]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for row in dmrs.to_dict("records"):
        genes, distance = nearest_genes_for_interval(
            genes_by_chrom.get(row["chrom"], []),
            int(row["start"]),
            int(row["end"]),
        )
        for gene in genes:
            annotated = dict(row)
            annotated.update(
                {
                    "nearest_gene_id": gene.gene_id,
                    "nearest_gene_symbol": gene.gene_symbol,
                    "nearest_gene_biotype": gene.gene_biotype,
                    "nearest_gene_tss": gene.tss,
                    "distance_to_nearest_tss": distance,
                }
            )
            rows.append(annotated)
    return pd.DataFrame(rows)


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def split_ids(value: object) -> list[str]:
    if value is None:
        return []
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return []
    return [item.strip() for item in re.split(r"[;,]", text) if item.strip()]


def read_project_dmrs(path: Path, source: str) -> pd.DataFrame:
    df = read_tsv(path)
    out = pd.DataFrame(
        {
            "source": source,
            "dmr_id": df["dmr_id"],
            "chrom": df["chrom"].map(normalize_chrom),
            "start": df["start"].astype(int),
            "end": df["end"].astype(int),
            "direction": df.get("direction", ""),
            "mean_delta_beta": df.get("mean_delta_beta", ""),
        }
    )
    return out


def interval_overlaps(query: pd.DataFrame, target: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    target_by_chrom = {chrom: sub.sort_values("start") for chrom, sub in target.groupby("chrom")}
    for q in query.to_dict("records"):
        sub = target_by_chrom.get(q["chrom"])
        if sub is None:
            continue
        starts = sub["start"].to_numpy()
        stop = bisect.bisect_right(starts, int(q["end"]))
        for t in sub.iloc[:stop].to_dict("records"):
            if int(t["end"]) < int(q["start"]):
                continue
            overlap_bp = min(int(q["end"]), int(t["end"])) - max(int(q["start"]), int(t["start"]))
            if overlap_bp > 0:
                rows.append(
                    {
                        "query_source": q["source"],
                        "query_dmr_id": q["dmr_id"],
                        "query_chrom": q["chrom"],
                        "query_start": q["start"],
                        "query_end": q["end"],
                        "target_source": t["source"],
                        "target_dmr_id": t["dmr_id"],
                        "target_chrom": t["chrom"],
                        "target_start": t["start"],
                        "target_end": t["end"],
                        "target_direction": t.get("direction", ""),
                        "target_mean_delta_beta": t.get("mean_delta_beta", ""),
                        "overlap_bp": overlap_bp,
                    }
                )
    return pd.DataFrame(rows)


def collapse_counts(df: pd.DataFrame, symbol_col: str) -> dict[str, int]:
    if df.empty:
        return {}
    return df.groupby(symbol_col)["dmr_id"].nunique().to_dict()


def semicolon(values: Iterable[object]) -> str:
    unique = sorted({str(value) for value in values if str(value)})
    return ";".join(unique)


def fisher_overlap(query: set[str], heyn: set[str], universe: set[str]) -> tuple[int, int, int, int, float, float]:
    query_input_size = len(set(query))
    query = set(query) & universe
    heyn = set(heyn) & universe
    overlap = query & heyn
    a = len(overlap)
    b = len(query - heyn)
    c = len(heyn - query)
    d = len(universe - query - heyn)
    odds, pvalue = fisher_exact([[a, b], [c, d]], alternative="greater")
    return query_input_size, len(query), len(heyn), a, float(odds), float(pvalue)


def bh_fdr(pvalues: list[float]) -> list[float]:
    n = len(pvalues)
    order = np.argsort(pvalues)
    out = [1.0] * n
    prev = 1.0
    for rank, idx in enumerate(order[::-1], start=1):
        original_rank = n - rank + 1
        value = min(prev, pvalues[idx] * n / original_rank)
        out[idx] = value
        prev = value
    return out


def build_final_candidate_table(
    final_df: pd.DataFrame,
    integration_df: pd.DataFrame,
    human_counts: dict[str, int],
    mouse_counts: dict[str, int],
    human_coord_overlaps: pd.DataFrame,
    mouse_coord_overlaps: pd.DataFrame,
) -> pd.DataFrame:
    merged = final_df.merge(
        integration_df[
            [
                "gene_symbol",
                "mouse_symbol",
                "human_dnmt3l_oe_dmr_ids",
                "mouse_dnmt3l_dmr_ids",
                "candidate_locus_mechanism_call",
                "methylation_dependence_inference",
            ]
        ],
        on="gene_symbol",
        how="left",
    )
    human_by_project = defaultdict(list)
    human_heyn_by_project = defaultdict(list)
    if not human_coord_overlaps.empty:
        for row in human_coord_overlaps.to_dict("records"):
            human_by_project[row["target_dmr_id"]].append(row["query_dmr_id"])
            human_heyn_by_project[row["target_dmr_id"]].append(row["query_dmr_id"])
    mouse_by_project = defaultdict(list)
    if not mouse_coord_overlaps.empty:
        for row in mouse_coord_overlaps.to_dict("records"):
            mouse_by_project[row["target_dmr_id"]].append(row["query_dmr_id"])

    rows: list[dict[str, object]] = []
    for row in merged.to_dict("records"):
        human_dmrs = split_ids(row.get("human_dnmt3l_oe_dmr_ids"))
        mouse_dmrs = split_ids(row.get("mouse_dnmt3l_dmr_ids"))
        human_overlap_heyn = sorted({item for dmr in human_dmrs for item in human_heyn_by_project.get(dmr, [])})
        mouse_overlap_heyn = sorted({item for dmr in mouse_dmrs for item in mouse_by_project.get(dmr, [])})
        human_gene_hit = int(human_counts.get(row["gene_symbol"], 0) > 0)
        mouse_symbol = row.get("mouse_symbol", "")
        mouse_gene_hit = int(bool(mouse_symbol) and mouse_counts.get(mouse_symbol, 0) > 0)
        coord_hit = int(bool(human_overlap_heyn) or bool(mouse_overlap_heyn))
        rows.append(
            {
                "gene_symbol": row["gene_symbol"],
                "mouse_symbol": mouse_symbol,
                "tier": row["tier"],
                "candidate_group": row["candidate_group"],
                "mechanism_class": row["mechanism_class"],
                "candidate_locus_mechanism_call": row.get("candidate_locus_mechanism_call", ""),
                "heyn_human_fibroblast_nearest_gene_hit": human_gene_hit,
                "heyn_human_fibroblast_nearest_dmr_count": human_counts.get(row["gene_symbol"], 0),
                "heyn_mouse_npc_nearest_gene_hit": mouse_gene_hit,
                "heyn_mouse_npc_nearest_dmr_count": mouse_counts.get(mouse_symbol, 0) if mouse_symbol else 0,
                "heyn_any_nearest_gene_hit": int(human_gene_hit or mouse_gene_hit),
                "human_dnmt3l_oe_dmr_ids": ";".join(human_dmrs),
                "human_dnmt3l_oe_dmr_count": len(human_dmrs),
                "heyn_human_coord_overlap_dmr_ids": ";".join(human_overlap_heyn),
                "heyn_human_coord_overlap_count": len(human_overlap_heyn),
                "mouse_dnmt3l_oe_dmr_ids": ";".join(mouse_dmrs),
                "mouse_dnmt3l_oe_dmr_count": len(mouse_dmrs),
                "heyn_mouse_coord_overlap_dmr_ids": ";".join(mouse_overlap_heyn),
                "heyn_mouse_coord_overlap_count": len(mouse_overlap_heyn),
                "heyn_any_coordinate_overlap": coord_hit,
                "interpretation": interpret_candidate(human_gene_hit, mouse_gene_hit, coord_hit, row["candidate_group"]),
            }
        )
    return pd.DataFrame(rows)


def interpret_candidate(human_hit: int, mouse_hit: int, coord_hit: int, group: str) -> str:
    if human_hit and mouse_hit:
        return "independent DNMT3A-GOF Polycomb-linked support in both human and mouse nearest-gene sets"
    if human_hit or mouse_hit:
        species = "human" if human_hit else "mouse NPC"
        return f"independent DNMT3A-GOF Polycomb-linked support in {species} nearest-gene set"
    if coord_hit:
        return "coordinate-level DMR overlap without nearest-gene set match"
    if "Methylation-heavy" in group:
        return "no direct Heyn hit; remains Polycomb/developmental by project-specific evidence"
    return "no direct Heyn hit"


def plot_summary(summary: pd.DataFrame, final_table: pd.DataFrame, figdir: Path) -> None:
    figdir.mkdir(parents=True, exist_ok=True)
    plot_df = summary.copy()
    plot_df["overlap_fraction"] = plot_df["overlap_count"] / plot_df["query_size"].replace(0, np.nan)
    labels = [
        "Final candidates\nhuman",
        "Final candidates\nmouse",
        "Human methylation-\nlinked genes",
        "Mouse ortholog\nmethylation-linked",
    ]
    fig, ax = plt.subplots(figsize=(7.4, 4.5))
    colors = ["#4f7cac", "#6d9f71", "#9b6aa0", "#c28b45"]
    bars = ax.bar(np.arange(len(plot_df)), plot_df["overlap_fraction"], color=colors, edgecolor="#333333", linewidth=0.8)
    for bar, (_, row) in zip(bars, plot_df.iterrows()):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.01,
            f"{int(row['overlap_count'])}/{int(row['query_size'])}",
            ha="center",
            va="bottom",
            fontsize=10,
        )
    ax.set_xticks(np.arange(len(plot_df)))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Overlap fraction with Heyn nearest-gene set", fontsize=10)
    ax.set_ylim(0, max(0.12, float(plot_df["overlap_fraction"].max()) * 1.35))
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#dddddd", linewidth=0.6)
    fig.tight_layout()
    fig.savefig(figdir / "heyn2019_dnmt3a_gof_overlap_summary.png", dpi=300)
    fig.savefig(figdir / "heyn2019_dnmt3a_gof_overlap_summary.pdf")
    plt.close(fig)

    columns = [
        ("Human nearest gene", "heyn_human_fibroblast_nearest_gene_hit"),
        ("Mouse NPC nearest gene", "heyn_mouse_npc_nearest_gene_hit"),
        ("Human DMR coord.", "heyn_human_coord_overlap_count"),
        ("Mouse DMR coord.", "heyn_mouse_coord_overlap_count"),
    ]
    heat = final_table[[col for _, col in columns]].copy()
    heat["heyn_human_coord_overlap_count"] = (heat["heyn_human_coord_overlap_count"].astype(float) > 0).astype(int)
    heat["heyn_mouse_coord_overlap_count"] = (heat["heyn_mouse_coord_overlap_count"].astype(float) > 0).astype(int)
    genes = final_table["gene_symbol"].tolist()
    fig_height = max(6.2, 0.24 * len(genes) + 1.4)
    fig, ax = plt.subplots(figsize=(6.7, fig_height))
    cmap = plt.matplotlib.colors.ListedColormap(["#f2f2ef", "#2f6f8f"])
    ax.imshow(heat.to_numpy(dtype=float), aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(range(len(columns)))
    ax.set_xticklabels([label for label, _ in columns], rotation=35, ha="right", fontsize=9)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=8.5, fontstyle="italic")
    ax.tick_params(length=0)
    for i in range(len(genes)):
        for j in range(len(columns)):
            value = int(heat.iloc[i, j])
            ax.text(j, i, "yes" if value else "", ha="center", va="center", fontsize=7.5, color="white")
    ax.set_title("Final Candidate Overlap With Heyn DNMT3A-GOF DMR Sets", fontsize=11, pad=10)
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.tight_layout()
    fig.savefig(figdir / "heyn2019_final_candidate_overlap_heatmap.png", dpi=300)
    fig.savefig(figdir / "heyn2019_final_candidate_overlap_heatmap.pdf")
    plt.close(fig)


def write_note(
    outdir: Path,
    summary: pd.DataFrame,
    final_table: pd.DataFrame,
    human_dmrs: pd.DataFrame,
    mouse_dmrs: pd.DataFrame,
) -> None:
    final_hits = final_table[final_table["heyn_any_nearest_gene_hit"].astype(int) == 1]
    coord_hits = final_table[final_table["heyn_any_coordinate_overlap"].astype(int) == 1]
    lines = [
        "# Heyn 2019 DNMT3A-GOF Polycomb Comparator",
        "",
        "## Source",
        "",
        "- Heyn et al., Nature Genetics 2019, DOI 10.1038/s41588-018-0274-x.",
        "- Supplementary Table 2: human patient fibroblast hypermethylated DMRs, hg19 coordinates.",
        "- Supplementary Table 5: mouse neural differentiation hypermethylated DMRs, mm10 coordinates.",
        "- Gene sets were reconstructed by nearest protein-coding TSS using the annotation releases described by the paper: Ensembl GRCh37 release 75 for human and GRCm38 release 91 for mouse.",
        "",
        "## Result",
        "",
        f"- Heyn human hyper-DMRs: {len(human_dmrs)} DMRs.",
        f"- Heyn mouse NPC hyper-DMRs: {len(mouse_dmrs)} DMRs.",
    ]
    for row in summary.to_dict("records"):
        genes = row["overlap_genes"] if row["overlap_genes"] else "none"
        mapped_note = (
            f"{int(row['query_size'])} testable genes"
            if int(row["query_input_size"]) == int(row["query_size"])
            else f"{int(row['query_size'])} of {int(row['query_input_size'])} query genes testable in the annotation universe"
        )
        lines.append(
            f"- {row['query_set']}: {row['overlap_count']}/{row['query_size']} ({mapped_note}) overlap with {row['heyn_set']} "
            f"(Fisher P={row['fisher_pvalue']:.3g}, FDR={row['fisher_fdr']:.3g}); genes: {genes}."
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This is useful as an independent DNMT3A/Polycomb comparator rather than primary DNMT3L evidence. "
            "A positive overlap means a DNMT3L-linked candidate or methylation-linked gene also lies near regions that become hypermethylated when DNMT3A PWWP gain-of-function disrupts normal Polycomb-protected methylation boundaries. "
            "It supports the plausibility of a DNMT3A-mediated Polycomb/developmental methylation axis, but it does not prove that DNMT3L directly causes the same event in DS brain.",
            "",
            f"Final candidates with nearest-gene Heyn support: {semicolon(final_hits['gene_symbol']) if not final_hits.empty else 'none'}.",
            f"Final candidates with coordinate-level DMR overlap: {semicolon(coord_hits['gene_symbol']) if not coord_hits.empty else 'none'}.",
        ]
    )
    (outdir / "heyn2019_dnmt3a_gof_overlap.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    figdir = Path(args.figdir)
    suppdir = Path(args.suppdir)
    outdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)
    suppdir.mkdir(parents=True, exist_ok=True)

    human_genes_by_chrom, human_universe = read_gene_tss(Path(args.human_gtf))
    mouse_genes_by_chrom, mouse_universe = read_gene_tss(Path(args.mouse_gtf))

    heyn_human_dmrs = read_bed4(Path(args.heyn_human_dmrs), "Heyn2019_human_fibroblast_hyper_DMR")
    heyn_mouse_dmrs = read_bed4(Path(args.heyn_mouse_dmrs), "Heyn2019_mouse_NPC_hyper_DMR")
    heyn_human_annot = annotate_nearest_genes(heyn_human_dmrs, human_genes_by_chrom)
    heyn_mouse_annot = annotate_nearest_genes(heyn_mouse_dmrs, mouse_genes_by_chrom)

    heyn_human_annot.to_csv(outdir / "heyn2019_human_hyper_dmr_nearest_tss.tsv", sep="\t", index=False)
    heyn_mouse_annot.to_csv(outdir / "heyn2019_mouse_npc_hyper_dmr_nearest_tss.tsv", sep="\t", index=False)

    heyn_human_symbols = set(heyn_human_annot["nearest_gene_symbol"].dropna())
    heyn_mouse_symbols = set(heyn_mouse_annot["nearest_gene_symbol"].dropna())
    human_counts = collapse_counts(heyn_human_annot, "nearest_gene_symbol")
    mouse_counts = collapse_counts(heyn_mouse_annot, "nearest_gene_symbol")

    human_project_dmrs = read_project_dmrs(Path(args.human_dnmt3l_dmrs), "GSE157078_DNMT3L_OE_DMR")
    mouse_project_dmrs = read_project_dmrs(Path(args.mouse_dnmt3l_dmrs), "GSE271949_Dnmt3l_OE_DMR")
    human_coord = interval_overlaps(heyn_human_dmrs, human_project_dmrs)
    mouse_coord = interval_overlaps(heyn_mouse_dmrs, mouse_project_dmrs)
    human_coord.to_csv(outdir / "heyn2019_human_vs_gse157078_dmr_coordinate_overlaps.tsv", sep="\t", index=False)
    mouse_coord.to_csv(outdir / "heyn2019_mouse_vs_gse271949_dmr_coordinate_overlaps.tsv", sep="\t", index=False)

    final_df = read_tsv(Path(args.final_candidates))
    integration_df = read_tsv(Path(args.candidate_integration))
    methyl_df = read_tsv(Path(args.human_methylation_genes))
    methyl_df["human_methylation_support_score_num"] = pd.to_numeric(
        methyl_df["human_methylation_support_score"], errors="coerce"
    ).fillna(0)
    methyl_linked = methyl_df[
        (methyl_df["human_methylation_support_score_num"] > 0) & (methyl_df["gene_type"] == "protein_coding")
    ].copy()
    mouse_methyl_linked = methyl_linked[
        (methyl_linked["mouse_candidate"] == "yes") & (methyl_linked["mouse_gene_symbol"] != "")
    ].copy()

    final_table = build_final_candidate_table(
        final_df,
        integration_df,
        human_counts,
        mouse_counts,
        human_coord,
        mouse_coord,
    )
    final_table.to_csv(outdir / "final_candidate_heyn2019_overlap.tsv", sep="\t", index=False)

    methyl_overlap = methyl_linked.copy()
    methyl_overlap["heyn_human_fibroblast_nearest_gene_hit"] = methyl_overlap["gene_symbol"].isin(
        heyn_human_symbols
    ).astype(int)
    methyl_overlap["heyn_human_fibroblast_nearest_dmr_count"] = methyl_overlap["gene_symbol"].map(
        human_counts
    ).fillna(0).astype(int)
    methyl_overlap["heyn_mouse_npc_nearest_gene_hit"] = methyl_overlap["mouse_gene_symbol"].isin(
        heyn_mouse_symbols
    ).astype(int)
    methyl_overlap["heyn_mouse_npc_nearest_dmr_count"] = methyl_overlap["mouse_gene_symbol"].map(
        mouse_counts
    ).fillna(0).astype(int)
    keep_cols = [
        "gene_symbol",
        "gene_type",
        "human_methylation_support_score",
        "combined_methylation_expression_score",
        "n_human_triple_phase3_any_ds_probes",
        "mouse_candidate",
        "mouse_gene_symbol",
        "mouse_evidence_score",
        "heyn_human_fibroblast_nearest_gene_hit",
        "heyn_human_fibroblast_nearest_dmr_count",
        "heyn_mouse_npc_nearest_gene_hit",
        "heyn_mouse_npc_nearest_dmr_count",
    ]
    methyl_overlap[keep_cols].to_csv(
        outdir / "methylation_linked_gene_heyn2019_overlap.tsv", sep="\t", index=False
    )

    summary_rows: list[dict[str, object]] = []
    comparisons = [
        (
            "final_candidates_human_symbols",
            set(final_df["gene_symbol"]),
            "Heyn2019_human_fibroblast_hyper_DMR_nearest_protein_coding_genes",
            heyn_human_symbols,
            human_universe,
        ),
        (
            "final_candidates_mouse_symbols",
            set(integration_df["mouse_symbol"]) - {""},
            "Heyn2019_mouse_NPC_hyper_DMR_nearest_protein_coding_genes",
            heyn_mouse_symbols,
            mouse_universe,
        ),
        (
            "human_methylation_linked_protein_coding_genes_score_gt0",
            set(methyl_linked["gene_symbol"]),
            "Heyn2019_human_fibroblast_hyper_DMR_nearest_protein_coding_genes",
            heyn_human_symbols,
            human_universe,
        ),
        (
            "mouse_candidate_orthologs_from_human_methylation_linked_genes",
            set(mouse_methyl_linked["mouse_gene_symbol"]),
            "Heyn2019_mouse_NPC_hyper_DMR_nearest_protein_coding_genes",
            heyn_mouse_symbols,
            mouse_universe,
        ),
    ]
    pvalues: list[float] = []
    for query_name, query, heyn_name, heyn, universe in comparisons:
        query_input_size, query_size, heyn_size, overlap_count, odds, pvalue = fisher_overlap(query, heyn, universe)
        pvalues.append(pvalue)
        summary_rows.append(
            {
                "query_set": query_name,
                "heyn_set": heyn_name,
                "universe_size": len(universe),
                "query_input_size": query_input_size,
                "query_size": query_size,
                "heyn_size": heyn_size,
                "overlap_count": overlap_count,
                "overlap_genes": semicolon(set(query) & set(heyn) & set(universe)),
                "fisher_odds_ratio": odds,
                "fisher_pvalue": pvalue,
            }
        )
    fdrs = bh_fdr(pvalues)
    for row, fdr in zip(summary_rows, fdrs):
        row["fisher_fdr"] = fdr
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(outdir / "heyn2019_overlap_enrichment_summary.tsv", sep="\t", index=False)

    group_rows: list[dict[str, object]] = []
    for group, sub in final_table.groupby("candidate_group", sort=False):
        for source_label, hit_col in [
            ("human_fibroblast_nearest_gene", "heyn_human_fibroblast_nearest_gene_hit"),
            ("mouse_npc_nearest_gene", "heyn_mouse_npc_nearest_gene_hit"),
            ("any_nearest_gene", "heyn_any_nearest_gene_hit"),
        ]:
            hits = sub[sub[hit_col].astype(int) == 1]
            group_rows.append(
                {
                    "candidate_group": group,
                    "heyn_source": source_label,
                    "group_size": len(sub),
                    "hit_count": len(hits),
                    "hit_fraction": len(hits) / len(sub) if len(sub) else np.nan,
                    "hit_genes": semicolon(hits["gene_symbol"]),
                }
            )
    group_summary = pd.DataFrame(group_rows)
    group_summary.to_csv(outdir / "final_candidate_group_heyn2019_overlap_summary.tsv", sep="\t", index=False)

    source_summary = pd.DataFrame(
        [
            {
                "source": "Heyn2019_human_fibroblast_hyper_DMR",
                "genome_build": "hg19_GRCh37",
                "dmr_count": len(heyn_human_dmrs),
                "nearest_protein_coding_gene_count": len(heyn_human_symbols),
                "annotation": "Ensembl_GRCh37_release_75_nearest_TSS",
                "source_file": args.heyn_human_dmrs,
            },
            {
                "source": "Heyn2019_mouse_NPC_hyper_DMR",
                "genome_build": "mm10_GRCm38",
                "dmr_count": len(heyn_mouse_dmrs),
                "nearest_protein_coding_gene_count": len(heyn_mouse_symbols),
                "annotation": "Ensembl_GRCm38_release_91_nearest_TSS",
                "source_file": args.heyn_mouse_dmrs,
            },
        ]
    )
    source_summary.to_csv(outdir / "heyn2019_gene_set_source_summary.tsv", sep="\t", index=False)

    plot_summary(summary, final_table, figdir)
    write_note(outdir, summary, final_table, heyn_human_dmrs, heyn_mouse_dmrs)

    # Paper-facing supplementary table files.
    source_summary.to_csv(suppdir / "ST19a_heyn2019_gene_set_source_summary.tsv", sep="\t", index=False)
    final_table.to_csv(suppdir / "ST19b_final_candidate_heyn2019_overlap.tsv", sep="\t", index=False)
    methyl_overlap[keep_cols].to_csv(
        suppdir / "ST19c_methylation_linked_gene_heyn2019_overlap.tsv", sep="\t", index=False
    )
    summary.to_csv(suppdir / "ST19d_heyn2019_overlap_enrichment_summary.tsv", sep="\t", index=False)
    group_summary.to_csv(suppdir / "ST19e_final_candidate_group_heyn2019_overlap_summary.tsv", sep="\t", index=False)
    human_coord.to_csv(suppdir / "ST19f_heyn2019_human_dmr_coordinate_overlaps.tsv", sep="\t", index=False)
    mouse_coord.to_csv(suppdir / "ST19g_heyn2019_mouse_dmr_coordinate_overlaps.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
