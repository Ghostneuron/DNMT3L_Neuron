#!/usr/bin/env python3
"""Draw locus-level methylation and chromatin panels for primary candidates."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
from collections import defaultdict
from glob import glob
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


DEFAULT_MATRIX = "results/final_candidate_evidence/final_candidate_evidence_matrix_with_regulatory_context.tsv"
DEFAULT_HUMAN_LINKS = "results/human_methylation_convergence/human_methylation_probe_gene_links.tsv"
DEFAULT_MOUSE_LINKS = "results/integrated_candidate_genes/mouse_triple_hit_dmr_gene_links.tsv"
DEFAULT_FIGDIR = "figures/primary_candidate_locus_panels"
DEFAULT_OUTDIR = "results/final_candidate_evidence"
DEFAULT_NOTE = "notes/primary_candidate_locus_panels.md"
DEFAULT_GENES = "PLPP3,FOXP1,FGFR2,ZFHX3,SIAH3"


HUMAN_TRACKS = [
    ("OCR", ["data/geo/GSE149268/supplementary/GSE149268_annotation-ocr-hg38.bed.gz"], "#7a9e3f"),
    ("ATAC peaks", ["data/geo/GSE149268/supplementary/GSM*_atac-*-optimal_idr-hg38.narrowPeak.gz"], "#2f7fb8"),
    ("H3K27ac", ["data/geo/GSE149268/supplementary/GSM*_chip-*_H3K27ac-hg38.bed.gz"], "#d9822b"),
    ("H3K4me1", ["data/geo/GSE149268/supplementary/GSM*_chip-*_H3K4me1-hg38.bed.gz"], "#c2a13a"),
    ("H3K27me3", ["data/geo/GSE149268/supplementary/GSM*_chip-*_H3K27me3-hg38.bed.gz"], "#7e5aa6"),
]


MOUSE_TRACKS = [
    ("ATAC peaks", ["data/geo/GSE172703/supplementary/*peaks_mm10.bed.gz"], "#2f7fb8"),
    ("H3K27ac", ["data/geo/GSE201400/supplementary/GSE201400_*_H3K27AC.bed.gz"], "#d9822b"),
    ("H3K4me3", ["data/geo/GSE201487/supplementary/GSE201487_*_H3K4ME3.bed.gz"], "#c2a13a"),
    ("H3K27me3", ["data/geo/GSE201488/supplementary/GSE201488_*_H3K27ME3.bed.gz"], "#7e5aa6"),
]


PROBE_COLORS = {
    "strict": "#b2182b",
    "any_ds": "#ef8a62",
    "dnmt3l_only": "#6b6b6b",
}

DMR_COLORS = {
    "hyper": "#b2182b",
    "hypo": "#2166ac",
    "unknown": "#6b6b6b",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", default=DEFAULT_MATRIX)
    parser.add_argument("--human-links", default=DEFAULT_HUMAN_LINKS)
    parser.add_argument("--mouse-links", default=DEFAULT_MOUSE_LINKS)
    parser.add_argument("--figdir", default=DEFAULT_FIGDIR)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--note", default=DEFAULT_NOTE)
    parser.add_argument("--genes", default=DEFAULT_GENES)
    parser.add_argument("--padding-bp", type=int, default=100000)
    return parser.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="", errors="replace")
    return path.open(newline="", errors="replace")


def clean_key(value: str | None) -> str:
    return (value or "").strip().replace("\r", "")


def is_yes(value: str | None) -> bool:
    return (value or "").strip().lower() == "yes"


def safe_float(value: object) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return math.nan
    return number if math.isfinite(number) else math.nan


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        reader.fieldnames = [clean_key(field) for field in reader.fieldnames or []]
        return [{clean_key(k): (v or "").strip() for k, v in row.items()} for row in reader]


def parse_bed_line(line: str) -> tuple[str, int, int] | None:
    if not line.strip() or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
        return None
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 3:
        return None
    try:
        start = int(float(parts[1]))
        end = int(float(parts[2]))
    except ValueError:
        return None
    if end <= start:
        return None
    return parts[0], start, end


def read_human_links(path: Path, genes: set[str]) -> dict[str, list[dict[str, object]]]:
    by_gene: dict[str, dict[str, dict[str, object]]] = defaultdict(dict)
    for row in read_tsv(path):
        gene = row.get("gene_symbol", "").upper()
        if gene not in genes:
            continue
        probe_id = row.get("probe_id", "")
        if not probe_id:
            continue
        try:
            start = int(float(row["probe_start"]))
            end = int(float(row["probe_end"]))
            gene_start = int(float(row["gene_start"]))
            gene_end = int(float(row["gene_end"]))
            gene_tss = int(float(row["gene_tss"]))
        except (KeyError, ValueError):
            continue
        if end <= start:
            end = start + 1
        strict = is_yes(row.get("human_triple_phase3_fetal_mean_positive"))
        any_ds = is_yes(row.get("human_triple_phase3_any_ds"))
        existing = by_gene[gene].get(probe_id)
        if existing is None:
            by_gene[gene][probe_id] = {
                "id": probe_id,
                "chrom": row["hg38_chrom"],
                "start": start,
                "end": end,
                "center": (start + end) / 2,
                "direction": row.get("direction", "unknown"),
                "gse157078_delta_beta": safe_float(row.get("gse157078_delta_beta")),
                "shsy_phase3_signed_delta": safe_float(row.get("shsy_phase3_signed_delta")),
                "ds_mean_signed_delta": safe_float(row.get("ds_mean_signed_delta")),
                "strict": strict,
                "any_ds": any_ds,
                "relationship": row.get("relationship", ""),
                "gene_chrom": row.get("gene_chrom", ""),
                "gene_start": gene_start,
                "gene_end": gene_end,
                "gene_strand": row.get("gene_strand", ""),
                "gene_tss": gene_tss,
            }
        else:
            existing["strict"] = bool(existing["strict"] or strict)
            existing["any_ds"] = bool(existing["any_ds"] or any_ds)
    return {gene: list(values.values()) for gene, values in by_gene.items()}


def read_mouse_links(path: Path, mouse_symbols: set[str]) -> dict[str, list[dict[str, object]]]:
    by_gene: dict[str, dict[str, dict[str, object]]] = defaultdict(dict)
    for row in read_tsv(path):
        gene = row.get("gene_symbol", "")
        if gene.upper() not in mouse_symbols:
            continue
        dmr_id = row.get("dmr_id", "")
        if not dmr_id:
            continue
        try:
            start = int(float(row["dmr_start"]))
            end = int(float(row["dmr_end"]))
            gene_start = int(float(row["gene_start"]))
            gene_end = int(float(row["gene_end"]))
            gene_tss = int(float(row["gene_tss"]))
        except (KeyError, ValueError):
            continue
        by_gene[gene.upper()][dmr_id] = {
            "id": dmr_id,
            "chrom": row["chrom"],
            "start": start,
            "end": end,
            "center": (start + end) / 2,
            "direction": row.get("dmr_direction", "unknown"),
            "delta_beta": safe_float(row.get("dmr_mean_delta_beta")),
            "relationship": row.get("relationship", ""),
            "gene_chrom": row.get("gene_chrom", ""),
            "gene_start": gene_start,
            "gene_end": gene_end,
            "gene_strand": row.get("gene_strand", ""),
            "gene_tss": gene_tss,
        }
    return {gene: list(values.values()) for gene, values in by_gene.items()}


def paths_from_globs(patterns: list[str]) -> list[Path]:
    paths: list[Path] = []
    for pattern in patterns:
        paths.extend(Path(path) for path in sorted(glob(pattern)))
    return paths


def collect_overlaps(patterns: list[str], chrom: str, start: int, end: int) -> list[tuple[int, int]]:
    intervals: list[tuple[int, int]] = []
    for path in paths_from_globs(patterns):
        with open_text(path) as handle:
            for line in handle:
                parsed = parse_bed_line(line)
                if parsed is None:
                    continue
                item_chrom, item_start, item_end = parsed
                if item_chrom != chrom:
                    continue
                if item_end <= start or item_start >= end:
                    continue
                intervals.append((max(start, item_start), min(end, item_end)))
    return merge_intervals(intervals)


def merge_intervals(intervals: list[tuple[int, int]], gap: int = 50) -> list[tuple[int, int]]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end + gap:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


def locus_bounds(gene_start: int, gene_end: int, intervals: list[dict[str, object]], padding: int) -> tuple[int, int]:
    starts = [gene_start]
    ends = [gene_end]
    for interval in intervals:
        starts.append(int(interval["start"]))
        ends.append(int(interval["end"]))
    raw_start = min(starts)
    raw_end = max(ends)
    width = raw_end - raw_start
    dynamic_padding = max(padding, min(250000, int(width * 0.08)))
    return max(0, raw_start - dynamic_padding), raw_end + dynamic_padding


def mb(value: float) -> float:
    return value / 1_000_000


def draw_segments(ax: plt.Axes, intervals: list[tuple[int, int]], y: float, color: str, linewidth: float = 4.0) -> None:
    for start, end in intervals:
        ax.plot([mb(start), mb(end)], [y, y], color=color, linewidth=linewidth, solid_capstyle="butt")


def draw_gene(ax: plt.Axes, start: int, end: int, tss: int, strand: str, y: float, label: str) -> None:
    ax.plot([mb(start), mb(end)], [y, y], color="#111111", linewidth=5, solid_capstyle="butt")
    marker = ">" if strand == "+" else "<"
    ax.scatter([mb(tss)], [y], marker=marker, s=70, color="#111111", zorder=5)
    ax.text(mb(start), y + 0.28, label, ha="left", va="bottom", fontsize=8, color="#111111")


def probe_kind(probe: dict[str, object]) -> str:
    if probe.get("strict"):
        return "strict"
    if probe.get("any_ds"):
        return "any_ds"
    return "dnmt3l_only"


def draw_human_panel(
    ax: plt.Axes,
    gene: str,
    probes: list[dict[str, object]],
    tracks: list[tuple[str, list[tuple[int, int]], str]],
    x_start: int,
    x_end: int,
) -> None:
    gene_start = int(probes[0]["gene_start"])
    gene_end = int(probes[0]["gene_end"])
    gene_tss = int(probes[0]["gene_tss"])
    gene_strand = str(probes[0]["gene_strand"])
    y_positions = {
        "Gene": 6.4,
        "DNMT3L CpGs": 5.4,
        "OCR": 4.3,
        "ATAC peaks": 3.4,
        "H3K27ac": 2.5,
        "H3K4me1": 1.6,
        "H3K27me3": 0.7,
    }
    draw_gene(ax, gene_start, gene_end, gene_tss, gene_strand, y_positions["Gene"], gene)

    probe_order = {"dnmt3l_only": 0, "any_ds": 1, "strict": 2}
    for probe in sorted(probes, key=lambda item: probe_order[probe_kind(item)]):
        x = mb(float(probe["center"]))
        kind = probe_kind(probe)
        ax.vlines(x, y_positions["DNMT3L CpGs"] - 0.25, y_positions["DNMT3L CpGs"] + 0.25, color=PROBE_COLORS[kind], linewidth=1.4)
        marker = "*" if kind == "strict" else "o"
        size = 75 if kind == "strict" else 30
        ax.scatter([x], [y_positions["DNMT3L CpGs"]], marker=marker, s=size, color=PROBE_COLORS[kind], edgecolor="white", linewidth=0.4, zorder=5)

    for label, intervals, color in tracks:
        draw_segments(ax, intervals, y_positions[label], color)

    ax.set_yticks([y_positions[label] for label in y_positions])
    ax.set_yticklabels(list(y_positions), fontsize=8)
    ax.set_xlim(mb(x_start), mb(x_end))
    ax.set_ylim(0.1, 7.1)
    ax.set_title(f"Human hg38: {probes[0]['chrom']}:{x_start:,}-{x_end:,}", fontsize=9, loc="left")
    ax.grid(axis="x", color="#dddddd", linewidth=0.6)
    ax.tick_params(axis="x", labelsize=8)
    ax.tick_params(axis="y", length=0)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)


def draw_mouse_panel(
    ax: plt.Axes,
    gene: str,
    dmrs: list[dict[str, object]],
    tracks: list[tuple[str, list[tuple[int, int]], str]],
    x_start: int,
    x_end: int,
) -> None:
    gene_start = int(dmrs[0]["gene_start"])
    gene_end = int(dmrs[0]["gene_end"])
    gene_tss = int(dmrs[0]["gene_tss"])
    gene_strand = str(dmrs[0]["gene_strand"])
    y_positions = {
        "Gene": 5.6,
        "Dnmt3l DMRs": 4.6,
        "ATAC peaks": 3.5,
        "H3K27ac": 2.6,
        "H3K4me3": 1.7,
        "H3K27me3": 0.8,
    }
    draw_gene(ax, gene_start, gene_end, gene_tss, gene_strand, y_positions["Gene"], gene)

    for dmr in dmrs:
        direction = str(dmr.get("direction", "unknown"))
        color = DMR_COLORS.get(direction, DMR_COLORS["unknown"])
        ax.plot([mb(int(dmr["start"])), mb(int(dmr["end"]))], [y_positions["Dnmt3l DMRs"], y_positions["Dnmt3l DMRs"]], color=color, linewidth=6, solid_capstyle="butt")
        ax.scatter([mb(float(dmr["center"]))], [y_positions["Dnmt3l DMRs"]], s=28, color=color, edgecolor="white", linewidth=0.4, zorder=5)

    for label, intervals, color in tracks:
        draw_segments(ax, intervals, y_positions[label], color)

    ax.set_yticks([y_positions[label] for label in y_positions])
    ax.set_yticklabels(list(y_positions), fontsize=8)
    ax.set_xlim(mb(x_start), mb(x_end))
    ax.set_ylim(0.2, 6.2)
    ax.set_xlabel("Genomic coordinate (Mb)", fontsize=9)
    ax.set_title(f"Mouse mm10: {dmrs[0]['chrom']}:{x_start:,}-{x_end:,}", fontsize=9, loc="left")
    ax.grid(axis="x", color="#dddddd", linewidth=0.6)
    ax.tick_params(axis="x", labelsize=8)
    ax.tick_params(axis="y", length=0)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)


def expression_summary(row: dict[str, str]) -> str:
    parts = [
        f"human expr {row.get('human_expression_support_level', '')}",
        f"mouse OE {row.get('mouse_gse268951_direction', '')} ({row.get('mouse_gse268951_log2fc', '')})",
        f"Dp10 {row.get('mouse_dp10_direction', '')} ({row.get('mouse_dp10_delta', '')})",
    ]
    if row.get("best_bulk_same_direction_context"):
        parts.append(f"best bulk {row['best_bulk_same_direction_context']} {row.get('best_bulk_same_direction_delta', '')}")
    if row.get("best_celltype_same_direction_label"):
        parts.append(f"best cell {row['best_celltype_same_direction_label']} {row.get('best_celltype_same_direction_delta', '')}")
    return " | ".join(parts)


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def plot_gene(
    row: dict[str, str],
    human_probes: list[dict[str, object]],
    mouse_dmrs: list[dict[str, object]],
    figdir: Path,
    padding: int,
) -> dict[str, str]:
    human_start, human_end = locus_bounds(
        int(human_probes[0]["gene_start"]),
        int(human_probes[0]["gene_end"]),
        human_probes,
        padding,
    )
    mouse_start, mouse_end = locus_bounds(
        int(mouse_dmrs[0]["gene_start"]),
        int(mouse_dmrs[0]["gene_end"]),
        mouse_dmrs,
        padding,
    )

    human_tracks = [
        (label, collect_overlaps(patterns, str(human_probes[0]["chrom"]), human_start, human_end), color)
        for label, patterns, color in HUMAN_TRACKS
    ]
    mouse_tracks = [
        (label, collect_overlaps(patterns, str(mouse_dmrs[0]["chrom"]), mouse_start, mouse_end), color)
        for label, patterns, color in MOUSE_TRACKS
    ]

    fig, axes = plt.subplots(2, 1, figsize=(12.4, 6.8), sharex=False)
    draw_human_panel(axes[0], row["gene_symbol"], human_probes, human_tracks, human_start, human_end)
    draw_mouse_panel(axes[1], row["mouse_symbol"], mouse_dmrs, mouse_tracks, mouse_start, mouse_end)

    fig.suptitle(f"{row['gene_symbol']} / {row['mouse_symbol']} locus evidence", fontsize=13, y=0.98)
    fig.text(0.02, 0.935, expression_summary(row), fontsize=8.4, ha="left", va="top")

    legend_handles = [
        Line2D([0], [0], marker="*", color="none", markerfacecolor=PROBE_COLORS["strict"], markeredgecolor="white", markersize=9, label="Strict human triple CpG"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=PROBE_COLORS["any_ds"], markeredgecolor="white", markersize=7, label="Human DNMT3L + DS CpG"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=PROBE_COLORS["dnmt3l_only"], markeredgecolor="white", markersize=7, label="Human DNMT3L CpG"),
        Line2D([0, 1], [0, 0], marker="o", color=DMR_COLORS["hyper"], markerfacecolor=DMR_COLORS["hyper"], markeredgecolor="white", linewidth=4, markersize=6, label="Mouse hyper-DMR"),
        Line2D([0, 1], [0, 0], marker="o", color=DMR_COLORS["hypo"], markerfacecolor=DMR_COLORS["hypo"], markeredgecolor="white", linewidth=4, markersize=6, label="Mouse hypo-DMR"),
        Patch(color="#7a9e3f", label="OCR"),
        Patch(color="#2f7fb8", label="ATAC"),
        Patch(color="#d9822b", label="H3K27ac"),
        Patch(color="#c2a13a", label="H3K4me1/3"),
        Patch(color="#7e5aa6", label="H3K27me3"),
    ]
    fig.legend(legend_handles, [handle.get_label() for handle in legend_handles], loc="lower center", ncols=5, frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0.09, 1, 0.93))

    output_base = figdir / f"{safe_name(row['gene_symbol'])}_locus_evidence"
    fig.savefig(output_base.with_suffix(".png"), dpi=300)
    fig.savefig(output_base.with_suffix(".pdf"))
    plt.close(fig)

    human_track_counts = {f"human_{label.lower().replace(' ', '_')}_intervals": str(len(intervals)) for label, intervals, _color in human_tracks}
    mouse_track_counts = {f"mouse_{label.lower().replace(' ', '_')}_intervals": str(len(intervals)) for label, intervals, _color in mouse_tracks}
    return {
        "gene_symbol": row["gene_symbol"],
        "mouse_symbol": row["mouse_symbol"],
        "human_locus": f"{human_probes[0]['chrom']}:{human_start}-{human_end}",
        "mouse_locus": f"{mouse_dmrs[0]['chrom']}:{mouse_start}-{mouse_end}",
        "human_probe_count": str(len(human_probes)),
        "human_strict_probe_count": str(sum(1 for probe in human_probes if probe.get("strict"))),
        "mouse_dmr_count": str(len(mouse_dmrs)),
        "png": str(output_base.with_suffix(".png")),
        "pdf": str(output_base.with_suffix(".pdf")),
        **human_track_counts,
        **mouse_track_counts,
    }


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_note(path: Path, summary_rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Primary candidate locus panels",
        "",
        "This pass generates locus-level panels for the five primary methylation-convergent candidates. Each panel shows human hg38 DNMT3L-sensitive CpGs and mouse mm10 Dnmt3l-sensitive DMRs separately, alongside fetal-brain ATAC and histone-mark annotations.",
        "",
        "Important interpretation: these figures show that candidate methylation intervals sit in plausible fetal-brain regulatory chromatin. They do not by themselves prove that methylation causally mediates the expression change.",
        "",
        "## Panel summary",
        "",
        "| Gene | Human locus | Human probes | Strict probes | Mouse locus | Mouse DMRs | Figure |",
        "|---|---|---:|---:|---|---:|---|",
    ]
    for row in summary_rows:
        lines.append(
            f"| {row['gene_symbol']} | {row['human_locus']} | {row['human_probe_count']} | "
            f"{row['human_strict_probe_count']} | {row['mouse_locus']} | {row['mouse_dmr_count']} | "
            f"`{row['png']}` |"
        )
    lines.extend(
        [
            "",
            "## Manuscript use",
            "",
            "These panels are best used as locus-level support for the overview candidate heatmap. A cautious results sentence would be: locus-level inspection showed that primary DNMT3L-linked candidates contain human DNMT3L-sensitive/DS-concordant CpGs and mouse Dnmt3l-sensitive DMRs within fetal-brain regulatory chromatin, including ATAC-accessible, H3K27ac-marked, or H3K27me3-marked regions.",
            "",
            "The strongest interpretation remains regulatory convergence, not confirmed causal methylation mediation.",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    genes = [gene.strip().upper() for gene in args.genes.split(",") if gene.strip()]
    gene_set = set(genes)
    matrix_rows = [row for row in read_tsv(Path(args.matrix)) if row.get("gene_symbol", "").upper() in gene_set]
    row_by_gene = {row["gene_symbol"].upper(): row for row in matrix_rows}
    mouse_symbols = {row["mouse_symbol"].upper() for row in matrix_rows}
    human_links = read_human_links(Path(args.human_links), gene_set)
    mouse_links = read_mouse_links(Path(args.mouse_links), mouse_symbols)

    figdir = Path(args.figdir)
    figdir.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict[str, str]] = []
    for gene in genes:
        row = row_by_gene.get(gene)
        if row is None:
            print(f"Skipping {gene}: not in matrix")
            continue
        probes = human_links.get(gene, [])
        dmrs = mouse_links.get(row["mouse_symbol"].upper(), [])
        if not probes or not dmrs:
            print(f"Skipping {gene}: probes={len(probes)} dmrs={len(dmrs)}")
            continue
        summary_rows.append(plot_gene(row, probes, dmrs, figdir, args.padding_bp))

    outdir = Path(args.outdir)
    write_tsv(outdir / "primary_candidate_locus_panel_summary.tsv", summary_rows)
    write_note(Path(args.note), summary_rows)
    print(f"Wrote {len(summary_rows)} locus panels to {figdir}")
    print(f"Wrote {outdir / 'primary_candidate_locus_panel_summary.tsv'}")
    print(f"Wrote {args.note}")


if __name__ == "__main__":
    main()
