#!/usr/bin/env python3
"""Count target-gene kmer evidence in GSE31454 Tc1/Tc0 brain bowtie files.

GSE31454 provides collapsed read alignments rather than a modern gene-count
matrix. The first field encodes a read-sequence count, e.g. ``1_3:1773:50``.
This script uses the read sequence column directly and counts conservative
target-gene kmers from GENCODE transcript FASTAs.

This is a first-pass target screen, not a replacement for full re-alignment and
gene-level quantification.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt


TARGETS = {
    "mouse": {"Dnmt3l", "Dnmt3a", "Dnmt3b", "Stat1", "Stat3"},
    "human": {"DNMT3L", "DNMT3A", "DNMT3B", "STAT1", "STAT3"},
}

SAMPLE_PATTERNS = [
    ("Tc0_mouse_brain", re.compile(r"Tc0_Brain")),
    ("Tc1_mouse_brain", re.compile(r"Tc1_Brain")),
    ("Human_brain", re.compile(r"Human_Brain")),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", default="data/geo/GSE31454/supplementary")
    parser.add_argument("--mouse-transcripts", default="data/reference/gencode/gencode.vM32.transcripts.fa.gz")
    parser.add_argument("--human-transcripts", default="data/reference/gencode/gencode.v39.transcripts.fa.gz")
    parser.add_argument("--outdir", default="results/gse31454_tc1_target_kmers")
    parser.add_argument("--figdir", default="figures/gse31454_tc1_target_kmers")
    parser.add_argument("--k", type=int, default=31)
    parser.add_argument("--max-base-fraction", type=float, default=0.80)
    return parser.parse_args()


def fasta_records(path: Path):
    name = None
    seq_parts = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq_parts).upper()
                name = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if name is not None:
        yield name, "".join(seq_parts).upper()


def gene_name_from_header(header: str) -> str:
    fields = header.split("|")
    return fields[5] if len(fields) > 5 else ""


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1].upper()


def is_low_complexity(kmer: str, max_base_fraction: float) -> bool:
    if "N" in kmer:
        return True
    length = len(kmer)
    return max(kmer.count(base) for base in "ACGT") / length > max_base_fraction


def iter_kmers(seq: str, k: int, max_base_fraction: float):
    if len(seq) < k:
        return
    for idx in range(0, len(seq) - k + 1):
        kmer = seq[idx : idx + k]
        if not is_low_complexity(kmer, max_base_fraction):
            yield kmer


def build_target_kmers(args: argparse.Namespace) -> tuple[dict[str, str], dict[str, int]]:
    raw_map: dict[str, str | None] = {}
    gene_kmer_counts: dict[str, set[str]] = defaultdict(set)
    for species, path in (
        ("mouse", Path(args.mouse_transcripts)),
        ("human", Path(args.human_transcripts)),
    ):
        for header, seq in fasta_records(path):
            gene = gene_name_from_header(header)
            if gene not in TARGETS[species]:
                continue
            key = f"{species}:{gene}"
            for source_seq in (seq, revcomp(seq)):
                for kmer in iter_kmers(source_seq, args.k, args.max_base_fraction):
                    current = raw_map.get(kmer)
                    if current is None and kmer in raw_map:
                        continue
                    if current is not None and current != key:
                        raw_map[kmer] = None
                        continue
                    raw_map[kmer] = key
                    gene_kmer_counts[key].add(kmer)

    unique_map = {kmer: key for kmer, key in raw_map.items() if key is not None}
    unique_counts = {key: 0 for species in TARGETS for key in [f"{species}:{gene}" for gene in TARGETS[species]]}
    for key in unique_map.values():
        unique_counts[key] = unique_counts.get(key, 0) + 1
    return unique_map, unique_counts


def sample_label(path: Path) -> str:
    for label, pattern in SAMPLE_PATTERNS:
        if pattern.search(path.name):
            return label
    return "unknown"


def read_weight(read_id: str) -> int:
    parts = read_id.rsplit(":", 2)
    if len(parts) >= 2 and parts[-2].isdigit():
        return int(parts[-2])
    return 1


def scan_sequence(seq: str, k: int, kmer_to_gene: dict[str, str]) -> set[str]:
    genes = set()
    seq = seq.upper()
    if len(seq) < k:
        return genes
    for idx in range(0, len(seq) - k + 1):
        gene = kmer_to_gene.get(seq[idx : idx + k])
        if gene:
            genes.add(gene)
            if len(genes) > 1:
                return genes
    return genes


def count_file(path: Path, kmer_to_gene: dict[str, str], k: int) -> tuple[dict[str, int], int, int, int]:
    counts: dict[str, int] = defaultdict(int)
    total_weighted_reads = 0
    matched_weighted_reads = 0
    ambiguous_weighted_reads = 0
    last_read_id = None
    with gzip.open(path, "rt", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            read_id = fields[0]
            if read_id == last_read_id:
                continue
            last_read_id = read_id
            weight = read_weight(read_id)
            total_weighted_reads += weight
            genes = scan_sequence(fields[4], k, kmer_to_gene)
            if len(genes) == 1:
                gene = next(iter(genes))
                counts[gene] += weight
                matched_weighted_reads += weight
            elif len(genes) > 1:
                ambiguous_weighted_reads += weight
    return counts, total_weighted_reads, matched_weighted_reads, ambiguous_weighted_reads


def fmt(value: float) -> str:
    if not math.isfinite(value):
        return ""
    return f"{value:.6g}"


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def plot_summary(summary_rows: list[dict[str, str]], figdir: Path) -> None:
    figdir.mkdir(parents=True, exist_ok=True)
    ordered_samples = ["Tc0_mouse_brain", "Tc1_mouse_brain", "Human_brain"]
    genes = sorted({row["gene_key"] for row in summary_rows})
    values = {(row["sample"], row["gene_key"]): float(row["cpm"]) for row in summary_rows}

    fig, ax = plt.subplots(figsize=(11, 5))
    x_positions = []
    labels = []
    heights = []
    colors = []
    palette = {
        "Tc0_mouse_brain": "#4C78A8",
        "Tc1_mouse_brain": "#F58518",
        "Human_brain": "#54A24B",
    }
    x = 0
    for gene in genes:
        for sample in ordered_samples:
            x_positions.append(x)
            labels.append(f"{gene}\n{sample.replace('_brain', '').replace('_', ' ')}")
            heights.append(values.get((sample, gene), 0.0))
            colors.append(palette[sample])
            x += 1
        x += 0.8
    ax.bar(x_positions, heights, color=colors, width=0.78)
    ax.set_ylabel("Target-kmer CPM")
    ax.set_yscale("symlog", linthresh=0.1)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(labels, rotation=70, ha="right", fontsize=8)
    ax.set_title("GSE31454 Tc1/Tc0 Brain Target-Kmer Signal")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(figdir / "gse31454_tc1_target_kmer_cpm.png", dpi=200)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    figdir = Path(args.figdir)
    kmer_to_gene, unique_kmer_counts = build_target_kmers(args)
    kmer_rows = [
        {"gene_key": key, "unique_kmers": str(count)}
        for key, count in sorted(unique_kmer_counts.items())
    ]
    write_tsv(outdir / "target_unique_kmer_counts.tsv", kmer_rows, ["gene_key", "unique_kmers"])

    files = sorted(Path(args.input_dir).glob("*Brain*.bowtie.gz"))
    file_rows = []
    combined_counts: dict[tuple[str, str], int] = defaultdict(int)
    combined_totals: dict[str, int] = defaultdict(int)
    combined_matched: dict[str, int] = defaultdict(int)
    combined_ambiguous: dict[str, int] = defaultdict(int)

    for path in files:
        label = sample_label(path)
        counts, total, matched, ambiguous = count_file(path, kmer_to_gene, args.k)
        combined_totals[label] += total
        combined_matched[label] += matched
        combined_ambiguous[label] += ambiguous
        for gene_key in sorted(unique_kmer_counts):
            value = counts.get(gene_key, 0)
            combined_counts[(label, gene_key)] += value
            cpm = value / total * 1_000_000 if total else float("nan")
            file_rows.append(
                {
                    "file": path.name,
                    "sample": label,
                    "gene_key": gene_key,
                    "weighted_target_reads": str(value),
                    "total_weighted_reads": str(total),
                    "matched_weighted_reads": str(matched),
                    "ambiguous_weighted_reads": str(ambiguous),
                    "cpm": fmt(cpm),
                }
            )

    write_tsv(
        outdir / "gse31454_file_target_kmer_counts.tsv",
        file_rows,
        [
            "file",
            "sample",
            "gene_key",
            "weighted_target_reads",
            "total_weighted_reads",
            "matched_weighted_reads",
            "ambiguous_weighted_reads",
            "cpm",
        ],
    )

    summary_rows = []
    for sample in sorted(combined_totals):
        total = combined_totals[sample]
        for gene_key in sorted(unique_kmer_counts):
            value = combined_counts[(sample, gene_key)]
            cpm = value / total * 1_000_000 if total else float("nan")
            summary_rows.append(
                {
                    "sample": sample,
                    "gene_key": gene_key,
                    "weighted_target_reads": str(value),
                    "total_weighted_reads": str(total),
                    "matched_weighted_reads": str(combined_matched[sample]),
                    "ambiguous_weighted_reads": str(combined_ambiguous[sample]),
                    "cpm": fmt(cpm),
                }
            )
    write_tsv(
        outdir / "gse31454_sample_target_kmer_counts.tsv",
        summary_rows,
        [
            "sample",
            "gene_key",
            "weighted_target_reads",
            "total_weighted_reads",
            "matched_weighted_reads",
            "ambiguous_weighted_reads",
            "cpm",
        ],
    )

    comparisons = []
    by_sample_gene = {(row["sample"], row["gene_key"]): float(row["cpm"]) for row in summary_rows}
    for gene_key in sorted(unique_kmer_counts):
        tc0 = by_sample_gene.get(("Tc0_mouse_brain", gene_key), 0.0)
        tc1 = by_sample_gene.get(("Tc1_mouse_brain", gene_key), 0.0)
        human = by_sample_gene.get(("Human_brain", gene_key), 0.0)
        comparisons.append(
            {
                "gene_key": gene_key,
                "tc0_mouse_brain_cpm": fmt(tc0),
                "tc1_mouse_brain_cpm": fmt(tc1),
                "human_brain_cpm": fmt(human),
                "tc1_vs_tc0_delta_log2_cpm_plus_0.1": fmt(math.log2(tc1 + 0.1) - math.log2(tc0 + 0.1)),
            }
        )
    write_tsv(
        outdir / "gse31454_tc1_vs_tc0_target_kmer_comparison.tsv",
        comparisons,
        [
            "gene_key",
            "tc0_mouse_brain_cpm",
            "tc1_mouse_brain_cpm",
            "human_brain_cpm",
            "tc1_vs_tc0_delta_log2_cpm_plus_0.1",
        ],
    )

    plot_summary(summary_rows, figdir)
    print(f"Built {len(kmer_to_gene):,} target-unique kmers")
    print(f"Processed {len(files)} GSE31454 brain bowtie files")
    print(f"Wrote outputs to {outdir} and {figdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
