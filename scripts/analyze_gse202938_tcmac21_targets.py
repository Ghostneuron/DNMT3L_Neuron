#!/usr/bin/env python3
"""First-pass target-gene analysis for GSE202938 TcMAC21 cerebellum RNA-seq.

The submitter-supplied normalized-count matrix uses identifiers like
``chr1.10263`` instead of gene symbols. This script first tests whether those
identifiers can be treated as coordinate-like bins on the custom reference
described by the paper (mouse mm39 plus human chr21). When that mapping is too
uncertain, it writes an audit table and stops short of reporting gene-level
effects.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
from dataclasses import dataclass
from pathlib import Path


MOUSE_TARGETS = {"Dnmt3l", "Dnmt3a", "Dnmt3b", "Stat1", "Stat3"}
HUMAN_TARGETS = {"DNMT3L", "DNMT3A", "DNMT3B", "STAT1", "STAT3"}

SAMPLE_GROUPS = {
    "TcMAC21_euploid": ["Eu1.so.bam", "Eu2.so.bam", "Eu3.so.bam", "Eu8.so.bam"],
    "TcMAC21": ["Mac1.so.bam", "Mac2.so.bam", "Mac3.so.bam", "Mac4.so.bam"],
    "Ts65Dn_euploid": ["Eu4.so.bam", "Eu5.so.bam", "Eu6.so.bam", "Eu7.so.bam"],
    "Ts65Dn": ["Ts1.so.bam", "Ts2.so.bam", "Ts3.so.bam", "Ts4.so.bam"],
}


@dataclass
class Gene:
    species: str
    gene_name: str
    chrom: str
    start: int
    end: int
    strand: str
    gene_id: str
    gene_type: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--counts",
        default="data/geo/GSE202938/supplementary/GSE202938_normalized_counts.all.csv.gz",
        help="GSE202938 submitter normalized count matrix.",
    )
    parser.add_argument(
        "--mouse-gtf",
        default="data/reference/gencode/gencode.vM32.annotation.gtf.gz",
        help="Mouse GENCODE GTF, GRCm39/mm39 coordinates.",
    )
    parser.add_argument(
        "--human-gtf",
        default="data/reference/gencode/gencode.v39.annotation.gtf.gz",
        help="Human GENCODE GTF, GRCh38 coordinates.",
    )
    parser.add_argument("--mouse-chrom-sizes", default="data/reference/gencode/mm39.chrom.sizes")
    parser.add_argument("--human-chrom-sizes", default="data/reference/gencode/hg38.chrom.sizes")
    parser.add_argument("--outdir", default="results/gse202938_tcmac21_targets")
    parser.add_argument(
        "--bin-size",
        type=float,
        default=5500.0,
        help="Approximate coordinate-bin size to test for chrN.number IDs.",
    )
    parser.add_argument(
        "--window",
        type=int,
        default=20,
        help="Number of row-ID bins around the predicted bin to inspect.",
    )
    return parser.parse_args()


def parse_attrs(value: str) -> dict[str, str]:
    return dict(re.findall(r'([A-Za-z_]+) "([^"]+)"', value))


def read_target_genes(mouse_gtf: Path, human_gtf: Path) -> list[Gene]:
    genes: list[Gene] = []
    for species, path, targets in (
        ("mouse", mouse_gtf, MOUSE_TARGETS),
        ("human", human_gtf, HUMAN_TARGETS),
    ):
        with gzip.open(path, "rt") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[2] != "gene":
                    continue
                attrs = parse_attrs(fields[8])
                gene_name = attrs.get("gene_name", "")
                if gene_name not in targets:
                    continue
                genes.append(
                    Gene(
                        species=species,
                        gene_name=gene_name,
                        chrom=fields[0],
                        start=int(fields[3]),
                        end=int(fields[4]),
                        strand=fields[6],
                        gene_id=attrs.get("gene_id", ""),
                        gene_type=attrs.get("gene_type", ""),
                    )
                )
    return genes


def read_all_genes(mouse_gtf: Path, human_gtf: Path) -> dict[tuple[str, str], list[Gene]]:
    genes_by_chrom: dict[tuple[str, str], list[Gene]] = {}
    for species, path in (("mouse", mouse_gtf), ("human", human_gtf)):
        with gzip.open(path, "rt") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[2] != "gene":
                    continue
                attrs = parse_attrs(fields[8])
                gene = Gene(
                    species=species,
                    gene_name=attrs.get("gene_name", ""),
                    chrom=fields[0],
                    start=int(fields[3]),
                    end=int(fields[4]),
                    strand=fields[6],
                    gene_id=attrs.get("gene_id", ""),
                    gene_type=attrs.get("gene_type", ""),
                )
                genes_by_chrom.setdefault((species, gene.chrom), []).append(gene)
    return genes_by_chrom


def read_counts(path: Path) -> tuple[list[str], dict[str, dict[str, float]]]:
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle)
        sample_columns = [name for name in reader.fieldnames or [] if name]
        rows = {}
        for row in reader:
            row_id = row.pop("") if "" in row else row.pop(None, "")
            rows[row_id] = {sample: safe_float(row.get(sample, "")) for sample in sample_columns}
    return sample_columns, rows


def read_chrom_sizes(path: Path) -> dict[str, int]:
    sizes = {}
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                sizes[fields[0]] = int(fields[1])
    return sizes


def custom_reference_offsets(
    mouse_sizes: dict[str, int],
    human_sizes: dict[str, int],
    bin_size: float,
) -> dict[tuple[str, str], float]:
    """Return bin offsets for the paper's custom mouse-plus-HSA21 reference.

    GSE202938 methods describe STAR alignment to mm39 with human chromosome 21
    added. The submitter row IDs are consistent with bins along a concatenated
    reference where human chr21 comes first, followed by primary mouse
    chromosomes in lexicographic order.
    """
    offsets: dict[tuple[str, str], float] = {("human", "chr21"): 0.0}
    cumulative_bp = human_sizes["chr21"]
    primary_mouse = [
        chrom
        for chrom in mouse_sizes
        if re.fullmatch(r"chr(?:[0-9]+|X|Y|M)", chrom) and chrom != "chrM"
    ]
    for chrom in sorted(primary_mouse):
        offsets[("mouse", chrom)] = cumulative_bp / bin_size
        cumulative_bp += mouse_sizes[chrom]
    return offsets


def safe_float(value: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return number if math.isfinite(number) else float("nan")


def row_numeric_ids(rows: dict[str, dict[str, float]]) -> dict[str, set[int]]:
    by_chrom: dict[str, set[int]] = {}
    for row_id in rows:
        match = re.fullmatch(r"(chr[^.]+)\.(\d+)", row_id)
        if not match:
            continue
        chrom, number = match.groups()
        by_chrom.setdefault(chrom, set()).add(int(number))
    return by_chrom


def fmt(value: float) -> str:
    if not math.isfinite(value):
        return ""
    return f"{value:.6g}"


def mean(values: list[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    return sum(finite) / len(finite) if finite else float("nan")


def summarize_candidate_row(row_values: dict[str, float], pseudocount: float = 0.1) -> dict[str, str]:
    group_means = {}
    for group, samples in SAMPLE_GROUPS.items():
        group_means[group] = mean([row_values[sample] for sample in samples if sample in row_values])
    tc_delta = math.log2(group_means["TcMAC21"] + pseudocount) - math.log2(
        group_means["TcMAC21_euploid"] + pseudocount
    )
    ts_delta = math.log2(group_means["Ts65Dn"] + pseudocount) - math.log2(
        group_means["Ts65Dn_euploid"] + pseudocount
    )
    return {
        "tcmac21_euploid_mean": fmt(group_means["TcMAC21_euploid"]),
        "tcmac21_mean": fmt(group_means["TcMAC21"]),
        "tcmac21_delta_log2": fmt(tc_delta),
        "ts65dn_euploid_mean": fmt(group_means["Ts65Dn_euploid"]),
        "ts65dn_mean": fmt(group_means["Ts65Dn"]),
        "ts65dn_delta_log2": fmt(ts_delta),
    }


def interval_distance(position: float, gene: Gene) -> float:
    if gene.start <= position <= gene.end:
        return 0.0
    return min(abs(position - gene.start), abs(position - gene.end))


def nearest_gene(position: float, genes: list[Gene]) -> tuple[Gene | None, float]:
    if not genes:
        return None, float("nan")
    gene = min(genes, key=lambda candidate: interval_distance(position, candidate))
    return gene, interval_distance(position, gene)


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    mouse_gtf = Path(args.mouse_gtf)
    human_gtf = Path(args.human_gtf)
    genes = read_target_genes(mouse_gtf, human_gtf)
    genes_by_chrom = read_all_genes(mouse_gtf, human_gtf)
    sample_columns, rows = read_counts(Path(args.counts))
    numeric_by_chrom = row_numeric_ids(rows)
    offsets = custom_reference_offsets(
        read_chrom_sizes(Path(args.mouse_chrom_sizes)),
        read_chrom_sizes(Path(args.human_chrom_sizes)),
        args.bin_size,
    )

    gene_rows = []
    candidate_rows = []
    for gene in genes:
        # Human chr21 in the custom reference is called chr21. Other human
        # targets are not expected to appear in this GSE202938 matrix because
        # the added human sequence was chr21 only.
        matrix_chrom = gene.chrom
        if gene.species == "human" and gene.chrom != "chr21":
            expected = "not_in_custom_reference"
            predicted_bin = ""
            nearest = ""
            distance = ""
        else:
            values = numeric_by_chrom.get(matrix_chrom, set())
            offset = offsets.get((gene.species, matrix_chrom))
            predicted_number = int(round((offset or 0.0) + gene.start / args.bin_size))
            expected = "global_coordinate_bin_audit" if offset is not None else "missing_reference_offset"
            nearest_number = min(values, key=lambda item: abs(item - predicted_number)) if values else None
            predicted_bin = str(predicted_number)
            nearest = f"{matrix_chrom}.{nearest_number}" if nearest_number is not None else ""
            distance = str(abs(nearest_number - predicted_number)) if nearest_number is not None else ""

            for number in sorted(n for n in values if abs(n - predicted_number) <= args.window):
                row_id = f"{matrix_chrom}.{number}"
                inferred_position = (number - (offset or 0.0)) * args.bin_size
                nearby_gene, nearby_distance = nearest_gene(
                    inferred_position,
                    genes_by_chrom.get((gene.species, matrix_chrom), []),
                )
                candidate_rows.append(
                    {
                        "species": gene.species,
                        "gene_symbol": gene.gene_name,
                        "gene_id": gene.gene_id,
                        "matrix_row_id": row_id,
                        "bin_distance_from_gene_start": str(number - predicted_number),
                        "inferred_position_bp": str(int(round(inferred_position))),
                        "nearest_gencode_gene": nearby_gene.gene_name if nearby_gene else "",
                        "nearest_gencode_gene_id": nearby_gene.gene_id if nearby_gene else "",
                        "nearest_gene_distance_bp": fmt(nearby_distance),
                        "nearest_is_target_gene": str(bool(nearby_gene and nearby_gene.gene_name == gene.gene_name)),
                        **summarize_candidate_row(rows[row_id]),
                    }
                )

        gene_rows.append(
            {
                "species": gene.species,
                "gene_symbol": gene.gene_name,
                "gene_id": gene.gene_id,
                "gene_type": gene.gene_type,
                "chrom": gene.chrom,
                "start": str(gene.start),
                "end": str(gene.end),
                "strand": gene.strand,
                "mapping_status": expected,
                "predicted_local_bin": predicted_bin,
                "nearest_matrix_row_id": nearest,
                "nearest_bin_distance": distance,
            }
        )

    write_tsv(
        outdir / "gse202938_target_gene_coordinate_audit.tsv",
        gene_rows,
        [
            "species",
            "gene_symbol",
            "gene_id",
            "gene_type",
            "chrom",
            "start",
            "end",
            "strand",
            "mapping_status",
            "predicted_local_bin",
            "nearest_matrix_row_id",
            "nearest_bin_distance",
        ],
    )
    write_tsv(
        outdir / "gse202938_candidate_bin_expression.tsv",
        candidate_rows,
        [
            "species",
            "gene_symbol",
            "gene_id",
            "matrix_row_id",
            "bin_distance_from_gene_start",
            "inferred_position_bp",
            "nearest_gencode_gene",
            "nearest_gencode_gene_id",
            "nearest_gene_distance_bp",
            "nearest_is_target_gene",
            "tcmac21_euploid_mean",
            "tcmac21_mean",
            "tcmac21_delta_log2",
            "ts65dn_euploid_mean",
            "ts65dn_mean",
            "ts65dn_delta_log2",
        ],
    )
    print(f"Wrote coordinate audit for {len(gene_rows)} target genes to {outdir}")
    print(f"Wrote {len(candidate_rows)} nearby candidate row summaries")
    print("No gene-level target-effect table was written because the matrix lacks a trusted gene-ID map.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
