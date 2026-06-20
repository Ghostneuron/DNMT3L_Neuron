#!/usr/bin/env python3
"""Overlap Dnmt3l-sensitive mouse methylation loci with GSE141587 references."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


DEFAULT_QUERY_DMRS = "results/mouse_candidate_dmrs_sex_consistent/GSE271949_candidate_dmrs.tsv"
DEFAULT_QUERY_DMCS = "results/mouse_candidate_dmrs_sex_consistent/GSE271949_candidate_dmcs.tsv"
DEFAULT_OUTDIR = "results/mouse_dnmt3l_dnmt3a_neuron_overlap"
TARGET_GENES = {"Dnmt3l", "Dnmt3a", "Dnmt3b", "Stat1", "Stat3"}


REFERENCE_SETS = [
    {
        "ref_set": "dnmt3a_p39_cko_hypo_cg_dmrs",
        "label": "P39 Dnmt3a cKO hypo CG-DMRs; DNMT3A-required neuronal mCG",
        "path": "data/geo/GSE141587/supplementary/GSE141587_MethylCseq_P39_Dnmt3aCKO_vs_Control_CG_DMR_delta0.3_cKO_hypo.bed.gz",
    },
    {
        "ref_set": "control_postnatal_mcg_gain_p0_hypo",
        "label": "Control P0-vs-P39 CG-DMRs where P0 is hypo; postnatal mCG gain",
        "path": "data/geo/GSE141587/supplementary/GSE141587_MethylCseq_P0_Control_vs_P39_Control_CG_DMR_delta0.3_P0_Control_hypo.bed.gz",
    },
    {
        "ref_set": "control_postnatal_mcg_loss_p39_hypo",
        "label": "Control P0-vs-P39 CG-DMRs where P39 is hypo; postnatal mCG loss",
        "path": "data/geo/GSE141587/supplementary/GSE141587_MethylCseq_P0_Control_vs_P39_Control_CG_DMR_delta0.3_P39_Control_hypo.bed.gz",
    },
    {
        "ref_set": "h3k27me3_p39_cko_up",
        "label": "P39 Dnmt3a cKO H3K27me3-up regions",
        "path": "data/geo/GSE141587/supplementary/GSE141587_ChIPseq_differential_H3K27me3_P39_Dnmt3aCKO_vs_P39_Control_up.bed.gz",
    },
    {
        "ref_set": "dmv_p0_control",
        "label": "P0 control DNA methylation valleys",
        "path": "data/geo/GSE141587/supplementary/GSE141587_MethylCseq_P0_Control_DMVs.bed.gz",
    },
    {
        "ref_set": "dmv_p39_control",
        "label": "P39 control DNA methylation valleys",
        "path": "data/geo/GSE141587/supplementary/GSE141587_MethylCseq_P39_Control_DMVs.bed.gz",
    },
    {
        "ref_set": "dmv_p39_dnmt3a_cko",
        "label": "P39 Dnmt3a cKO DNA methylation valleys",
        "path": "data/geo/GSE141587/supplementary/GSE141587_MethylCseq_P39_Dnmt3aCKO_DMVs.bed.gz",
    },
]


@dataclass(frozen=True)
class Interval:
    chrom: str
    start: int
    end: int
    item_id: str
    direction: str = ""
    delta_beta: float = float("nan")
    n_sites: str = ""
    value_a: str = ""
    value_b: str = ""
    label: str = ""

    @property
    def width(self) -> int:
        return max(0, self.end - self.start)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query-dmrs", default=DEFAULT_QUERY_DMRS)
    parser.add_argument("--query-dmcs", default=DEFAULT_QUERY_DMCS)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument(
        "--rna-de",
        default="data/geo/GSE141587/supplementary/GSE141587_RNAseq_P39_DEgenes_edgeR_qlfTest_Dnmt3aCKO_vs_Control.txt.gz",
    )
    parser.add_argument(
        "--rna-tpm",
        default="data/geo/GSE141587/supplementary/GSE141587_RNAseq_all_samples_genes_TPM.txt.gz",
    )
    return parser.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="", errors="replace")
    return path.open(newline="", errors="replace")


def safe_float(value: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return number if math.isfinite(number) else float("nan")


def fmt(value: float) -> str:
    return f"{value:.6g}" if math.isfinite(value) else ""


def mean(values: list[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    return sum(finite) / len(finite) if finite else float("nan")


def read_query_dmrs(path: Path) -> list[Interval]:
    intervals = []
    with open_text(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            intervals.append(
                Interval(
                    chrom=row["chrom"],
                    start=int(row["start"]),
                    end=int(row["end"]),
                    item_id=row["dmr_id"],
                    direction=row["direction"],
                    delta_beta=safe_float(row["mean_delta_beta"]),
                    n_sites=row["n_cpgs"],
                )
            )
    return sort_intervals(intervals)


def read_query_dmcs(path: Path) -> list[Interval]:
    intervals = []
    with open_text(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            intervals.append(
                Interval(
                    chrom=row["chrom"],
                    start=int(row["start"]),
                    end=int(row["end"]),
                    item_id=row["locus_id"],
                    direction=row["direction"],
                    delta_beta=safe_float(row["delta_beta"]),
                )
            )
    return sort_intervals(intervals)


def read_bed(path: Path, label: str) -> list[Interval]:
    intervals = []
    with open_text(path) as handle:
        for index, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            name = parts[3] if len(parts) > 3 else f"region_{index}"
            tokens = name.split("|")
            item_id = tokens[0] if tokens else name
            n_sites = tokens[1] if len(tokens) > 1 else ""
            value_a = tokens[2] if len(tokens) > 2 else ""
            value_b = tokens[3] if len(tokens) > 3 else ""
            intervals.append(
                Interval(
                    chrom=chrom,
                    start=start,
                    end=end,
                    item_id=item_id,
                    n_sites=n_sites,
                    value_a=value_a,
                    value_b=value_b,
                    label=label,
                )
            )
    return sort_intervals(intervals)


def sort_intervals(intervals: list[Interval]) -> list[Interval]:
    return sorted(intervals, key=lambda item: (chrom_key(item.chrom), item.start, item.end, item.item_id))


def chrom_key(chrom: str) -> tuple[int, int | str]:
    name = chrom.removeprefix("chr")
    if name.isdigit():
        return (0, int(name))
    special = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    return (0, special[name]) if name in special else (1, name)


def by_chrom(intervals: list[Interval]) -> dict[str, list[Interval]]:
    grouped: dict[str, list[Interval]] = defaultdict(list)
    for interval in intervals:
        grouped[interval.chrom].append(interval)
    return grouped


def overlap_len(a: Interval, b: Interval) -> int:
    return max(0, min(a.end, b.end) - max(a.start, b.start))


def interval_overlaps(queries: list[Interval], refs: list[Interval]) -> list[tuple[Interval, Interval, int]]:
    ref_by_chrom = by_chrom(refs)
    pairs = []
    for chrom, qlist in by_chrom(queries).items():
        rlist = ref_by_chrom.get(chrom, [])
        cursor = 0
        for query in qlist:
            while cursor < len(rlist) and rlist[cursor].end <= query.start:
                cursor += 1
            scan = cursor
            while scan < len(rlist) and rlist[scan].start < query.end:
                ref = rlist[scan]
                length = overlap_len(query, ref)
                if length > 0:
                    pairs.append((query, ref, length))
                scan += 1
    return pairs


def subset_queries(queries: list[Interval]) -> dict[str, list[Interval]]:
    return {
        "all": queries,
        "hyper": [query for query in queries if query.direction == "hyper"],
        "hypo": [query for query in queries if query.direction == "hypo"],
    }


def summarize_overlap(
    query_level: str,
    subset_name: str,
    queries: list[Interval],
    ref_set: str,
    ref_label: str,
    refs: list[Interval],
    pairs: list[tuple[Interval, Interval, int]],
) -> dict[str, str]:
    overlapped_ids = {query.item_id for query, _, _ in pairs}
    overlap_delta = [query.delta_beta for query in queries if query.item_id in overlapped_ids]
    return {
        "query_level": query_level,
        "query_subset": subset_name,
        "ref_set": ref_set,
        "ref_label": ref_label,
        "query_regions": str(len(queries)),
        "query_bp": str(sum(query.width for query in queries)),
        "ref_regions": str(len(refs)),
        "ref_bp": str(sum(ref.width for ref in refs)),
        "query_with_overlap": str(len(overlapped_ids)),
        "query_overlap_fraction": fmt(len(overlapped_ids) / len(queries)) if queries else "",
        "overlap_pairs": str(len(pairs)),
        "overlap_bp_sum": str(sum(length for _, _, length in pairs)),
        "mean_query_delta_beta_all": fmt(mean([query.delta_beta for query in queries])),
        "mean_query_delta_beta_overlapped": fmt(mean(overlap_delta)),
    }


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_pair_rows(
    query_level: str,
    ref_set: str,
    pairs: list[tuple[Interval, Interval, int]],
) -> list[dict[str, str]]:
    rows = []
    for query, ref, length in pairs:
        rows.append(
            {
                "query_level": query_level,
                "ref_set": ref_set,
                "query_id": query.item_id,
                "query_chrom": query.chrom,
                "query_start": str(query.start),
                "query_end": str(query.end),
                "query_direction": query.direction,
                "query_delta_beta": fmt(query.delta_beta),
                "query_n_sites": query.n_sites,
                "ref_id": ref.item_id,
                "ref_chrom": ref.chrom,
                "ref_start": str(ref.start),
                "ref_end": str(ref.end),
                "ref_n_sites_or_score": ref.n_sites,
                "ref_value_a": ref.value_a,
                "ref_value_b": ref.value_b,
                "overlap_bp": str(length),
            }
        )
    return rows


def build_flag_rows(dmrs: list[Interval], pair_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    by_query_ref: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    ref_sets = [ref["ref_set"] for ref in REFERENCE_SETS]
    for row in pair_rows:
        if row["query_level"] == "DMR":
            by_query_ref[(row["query_id"], row["ref_set"])].append(row)

    rows = []
    for dmr in dmrs:
        row = {
            "dmr_id": dmr.item_id,
            "chrom": dmr.chrom,
            "start": str(dmr.start),
            "end": str(dmr.end),
            "direction": dmr.direction,
            "n_cpgs": dmr.n_sites,
            "mean_delta_beta": fmt(dmr.delta_beta),
        }
        for ref_set in ref_sets:
            overlaps = by_query_ref.get((dmr.item_id, ref_set), [])
            row[f"{ref_set}_overlap"] = "yes" if overlaps else "no"
            row[f"{ref_set}_overlap_bp"] = str(sum(int(hit["overlap_bp"]) for hit in overlaps)) if overlaps else "0"
            row[f"{ref_set}_n_regions"] = str(len({hit["ref_id"] for hit in overlaps}))
        rows.append(row)
    return rows


def target_gene_expression_rows(de_path: Path, tpm_path: Path) -> list[dict[str, str]]:
    de_rows: dict[str, dict[str, str]] = {}
    with open_text(de_path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("geneName") in TARGET_GENES:
                de_rows[row["geneName"]] = row

    tpm_rows: dict[str, dict[str, str]] = {}
    with open_text(tpm_path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("geneName") in TARGET_GENES:
                tpm_rows[row["geneName"]] = row

    rows = []
    for gene in sorted(TARGET_GENES):
        de = de_rows.get(gene, {})
        tpm = tpm_rows.get(gene, {})
        control_values = [
            safe_float(tpm.get("P39_Control_rep1", "")),
            safe_float(tpm.get("P39_Control_rep2", "")),
        ]
        cko_values = [
            safe_float(tpm.get("P39_Dnmt3aKO_rep1", "")),
            safe_float(tpm.get("P39_Dnmt3aKO_rep2", "")),
        ]
        rows.append(
            {
                "gene": gene,
                "gene_id": tpm.get("geneID") or de.get("geneID", ""),
                "gene_type": tpm.get("geneType") or de.get("geneType", ""),
                "control_mean_tpm": fmt(mean(control_values)),
                "dnmt3a_cko_mean_tpm": fmt(mean(cko_values)),
                "logFC_cko_vs_control": de.get("logFC", ""),
                "FDR": de.get("FDR", ""),
                "in_edgeR_filtered_table": "yes" if de else "no",
                "TPM_P39_Control_rep1": tpm.get("P39_Control_rep1", ""),
                "TPM_P39_Control_rep2": tpm.get("P39_Control_rep2", ""),
                "TPM_P39_Dnmt3aKO_rep1": tpm.get("P39_Dnmt3aKO_rep1", ""),
                "TPM_P39_Dnmt3aKO_rep2": tpm.get("P39_Dnmt3aKO_rep2", ""),
            }
        )
    return rows


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    dmrs = read_query_dmrs(Path(args.query_dmrs))
    dmcs = read_query_dmcs(Path(args.query_dmcs))

    summary_rows = []
    pair_rows = []
    for ref_spec in REFERENCE_SETS:
        path = Path(ref_spec["path"])
        if not path.exists():
            continue
        refs = read_bed(path, label=ref_spec["label"])
        for query_level, queries in [("DMR", dmrs), ("CpG", dmcs)]:
            for subset_name, subset in subset_queries(queries).items():
                pairs = interval_overlaps(subset, refs)
                summary_rows.append(
                    summarize_overlap(
                        query_level=query_level,
                        subset_name=subset_name,
                        queries=subset,
                        ref_set=ref_spec["ref_set"],
                        ref_label=ref_spec["label"],
                        refs=refs,
                        pairs=pairs,
                    )
                )
                if subset_name == "all":
                    pair_rows.extend(build_pair_rows(query_level, ref_spec["ref_set"], pairs))

    summary_fields = [
        "query_level",
        "query_subset",
        "ref_set",
        "ref_label",
        "query_regions",
        "query_bp",
        "ref_regions",
        "ref_bp",
        "query_with_overlap",
        "query_overlap_fraction",
        "overlap_pairs",
        "overlap_bp_sum",
        "mean_query_delta_beta_all",
        "mean_query_delta_beta_overlapped",
    ]
    pair_fields = [
        "query_level",
        "ref_set",
        "query_id",
        "query_chrom",
        "query_start",
        "query_end",
        "query_direction",
        "query_delta_beta",
        "query_n_sites",
        "ref_id",
        "ref_chrom",
        "ref_start",
        "ref_end",
        "ref_n_sites_or_score",
        "ref_value_a",
        "ref_value_b",
        "overlap_bp",
    ]
    flag_rows = build_flag_rows(dmrs, pair_rows)
    flag_fields = list(flag_rows[0].keys()) if flag_rows else []
    priority_rows = [
        row
        for row in flag_rows
        if row["dnmt3a_p39_cko_hypo_cg_dmrs_overlap"] == "yes"
        and row["control_postnatal_mcg_gain_p0_hypo_overlap"] == "yes"
        and row["h3k27me3_p39_cko_up_overlap"] == "yes"
    ]
    priority_rows.sort(key=lambda row: abs(safe_float(row["mean_delta_beta"])), reverse=True)
    target_rows = target_gene_expression_rows(Path(args.rna_de), Path(args.rna_tpm))
    target_fields = [
        "gene",
        "gene_id",
        "gene_type",
        "control_mean_tpm",
        "dnmt3a_cko_mean_tpm",
        "logFC_cko_vs_control",
        "FDR",
        "in_edgeR_filtered_table",
        "TPM_P39_Control_rep1",
        "TPM_P39_Control_rep2",
        "TPM_P39_Dnmt3aKO_rep1",
        "TPM_P39_Dnmt3aKO_rep2",
    ]

    write_tsv(outdir / "gse271949_gse141587_overlap_summary.tsv", summary_rows, summary_fields)
    write_tsv(outdir / "gse271949_gse141587_overlap_pairs.tsv", pair_rows, pair_fields)
    write_tsv(outdir / "gse271949_dmr_reference_flags.tsv", flag_rows, flag_fields)
    write_tsv(outdir / "gse271949_dmr_priority_triple_hits.tsv", priority_rows, flag_fields)
    write_tsv(outdir / "gse141587_target_gene_expression.tsv", target_rows, target_fields)

    print(f"Wrote {len(summary_rows)} overlap summary rows")
    print(f"Wrote {len(pair_rows)} overlap pair rows")
    print(f"Wrote {len(flag_rows)} DMR flag rows")
    print(f"Wrote {len(priority_rows)} triple-hit priority DMR rows")
    print(f"Wrote {len(target_rows)} target gene expression rows")


if __name__ == "__main__":
    main()
