#!/usr/bin/env python3
"""Run matched-region permutation tests for candidate regulatory overlap."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from bisect import bisect_left
from collections import defaultdict
from glob import glob
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DEFAULT_MATRIX = "results/final_candidate_evidence/final_candidate_evidence_matrix_with_regulatory_context.tsv"
DEFAULT_HUMAN_LINKS = "results/human_methylation_convergence/human_methylation_probe_gene_links.tsv"
DEFAULT_MOUSE_LINKS = "results/integrated_candidate_genes/mouse_triple_hit_dmr_gene_links.tsv"
DEFAULT_OUTDIR = "results/mechanism_classification"
DEFAULT_FIGDIR = "figures/mechanism_classification"
DEFAULT_NOTE = "notes/matched_regulatory_permutation_first_pass.md"


HUMAN_TRACKS = {
    "human_any_regulatory_probe": [
        "data/geo/GSE149268/supplementary/GSE149268_annotation-ocr-hg38.bed.gz",
        "data/geo/GSE149268/supplementary/GSE149268_annotation-pre-hg38.bed.gz",
        "data/geo/GSE149268/supplementary/GSM*_atac-*-optimal_idr-hg38.narrowPeak.gz",
        "data/geo/GSE149268/supplementary/GSM*_chip-*_H3K27ac-hg38.bed.gz",
        "data/geo/GSE149268/supplementary/GSM*_chip-*_H3K4me1-hg38.bed.gz",
        "data/geo/GSE149268/supplementary/GSM*_chip-*_H3K27me3-hg38.bed.gz",
    ],
    "human_atac_or_ocr_probe": [
        "data/geo/GSE149268/supplementary/GSE149268_annotation-ocr-hg38.bed.gz",
        "data/geo/GSE149268/supplementary/GSM*_atac-*-optimal_idr-hg38.narrowPeak.gz",
    ],
    "human_h3k27ac_probe": [
        "data/geo/GSE149268/supplementary/GSM*_chip-*_H3K27ac-hg38.bed.gz",
    ],
    "human_h3k27me3_probe": [
        "data/geo/GSE149268/supplementary/GSM*_chip-*_H3K27me3-hg38.bed.gz",
    ],
}


MOUSE_TRACKS = {
    "mouse_any_regulatory_dmr": [
        "data/geo/GSE172703/supplementary/*peaks_mm10.bed.gz",
        "data/geo/GSE201400/supplementary/GSE201400_*_H3K27AC.bed.gz",
        "data/geo/GSE201487/supplementary/GSE201487_*_H3K4ME3.bed.gz",
        "data/geo/GSE201488/supplementary/GSE201488_*_H3K27ME3.bed.gz",
    ],
    "mouse_atac_dmr": [
        "data/geo/GSE172703/supplementary/*peaks_mm10.bed.gz",
    ],
    "mouse_h3k27ac_dmr": [
        "data/geo/GSE201400/supplementary/GSE201400_*_H3K27AC.bed.gz",
    ],
    "mouse_h3k27me3_dmr": [
        "data/geo/GSE201488/supplementary/GSE201488_*_H3K27ME3.bed.gz",
    ],
}


OUTPUT_FIELDS = [
    "candidate_set",
    "assay",
    "species",
    "background",
    "match_strategy",
    "candidate_records",
    "candidate_successes",
    "candidate_success_fraction",
    "null_mean_successes",
    "null_sd_successes",
    "null_mean_fraction",
    "fold_over_null",
    "empirical_p_greater",
    "empirical_p_less",
    "z_score",
    "n_permutations",
    "candidate_success_genes",
    "candidate_records_by_stratum",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", default=DEFAULT_MATRIX)
    parser.add_argument("--human-links", default=DEFAULT_HUMAN_LINKS)
    parser.add_argument("--mouse-links", default=DEFAULT_MOUSE_LINKS)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--figdir", default=DEFAULT_FIGDIR)
    parser.add_argument("--note", default=DEFAULT_NOTE)
    parser.add_argument("--permutations", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20260609)
    return parser.parse_args()


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


def fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}" if math.isfinite(value) else ""
    return str(value)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        reader.fieldnames = [clean_key(field) for field in reader.fieldnames or []]
        return [{clean_key(k): (v or "").strip() for k, v in row.items()} for row in reader]


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="", errors="replace")
    return path.open(newline="", errors="replace")


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


def paths_from_globs(patterns: list[str]) -> list[Path]:
    paths: list[Path] = []
    for pattern in patterns:
        paths.extend(Path(path) for path in sorted(glob(pattern)))
    return paths


def load_annotation_group(patterns: list[str]) -> dict[str, tuple[list[int], list[int]]]:
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for path in paths_from_globs(patterns):
        with open_text(path) as handle:
            for line in handle:
                parsed = parse_bed_line(line)
                if parsed is None:
                    continue
                chrom, start, end = parsed
                intervals[chrom].append((start, end))
    indexed: dict[str, tuple[list[int], list[int]]] = {}
    for chrom, values in intervals.items():
        values.sort()
        starts = [start for start, _end in values]
        prefix_max_ends: list[int] = []
        current = 0
        for _start, end in values:
            current = max(current, end)
            prefix_max_ends.append(current)
        indexed[chrom] = (starts, prefix_max_ends)
    return indexed


def has_overlap(index: dict[str, tuple[list[int], list[int]]], chrom: str, start: int, end: int) -> bool:
    starts, prefix_max_ends = index.get(chrom, ([], []))
    if not starts:
        return False
    pos = bisect_left(starts, end) - 1
    return pos >= 0 and prefix_max_ends[pos] > start


def load_indexes(tracks: dict[str, list[str]]) -> dict[str, dict[str, tuple[list[int], list[int]]]]:
    return {name: load_annotation_group(patterns) for name, patterns in tracks.items()}


def relationship_group(row: dict[str, str]) -> str:
    if is_yes(row.get("promoter_2kb")):
        return "promoter_2kb"
    if is_yes(row.get("gene_body")):
        return "gene_body"
    if is_yes(row.get("nearest_100kb")):
        return "nearest_100kb"
    return row.get("relationship", "other") or "other"


def build_human_records(
    path: Path,
    indexes: dict[str, dict[str, tuple[list[int], list[int]]]],
) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    seen: set[tuple[str, str]] = set()
    for row in read_tsv(path):
        if row.get("gene_type") != "protein_coding":
            continue
        symbol = row.get("gene_symbol", "").upper()
        probe_id = row.get("probe_id", "")
        chrom = row.get("hg38_chrom", "")
        if not symbol or not probe_id or not chrom:
            continue
        key = (symbol, probe_id)
        if key in seen:
            continue
        seen.add(key)
        try:
            start = int(float(row["probe_start"]))
            end = int(float(row["probe_end"]))
        except (KeyError, ValueError):
            continue
        if end <= start:
            end = start + 1
        records.append(
            {
                "gene_symbol": symbol,
                "id": probe_id,
                "chrom": chrom,
                "start": start,
                "end": end,
                "strict": is_yes(row.get("human_triple_phase3_fetal_mean_positive")),
                "direction": row.get("direction", "unknown") or "unknown",
                "stratum": f"{relationship_group(row)}|{row.get('direction', 'unknown') or 'unknown'}",
                "success": {
                    name: has_overlap(index, chrom, start, end)
                    for name, index in indexes.items()
                },
            }
        )
    return records


def build_mouse_records(
    path: Path,
    indexes: dict[str, dict[str, tuple[list[int], list[int]]]],
) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    seen: set[tuple[str, str]] = set()
    for row in read_tsv(path):
        symbol = row.get("gene_symbol", "").upper()
        dmr_id = row.get("dmr_id", "")
        chrom = row.get("chrom", "")
        if not symbol or not dmr_id or not chrom:
            continue
        key = (symbol, dmr_id)
        if key in seen:
            continue
        seen.add(key)
        try:
            start = int(float(row["dmr_start"]))
            end = int(float(row["dmr_end"]))
        except (KeyError, ValueError):
            continue
        direction = row.get("dmr_direction", "unknown") or "unknown"
        records.append(
            {
                "gene_symbol": symbol,
                "id": dmr_id,
                "chrom": chrom,
                "start": start,
                "end": end,
                "direction": direction,
                "stratum": f"{relationship_group(row)}|{direction}",
                "success": {
                    name: has_overlap(index, chrom, start, end)
                    for name, index in indexes.items()
                },
            }
        )
    return records


def candidate_sets(matrix_rows: list[dict[str, str]]) -> dict[str, dict[str, set[str]]]:
    def symbols(rows: list[dict[str, str]]) -> dict[str, set[str]]:
        return {
            "human": {row["gene_symbol"].upper() for row in rows},
            "mouse": {row["mouse_symbol"].upper() for row in rows if row.get("mouse_symbol")},
        }

    return {
        "primary_5": symbols([row for row in matrix_rows if row["tier"] == "1_primary"]),
        "final_24": symbols(matrix_rows),
        "methylation_supported_final": symbols(
            [row for row in matrix_rows if row["tier"] in {"1_primary", "2_secondary_balanced", "3_methylation_heavy"}]
        ),
        "expression_pathway_final": symbols([row for row in matrix_rows if row["tier"] == "4_expression_pathway"]),
    }


def stratified_null_counts(
    background_by_stratum: dict[str, np.ndarray],
    candidate_counts_by_stratum: dict[str, int],
    rng: np.random.Generator,
    n_permutations: int,
) -> np.ndarray:
    null = np.zeros(n_permutations, dtype=float)
    for stratum, n in candidate_counts_by_stratum.items():
        if n <= 0:
            continue
        values = background_by_stratum.get(stratum)
        if values is None or values.size == 0:
            continue
        replace = values.size < n
        for idx in range(n_permutations):
            sampled = rng.choice(values, size=n, replace=replace)
            null[idx] += float(sampled.sum())
    return null


def run_one_test(
    candidate_set: str,
    assay: str,
    species: str,
    background_name: str,
    records: list[dict[str, object]],
    candidate_symbols: set[str],
    strict_only: bool,
    rng: np.random.Generator,
    n_permutations: int,
) -> dict[str, object]:
    background = [
        record for record in records
        if (not strict_only or bool(record.get("strict")))
    ]
    candidates = [
        record for record in background
        if str(record["gene_symbol"]).upper() in candidate_symbols
    ]
    if not candidates:
        return {
            "candidate_set": candidate_set,
            "assay": assay,
            "species": species,
            "background": background_name,
            "match_strategy": "relationship_direction_stratified_link_level",
            "candidate_records": 0,
            "n_permutations": n_permutations,
        }

    observed_successes = sum(1 for record in candidates if bool(record["success"].get(assay)))
    candidate_counts_by_stratum: dict[str, int] = defaultdict(int)
    for record in candidates:
        candidate_counts_by_stratum[str(record["stratum"])] += 1

    background_by_stratum: dict[str, list[int]] = defaultdict(list)
    for record in background:
        background_by_stratum[str(record["stratum"])].append(1 if bool(record["success"].get(assay)) else 0)
    background_arrays = {key: np.array(values, dtype=int) for key, values in background_by_stratum.items()}

    null = stratified_null_counts(background_arrays, candidate_counts_by_stratum, rng, n_permutations)
    null_mean = float(null.mean())
    null_sd = float(null.std(ddof=1)) if n_permutations > 1 else math.nan
    observed_fraction = observed_successes / len(candidates)
    null_fraction = null_mean / len(candidates)
    fold = observed_fraction / null_fraction if null_fraction > 0 else math.nan
    p_greater = (float(np.sum(null >= observed_successes)) + 1.0) / (n_permutations + 1.0)
    p_less = (float(np.sum(null <= observed_successes)) + 1.0) / (n_permutations + 1.0)
    z = (observed_successes - null_mean) / null_sd if null_sd and math.isfinite(null_sd) and null_sd > 0 else math.nan
    success_genes = sorted({str(record["gene_symbol"]) for record in candidates if bool(record["success"].get(assay))})
    stratum_summary = ";".join(f"{stratum}:{count}" for stratum, count in sorted(candidate_counts_by_stratum.items()))
    return {
        "candidate_set": candidate_set,
        "assay": assay,
        "species": species,
        "background": background_name,
        "match_strategy": "relationship_direction_stratified_link_level",
        "candidate_records": len(candidates),
        "candidate_successes": observed_successes,
        "candidate_success_fraction": observed_fraction,
        "null_mean_successes": null_mean,
        "null_sd_successes": null_sd,
        "null_mean_fraction": null_fraction,
        "fold_over_null": fold,
        "empirical_p_greater": p_greater,
        "empirical_p_less": p_less,
        "z_score": z,
        "n_permutations": n_permutations,
        "candidate_success_genes": ",".join(success_genes),
        "candidate_records_by_stratum": stratum_summary,
    }


def run_tests(
    matrix_rows: list[dict[str, str]],
    human_records: list[dict[str, object]],
    mouse_records: list[dict[str, object]],
    n_permutations: int,
    seed: int,
) -> list[dict[str, object]]:
    rng = np.random.default_rng(seed)
    sets = candidate_sets(matrix_rows)
    rows: list[dict[str, object]] = []
    for set_name, set_symbols in sets.items():
        for assay in HUMAN_TRACKS:
            rows.append(
                run_one_test(
                    candidate_set=set_name,
                    assay=assay,
                    species="human",
                    background_name="human_dnmt3l_probe_gene_links",
                    records=human_records,
                    candidate_symbols=set_symbols["human"],
                    strict_only=False,
                    rng=rng,
                    n_permutations=n_permutations,
                )
            )
        rows.append(
            run_one_test(
                candidate_set=set_name,
                assay="human_strict_any_regulatory_probe",
                species="human",
                background_name="human_strict_triple_probe_gene_links",
                records=[
                    {**record, "success": {"human_strict_any_regulatory_probe": any(record["success"].values())}}
                    for record in human_records
                ],
                candidate_symbols=set_symbols["human"],
                strict_only=True,
                rng=rng,
                n_permutations=n_permutations,
            )
        )
        for assay in MOUSE_TRACKS:
            rows.append(
                run_one_test(
                    candidate_set=set_name,
                    assay=assay,
                    species="mouse",
                    background_name="mouse_triple_hit_dmr_gene_links",
                    records=mouse_records,
                    candidate_symbols=set_symbols["mouse"],
                    strict_only=False,
                    rng=rng,
                    n_permutations=n_permutations,
                )
            )
    return rows


def plot_results(rows: list[dict[str, object]], figdir: Path) -> None:
    figdir.mkdir(parents=True, exist_ok=True)
    assays = [
        "human_strict_any_regulatory_probe",
        "human_atac_or_ocr_probe",
        "human_h3k27ac_probe",
        "human_h3k27me3_probe",
        "mouse_any_regulatory_dmr",
        "mouse_atac_dmr",
        "mouse_h3k27ac_dmr",
        "mouse_h3k27me3_dmr",
    ]
    set_names = ["primary_5", "final_24", "methylation_supported_final", "expression_pathway_final"]
    by_key = {(row["candidate_set"], row["assay"]): row for row in rows}
    values = np.array(
        [[safe_float(by_key.get((set_name, assay), {}).get("fold_over_null")) for assay in assays] for set_name in set_names],
        dtype=float,
    )
    values = np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)
    x = np.arange(len(assays))
    width = 0.19
    colors = ["#2f6f8f", "#6f9e5a", "#b68a3a", "#9a5c83"]
    fig, ax = plt.subplots(figsize=(12.2, 5.0))
    for idx, set_name in enumerate(set_names):
        ax.bar(x + (idx - 1.5) * width, values[idx], width=width, color=colors[idx], label=set_name)
    ax.axhline(1, color="#333333", linestyle="--", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([assay.replace("_", "\n") for assay in assays], fontsize=7.5)
    ax.set_ylabel("Observed / matched-null mean")
    ax.set_title("Matched-Region Regulatory Overlap Permutation Tests", fontsize=12, pad=10)
    ax.legend(frameon=False, fontsize=8, ncols=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#dddddd", linewidth=0.6)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(figdir / "matched_regulatory_permutation_fold_over_null.png", dpi=300)
    fig.savefig(figdir / "matched_regulatory_permutation_fold_over_null.pdf")
    plt.close(fig)


def write_note(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    key_rows = [
        row for row in rows
        if row["candidate_set"] in {"primary_5", "final_24", "methylation_supported_final", "expression_pathway_final"}
        and row["assay"] in {
            "human_strict_any_regulatory_probe",
            "human_atac_or_ocr_probe",
            "human_h3k27ac_probe",
            "human_h3k27me3_probe",
            "mouse_any_regulatory_dmr",
            "mouse_atac_dmr",
            "mouse_h3k27ac_dmr",
            "mouse_h3k27me3_dmr",
        }
    ]
    lines = [
        "# Matched regulatory permutation tests",
        "",
        "This pass uses relationship/direction-stratified link-level permutation tests. Human CpG-gene links are matched by gene relationship and DNMT3L methylation direction. Mouse DMR-gene links are matched by gene relationship and DMR direction. The statistic is the number of candidate-linked methylation records overlapping a fetal-brain regulatory annotation.",
        "",
        "These are stricter than the gene-level hypergeometric tests, but still first-pass: they do not match CpG density, CpG-island state, gene length, local mappability, or chromosomal neighborhood.",
        "",
        "## Key permutation results",
        "",
        "| Set | Assay | Records | Observed successes | Null mean | Fold over null | Empirical p greater |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in key_rows:
        if int(float(row.get("candidate_records") or 0)) == 0:
            continue
        lines.append(
            f"| {row['candidate_set']} | {row['assay']} | {row['candidate_records']} | "
            f"{row['candidate_successes']} | {float(row['null_mean_successes']):.2f} | "
            f"{float(row['fold_over_null']):.3f} | {float(row['empirical_p_greater']):.3g} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The matched tests ask a different question from the locus panels. The locus panels show that the selected primary genes have regulatory-context support. The matched tests ask whether their regulatory overlap is higher than matched DNMT3L-sensitive background intervals.",
            "",
            "For manuscript wording, use these tests conservatively: they support regulatory annotation of the candidate loci, but strong enrichment claims should be reserved for assays/sets with empirical p values below the chosen threshold or followed up with a stricter region-matched permutation design.",
            "",
            "## Outputs",
            "",
            "- `results/mechanism_classification/matched_regulatory_permutation_tests.tsv`",
            "- `figures/mechanism_classification/matched_regulatory_permutation_fold_over_null.png`",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    matrix_rows = read_tsv(Path(args.matrix))
    print("Loading human regulatory annotations...")
    human_indexes = load_indexes(HUMAN_TRACKS)
    print("Loading mouse regulatory annotations...")
    mouse_indexes = load_indexes(MOUSE_TRACKS)
    print("Building human probe-gene records...")
    human_records = build_human_records(Path(args.human_links), human_indexes)
    print("Building mouse DMR-gene records...")
    mouse_records = build_mouse_records(Path(args.mouse_links), mouse_indexes)
    print(f"Running {args.permutations} permutations per test...")
    rows = run_tests(matrix_rows, human_records, mouse_records, args.permutations, args.seed)

    outdir = Path(args.outdir)
    write_tsv(outdir / "matched_regulatory_permutation_tests.tsv", rows, OUTPUT_FIELDS)
    plot_results(rows, Path(args.figdir))
    write_note(Path(args.note), rows)
    print(f"Wrote {len(rows)} matched permutation tests")
    print(f"Wrote {outdir / 'matched_regulatory_permutation_tests.tsv'}")
    print(f"Wrote {args.note}")


if __name__ == "__main__":
    main()
