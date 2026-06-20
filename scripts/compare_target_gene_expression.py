#!/usr/bin/env python3
"""Compare target-gene expression in the direct DNMT3L perturbation RNA-seq datasets."""

from __future__ import annotations

import argparse
import csv
import gzip
import itertools
import math
import re
from pathlib import Path


DATASETS = {
    "GSE157080": {
        "species": "human",
        "matrix": "data/geo/GSE157080/supplementary/GSE157080_all.fpkm_anno.txt.gz",
        "sample_columns": [
            {"sample": "Z1", "column": "Z1_FPKM", "condition": "control", "sex": "female"},
            {"sample": "Z2", "column": "Z2_FPKM", "condition": "control", "sex": "female"},
            {"sample": "Z3", "column": "Z3_FPKM", "condition": "control", "sex": "female"},
            {"sample": "D1", "column": "D1_FPKM", "condition": "DNMT3L_overexpression", "sex": "female"},
            {"sample": "D2", "column": "D2_FPKM", "condition": "DNMT3L_overexpression", "sex": "female"},
            {"sample": "D3", "column": "D3_FPKM", "condition": "DNMT3L_overexpression", "sex": "female"},
        ],
        "comparisons": [
            {
                "comparison": "human_all_control_vs_DNMT3L_overexpression",
                "control": ["Z1", "Z2", "Z3"],
                "case": ["D1", "D2", "D3"],
            },
        ],
    },
    "GSE268951": {
        "species": "mouse",
        "matrix": "data/geo/GSE268951/supplementary/GSE268951_all.fpkm_mRNA_anno.txt.gz",
        "sample_columns": [
            {"sample": "A-1", "column": "A-1_FPKM", "condition": "control", "sex": "male"},
            {"sample": "A-2", "column": "A-2_FPKM", "condition": "control", "sex": "male"},
            {"sample": "A-3", "column": "A-3_FPKM", "condition": "control", "sex": "male"},
            {"sample": "B-1", "column": "B-1_FPKM", "condition": "control", "sex": "female"},
            {"sample": "B-2", "column": "B-2_FPKM", "condition": "control", "sex": "female"},
            {"sample": "B-3", "column": "B-3_FPKM", "condition": "control", "sex": "female"},
            {"sample": "C-1", "column": "C-1_FPKM", "condition": "Dnmt3l_overexpression", "sex": "male"},
            {"sample": "C-2", "column": "C-2_FPKM", "condition": "Dnmt3l_overexpression", "sex": "male"},
            {"sample": "C-3", "column": "C-3_FPKM", "condition": "Dnmt3l_overexpression", "sex": "male"},
            {"sample": "D-1", "column": "D-1_FPKM", "condition": "Dnmt3l_overexpression", "sex": "female"},
            {"sample": "D-2", "column": "D-2_FPKM", "condition": "Dnmt3l_overexpression", "sex": "female"},
            {"sample": "D-3", "column": "D-3_FPKM", "condition": "Dnmt3l_overexpression", "sex": "female"},
        ],
        "comparisons": [
            {
                "comparison": "mouse_all_control_vs_Dnmt3l_overexpression",
                "control": ["A-1", "A-2", "A-3", "B-1", "B-2", "B-3"],
                "case": ["C-1", "C-2", "C-3", "D-1", "D-2", "D-3"],
            },
            {
                "comparison": "mouse_male_control_vs_Dnmt3l_overexpression",
                "control": ["A-1", "A-2", "A-3"],
                "case": ["C-1", "C-2", "C-3"],
            },
            {
                "comparison": "mouse_female_control_vs_Dnmt3l_overexpression",
                "control": ["B-1", "B-2", "B-3"],
                "case": ["D-1", "D-2", "D-3"],
            },
        ],
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", default="config/target_genes.tsv", help="Target gene TSV.")
    parser.add_argument("--outdir", default="results/first_pass_expression", help="Output directory.")
    parser.add_argument("--pseudocount", type=float, default=0.1, help="Pseudocount for log2 expression.")
    parser.add_argument("--datasets", nargs="*", default=sorted(DATASETS), help="Datasets to analyze.")
    return parser.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="", errors="replace")
    return path.open(newline="", errors="replace")


def normalize_token(value: str) -> str:
    value = value.strip().strip('"').strip("'")
    value = re.sub(r"\.\d+$", "", value)
    return value.upper()


def split_aliases(value: str) -> list[str]:
    return [token for token in re.split(r"[|;,]\s*|\s+", value.strip()) if token]


def read_targets(path: Path) -> tuple[dict[str, list[dict[str, str]]], list[dict[str, str]]]:
    alias_map: dict[str, list[dict[str, str]]] = {}
    targets: list[dict[str, str]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            targets.append(row)
            aliases = [
                row["gene_symbol"],
                *split_aliases(row.get("aliases", "")),
                *split_aliases(row.get("ensembl_hint", "")),
            ]
            for alias in aliases:
                alias_map.setdefault(normalize_token(alias), []).append(row)
    return alias_map, targets


def row_tokens(row: dict[str, str]) -> set[str]:
    tokens = set()
    for key in ("gene_id", "GeneID", "GeneSymbol", "Description"):
        value = row.get(key) or ""
        for token in re.split(r"[|;,()\[\]\s]+", value):
            if token:
                tokens.add(normalize_token(token))
    return tokens


def find_matching_target(row: dict[str, str], alias_map: dict[str, list[dict[str, str]]], species: str) -> list[str]:
    matches = []
    seen = set()
    for token in row_tokens(row):
        for target in alias_map.get(token, []):
            if target["species"] != species:
                continue
            gene = target["gene_symbol"]
            if gene not in seen:
                seen.add(gene)
                matches.append(gene)
    return matches


def read_target_rows(
    matrix: Path,
    species: str,
    alias_map: dict[str, list[dict[str, str]]],
    sample_columns: list[dict[str, str]],
) -> dict[str, dict[str, str]]:
    sample_column_names = [sample["column"] for sample in sample_columns]
    candidates: dict[str, list[dict[str, str]]] = {}
    with open_text(matrix) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row_number, row in enumerate(reader, start=2):
            matches = find_matching_target(row, alias_map, species)
            if not matches:
                continue
            row["_row_number"] = str(row_number)
            for gene in matches:
                candidates.setdefault(gene, []).append(row)

    selected = {}
    for gene, rows in candidates.items():
        selected[gene] = max(rows, key=lambda row: mean([safe_float(row[col]) for col in sample_column_names]))
    return selected


def safe_float(value: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return number if math.isfinite(number) else float("nan")


def mean(values: list[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    return sum(finite) / len(finite) if finite else float("nan")


def fmt(value: float) -> str:
    if not math.isfinite(value):
        return ""
    return f"{value:.6g}"


def log2_values(values: list[float], pseudocount: float) -> list[float]:
    return [math.log2(max(value, 0.0) + pseudocount) for value in values if math.isfinite(value)]


def exact_permutation_pvalue(control: list[float], case: list[float]) -> float:
    observed = mean(case) - mean(control)
    pooled = control + case
    case_size = len(case)
    if len(control) == 0 or len(case) == 0:
        return float("nan")

    total = 0
    extreme = 0
    all_indices = range(len(pooled))
    for case_indices in itertools.combinations(all_indices, case_size):
        case_set = set(case_indices)
        perm_case = [pooled[idx] for idx in all_indices if idx in case_set]
        perm_control = [pooled[idx] for idx in all_indices if idx not in case_set]
        diff = mean(perm_case) - mean(perm_control)
        total += 1
        if abs(diff) >= abs(observed) - 1e-12:
            extreme += 1
    return extreme / total if total else float("nan")


def direction(delta: float, threshold: float = 0.1) -> str:
    if not math.isfinite(delta) or abs(delta) < threshold:
        return "flat"
    return "up" if delta > 0 else "down"


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def analyze_dataset(
    dataset: str,
    config: dict,
    alias_map: dict[str, list[dict[str, str]]],
    pseudocount: float,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    matrix = Path(config["matrix"])
    sample_columns = config["sample_columns"]
    sample_meta = {sample["sample"]: sample for sample in sample_columns}
    selected_rows = read_target_rows(matrix, config["species"], alias_map, sample_columns)

    sample_value_rows = []
    comparison_rows = []
    for gene, row in selected_rows.items():
        values_by_sample = {}
        for sample in sample_columns:
            fpkm = safe_float(row[sample["column"]])
            values_by_sample[sample["sample"]] = fpkm
            sample_value_rows.append(
                {
                    "dataset": dataset,
                    "species": config["species"],
                    "gene_symbol": gene,
                    "row_number": row.get("_row_number", ""),
                    "gene_id": row.get("gene_id", ""),
                    "GeneID": row.get("GeneID", ""),
                    "GeneSymbol": row.get("GeneSymbol", ""),
                    "sample": sample["sample"],
                    "condition": sample["condition"],
                    "sex": sample["sex"],
                    "fpkm": fmt(fpkm),
                    "log2_fpkm_plus_pseudocount": fmt(math.log2(max(fpkm, 0.0) + pseudocount)),
                }
            )

        for comparison in config["comparisons"]:
            control_fpkm = [values_by_sample[sample] for sample in comparison["control"]]
            case_fpkm = [values_by_sample[sample] for sample in comparison["case"]]
            control_log2 = log2_values(control_fpkm, pseudocount)
            case_log2 = log2_values(case_fpkm, pseudocount)
            mean_control = mean(control_fpkm)
            mean_case = mean(case_fpkm)
            mean_control_log2 = mean(control_log2)
            mean_case_log2 = mean(case_log2)
            delta_log2 = mean_case_log2 - mean_control_log2
            log2fc_mean = math.log2((mean_case + pseudocount) / (mean_control + pseudocount))
            pvalue = exact_permutation_pvalue(control_log2, case_log2)
            comparison_rows.append(
                {
                    "dataset": dataset,
                    "comparison": comparison["comparison"],
                    "species": config["species"],
                    "gene_symbol": gene,
                    "row_number": row.get("_row_number", ""),
                    "gene_id": row.get("gene_id", ""),
                    "GeneID": row.get("GeneID", ""),
                    "GeneSymbol": row.get("GeneSymbol", ""),
                    "n_control": str(len(control_fpkm)),
                    "n_case": str(len(case_fpkm)),
                    "mean_control_fpkm": fmt(mean_control),
                    "mean_case_fpkm": fmt(mean_case),
                    "log2fc_mean_fpkm": fmt(log2fc_mean),
                    "mean_log2_control": fmt(mean_control_log2),
                    "mean_log2_case": fmt(mean_case_log2),
                    "delta_mean_log2_case_minus_control": fmt(delta_log2),
                    "direction": direction(delta_log2),
                    "exact_label_permutation_p": fmt(pvalue),
                    "control_samples": ",".join(comparison["control"]),
                    "case_samples": ",".join(comparison["case"]),
                    "control_fpkm": ",".join(fmt(value) for value in control_fpkm),
                    "case_fpkm": ",".join(fmt(value) for value in case_fpkm),
                }
            )

    return comparison_rows, sample_value_rows


def main() -> int:
    args = parse_args()
    alias_map, _targets = read_targets(Path(args.targets))

    all_comparisons = []
    all_sample_values = []
    for dataset in args.datasets:
        if dataset not in DATASETS:
            raise SystemExit(f"Unknown dataset: {dataset}")
        comparison_rows, sample_value_rows = analyze_dataset(
            dataset,
            DATASETS[dataset],
            alias_map,
            pseudocount=args.pseudocount,
        )
        all_comparisons.extend(comparison_rows)
        all_sample_values.extend(sample_value_rows)

    outdir = Path(args.outdir)
    write_tsv(
        outdir / "target_gene_group_comparison.tsv",
        sorted(all_comparisons, key=lambda row: (row["dataset"], row["gene_symbol"], row["comparison"])),
        [
            "dataset",
            "comparison",
            "species",
            "gene_symbol",
            "row_number",
            "gene_id",
            "GeneID",
            "GeneSymbol",
            "n_control",
            "n_case",
            "mean_control_fpkm",
            "mean_case_fpkm",
            "log2fc_mean_fpkm",
            "mean_log2_control",
            "mean_log2_case",
            "delta_mean_log2_case_minus_control",
            "direction",
            "exact_label_permutation_p",
            "control_samples",
            "case_samples",
            "control_fpkm",
            "case_fpkm",
        ],
    )
    write_tsv(
        outdir / "target_gene_sample_values.tsv",
        sorted(all_sample_values, key=lambda row: (row["dataset"], row["gene_symbol"], row["sample"])),
        [
            "dataset",
            "species",
            "gene_symbol",
            "row_number",
            "gene_id",
            "GeneID",
            "GeneSymbol",
            "sample",
            "condition",
            "sex",
            "fpkm",
            "log2_fpkm_plus_pseudocount",
        ],
    )
    print(f"Wrote {len(all_comparisons)} comparisons and {len(all_sample_values)} sample values to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
