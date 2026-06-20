#!/usr/bin/env python3
"""Separate 5mC and 5hmC components in GSE74486 OXBS cerebellum data."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path

from scipy import stats


CONTROL_SAMPLES = ["tb8486", "tb8487", "tb8488"]
DS_SAMPLES = ["tb8489", "tb8490", "tb8491"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--oxbs-matrix", default="data/geo/GSE74486/supplementary/GSE74486_DS_OXBS_processed_matrix_3.txt.gz")
    parser.add_argument("--dnmt3l-dmcs", default="results/human_candidate_dmrs/GSE157078_candidate_dmcs.tsv")
    parser.add_argument("--probe-support", default="results/human_methylation_convergence/human_methylation_probe_support.tsv")
    parser.add_argument("--probe-gene-links", default="results/human_methylation_convergence/human_methylation_probe_gene_links.tsv")
    parser.add_argument("--final-candidates", default="results/manuscript_summary/final_candidate_manuscript_summary.tsv")
    parser.add_argument("--outdir", default="results/mechanism_question")
    parser.add_argument("--note", default="notes/gse74486_oxbs_5mc_5hmc_first_pass.md")
    parser.add_argument("--detection-p", type=float, default=0.005)
    parser.add_argument("--min-group-samples", type=int, default=2)
    return parser.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="", errors="replace")
    return path.open(newline="", errors="replace")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", errors="replace") as handle:
        return [{k: (v or "").strip() for k, v in row.items()} for row in csv.DictReader(handle, delimiter="\t")]


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})


def safe_float(value: object) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return math.nan
    return number if math.isfinite(number) else math.nan


def mean(values: list[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    return sum(finite) / len(finite) if finite else math.nan


def median(values: list[float]) -> float:
    finite = sorted(value for value in values if math.isfinite(value))
    if not finite:
        return math.nan
    midpoint = len(finite) // 2
    if len(finite) % 2:
        return finite[midpoint]
    return (finite[midpoint - 1] + finite[midpoint]) / 2


def fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}" if math.isfinite(value) else ""
    return str(value)


def yes_no(value: bool) -> str:
    return "yes" if value else "no"


def safe_ttest(control: list[float], case: list[float]) -> float:
    control = [value for value in control if math.isfinite(value)]
    case = [value for value in case if math.isfinite(value)]
    if len(control) < 2 or len(case) < 2:
        return math.nan
    _stat, pvalue = stats.ttest_ind(case, control, equal_var=False)
    return pvalue if math.isfinite(pvalue) else math.nan


def read_dnmt3l_signature(path: Path) -> dict[str, dict[str, str]]:
    signature: dict[str, dict[str, str]] = {}
    with path.open(newline="", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            signature[row["probe_id"]] = row
    return signature


def read_probe_support(path: Path) -> dict[str, dict[str, str]]:
    return {row["probe_id"]: row for row in read_tsv(path)}


def read_final_candidates(path: Path) -> set[str]:
    return {row["gene_symbol"].upper() for row in read_tsv(path) if row.get("gene_symbol")}


def read_probe_gene_links(path: Path, final_genes: set[str]) -> dict[str, set[str]]:
    links: dict[str, set[str]] = defaultdict(set)
    for row in read_tsv(path):
        gene = row.get("gene_symbol", "").upper()
        probe = row.get("probe_id", "")
        if gene in final_genes and probe:
            links[gene].add(probe)
    return links


def read_oxbs_probe_values(args: argparse.Namespace, signature: dict[str, dict[str, str]], support: dict[str, dict[str, str]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with open_text(Path(args.oxbs_matrix)) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            probe = row.get("ID_REF", "")
            sig = signature.get(probe)
            if sig is None:
                continue
            control_total, control_mc, control_hmc, control_hmc_clipped = sample_components(row, CONTROL_SAMPLES, args.detection_p)
            ds_total, ds_mc, ds_hmc, ds_hmc_clipped = sample_components(row, DS_SAMPLES, args.detection_p)
            if len(control_total) < args.min_group_samples or len(ds_total) < args.min_group_samples:
                continue
            delta_total = mean(ds_total) - mean(control_total)
            delta_mc = mean(ds_mc) - mean(control_mc)
            delta_hmc = mean(ds_hmc) - mean(control_hmc)
            delta_hmc_clipped = mean(ds_hmc_clipped) - mean(control_hmc_clipped)
            expected = 1 if sig["direction"] == "hyper" else -1
            probe_support = support.get(probe, {})
            rows.append(
                {
                    "probe_id": probe,
                    "dnmt3l_direction": sig["direction"],
                    "dnmt3l_delta_beta": safe_float(sig["delta_beta"]),
                    "bs_total_control_mean": mean(control_total),
                    "bs_total_ds_mean": mean(ds_total),
                    "bs_total_ds_delta": delta_total,
                    "oxbs_5mc_control_mean": mean(control_mc),
                    "oxbs_5mc_ds_mean": mean(ds_mc),
                    "oxbs_5mc_ds_delta": delta_mc,
                    "inferred_5hmc_control_mean": mean(control_hmc),
                    "inferred_5hmc_ds_mean": mean(ds_hmc),
                    "inferred_5hmc_ds_delta": delta_hmc,
                    "inferred_5hmc_clipped_control_mean": mean(control_hmc_clipped),
                    "inferred_5hmc_clipped_ds_mean": mean(ds_hmc_clipped),
                    "inferred_5hmc_clipped_ds_delta": delta_hmc_clipped,
                    "signed_total_delta_by_dnmt3l": delta_total * expected,
                    "signed_5mc_delta_by_dnmt3l": delta_mc * expected,
                    "signed_5hmc_delta_by_dnmt3l": delta_hmc * expected,
                    "signed_5hmc_clipped_delta_by_dnmt3l": delta_hmc_clipped * expected,
                    "total_concordant": yes_no(delta_total * expected > 0),
                    "five_mc_concordant": yes_no(delta_mc * expected > 0),
                    "five_hmc_concordant": yes_no(delta_hmc * expected > 0),
                    "dominant_component": dominant_component(delta_mc, delta_hmc),
                    "control_n": len(control_total),
                    "ds_n": len(ds_total),
                    "bs_total_welch_pvalue": safe_ttest(control_total, ds_total),
                    "oxbs_5mc_welch_pvalue": safe_ttest(control_mc, ds_mc),
                    "inferred_5hmc_welch_pvalue": safe_ttest(control_hmc, ds_hmc),
                    "human_triple_phase3_fetal_mean_positive": probe_support.get("human_triple_phase3_fetal_mean_positive", ""),
                    "human_triple_phase3_any_ds": probe_support.get("human_triple_phase3_any_ds", ""),
                    "ds_mean_signed_delta_prior": probe_support.get("ds_mean_signed_delta", ""),
                    "ucsc_refgene_name": sig.get("ucsc_refgene_name", ""),
                    "ucsc_refgene_group": sig.get("ucsc_refgene_group", ""),
                    "relation_to_cpg_island": sig.get("relation_to_cpg_island", ""),
                }
            )
    return rows


def sample_components(row: dict[str, str], samples: list[str], detection_p: float) -> tuple[list[float], list[float], list[float], list[float]]:
    total_values: list[float] = []
    mc_values: list[float] = []
    hmc_values: list[float] = []
    hmc_clipped_values: list[float] = []
    for sample in samples:
        bs = safe_float(row.get(f"{sample}_BS"))
        oxbs = safe_float(row.get(f"{sample}_OXBS"))
        bs_p = safe_float(row.get(f"{sample}_BS_DetectionPval"))
        oxbs_p = safe_float(row.get(f"{sample}_OXBS_DetectionPval"))
        if not all(math.isfinite(value) for value in (bs, oxbs, bs_p, oxbs_p)):
            continue
        if bs_p > detection_p or oxbs_p > detection_p:
            continue
        total_values.append(bs)
        mc_values.append(oxbs)
        hmc = bs - oxbs
        hmc_values.append(hmc)
        hmc_clipped_values.append(max(0.0, hmc))
    return total_values, mc_values, hmc_values, hmc_clipped_values


def dominant_component(delta_mc: float, delta_hmc: float) -> str:
    if not math.isfinite(delta_mc) or not math.isfinite(delta_hmc):
        return "unresolved"
    if abs(delta_mc) < 0.005 and abs(delta_hmc) < 0.005:
        return "small_delta"
    if abs(delta_mc) >= 1.5 * abs(delta_hmc):
        return "5mC_dominant"
    if abs(delta_hmc) >= 1.5 * abs(delta_mc):
        return "5hmC_dominant"
    return "mixed_5mC_5hmC"


def summarize_group(label: str, rows: list[dict[str, object]]) -> dict[str, object]:
    signed_total = [safe_float(row["signed_total_delta_by_dnmt3l"]) for row in rows]
    signed_mc = [safe_float(row["signed_5mc_delta_by_dnmt3l"]) for row in rows]
    signed_hmc = [safe_float(row["signed_5hmc_delta_by_dnmt3l"]) for row in rows]
    signed_hmc_clipped = [safe_float(row["signed_5hmc_clipped_delta_by_dnmt3l"]) for row in rows]
    return {
        "group": label,
        "n_probes": len(rows),
        "n_hyper": sum(1 for row in rows if row["dnmt3l_direction"] == "hyper"),
        "n_hypo": sum(1 for row in rows if row["dnmt3l_direction"] == "hypo"),
        "mean_signed_total_delta": mean(signed_total),
        "mean_signed_5mc_delta": mean(signed_mc),
        "mean_signed_5hmc_delta": mean(signed_hmc),
        "mean_signed_5hmc_clipped_delta": mean(signed_hmc_clipped),
        "median_signed_total_delta": median(signed_total),
        "median_signed_5mc_delta": median(signed_mc),
        "median_signed_5hmc_delta": median(signed_hmc),
        "total_concordant_fraction": fraction_positive(signed_total),
        "five_mc_concordant_fraction": fraction_positive(signed_mc),
        "five_hmc_concordant_fraction": fraction_positive(signed_hmc),
        "n_5mc_dominant": sum(1 for row in rows if row["dominant_component"] == "5mC_dominant"),
        "n_5hmc_dominant": sum(1 for row in rows if row["dominant_component"] == "5hmC_dominant"),
        "n_mixed_5mc_5hmc": sum(1 for row in rows if row["dominant_component"] == "mixed_5mC_5hmC"),
        "n_small_delta": sum(1 for row in rows if row["dominant_component"] == "small_delta"),
    }


def fraction_positive(values: list[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    return sum(1 for value in finite if value > 0) / len(finite) if finite else math.nan


def build_group_summaries(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    groups = {
        "all_dnmt3l_sensitive_overlap": rows,
        "dnmt3l_hyper": [row for row in rows if row["dnmt3l_direction"] == "hyper"],
        "dnmt3l_hypo": [row for row in rows if row["dnmt3l_direction"] == "hypo"],
        "strict_human_triple_fetal_mean_positive": [
            row for row in rows if str(row.get("human_triple_phase3_fetal_mean_positive", "")).lower() == "yes"
        ],
        "human_triple_any_ds": [row for row in rows if str(row.get("human_triple_phase3_any_ds", "")).lower() == "yes"],
    }
    return [summarize_group(label, group_rows) for label, group_rows in groups.items() if group_rows]


def build_gene_summaries(
    rows: list[dict[str, object]],
    probe_links: dict[str, set[str]],
    final_rows: list[dict[str, str]],
) -> list[dict[str, object]]:
    by_probe = {str(row["probe_id"]): row for row in rows}
    summaries: list[dict[str, object]] = []
    for final in final_rows:
        gene = final.get("gene_symbol", "").upper()
        if not gene:
            continue
        probes = probe_links.get(gene, set())
        probe_rows = [by_probe[probe] for probe in sorted(probes) if probe in by_probe]
        summary = summarize_group(gene, probe_rows)
        summary.update(
            {
                "gene_symbol": gene,
                "tier": final.get("tier", ""),
                "candidate_group": final.get("candidate_group", ""),
                "mechanism_class": final.get("mechanism_class", ""),
                "n_linked_probes_tested": len(probe_rows),
                "linked_probe_ids": ";".join(sorted({str(row["probe_id"]) for row in probe_rows})),
                "gene_level_interpretation": interpret_gene(summary, had_linked_probes=bool(probes)),
            }
        )
        summaries.append(summary)
    return summaries


def interpret_gene(summary: dict[str, object], had_linked_probes: bool = True) -> str:
    mc = safe_float(summary.get("mean_signed_5mc_delta"))
    hmc = safe_float(summary.get("mean_signed_5hmc_delta"))
    total = safe_float(summary.get("mean_signed_total_delta"))
    if not math.isfinite(total):
        if had_linked_probes:
            return "not testable in this OXBS subset after detection-p filtering"
        return "not testable in this OXBS subset because no linked candidate probe was available"
    if total > 0 and mc > 0 and abs(mc) >= abs(hmc):
        return "DS cerebellum concordance is mainly 5mC-supported"
    if total > 0 and hmc > 0 and abs(hmc) > abs(mc):
        return "DS cerebellum concordance is mainly 5hmC-supported"
    if total > 0:
        return "DS cerebellum concordance is mixed or small"
    if mc > 0 and hmc < 0:
        return "5mC component is concordant but total signal is offset by 5hmC"
    return "no DS cerebellum OXBS concordance in this subset"


def write_note(path: Path, group_summaries: list[dict[str, object]], gene_summaries: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    primary = [row for row in gene_summaries if str(row.get("tier", "")).startswith("1_")]
    with path.open("w") as handle:
        handle.write("# GSE74486 OXBS 5mC/5hmC First Pass\n\n")
        handle.write("## Question\n\n")
        handle.write("Does adult DS cerebellum methylation concordance at DNMT3L-sensitive loci reflect true 5mC, inferred 5hmC, or both?\n\n")
        handle.write("## Method\n\n")
        handle.write("The processed `GSE74486_DS_OXBS_processed_matrix_3` matrix was used. BS beta values were treated as total modified cytosine signal, OXBS beta values as 5mC, and inferred 5hmC as BS minus OXBS. Three control adult cerebellar cortex samples (`tb8486-8488`) were compared with three DS adult cerebellar cortex samples (`tb8489-8491`). Probes were retained when both BS and OXBS detection p values were <= 0.005 in at least two controls and two DS samples.\n\n")
        handle.write("## Group Summary\n\n")
        handle.write("| Group | Probes | mean signed total | mean signed 5mC | mean signed 5hmC | total concordance | 5mC concordance | 5hmC concordance |\n")
        handle.write("|---|---:|---:|---:|---:|---:|---:|---:|\n")
        for row in group_summaries:
            handle.write(
                f"| {row['group']} | {row['n_probes']} | {fmt(row['mean_signed_total_delta'])} | {fmt(row['mean_signed_5mc_delta'])} | {fmt(row['mean_signed_5hmc_delta'])} | {fmt(row['total_concordant_fraction'])} | {fmt(row['five_mc_concordant_fraction'])} | {fmt(row['five_hmc_concordant_fraction'])} |\n"
            )
        handle.write("\n## Primary Candidate Summary\n\n")
        for row in primary:
            if safe_float(row.get("n_linked_probes_tested")) == 0:
                handle.write(f"- {row['gene_symbol']}: no usable linked OXBS probes; {row['gene_level_interpretation']}.\n")
            else:
                handle.write(
                    f"- {row['gene_symbol']}: {row['n_linked_probes_tested']} linked probes; mean signed total {fmt(row['mean_signed_total_delta'])}; 5mC {fmt(row['mean_signed_5mc_delta'])}; 5hmC {fmt(row['mean_signed_5hmc_delta'])}; {row['gene_level_interpretation']}.\n"
                )
        handle.write("\n## Interpretation\n\n")
        handle.write("This OXBS subset is adult cerebellar cortex only, so it should be treated as a mature-brain context rather than fetal neurodevelopment. A positive signed 5mC delta supports true methylation concordance with the DNMT3L-sensitive direction. A positive signed 5hmC delta suggests that hydroxymethylation contributes to the apparent BS methylation signal, which is especially important in brain tissue.\n")


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    signature = read_dnmt3l_signature(Path(args.dnmt3l_dmcs))
    support = read_probe_support(Path(args.probe_support))
    final_rows = read_tsv(Path(args.final_candidates))
    final_genes = read_final_candidates(Path(args.final_candidates))
    probe_links = read_probe_gene_links(Path(args.probe_gene_links), final_genes)

    probe_rows = read_oxbs_probe_values(args, signature, support)
    probe_fields = [
        "probe_id",
        "dnmt3l_direction",
        "dnmt3l_delta_beta",
        "bs_total_control_mean",
        "bs_total_ds_mean",
        "bs_total_ds_delta",
        "oxbs_5mc_control_mean",
        "oxbs_5mc_ds_mean",
        "oxbs_5mc_ds_delta",
        "inferred_5hmc_control_mean",
        "inferred_5hmc_ds_mean",
        "inferred_5hmc_ds_delta",
        "inferred_5hmc_clipped_control_mean",
        "inferred_5hmc_clipped_ds_mean",
        "inferred_5hmc_clipped_ds_delta",
        "signed_total_delta_by_dnmt3l",
        "signed_5mc_delta_by_dnmt3l",
        "signed_5hmc_delta_by_dnmt3l",
        "signed_5hmc_clipped_delta_by_dnmt3l",
        "total_concordant",
        "five_mc_concordant",
        "five_hmc_concordant",
        "dominant_component",
        "control_n",
        "ds_n",
        "bs_total_welch_pvalue",
        "oxbs_5mc_welch_pvalue",
        "inferred_5hmc_welch_pvalue",
        "human_triple_phase3_fetal_mean_positive",
        "human_triple_phase3_any_ds",
        "ds_mean_signed_delta_prior",
        "ucsc_refgene_name",
        "ucsc_refgene_group",
        "relation_to_cpg_island",
    ]
    write_tsv(outdir / "gse74486_oxbs_5mc_5hmc_probe_values.tsv", probe_rows, probe_fields)

    group_summaries = build_group_summaries(probe_rows)
    summary_fields = [
        "group",
        "n_probes",
        "n_hyper",
        "n_hypo",
        "mean_signed_total_delta",
        "mean_signed_5mc_delta",
        "mean_signed_5hmc_delta",
        "mean_signed_5hmc_clipped_delta",
        "median_signed_total_delta",
        "median_signed_5mc_delta",
        "median_signed_5hmc_delta",
        "total_concordant_fraction",
        "five_mc_concordant_fraction",
        "five_hmc_concordant_fraction",
        "n_5mc_dominant",
        "n_5hmc_dominant",
        "n_mixed_5mc_5hmc",
        "n_small_delta",
    ]
    write_tsv(outdir / "gse74486_oxbs_5mc_5hmc_group_summary.tsv", group_summaries, summary_fields)

    gene_summaries = build_gene_summaries(probe_rows, probe_links, final_rows)
    gene_fields = [
        "gene_symbol",
        "tier",
        "candidate_group",
        "mechanism_class",
        "n_linked_probes_tested",
        "n_hyper",
        "n_hypo",
        "mean_signed_total_delta",
        "mean_signed_5mc_delta",
        "mean_signed_5hmc_delta",
        "mean_signed_5hmc_clipped_delta",
        "median_signed_total_delta",
        "median_signed_5mc_delta",
        "median_signed_5hmc_delta",
        "total_concordant_fraction",
        "five_mc_concordant_fraction",
        "five_hmc_concordant_fraction",
        "n_5mc_dominant",
        "n_5hmc_dominant",
        "n_mixed_5mc_5hmc",
        "n_small_delta",
        "linked_probe_ids",
        "gene_level_interpretation",
    ]
    write_tsv(outdir / "gse74486_oxbs_5mc_5hmc_candidate_gene_summary.tsv", gene_summaries, gene_fields)
    write_note(Path(args.note), group_summaries, gene_summaries)


if __name__ == "__main__":
    main()
