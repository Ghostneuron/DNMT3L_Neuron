"""Reproduce the revision analyses from the archived processed results.

Usage: python reanalyze.py <DNMT3L_Neuron directory> [--outdir PATH]
Requires Python 3.10 or newer. Inputs are read only.
"""

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean, median


FETAL_CONTEXTS = (
    "GSE73747_fetal_all",
    "GSE73747_fetal_frontal",
    "GSE73747_fetal_temporal",
    "GSE74486_fetus_cerebrum_whole_cerebrum",
)
COMPARISON_CONTEXTS = FETAL_CONTEXTS + (
    "GSE74486_adult_frontal_cortex_fluorescence_activated_nuclei_sorted_neurons",
    "GSE74486_adult_frontal_cortex_fluorescence_activated_nuclei_sorted_glia",
)
RANK_VARIANTS = (
    "corrected", "equal_weight", "no_expression", "no_human_methylation",
    "no_mouse_methylation", "no_dp10", "no_distal_mouse",
)


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def numeric(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def integer(value):
    number = numeric(value)
    return int(number) if math.isfinite(number) else 0


def direction(effect):
    if effect >= 0.25:
        return "up"
    if effect <= -0.25:
        return "down"
    return "flat" if math.isfinite(effect) else "unavailable"


def pairing_summary(rows):
    counts = Counter((row["dnmt3l_direction"], row["ds_direction"]) for row in rows)
    n = len(rows)
    dh = sum(counts["hyper", sign] for sign in ("up", "down", "zero"))
    dp = counts["hyper", "up"] + counts["hypo", "up"]
    dn = counts["hyper", "down"] + counts["hypo", "down"]
    observed = counts["hyper", "up"] + counts["hypo", "down"]
    expected = (dh * dp + (n - dh) * dn) / n
    return n, observed, expected


def human_concordance(root, out):
    directory = root / "results/ds_methylation_dnmt3l_concordance"
    rows = read_tsv(directory / "ds_dnmt3l_candidate_probe_values.tsv")
    by_context = defaultdict(list)
    for row in rows:
        if row["context"] not in COMPARISON_CONTEXTS:
            continue
        effect = float(row["ds_delta_beta"])
        row["ds_direction"] = "up" if effect > 0 else "down" if effect < 0 else "zero"
        row["baseline_bin"] = min(9, int(float(row["ds_control_mean_beta"]) * 10))
        by_context[row["context"]].append(row)

    summaries = []
    for context in COMPARISON_CONTEXTS:
        group = by_context[context]
        assert len({r["probe_id"] for r in group}) == len(group)
        n, concordant, expected = pairing_summary(group)
        strata = defaultdict(list)
        for row in group:
            strata[row["chrom"], row["relation_to_cpg_island"] or "unannotated", row["baseline_bin"]].append(row)
        blocks = defaultdict(list)
        for (chrom, _, _), records in strata.items():
            blocks[chrom].append(pairing_summary(records))
        block_stats = {
            chrom: tuple(sum(values) for values in zip(*data))
            for chrom, data in blocks.items()
        }
        expected_matched = sum(stat[2] for stat in block_stats.values())
        omitted = [
            ((concordant - matches) - (expected_matched - exp)) / (n - size)
            for size, matches, exp in block_stats.values()
        ]
        hyper = [r for r in group if r["dnmt3l_direction"] == "hyper"]
        hypo = [r for r in group if r["dnmt3l_direction"] == "hypo"]
        summaries.append({
            "context": context, "n_cpg": n, "n_oe_hyper": len(hyper), "n_oe_hypo": len(hypo),
            "n_ds_hyper": sum(r["ds_direction"] == "up" for r in group),
            "n_ds_hypo": sum(r["ds_direction"] == "down" for r in group),
            "observed_agreement": concordant / n, "expected_sign_margins": expected / n,
            "expected_stratified": expected_matched / n,
            "excess_stratified_pp": 100 * (concordant - expected_matched) / n,
            "omit_one_chrom_min_pp": 100 * min(omitted),
            "omit_one_chrom_max_pp": 100 * max(omitted),
            "oe_hyper_agreement": mean(r["ds_direction"] == "up" for r in hyper),
            "oe_hypo_agreement": mean(r["ds_direction"] == "down" for r in hypo),
            "mean_signed_ds_delta_beta": mean(float(r["signed_ds_delta_by_dnmt3l_direction"]) for r in group),
        })
    write_tsv(out / "human_concordance.tsv", summaries)
    return summaries


def mouse_concordance(root, out):
    source = root / "results/dp10_methylation_dnmt3l_concordance/dp10_dp16_dnmt3l_cpg_values.tsv"
    models = defaultdict(dict)
    for row in read_tsv(source):
        key = row["chrom"], int(row["start"]), int(row["end"])
        assert key not in models[row["model"]]
        models[row["model"]][key] = row
    ten, sixteen = models["Dp10"], models["Dp16"]
    common = sorted(ten.keys() & sixteen.keys())
    pairs = []
    for key in common:
        a, b = ten[key], sixteen[key]
        assert a["dnmt3l_direction"] == b["dnmt3l_direction"]
        assert a["dnmt3l_delta_beta"] == b["dnmt3l_delta_beta"]
        pairs.append({
            "chr": key[0], "start": key[1], "end": key[2],
            "perturbation_direction": a["dnmt3l_direction"],
            "dp10_delta_beta": float(a["model_delta_beta"]),
            "dp16_delta_beta": float(b["model_delta_beta"]),
            "dp10_signed": float(a["signed_model_delta_by_dnmt3l_direction"]),
            "dp16_signed": float(b["signed_model_delta_by_dnmt3l_direction"]),
        })
    summaries = []
    for stratum in ("all", "hyper", "hypo"):
        group = [p for p in pairs if stratum == "all" or p["perturbation_direction"] == stratum]
        agree_10 = mean(p["dp10_signed"] > 0 for p in group)
        agree_16 = mean(p["dp16_signed"] > 0 for p in group)
        summaries.append({
            "dnmt3l_stratum": stratum, "shared_cpg": len(group),
            "dp10_agreement": agree_10, "dp16_agreement": agree_16,
            "dp10_minus_dp16_pp": 100 * (agree_10 - agree_16),
            "dp10_mean_signed_delta": mean(p["dp10_signed"] for p in group),
            "dp16_mean_signed_delta": mean(p["dp16_signed"] for p in group),
            "dp10_only_agrees": sum(p["dp10_signed"] > 0 >= p["dp16_signed"] for p in group),
            "dp16_only_agrees": sum(p["dp16_signed"] > 0 >= p["dp10_signed"] for p in group),
        })
    write_tsv(out / "mouse_shared_cpg_summary.tsv", summaries)
    write_tsv(out / "mouse_shared_cpg_values.tsv", pairs)
    return summaries


def rank_candidates(root, out):
    def archive(pattern):
        return read_tsv(next((root / "manuscript/Supplementary Tables").glob(pattern)))

    directory = root / "results"
    human = archive("ST4c_*.tsv")
    de = {r["gene_id"]: r for r in archive("ST17_*.tsv")}
    nominations = {r["gene_symbol"]: r for r in archive("ST6b_*.tsv")}
    mouse_file = directory / "integrated_candidate_genes/mouse_triple_hit_candidate_gene_ranking.tsv"
    mouse = {}
    for row in read_tsv(mouse_file):
        mouse.setdefault(row["gene_symbol"], row)
    comparison_files = directory / "cross_species_candidate_expression"
    bulk = defaultdict(list)
    cells = defaultdict(list)
    for row in read_tsv(comparison_files / "human_bulk_candidate_effects.tsv"):
        bulk[row["human_symbol"]].append(row)
    for row in read_tsv(comparison_files / "gse305153_candidate_celltype_effects.tsv"):
        if row["label_type"] == "author_cell_type":
            cells[row["human_symbol"]].append(row)

    results = []
    for row in human:
        symbol = row["gene_symbol"]
        curator = nominations.get(symbol, {})
        mouse_symbol = row["mouse_gene_symbol"] or curator.get("mouse_symbol", "")
        m = mouse.get(mouse_symbol, {})
        gene_result = de.get(m.get("gene_id", ""), {})
        fc = numeric(gene_result.get("deseq2_log2fc"))
        q = numeric(gene_result.get("deseq2_padj"))
        ref_dir = direction(fc)
        obs = m.get("gse268951_max_group_mean_fpkm", "")
        if not obs:
            obs = max(numeric(m.get("gse268951_control_mean_fpkm")), numeric(m.get("gse268951_dnmt3l_oe_mean_fpkm")))
        mp = integer(m.get("n_promoter_2kb_dmrs"))
        mb = integer(m.get("n_gene_body_dmrs"))
        mn = integer(m.get("n_nearest_100kb_dmrs"))
        direct_mouse_m = min(9, 3 * mp) + min(4, 2 * mb)
        mouse_m = direct_mouse_m + min(2, mn)
        same_sex = m.get("gse268951_sex_consistent_direction") == "yes"
        detection = numeric(obs) >= 1
        expr = sum((detection, abs(fc) >= 0.58, abs(fc) >= 1,
                    q <= 0.05 and abs(fc) >= 1, same_sex))
        dp = m.get("dp10_expression_direction", "")
        dp_effect = numeric(m.get("dp10_delta_log2_cpm_plus_1"))
        dp_score = int(abs(dp_effect) >= 0.25 or dp in ("up", "down"))
        dp_score += int(dp in ("up", "down") and dp == ref_dir)
        has_dir = ref_dir in ("up", "down")
        matches_bulk = sum(r["direction"] == ref_dir for r in bulk[symbol]) if has_dir else 0
        matches_cells = sum(r["direction"] == ref_dir for r in cells[symbol]) if has_dir else 0
        human_e = min(2, matches_bulk) + min(4, matches_cells) + int(matches_bulk > 0 and matches_cells > 0)
        human_m = integer(row["human_methylation_support_score"])
        calculated_human_m = sum(
            min(cap, multiplier * integer(row[field]))
            for field, cap, multiplier in [
                ("n_human_triple_phase3_fetal_mean_positive_probes", 10, 3),
                ("n_human_triple_phase3_fetal_ds_majority_probes", 5, 2),
                ("n_human_triple_phase3_any_ds_probes", 3, 1),
                ("n_promoter_triple_any_ds_probes", 3, 1),
                ("n_gene_body_triple_any_ds_probes", 2, 1),
            ]
        )
        assert human_m == calculated_human_m, symbol
        corrected = human_m + mouse_m + expr + dp_score + human_e
        results.append({
            "gene_id": row["gene_id"], "gene_symbol": symbol, "prior_tier": curator.get("tier", ""),
            "provisional_mouse_symbol": mouse_symbol, "mouse_gene_id": m.get("gene_id", ""),
            "mouse_promoter_dmrs": mp, "mouse_body_dmrs": mb, "mouse_nearest_dmrs": mn,
            "mouse_gene_log2fc": fc if math.isfinite(fc) else "",
            "mouse_gene_deseq2_fdr": q if math.isfinite(q) else "",
            "mouse_expression_direction": ref_dir,
            "mouse_fpkm_at_least_1": int(detection), "mouse_fpkm_same_sex_direction": int(same_sex),
            "dp10_expression_direction": dp, "n_bulk_direction_agreements": matches_bulk,
            "n_celltype_direction_agreements": matches_cells,
            "n_bulk_contexts_available": len(bulk[symbol]),
            "n_celltype_contexts_available": len(cells[symbol]),
            "human_fetal_mean_positive_cpgs": integer(row["n_human_triple_phase3_fetal_mean_positive_probes"]),
            "human_fetal_majority_cpgs": integer(row["n_human_triple_phase3_fetal_ds_majority_probes"]),
            "human_all_triple_cpgs": integer(row["n_human_triple_phase3_any_ds_probes"]),
            "human_promoter_triple_cpgs": integer(row["n_promoter_triple_any_ds_probes"]),
            "human_body_triple_cpgs": integer(row["n_gene_body_triple_any_ds_probes"]),
            "human_methylation_score": human_m, "mouse_methylation_score": mouse_m,
            "mouse_expression_score": expr, "dp10_expression_score": dp_score,
            "human_concordant_expression_score": human_e,
            "submitted_st4c_score": integer(row["combined_methylation_expression_score"]),
            "corrected": corrected,
            "equal_weight": round(human_m / 23 + mouse_m / 15 + expr / 5 + dp_score / 2 + human_e / 7, 9),
            "no_expression": human_m + mouse_m,
            "no_human_methylation": corrected - human_m,
            "no_mouse_methylation": corrected - mouse_m,
            "no_dp10": corrected - dp_score,
            "no_distal_mouse": corrected - mouse_m + direct_mouse_m,
        })

    assert len({r["gene_id"] for r in results}) == len(results)
    missing = set(nominations) - {r["gene_symbol"] for r in results}
    top = {}
    for scenario in ("submitted_st4c_score", *RANK_VARIANTS):
        ranked = sorted(results, key=lambda r: (-r[scenario], r["gene_symbol"], r["gene_id"]))
        for index, row in enumerate(ranked, 1):
            row[f"{scenario}_rank"] = index
        top[scenario] = {row["gene_id"] for row in ranked[:24]}
    primary = {k for k, r in nominations.items() if r["tier"] == "1_primary"}
    summary = [
        {
            "scenario": scenario, "n_ranked_gene_ids": len(results),
            "top24_overlap_with_corrected": len(top["corrected"] & top[scenario]),
            "original_primary_in_top24": sum(
                r["gene_symbol"] in primary and r["gene_id"] in top[scenario] for r in results
            ),
        }
        for scenario in top
    ]
    nominated = sorted(
        (r for r in results if r["gene_symbol"] in nominations),
        key=lambda r: (r["prior_tier"], r["corrected_rank"]),
    )
    audited_nominations = [
        {**row, "discovery_universe_status": "linked and ranked"} for row in nominated
    ]
    for symbol in sorted(missing):
        row = dict.fromkeys(nominated[0], "")
        row.update(gene_symbol=symbol, prior_tier=nominations[symbol]["tier"],
                   provisional_mouse_symbol=nominations[symbol]["mouse_symbol"],
                   discovery_universe_status="outside the human probe-linked discovery universe")
        audited_nominations.append(row)
    assert len(audited_nominations) == len(nominations)
    write_tsv(out / "rank_sensitivity.tsv", summary)
    write_tsv(out / "all_discovery_ranks.tsv", sorted(results, key=lambda r: r["corrected_rank"]))
    write_tsv(out / "curated_candidate_score_audit.tsv", audited_nominations)
    return {"sensitivity": summary, "missing": sorted(missing), "primary": [
        {k: r[k] for k in ("gene_symbol", "corrected", "mouse_gene_log2fc", "mouse_gene_deseq2_fdr",
                            "human_methylation_score", "mouse_methylation_score", "mouse_expression_score",
                            "dp10_expression_score", "human_concordant_expression_score",
                            "submitted_st4c_score_rank", *(f"{s}_rank" for s in RANK_VARIANTS))}
        for r in nominated if r["gene_symbol"] in primary
    ]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--outdir", type=Path, default=Path("r1_results"))
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    results = {
        "human": human_concordance(args.project_root, args.outdir),
        "mouse": mouse_concordance(args.project_root, args.outdir),
        "mouse_unpaired": read_tsv(
            args.project_root / "results/dp10_methylation_dnmt3l_concordance/dp10_dp16_dnmt3l_cpg_concordance_summary.tsv"
        ),
        "ranks": rank_candidates(args.project_root, args.outdir),
    }
    (args.outdir / "results.json").write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    main()
