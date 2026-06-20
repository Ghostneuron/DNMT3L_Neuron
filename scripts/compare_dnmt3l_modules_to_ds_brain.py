#!/usr/bin/env python3
"""Compare Dnmt3l OE module signatures with DS brain expression datasets."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


DEFAULT_DNMT3L_EXPR = "results/neurogenesis_mechanism/gse268951_rederived_gene_expression.tsv"
DEFAULT_MEMBERSHIP = "results/neurogenesis_mechanism/developmental_module_gene_membership.tsv"
DEFAULT_DP10 = "results/gse214273_dp10_scrna/gse214273_gene_effect_summary.tsv"
DEFAULT_GSE301886_COUNTS = "data/geo/GSE301886/supplementary/GSE301886_T21_raw_counts.txt.gz"
DEFAULT_GSE301886_META = "data/geo/GSE301886/metadata/sample_metadata.tsv"
DEFAULT_GSE305139_COUNTS = "data/geo/GSE305139/supplementary/GSE305139_raw_counts_tissue_bulk.tsv.gz"
DEFAULT_GSE305139_META = "data/geo/GSE305139/metadata/sample_metadata.tsv"
DEFAULT_HUMAN_DS_BULK_COUNT_DE = "results/human_ds_bulk_count_de/human_ds_bulk_count_de_all_contexts.tsv"
DEFAULT_OUTDIR = "results/neurogenesis_mechanism"
DEFAULT_FIGDIR = "figures/neurogenesis_mechanism"
DEFAULT_NOTE = "notes/dnmt3l_module_ds_brain_concordance.md"

MODULES = [
    "glial_maturation_myelination",
    "neurodevelopment_signaling",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dnmt3l-expr", default=DEFAULT_DNMT3L_EXPR)
    parser.add_argument("--membership", default=DEFAULT_MEMBERSHIP)
    parser.add_argument("--dp10", default=DEFAULT_DP10)
    parser.add_argument("--gse301886-counts", default=DEFAULT_GSE301886_COUNTS)
    parser.add_argument("--gse301886-meta", default=DEFAULT_GSE301886_META)
    parser.add_argument("--gse305139-counts", default=DEFAULT_GSE305139_COUNTS)
    parser.add_argument("--gse305139-meta", default=DEFAULT_GSE305139_META)
    parser.add_argument("--human-ds-bulk-count-de", default=DEFAULT_HUMAN_DS_BULK_COUNT_DE)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--figdir", default=DEFAULT_FIGDIR)
    parser.add_argument("--note", default=DEFAULT_NOTE)
    parser.add_argument("--min-expression", type=float, default=0.5)
    parser.add_argument("--dnmt3l-change-threshold", type=float, default=0.25)
    parser.add_argument("--ds-change-threshold", type=float, default=0.05)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", errors="replace") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


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


def heatmap_text_color(value: float, vmin: float, vmax: float, cmap_name: str = "RdBu_r") -> str:
    if not math.isfinite(value):
        return "black"
    fraction = min(1.0, max(0.0, (value - vmin) / (vmax - vmin)))
    red, green, blue, _alpha = plt.get_cmap(cmap_name)(fraction)
    luminance = 0.2126 * red + 0.7152 * green + 0.0722 * blue
    return "white" if luminance < 0.45 else "black"


def fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}" if math.isfinite(value) else ""
    if isinstance(value, int):
        return str(value)
    return str(value)


def mean(values: list[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    return sum(finite) / len(finite) if finite else math.nan


def median(values: list[float]) -> float:
    finite = sorted(value for value in values if math.isfinite(value))
    if not finite:
        return math.nan
    mid = len(finite) // 2
    return finite[mid] if len(finite) % 2 else (finite[mid - 1] + finite[mid]) / 2


def clean_ensembl(value: str) -> str:
    return value.split(".", 1)[0]


def parse_characteristics(text: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for part in text.split(";"):
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        values[key.strip().lower()] = value.strip()
    return values


def normalize_group(value: str) -> str:
    norm = value.strip().lower()
    if norm in {"t21", "trisomy 21", "ds", "down syndrome"}:
        return "DS"
    if norm in {"control", "con", "disomic", "euploid"}:
        return "Control"
    return value.strip()


def load_dnmt3l_signature(path: Path, min_expression: float) -> dict[str, dict[str, object]]:
    rows: dict[str, dict[str, object]] = {}
    for row in read_tsv(path):
        symbol = row.get("gene_symbol", "")
        if not symbol:
            continue
        control = safe_float(row.get("gse268951_control_mean_fpkm"))
        oe = safe_float(row.get("gse268951_dnmt3l_oe_mean_fpkm"))
        max_expr = max(control, oe)
        logfc = safe_float(row.get("gse268951_log2fc_oe_vs_control"))
        if max_expr < min_expression or not math.isfinite(logfc):
            continue
        rows[symbol] = {
            "mouse_symbol": symbol,
            "human_symbol": symbol.upper(),
            "dnmt3l_log2fc": logfc,
            "dnmt3l_abs_log2fc": abs(logfc),
            "dnmt3l_direction": direction(logfc, 0.0),
        }
    return rows


def load_modules(path: Path) -> dict[str, set[str]]:
    modules: dict[str, set[str]] = defaultdict(set)
    for row in read_tsv(path):
        module = row.get("module", "")
        symbol = row.get("gene_symbol", "")
        if module and symbol:
            modules[module].add(symbol)
    return modules


def direction(value: float, threshold: float) -> str:
    if not math.isfinite(value) or abs(value) <= threshold:
        return "flat"
    return "up" if value > 0 else "down"


def sign(value: float) -> int:
    if not math.isfinite(value) or value == 0:
        return 0
    return 1 if value > 0 else -1


def load_dp10_effects(path: Path) -> dict[str, dict[str, object]]:
    effects: dict[str, dict[str, object]] = {}
    for row in read_tsv(path):
        symbol = row.get("gene_symbol", "")
        delta = safe_float(row.get("delta_log2_cpm_plus_1"))
        if symbol and math.isfinite(delta):
            effects[symbol] = {
                "dataset": "GSE214273",
                "context": "Dp10_hippocampus_pseudobulk",
                "species": "mouse",
                "gene_symbol": symbol,
                "ds_delta": delta,
            }
    return effects


def gse305139_symbol_map(counts_path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with gzip.open(counts_path, "rt", newline="", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene_id = clean_ensembl(row.get("gene_id", ""))
            gene_name = row.get("gene_name", "")
            if gene_id and gene_name:
                mapping[gene_id] = gene_name
    return mapping


def load_gse301886_metadata(path: Path) -> dict[str, dict[str, str]]:
    metadata: dict[str, dict[str, str]] = {}
    for row in read_tsv(path):
        chars = parse_characteristics(row.get("characteristics", ""))
        label = row.get("title", "").split(",")[-1].strip()
        if not label:
            continue
        metadata[label] = {
            "group": normalize_group(chars.get("karyotype", "")),
            "region": chars.get("tissue", row.get("source_name", "")),
            "age": chars.get("gestational age", ""),
        }
    return metadata


def load_gse305139_metadata(path: Path) -> dict[str, dict[str, str]]:
    metadata: dict[str, dict[str, str]] = {}
    for row in read_tsv(path):
        chars = parse_characteristics(row.get("characteristics", ""))
        label = chars.get("library igf_id", "")
        if not label:
            continue
        metadata[label] = {
            "group": normalize_group(chars.get("group", chars.get("genotype", ""))),
            "region": "fetal_brain_final",
            "age": chars.get("dev stage_pcw", ""),
            "final_analysis": chars.get("in final_analysis", ""),
        }
    return metadata


def sample_contexts(sample_labels: list[str], metadata: dict[str, dict[str, str]], dataset: str) -> dict[str, list[str]]:
    contexts: dict[str, list[str]] = defaultdict(list)
    for sample in sample_labels:
        meta = metadata.get(sample, {})
        if meta.get("group") not in {"DS", "Control"}:
            continue
        if meta.get("final_analysis", "").upper() == "FALSE":
            continue
        region = meta.get("region", "")
        if dataset == "GSE301886":
            context = region or "all"
        else:
            context = "fetal_brain_final"
        contexts[context].append(sample)
    return contexts


def count_matrix_effects(
    dataset: str,
    counts_path: Path,
    metadata_path: Path,
    id_to_symbol: dict[str, str] | None = None,
) -> dict[tuple[str, str], dict[str, dict[str, object]]]:
    if dataset == "GSE301886":
        metadata = load_gse301886_metadata(metadata_path)
    elif dataset == "GSE305139":
        metadata = load_gse305139_metadata(metadata_path)
    else:
        raise ValueError(dataset)

    with gzip.open(counts_path, "rt", newline="", errors="replace") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        has_gene_name = "gene_name" in header
        first_sample = 2 if has_gene_name else 1
        sample_labels = header[first_sample:]
        contexts = sample_contexts(sample_labels, metadata, dataset)
        library_sizes = {sample: 0.0 for sample in sample_labels}
        raw_rows: list[tuple[str, list[float]]] = []
        for row in reader:
            if not row:
                continue
            gene_id = clean_ensembl(row[0])
            symbol = row[1] if has_gene_name else (id_to_symbol or {}).get(gene_id, "")
            if not symbol:
                continue
            counts = [safe_float(value) for value in row[first_sample:]]
            counts = [0.0 if not math.isfinite(value) else value for value in counts]
            for sample, count in zip(sample_labels, counts, strict=False):
                library_sizes[sample] += count
            raw_rows.append((symbol, counts))

    effects_by_context: dict[tuple[str, str], dict[str, dict[str, object]]] = {}
    sample_index = {sample: idx for idx, sample in enumerate(sample_labels)}
    for context, samples in contexts.items():
        ds_samples = [sample for sample in samples if metadata.get(sample, {}).get("group") == "DS"]
        control_samples = [sample for sample in samples if metadata.get(sample, {}).get("group") == "Control"]
        if len(ds_samples) < 2 or len(control_samples) < 2:
            continue
        context_effects: dict[str, dict[str, object]] = {}
        for symbol, counts in raw_rows:
            values_ds = []
            values_control = []
            for sample in ds_samples:
                idx = sample_index[sample]
                cpm = counts[idx] / library_sizes[sample] * 1_000_000 if library_sizes[sample] else math.nan
                values_ds.append(math.log2(cpm + 1.0) if math.isfinite(cpm) else math.nan)
            for sample in control_samples:
                idx = sample_index[sample]
                cpm = counts[idx] / library_sizes[sample] * 1_000_000 if library_sizes[sample] else math.nan
                values_control.append(math.log2(cpm + 1.0) if math.isfinite(cpm) else math.nan)
            delta = mean(values_ds) - mean(values_control)
            if math.isfinite(delta):
                context_effects[symbol.upper()] = {
                    "dataset": dataset,
                    "context": context,
                    "species": "human",
                    "gene_symbol": symbol.upper(),
                    "ds_delta": delta,
                    "n_ds": len(ds_samples),
                    "n_control": len(control_samples),
                }
        effects_by_context[(dataset, context)] = context_effects
    return effects_by_context


def load_human_count_de_effects(path: Path) -> dict[tuple[str, str], dict[str, dict[str, object]]]:
    effects_by_context: dict[tuple[str, str], dict[str, dict[str, object]]] = defaultdict(dict)
    if not path.exists():
        return effects_by_context
    for row in read_tsv(path):
        dataset = row.get("dataset", "")
        context = row.get("context", "")
        symbol = row.get("gene_symbol", "").upper()
        effect = safe_float(row.get("effect_log2fc_for_claim"))
        if dataset not in {"GSE301886", "GSE305139"} or not context or not symbol or not math.isfinite(effect):
            continue
        max_expr = max(safe_float(row.get("mean_ds_log2_cpm_plus_1")), safe_float(row.get("mean_control_log2_cpm_plus_1")))
        current = effects_by_context[(dataset, context)].get(symbol)
        current_max = safe_float(current.get("max_mean_log2_cpm", "")) if current else -math.inf
        if current is not None and current_max >= max_expr:
            continue
        effects_by_context[(dataset, context)][symbol] = {
            "dataset": dataset,
            "context": context,
            "species": "human",
            "gene_symbol": symbol,
            "ds_delta": effect,
            "max_mean_log2_cpm": max_expr,
            "descriptive_delta_log2_cpm": safe_float(row.get("delta_ds_minus_control_log2_cpm_plus_1")),
            "effect_source": row.get("effect_source", ""),
            "edgeR_logFC": safe_float(row.get("logFC")),
            "edgeR_PValue": safe_float(row.get("PValue")),
            "edgeR_FDR": safe_float(row.get("FDR")),
            "tested_by_edgeR": row.get("tested_by_edgeR", ""),
            "n_ds": safe_float(row.get("n_ds")),
            "n_control": safe_float(row.get("n_control")),
        }
    return effects_by_context


def paired_rows_for_context(
    signature: dict[str, dict[str, object]],
    ds_effects: dict[str, dict[str, object]],
    module_genes: set[str],
    species: str,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    paired: list[dict[str, object]] = []
    background: list[dict[str, object]] = []
    for mouse_symbol, sig in signature.items():
        key = mouse_symbol if species == "mouse" else str(sig["human_symbol"])
        ds = ds_effects.get(key)
        if ds is None:
            continue
        dnmt3l = safe_float(sig["dnmt3l_log2fc"])
        ds_delta = safe_float(ds["ds_delta"])
        if not math.isfinite(dnmt3l) or not math.isfinite(ds_delta):
            continue
        row = {
            "mouse_symbol": mouse_symbol,
            "matched_symbol": key,
            "dnmt3l_log2fc": dnmt3l,
            "ds_delta": ds_delta,
            "signed_ds_delta_by_dnmt3l_direction": ds_delta * sign(dnmt3l),
            "same_direction": "yes" if sign(dnmt3l) == sign(ds_delta) and sign(dnmt3l) != 0 else "no",
            "in_module": "yes" if mouse_symbol in module_genes else "no",
            "ds_effect_source": ds.get("effect_source", ""),
            "ds_edgeR_FDR": ds.get("edgeR_FDR", ""),
            "ds_tested_by_edgeR": ds.get("tested_by_edgeR", ""),
        }
        if mouse_symbol in module_genes:
            paired.append(row)
        else:
            background.append(row)
    return paired, background


def summarize_match(
    dataset: str,
    context: str,
    species: str,
    module: str,
    paired: list[dict[str, object]],
    background: list[dict[str, object]],
    dnmt3l_threshold: float,
    ds_threshold: float,
) -> dict[str, object]:
    n = len(paired)
    bg_n = len(background)
    dnmt3l = [safe_float(row["dnmt3l_log2fc"]) for row in paired]
    ds = [safe_float(row["ds_delta"]) for row in paired]
    rho, rho_p = (math.nan, math.nan)
    if n >= 5 and len(set(round(x, 12) for x in dnmt3l)) > 1 and len(set(round(y, 12) for y in ds)) > 1:
        rho, rho_p = stats.spearmanr(dnmt3l, ds)
    score = [safe_float(row["signed_ds_delta_by_dnmt3l_direction"]) for row in paired]
    bg_score = [safe_float(row["signed_ds_delta_by_dnmt3l_direction"]) for row in background]
    mw_p = math.nan
    if len(score) >= 3 and len(bg_score) >= 10:
        mw_p = float(stats.mannwhitneyu(score, bg_score, alternative="greater").pvalue)

    changed = [
        row
        for row in paired
        if abs(safe_float(row["dnmt3l_log2fc"])) >= dnmt3l_threshold
        and abs(safe_float(row["ds_delta"])) >= ds_threshold
    ]
    same_changed = [row for row in changed if row["same_direction"] == "yes"]
    background_changed = [
        row
        for row in background
        if abs(safe_float(row["dnmt3l_log2fc"])) >= dnmt3l_threshold
        and abs(safe_float(row["ds_delta"])) >= ds_threshold
    ]
    background_same = [row for row in background_changed if row["same_direction"] == "yes"]
    odds_ratio, fisher_p = (math.nan, math.nan)
    if changed and background_changed:
        table = [
            [len(same_changed), len(changed) - len(same_changed)],
            [len(background_same), len(background_changed) - len(background_same)],
        ]
        odds_ratio, fisher_p = stats.fisher_exact(table, alternative="greater")

    top_same = sorted(
        [row for row in changed if row["same_direction"] == "yes"],
        key=lambda row: -abs(safe_float(row["dnmt3l_log2fc"]) * safe_float(row["ds_delta"])),
    )[:20]
    return {
        "dataset": dataset,
        "context": context,
        "species": species,
        "module": module,
        "module_genes_matched": n,
        "background_genes_matched": bg_n,
        "module_mean_dnmt3l_log2fc": mean(dnmt3l),
        "module_median_dnmt3l_log2fc": median(dnmt3l),
        "module_mean_ds_delta": mean(ds),
        "module_median_ds_delta": median(ds),
        "spearman_rho_dnmt3l_vs_ds": float(rho) if math.isfinite(float(rho)) else math.nan,
        "spearman_p": float(rho_p) if math.isfinite(float(rho_p)) else math.nan,
        "mean_directional_score": mean(score),
        "background_mean_directional_score": mean(bg_score),
        "median_directional_score": median(score),
        "background_median_directional_score": median(bg_score),
        "mannwhitney_p_directional_score_gt_background": mw_p,
        "changed_gene_count": len(changed),
        "same_direction_changed_genes": len(same_changed),
        "same_direction_changed_fraction": len(same_changed) / len(changed) if changed else math.nan,
        "background_changed_gene_count": len(background_changed),
        "background_same_direction_fraction": len(background_same) / len(background_changed) if background_changed else math.nan,
        "same_direction_enrichment_odds_ratio": odds_ratio,
        "same_direction_enrichment_fisher_p": fisher_p,
        "top_same_direction_genes": ",".join(str(row["mouse_symbol"]) for row in top_same),
    }


def make_gene_rows(
    dataset: str,
    context: str,
    species: str,
    module: str,
    paired: list[dict[str, object]],
    dnmt3l_threshold: float,
    ds_threshold: float,
) -> list[dict[str, object]]:
    rows = []
    for row in paired:
        dnmt3l = safe_float(row["dnmt3l_log2fc"])
        ds_delta = safe_float(row["ds_delta"])
        rows.append(
            {
                "dataset": dataset,
                "context": context,
                "species": species,
                "module": module,
                "mouse_symbol": row["mouse_symbol"],
                "matched_symbol": row["matched_symbol"],
                "dnmt3l_log2fc": dnmt3l,
                "ds_delta": ds_delta,
                "signed_ds_delta_by_dnmt3l_direction": row["signed_ds_delta_by_dnmt3l_direction"],
                "dnmt3l_direction": direction(dnmt3l, dnmt3l_threshold),
                "ds_direction": direction(ds_delta, ds_threshold),
                "same_direction": row["same_direction"],
                "ds_effect_source": row.get("ds_effect_source", ""),
                "ds_edgeR_FDR": row.get("ds_edgeR_FDR", ""),
                "ds_tested_by_edgeR": row.get("ds_tested_by_edgeR", ""),
                "passes_change_thresholds": (
                    "yes" if abs(dnmt3l) >= dnmt3l_threshold and abs(ds_delta) >= ds_threshold else "no"
                ),
            }
        )
    rows.sort(
        key=lambda item: (
            item["passes_change_thresholds"] != "yes",
            item["same_direction"] != "yes",
            -abs(safe_float(item["signed_ds_delta_by_dnmt3l_direction"])),
            str(item["mouse_symbol"]),
        )
    )
    return rows


def make_heatmap(summary_rows: list[dict[str, object]], figdir: Path) -> None:
    figdir.mkdir(parents=True, exist_ok=True)
    modules = MODULES
    contexts = []
    for row in summary_rows:
        label = f"{row['dataset']}:{row['context']}"
        if label not in contexts:
            contexts.append(label)
    values = np.full((len(modules), len(contexts)), np.nan)
    for row in summary_rows:
        module_idx = modules.index(str(row["module"]))
        context_idx = contexts.index(f"{row['dataset']}:{row['context']}")
        values[module_idx, context_idx] = safe_float(row["mean_directional_score"])

    vmin, vmax = -0.08, 0.08
    fig, ax = plt.subplots(figsize=(max(7.4, len(contexts) * 1.55), 4.5), constrained_layout=True)
    im = ax.imshow(values, cmap="RdBu_r", vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(contexts)))
    ax.set_xticklabels(contexts, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(range(len(modules)))
    ax.set_yticklabels(["Glial/myelination", "Neurodevelopment signaling"], fontsize=9)
    ax.set_title("DS brain concordance with Dnmt3l OE module direction")
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            value = values[i, j]
            if math.isfinite(value):
                ax.text(
                    j,
                    i,
                    f"{value:.3f}",
                    ha="center",
                    va="center",
                    fontsize=10.5,
                    color=heatmap_text_color(value, vmin, vmax),
                )
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.025)
    cbar.set_label("Mean signed DS delta\nby Dnmt3l OE direction", fontsize=9, labelpad=10)
    fig.savefig(figdir / "dnmt3l_module_ds_brain_directional_score_heatmap.png", dpi=260, bbox_inches="tight")
    fig.savefig(figdir / "dnmt3l_module_ds_brain_directional_score_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)


def write_note(path: Path, summary_rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Dnmt3l OE Module Concordance With DS Brain",
        "",
        "Date: 2026-06-09",
        "",
        "## Question",
        "",
        "Do Dnmt3l-overexpression module shifts match DS brain expression changes?",
        "",
        "## Method",
        "",
        "For each module gene, the Dnmt3l OE mouse log2FC was compared with DS-control expression deltas from Dp10 mouse hippocampus and human DS fetal brain bulk datasets. A positive directional score means the DS delta moves in the same direction as the Dnmt3l OE effect.",
        "",
        "## Key Results",
        "",
    ]
    for row in summary_rows:
        lines.append(
            "- "
            f"{row['module']} in {row['dataset']}:{row['context']}: "
            f"matched genes {row['module_genes_matched']}, "
            f"mean directional score {fmt(row['mean_directional_score'])}, "
            f"Spearman rho {fmt(row['spearman_rho_dnmt3l_vs_ds'])}, "
            f"same-direction changed fraction {fmt(row['same_direction_changed_fraction'])}."
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "Concordance should be interpreted as module-level support, not proof that the same causal mechanism operates in every dataset. The mouse Dp10 comparison is species- and dosage-relevant but is hippocampal single-cell pseudobulk with very small n. The human bulk comparisons use edgeR count-based DS-control logFC values with age and sex covariates where the design is estimable, but they still do not fully resolve cell-composition effects.",
            "",
            "## Outputs",
            "",
            "- `results/neurogenesis_mechanism/dnmt3l_module_ds_brain_concordance_summary.tsv`",
            "- `results/neurogenesis_mechanism/dnmt3l_module_ds_brain_gene_concordance.tsv`",
            "- `figures/neurogenesis_mechanism/dnmt3l_module_ds_brain_directional_score_heatmap.png`",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    figdir = Path(args.figdir)
    signature = load_dnmt3l_signature(Path(args.dnmt3l_expr), args.min_expression)
    modules = load_modules(Path(args.membership))

    context_effects: dict[tuple[str, str, str], dict[str, dict[str, object]]] = {}
    context_effects[("GSE214273", "Dp10_hippocampus_pseudobulk", "mouse")] = load_dp10_effects(Path(args.dp10))

    human_count_de = load_human_count_de_effects(Path(args.human_ds_bulk_count_de))
    if human_count_de:
        for (dataset, context), effects in human_count_de.items():
            context_effects[(dataset, context, "human")] = effects
    else:
        id_to_symbol = gse305139_symbol_map(Path(args.gse305139_counts))
        for (dataset, context), effects in count_matrix_effects(
            "GSE301886",
            Path(args.gse301886_counts),
            Path(args.gse301886_meta),
            id_to_symbol=id_to_symbol,
        ).items():
            context_effects[(dataset, context, "human")] = effects
        for (dataset, context), effects in count_matrix_effects(
            "GSE305139",
            Path(args.gse305139_counts),
            Path(args.gse305139_meta),
        ).items():
            context_effects[(dataset, context, "human")] = effects

    summary_rows: list[dict[str, object]] = []
    gene_rows: list[dict[str, object]] = []
    for (dataset, context, species), effects in sorted(context_effects.items()):
        for module in MODULES:
            module_genes = modules.get(module, set())
            paired, background = paired_rows_for_context(signature, effects, module_genes, species)
            summary_rows.append(
                summarize_match(
                    dataset,
                    context,
                    species,
                    module,
                    paired,
                    background,
                    args.dnmt3l_change_threshold,
                    args.ds_change_threshold,
                )
            )
            gene_rows.extend(
                make_gene_rows(
                    dataset,
                    context,
                    species,
                    module,
                    paired,
                    args.dnmt3l_change_threshold,
                    args.ds_change_threshold,
                )
            )

    summary_fields = [
        "dataset",
        "context",
        "species",
        "module",
        "module_genes_matched",
        "background_genes_matched",
        "module_mean_dnmt3l_log2fc",
        "module_median_dnmt3l_log2fc",
        "module_mean_ds_delta",
        "module_median_ds_delta",
        "spearman_rho_dnmt3l_vs_ds",
        "spearman_p",
        "mean_directional_score",
        "background_mean_directional_score",
        "median_directional_score",
        "background_median_directional_score",
        "mannwhitney_p_directional_score_gt_background",
        "changed_gene_count",
        "same_direction_changed_genes",
        "same_direction_changed_fraction",
        "background_changed_gene_count",
        "background_same_direction_fraction",
        "same_direction_enrichment_odds_ratio",
        "same_direction_enrichment_fisher_p",
        "top_same_direction_genes",
    ]
    gene_fields = [
        "dataset",
        "context",
        "species",
        "module",
        "mouse_symbol",
        "matched_symbol",
        "dnmt3l_log2fc",
        "ds_delta",
        "signed_ds_delta_by_dnmt3l_direction",
        "dnmt3l_direction",
        "ds_direction",
        "same_direction",
        "ds_effect_source",
        "ds_edgeR_FDR",
        "ds_tested_by_edgeR",
        "passes_change_thresholds",
    ]
    write_tsv(outdir / "dnmt3l_module_ds_brain_concordance_summary.tsv", summary_rows, summary_fields)
    write_tsv(outdir / "dnmt3l_module_ds_brain_gene_concordance.tsv", gene_rows, gene_fields)
    make_heatmap(summary_rows, figdir)
    write_note(Path(args.note), summary_rows)
    print(f"Wrote {outdir / 'dnmt3l_module_ds_brain_concordance_summary.tsv'}")
    print(f"Wrote {outdir / 'dnmt3l_module_ds_brain_gene_concordance.tsv'}")
    print(f"Wrote {figdir / 'dnmt3l_module_ds_brain_directional_score_heatmap.png'}")
    print(f"Wrote {args.note}")


if __name__ == "__main__":
    main()
