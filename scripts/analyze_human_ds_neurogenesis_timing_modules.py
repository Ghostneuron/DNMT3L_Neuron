#!/usr/bin/env python3
"""Test neurogenesis and maturation module shifts in human DS brain samples."""

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


DEFAULT_MEMBERSHIP = "results/neurogenesis_mechanism/developmental_module_gene_membership.tsv"
DEFAULT_GSE301886_COUNTS = "data/geo/GSE301886/supplementary/GSE301886_T21_raw_counts.txt.gz"
DEFAULT_GSE301886_META = "data/geo/GSE301886/metadata/sample_metadata.tsv"
DEFAULT_GSE305139_COUNTS = "data/geo/GSE305139/supplementary/GSE305139_raw_counts_tissue_bulk.tsv.gz"
DEFAULT_GSE305139_META = "data/geo/GSE305139/metadata/sample_metadata.tsv"
DEFAULT_HUMAN_DS_BULK_COUNT_DE = "results/human_ds_bulk_count_de/human_ds_bulk_count_de_all_contexts.tsv"
DEFAULT_OUTDIR = "results/neurogenesis_mechanism"
DEFAULT_FIGDIR = "figures/neurogenesis_mechanism"
DEFAULT_NOTE = "notes/human_ds_neurogenesis_timing_modules.md"

MODULES = [
    "cell_cycle_progenitor_proliferation",
    "early_neurogenesis_neuron_differentiation",
    "cortical_patterning_region_identity",
    "axon_dendrite_neuron_projection",
    "synapse_formation_maturation",
]

MODULE_LABELS = {
    "cell_cycle_progenitor_proliferation": "Cell cycle/progenitor",
    "early_neurogenesis_neuron_differentiation": "Early neurogenesis",
    "cortical_patterning_region_identity": "Cortical patterning",
    "axon_dendrite_neuron_projection": "Axon/dendrite projection",
    "synapse_formation_maturation": "Synapse maturation",
}

TIMING_MARKERS = [
    "PAX6",
    "SOX2",
    "NES",
    "MKI67",
    "EOMES",
    "ASCL1",
    "HES1",
    "HES5",
    "NEUROD1",
    "NEUROD2",
    "DCX",
    "TUBB3",
    "RBFOX3",
    "SATB2",
    "BCL11B",
    "TLE4",
    "CUX2",
    "RELN",
    "MAP2",
    "SYN1",
    "SNAP25",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--membership", default=DEFAULT_MEMBERSHIP)
    parser.add_argument("--gse301886-counts", default=DEFAULT_GSE301886_COUNTS)
    parser.add_argument("--gse301886-meta", default=DEFAULT_GSE301886_META)
    parser.add_argument("--gse305139-counts", default=DEFAULT_GSE305139_COUNTS)
    parser.add_argument("--gse305139-meta", default=DEFAULT_GSE305139_META)
    parser.add_argument("--human-ds-bulk-count-de", default=DEFAULT_HUMAN_DS_BULK_COUNT_DE)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--figdir", default=DEFAULT_FIGDIR)
    parser.add_argument("--note", default=DEFAULT_NOTE)
    parser.add_argument("--min-mean-log2cpm", type=float, default=0.5)
    parser.add_argument("--delta-threshold", type=float, default=0.05)
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


def bh_fdr(pvalues: list[float]) -> list[float]:
    n = len(pvalues)
    order = sorted(range(n), key=lambda idx: math.inf if not math.isfinite(pvalues[idx]) else pvalues[idx])
    qvalues = [math.nan] * n
    running = 1.0
    for rank, idx in reversed(list(enumerate(order, start=1))):
        pvalue = pvalues[idx]
        if not math.isfinite(pvalue):
            continue
        running = min(running, pvalue * n / rank)
        qvalues[idx] = min(running, 1.0)
    return qvalues


def load_modules(path: Path) -> dict[str, set[str]]:
    modules: dict[str, set[str]] = defaultdict(set)
    for row in read_tsv(path):
        module = row.get("module", "")
        symbol = row.get("gene_symbol", "")
        if module in MODULES and symbol:
            modules[module].add(symbol.upper())
    return modules


def gse305139_symbol_map(counts_path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with gzip.open(counts_path, "rt", newline="", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene_id = clean_ensembl(row.get("gene_id", ""))
            gene_name = row.get("gene_name", "")
            if gene_id and gene_name:
                mapping[gene_id] = gene_name.upper()
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
        context = meta.get("region", "") if dataset == "GSE301886" else "fetal_brain_final"
        contexts[context or "all"].append(sample)
    return contexts


def count_matrix_effects(
    dataset: str,
    counts_path: Path,
    metadata_path: Path,
    id_to_symbol: dict[str, str] | None = None,
) -> dict[tuple[str, str], dict[str, dict[str, object]]]:
    metadata = load_gse301886_metadata(metadata_path) if dataset == "GSE301886" else load_gse305139_metadata(metadata_path)
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
            symbol = row[1].upper() if has_gene_name else (id_to_symbol or {}).get(gene_id, "")
            if not symbol:
                continue
            counts = [safe_float(value) for value in row[first_sample:]]
            counts = [0.0 if not math.isfinite(value) else value for value in counts]
            for sample, count in zip(sample_labels, counts, strict=False):
                library_sizes[sample] += count
            raw_rows.append((symbol, counts))

    sample_index = {sample: idx for idx, sample in enumerate(sample_labels)}
    effects_by_context: dict[tuple[str, str], dict[str, dict[str, object]]] = {}
    for context, samples in contexts.items():
        ds_samples = [sample for sample in samples if metadata.get(sample, {}).get("group") == "DS"]
        control_samples = [sample for sample in samples if metadata.get(sample, {}).get("group") == "Control"]
        if len(ds_samples) < 2 or len(control_samples) < 2:
            continue
        effects: dict[str, dict[str, object]] = {}
        for symbol, counts in raw_rows:
            ds_values: list[float] = []
            control_values: list[float] = []
            for sample in ds_samples:
                idx = sample_index[sample]
                cpm = counts[idx] / library_sizes[sample] * 1_000_000 if library_sizes[sample] else math.nan
                ds_values.append(math.log2(cpm + 1.0) if math.isfinite(cpm) else math.nan)
            for sample in control_samples:
                idx = sample_index[sample]
                cpm = counts[idx] / library_sizes[sample] * 1_000_000 if library_sizes[sample] else math.nan
                control_values.append(math.log2(cpm + 1.0) if math.isfinite(cpm) else math.nan)
            ds_mean = mean(ds_values)
            control_mean = mean(control_values)
            max_mean = max(ds_mean, control_mean)
            delta = ds_mean - control_mean
            if math.isfinite(delta):
                effects[symbol] = {
                    "dataset": dataset,
                    "context": context,
                    "gene_symbol": symbol,
                    "mean_log2_cpm_ds": ds_mean,
                    "mean_log2_cpm_control": control_mean,
                    "max_mean_log2_cpm": max_mean,
                    "delta_ds_minus_control": delta,
                    "n_ds": len(ds_samples),
                    "n_control": len(control_samples),
                }
        effects_by_context[(dataset, context)] = effects
    return effects_by_context


def load_count_de_effects(path: Path) -> dict[tuple[str, str], dict[str, dict[str, object]]]:
    effects_by_context: dict[tuple[str, str], dict[str, dict[str, object]]] = defaultdict(dict)
    if not path.exists():
        return effects_by_context
    for row in read_tsv(path):
        dataset = row.get("dataset", "")
        context = row.get("context", "")
        symbol = row.get("gene_symbol", "").upper()
        if dataset not in {"GSE301886", "GSE305139"} or not context or not symbol:
            continue
        ds_mean = safe_float(row.get("mean_ds_log2_cpm_plus_1"))
        control_mean = safe_float(row.get("mean_control_log2_cpm_plus_1"))
        effect = safe_float(row.get("effect_log2fc_for_claim"))
        if not math.isfinite(effect):
            continue
        current = effects_by_context[(dataset, context)].get(symbol)
        max_expr = max(ds_mean, control_mean)
        current_max = safe_float(current.get("max_mean_log2_cpm", "")) if current else -math.inf
        if current is not None and current_max >= max_expr:
            continue
        effects_by_context[(dataset, context)][symbol] = {
            "dataset": dataset,
            "context": context,
            "gene_symbol": symbol,
            "mean_log2_cpm_ds": ds_mean,
            "mean_log2_cpm_control": control_mean,
            "max_mean_log2_cpm": max_expr,
            "delta_ds_minus_control": effect,
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


def module_result(
    dataset: str,
    context: str,
    module: str,
    genes: set[str],
    effects: dict[str, dict[str, object]],
    min_expression: float,
    delta_threshold: float,
) -> tuple[dict[str, object] | None, list[dict[str, object]]]:
    universe = {
        symbol
        for symbol, row in effects.items()
        if safe_float(row.get("max_mean_log2_cpm")) >= min_expression
        and math.isfinite(safe_float(row.get("delta_ds_minus_control")))
    }
    module_genes = genes & universe
    background = universe - module_genes
    if len(module_genes) < 10 or len(background) < 100:
        return None, []

    module_values = [safe_float(effects[symbol]["delta_ds_minus_control"]) for symbol in module_genes]
    background_values = [safe_float(effects[symbol]["delta_ds_minus_control"]) for symbol in background]
    u_test = stats.mannwhitneyu(module_values, background_values, alternative="two-sided")
    up_test = stats.mannwhitneyu(module_values, background_values, alternative="greater")
    down_test = stats.mannwhitneyu(module_values, background_values, alternative="less")
    n1 = len(module_values)
    n2 = len(background_values)
    rank_effect = (2 * float(u_test.statistic) / (n1 * n2)) - 1
    up_genes = [symbol for symbol in module_genes if safe_float(effects[symbol]["delta_ds_minus_control"]) >= delta_threshold]
    down_genes = [symbol for symbol in module_genes if safe_float(effects[symbol]["delta_ds_minus_control"]) <= -delta_threshold]

    gene_rows = []
    for symbol in sorted(module_genes):
        row = effects[symbol]
        delta = safe_float(row["delta_ds_minus_control"])
        gene_rows.append(
            {
                "dataset": dataset,
                "context": context,
                "module": module,
                "gene_symbol": symbol,
                "mean_log2_cpm_ds": row["mean_log2_cpm_ds"],
                "mean_log2_cpm_control": row["mean_log2_cpm_control"],
                "delta_ds_minus_control": delta,
                "descriptive_delta_log2_cpm": row.get("descriptive_delta_log2_cpm", ""),
                "effect_source": row.get("effect_source", ""),
                "edgeR_logFC": row.get("edgeR_logFC", ""),
                "edgeR_PValue": row.get("edgeR_PValue", ""),
                "edgeR_FDR": row.get("edgeR_FDR", ""),
                "tested_by_edgeR": row.get("tested_by_edgeR", ""),
                "direction": "up" if delta >= delta_threshold else ("down" if delta <= -delta_threshold else "flat"),
            }
        )
    gene_rows.sort(key=lambda row: -abs(safe_float(row["delta_ds_minus_control"])))
    result = {
        "dataset": dataset,
        "context": context,
        "module": module,
        "module_label": MODULE_LABELS.get(module, module),
        "module_gene_count": len(module_values),
        "background_gene_count": len(background_values),
        "module_mean_delta": mean(module_values),
        "background_mean_delta": mean(background_values),
        "delta_mean_vs_background": mean(module_values) - mean(background_values),
        "module_median_delta": median(module_values),
        "background_median_delta": median(background_values),
        "delta_median_vs_background": median(module_values) - median(background_values),
        "rank_biserial_effect": rank_effect,
        "mannwhitney_p_two_sided": float(u_test.pvalue),
        "mannwhitney_p_up": float(up_test.pvalue),
        "mannwhitney_p_down": float(down_test.pvalue),
        "module_up_genes": len(up_genes),
        "module_down_genes": len(down_genes),
        "module_up_fraction": len(up_genes) / len(module_values),
        "module_down_fraction": len(down_genes) / len(module_values),
        "leading_up_genes": ",".join(symbol for symbol, _delta in top_genes(module_genes, effects, up=True)),
        "leading_down_genes": ",".join(symbol for symbol, _delta in top_genes(module_genes, effects, up=False)),
    }
    return result, gene_rows


def top_genes(genes: set[str], effects: dict[str, dict[str, object]], up: bool, n: int = 15) -> list[tuple[str, float]]:
    selected = [(symbol, safe_float(effects[symbol]["delta_ds_minus_control"])) for symbol in genes]
    selected.sort(key=lambda item: item[1], reverse=up)
    return selected[:n]


def timing_call(results_by_module: dict[str, dict[str, object]]) -> tuple[str, str]:
    early = safe_float(results_by_module.get("early_neurogenesis_neuron_differentiation", {}).get("delta_median_vs_background"))
    progenitor = safe_float(results_by_module.get("cell_cycle_progenitor_proliferation", {}).get("delta_median_vs_background"))
    axon = safe_float(results_by_module.get("axon_dendrite_neuron_projection", {}).get("delta_median_vs_background"))
    synapse = safe_float(results_by_module.get("synapse_formation_maturation", {}).get("delta_median_vs_background"))
    if early > 0 and progenitor > 0 and (synapse <= 0 or axon <= 0):
        return "delayed/immature", "early/progenitor programs are elevated while maturation modules are not comparably elevated"
    if early > 0 and (axon > 0 or synapse > 0):
        return "mixed accelerated differentiation", "early neurogenesis and maturation/projection programs both shift upward"
    if early < 0 and (axon > 0 or synapse > 0):
        return "premature maturation", "early neurogenesis is reduced while later neuronal maturation is elevated"
    if early < 0 and synapse < 0 and axon < 0:
        return "broad suppression/delay", "early and later neuronal modules are both reduced"
    return "mixed/unclear", "module directions do not form a simple delay-versus-premature pattern"


def marker_rows(
    dataset: str,
    context: str,
    effects: dict[str, dict[str, object]],
    delta_threshold: float,
) -> list[dict[str, object]]:
    rows = []
    for symbol in TIMING_MARKERS:
        row = effects.get(symbol)
        if not row:
            rows.append({"dataset": dataset, "context": context, "gene_symbol": symbol, "status": "not_found"})
            continue
        delta = safe_float(row["delta_ds_minus_control"])
        rows.append(
            {
                "dataset": dataset,
                "context": context,
                "gene_symbol": symbol,
                "status": "found",
                "mean_log2_cpm_ds": row["mean_log2_cpm_ds"],
                "mean_log2_cpm_control": row["mean_log2_cpm_control"],
                "delta_ds_minus_control": delta,
                "descriptive_delta_log2_cpm": row.get("descriptive_delta_log2_cpm", ""),
                "effect_source": row.get("effect_source", ""),
                "edgeR_logFC": row.get("edgeR_logFC", ""),
                "edgeR_PValue": row.get("edgeR_PValue", ""),
                "edgeR_FDR": row.get("edgeR_FDR", ""),
                "tested_by_edgeR": row.get("tested_by_edgeR", ""),
                "direction": "up" if delta >= delta_threshold else ("down" if delta <= -delta_threshold else "flat"),
            }
        )
    return rows


def make_heatmap(summary_rows: list[dict[str, object]], figdir: Path) -> None:
    figdir.mkdir(parents=True, exist_ok=True)
    contexts = []
    for row in summary_rows:
        label = f"{row['dataset']}:{row['context']}"
        if label not in contexts:
            contexts.append(label)
    values = np.full((len(MODULES), len(contexts)), np.nan)
    for row in summary_rows:
        i = MODULES.index(str(row["module"]))
        j = contexts.index(f"{row['dataset']}:{row['context']}")
        values[i, j] = safe_float(row["delta_median_vs_background"])
    vmin, vmax = -0.08, 0.08
    fig, ax = plt.subplots(figsize=(max(7.8, len(contexts) * 1.65), 5.2), constrained_layout=True)
    im = ax.imshow(values, cmap="RdBu_r", vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(contexts)))
    ax.set_xticklabels(contexts, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(range(len(MODULES)))
    ax.set_yticklabels([MODULE_LABELS[module] for module in MODULES], fontsize=8.5)
    ax.set_title("Human DS neurogenesis/maturation module shifts")
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
    cbar.set_label("Module median DS delta\nminus background median", fontsize=9, labelpad=10)
    fig.savefig(figdir / "human_ds_neurogenesis_timing_module_heatmap.png", dpi=260, bbox_inches="tight")
    fig.savefig(figdir / "human_ds_neurogenesis_timing_module_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)


def write_note(path: Path, timing_rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Human DS Neurogenesis Timing Modules",
        "",
        "Date: 2026-06-09",
        "",
        "## Question",
        "",
        "Are neurogenesis-related programs in human DS brain delayed or premature?",
        "",
        "## Interpretation Rule",
        "",
        "A delayed/immature pattern means early neurogenesis or progenitor modules remain elevated relative to background, especially if later axon/synapse maturation is not equally elevated. A premature maturation pattern means early neurogenesis is reduced while later neuronal maturation modules are elevated.",
        "",
        "## Context Calls",
        "",
    ]
    for row in timing_rows:
        lines.append(f"- {row['dataset']}:{row['context']}: {row['timing_call']} ({row['reason']}).")
    lines.extend(
        [
            "",
            "## Caveats",
            "",
            "These are module calls from bulk RNA-seq DS-control edgeR logFC values. They do not fully separate cell-composition effects from within-cell transcriptional timing, and the mouse-derived module symbols were mapped to human mostly by conserved gene symbols.",
            "",
            "## Outputs",
            "",
            "- `results/neurogenesis_mechanism/human_ds_neurogenesis_timing_module_summary.tsv`",
            "- `results/neurogenesis_mechanism/human_ds_neurogenesis_timing_gene_effects.tsv`",
            "- `results/neurogenesis_mechanism/human_ds_neurogenesis_timing_marker_effects.tsv`",
            "- `figures/neurogenesis_mechanism/human_ds_neurogenesis_timing_module_heatmap.png`",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    modules = load_modules(Path(args.membership))
    all_effects: dict[tuple[str, str], dict[str, dict[str, object]]] = {}
    all_effects.update(load_count_de_effects(Path(args.human_ds_bulk_count_de)))
    if not all_effects:
        id_to_symbol = gse305139_symbol_map(Path(args.gse305139_counts))
        all_effects.update(
            count_matrix_effects(
                "GSE301886",
                Path(args.gse301886_counts),
                Path(args.gse301886_meta),
                id_to_symbol=id_to_symbol,
            )
        )
        all_effects.update(
            count_matrix_effects(
                "GSE305139",
                Path(args.gse305139_counts),
                Path(args.gse305139_meta),
            )
        )

    summary_rows: list[dict[str, object]] = []
    gene_rows: list[dict[str, object]] = []
    marker_effect_rows: list[dict[str, object]] = []
    timing_rows: list[dict[str, object]] = []

    for (dataset, context), effects in sorted(all_effects.items()):
        context_results: dict[str, dict[str, object]] = {}
        for module in MODULES:
            result, genes = module_result(
                dataset,
                context,
                module,
                modules[module],
                effects,
                args.min_mean_log2cpm,
                args.delta_threshold,
            )
            if result:
                context_results[module] = result
                summary_rows.append(result)
                gene_rows.extend(genes)
        call, reason = timing_call(context_results)
        timing_rows.append({"dataset": dataset, "context": context, "timing_call": call, "reason": reason})
        marker_effect_rows.extend(marker_rows(dataset, context, effects, args.delta_threshold))

    fdrs = bh_fdr([safe_float(row["mannwhitney_p_two_sided"]) for row in summary_rows])
    for row, fdr in zip(summary_rows, fdrs):
        row["mannwhitney_fdr_two_sided"] = fdr

    summary_fields = [
        "dataset",
        "context",
        "module",
        "module_label",
        "module_gene_count",
        "background_gene_count",
        "module_mean_delta",
        "background_mean_delta",
        "delta_mean_vs_background",
        "module_median_delta",
        "background_median_delta",
        "delta_median_vs_background",
        "rank_biserial_effect",
        "mannwhitney_p_two_sided",
        "mannwhitney_fdr_two_sided",
        "mannwhitney_p_up",
        "mannwhitney_p_down",
        "module_up_genes",
        "module_down_genes",
        "module_up_fraction",
        "module_down_fraction",
        "leading_up_genes",
        "leading_down_genes",
    ]
    gene_fields = [
        "dataset",
        "context",
        "module",
        "gene_symbol",
        "mean_log2_cpm_ds",
        "mean_log2_cpm_control",
        "delta_ds_minus_control",
        "descriptive_delta_log2_cpm",
        "effect_source",
        "edgeR_logFC",
        "edgeR_PValue",
        "edgeR_FDR",
        "tested_by_edgeR",
        "direction",
    ]
    marker_fields = [
        "dataset",
        "context",
        "gene_symbol",
        "status",
        "mean_log2_cpm_ds",
        "mean_log2_cpm_control",
        "delta_ds_minus_control",
        "descriptive_delta_log2_cpm",
        "effect_source",
        "edgeR_logFC",
        "edgeR_PValue",
        "edgeR_FDR",
        "tested_by_edgeR",
        "direction",
    ]
    timing_fields = ["dataset", "context", "timing_call", "reason"]
    outdir = Path(args.outdir)
    write_tsv(outdir / "human_ds_neurogenesis_timing_module_summary.tsv", summary_rows, summary_fields)
    write_tsv(outdir / "human_ds_neurogenesis_timing_gene_effects.tsv", gene_rows, gene_fields)
    write_tsv(outdir / "human_ds_neurogenesis_timing_marker_effects.tsv", marker_effect_rows, marker_fields)
    write_tsv(outdir / "human_ds_neurogenesis_timing_context_calls.tsv", timing_rows, timing_fields)
    make_heatmap(summary_rows, Path(args.figdir))
    write_note(Path(args.note), timing_rows)
    print(f"Wrote {outdir / 'human_ds_neurogenesis_timing_module_summary.tsv'}")
    print(f"Wrote {outdir / 'human_ds_neurogenesis_timing_gene_effects.tsv'}")
    print(f"Wrote {outdir / 'human_ds_neurogenesis_timing_marker_effects.tsv'}")
    print(f"Wrote {outdir / 'human_ds_neurogenesis_timing_context_calls.tsv'}")
    print(f"Wrote {args.figdir}/human_ds_neurogenesis_timing_module_heatmap.png")
    print(f"Wrote {args.note}")


if __name__ == "__main__":
    main()
