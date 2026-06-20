#!/usr/bin/env python3
"""Rank-based module tests for neurogenesis and synapse-window gene sets."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


DEFAULT_EXPR = "results/neurogenesis_mechanism/gse268951_rederived_gene_expression.tsv"
DEFAULT_OVERLAP = "results/neurogenesis_mechanism/dnmt3l_vs_dnmt3a_transcriptomic_overlap.tsv"
DEFAULT_FINAL = "results/neurogenesis_mechanism/final_candidate_neurogenesis_annotation.tsv"
DEFAULT_OUTDIR = "results/neurogenesis_mechanism"
DEFAULT_FIGDIR = "figures/neurogenesis_mechanism"
DEFAULT_NOTE = "notes/developmental_module_shift_tests.md"


MODULE_KEYWORDS = {
    "early_neurogenesis_neuron_differentiation": [
        "neurogenesis",
        "generation of neurons",
        "neuron differentiation",
        "neuronal differentiation",
        "neuron development",
        "neural precursor",
        "neural progenitor",
        "neural stem",
    ],
    "cortical_patterning_region_identity": [
        "forebrain development",
        "telencephalon",
        "cerebral cortex",
        "cortical",
        "diencephalon",
        "midbrain development",
        "pattern specification",
        "regionalization",
    ],
    "axon_dendrite_neuron_projection": [
        "axon",
        "dendrite",
        "neuron projection",
        "neurite",
        "growth cone",
        "projection neuron",
    ],
    "synapse_formation_maturation": [
        "synapse assembly",
        "synapse organization",
        "synapse maturation",
        "synaptic plasticity",
        "synaptic transmission",
        "postsynaptic",
        "presynaptic",
        "synaptic vesicle",
        "dendritic spine",
        "synaptogenesis",
        "synapse structure",
    ],
    "glial_maturation_myelination": [
        "gliogenesis",
        "glial",
        "astrocyte",
        "oligodendrocyte",
        "myelination",
        "myelin",
    ],
    "neurodevelopment_signaling": [
        "wnt",
        "notch",
        "bmp",
        "fibroblast growth factor",
        "fgf",
        "retinoic acid",
        "tgf",
        "mapk",
        "pi3k",
        "akt",
    ],
    "cell_cycle_progenitor_proliferation": [
        "cell cycle",
        "cell division",
        "mitotic",
        "proliferation",
        "stem cell proliferation",
        "neural precursor cell proliferation",
        "neural progenitor cell proliferation",
    ],
}


CONTRASTS = [
    {
        "contrast": "Dnmt3l_OE_vs_control_P28_medial_frontal_cortex",
        "score_field": "gse268951_log2fc_oe_vs_control",
        "source": "GSE268951",
    },
    {
        "contrast": "Dnmt3a_cKO_vs_control_P39_excitatory_neurons",
        "score_field": "dnmt3a_logfc_cko_vs_control",
        "source": "GSE141587",
    },
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expr", default=DEFAULT_EXPR)
    parser.add_argument("--overlap", default=DEFAULT_OVERLAP)
    parser.add_argument("--final", default=DEFAULT_FINAL)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--figdir", default=DEFAULT_FIGDIR)
    parser.add_argument("--note", default=DEFAULT_NOTE)
    parser.add_argument("--min-module-size", type=int, default=10)
    return parser.parse_args()


def clean_key(value: str | None) -> str:
    return (value or "").strip().replace("\r", "")


def safe_float(value: object) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return math.nan
    return number if math.isfinite(number) else math.nan


def fmt(value: object, digits: int = 6) -> str:
    number = safe_float(value)
    if not math.isfinite(number):
        return ""
    return f"{number:.{digits}g}"


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
        writer.writerows(rows)


def bh_fdr(pvalues: list[float]) -> list[float]:
    n = len(pvalues)
    order = sorted(range(n), key=lambda i: math.inf if not math.isfinite(pvalues[i]) else pvalues[i])
    qvalues = [math.nan] * n
    running = 1.0
    for rank, idx in reversed(list(enumerate(order, start=1))):
        pvalue = pvalues[idx]
        if not math.isfinite(pvalue):
            continue
        running = min(running, pvalue * n / rank)
        qvalues[idx] = min(running, 1.0)
    return qvalues


def keyword_hit(text: str, keywords: list[str]) -> bool:
    lowered = text.lower()
    return any(keyword in lowered for keyword in keywords)


def load_expression_rows(expr_path: Path, overlap_path: Path) -> dict[str, dict[str, object]]:
    overlap_by_symbol = {row["gene_symbol"]: row for row in read_tsv(overlap_path)}
    rows: dict[str, dict[str, object]] = {}
    for row in read_tsv(expr_path):
        symbol = row["gene_symbol"]
        overlap = overlap_by_symbol.get(symbol, {})
        mean_fpkm = max(
            safe_float(row.get("gse268951_control_mean_fpkm")),
            safe_float(row.get("gse268951_dnmt3l_oe_mean_fpkm")),
        )
        record = {
            "gene_symbol": symbol,
            "gene_id": row.get("gene_id", ""),
            "description": row.get("description", ""),
            "go_bp_terms": row.get("go_bp_terms", ""),
            "mean_fpkm_max": mean_fpkm,
            "gse268951_log2fc_oe_vs_control": safe_float(row.get("gse268951_log2fc_oe_vs_control")),
            "dnmt3a_logfc_cko_vs_control": safe_float(overlap.get("dnmt3a_logfc_cko_vs_control")),
            "dnmt3a_fdr": safe_float(overlap.get("dnmt3a_fdr")),
        }
        rows[symbol] = record
    return rows


def build_modules(rows: dict[str, dict[str, object]], final_rows: list[dict[str, str]]) -> dict[str, set[str]]:
    modules: dict[str, set[str]] = defaultdict(set)
    for symbol, row in rows.items():
        text = f"{row.get('go_bp_terms', '')}; {row.get('description', '')}"
        for module, keywords in MODULE_KEYWORDS.items():
            if keyword_hit(text, keywords):
                modules[module].add(symbol)

    final_all = {row["mouse_symbol"] for row in final_rows if row.get("mouse_symbol")}
    final_neurodevelopment = {
        row["mouse_symbol"]
        for row in final_rows
        if row.get("mouse_symbol") and row.get("neurogenesis_related") == "yes"
    }
    final_synapse_window = {
        row["mouse_symbol"]
        for row in final_rows
        if row.get("mouse_symbol")
        and (
            row.get("neurite_synapse_migration") == "yes"
            or row.get("neurogenesis_signaling") == "yes"
            or row.get("core_neurogenesis") == "yes"
        )
    }
    final_primary = {
        row["mouse_symbol"]
        for row in final_rows
        if row.get("mouse_symbol") and row.get("tier") == "1_primary"
    }
    modules["final24_all_candidates"] = final_all
    modules["final16_neurodevelopment_candidates"] = final_neurodevelopment
    modules["final_synapse_window_candidates"] = final_synapse_window
    modules["primary5_candidates"] = final_primary
    return dict(modules)


def module_test(
    rows: dict[str, dict[str, object]],
    module_name: str,
    genes: set[str],
    contrast: dict[str, str],
    min_module_size: int,
) -> dict[str, object] | None:
    field = contrast["score_field"]
    universe: list[tuple[str, float]] = []
    for symbol, row in rows.items():
        score = safe_float(row.get(field))
        if not math.isfinite(score):
            continue
        if field == "gse268951_log2fc_oe_vs_control" and safe_float(row.get("mean_fpkm_max")) < 0.5:
            continue
        universe.append((symbol, score))
    module_values = [score for symbol, score in universe if symbol in genes]
    background_values = [score for symbol, score in universe if symbol not in genes]
    effective_min = 3 if module_name.startswith("final") or module_name.startswith("primary") else min_module_size
    if len(module_values) < effective_min or len(background_values) < min_module_size:
        return None

    u_test = stats.mannwhitneyu(module_values, background_values, alternative="two-sided")
    up_test = stats.mannwhitneyu(module_values, background_values, alternative="greater")
    down_test = stats.mannwhitneyu(module_values, background_values, alternative="less")
    module_mean = float(np.mean(module_values))
    background_mean = float(np.mean(background_values))
    module_median = float(np.median(module_values))
    background_median = float(np.median(background_values))
    direction = "up_shift" if module_median > background_median else ("down_shift" if module_median < background_median else "flat")
    # Rank-biserial effect from Mann-Whitney U: positive means module ranks higher.
    n1 = len(module_values)
    n2 = len(background_values)
    rank_biserial = (2 * float(u_test.statistic) / (n1 * n2)) - 1
    return {
        "contrast": contrast["contrast"],
        "source": contrast["source"],
        "module": module_name,
        "module_gene_count": n1,
        "background_gene_count": n2,
        "module_mean_logfc": fmt(module_mean),
        "background_mean_logfc": fmt(background_mean),
        "delta_mean_logfc": fmt(module_mean - background_mean),
        "module_median_logfc": fmt(module_median),
        "background_median_logfc": fmt(background_median),
        "delta_median_logfc": fmt(module_median - background_median),
        "rank_biserial_effect": fmt(rank_biserial),
        "shift_direction": direction,
        "mannwhitney_p_two_sided": fmt(float(u_test.pvalue)),
        "mannwhitney_p_up": fmt(float(up_test.pvalue)),
        "mannwhitney_p_down": fmt(float(down_test.pvalue)),
        "leading_up_genes": leading_genes(universe, genes, up=True),
        "leading_down_genes": leading_genes(universe, genes, up=False),
    }


def leading_genes(universe: list[tuple[str, float]], genes: set[str], up: bool, n: int = 12) -> str:
    selected = [(symbol, score) for symbol, score in universe if symbol in genes]
    selected.sort(key=lambda item: item[1], reverse=up)
    return ",".join(f"{symbol}:{score:.3g}" for symbol, score in selected[:n])


def run_tests(
    rows: dict[str, dict[str, object]],
    modules: dict[str, set[str]],
    min_module_size: int,
) -> list[dict[str, object]]:
    test_rows: list[dict[str, object]] = []
    for contrast in CONTRASTS:
        contrast_rows: list[dict[str, object]] = []
        for module, genes in modules.items():
            result = module_test(rows, module, genes, contrast, min_module_size)
            if result:
                contrast_rows.append(result)
        qvalues = bh_fdr([safe_float(row["mannwhitney_p_two_sided"]) for row in contrast_rows])
        for row, qvalue in zip(contrast_rows, qvalues):
            row["mannwhitney_fdr_two_sided"] = fmt(qvalue)
            test_rows.append(row)
    return test_rows


def membership_rows(modules: dict[str, set[str]], rows: dict[str, dict[str, object]]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for module, genes in sorted(modules.items()):
        for symbol in sorted(genes):
            row = rows.get(symbol, {})
            out.append(
                {
                    "module": module,
                    "gene_symbol": symbol,
                    "description": row.get("description", ""),
                    "gse268951_log2fc_oe_vs_control": fmt(row.get("gse268951_log2fc_oe_vs_control")),
                    "dnmt3a_logfc_cko_vs_control": fmt(row.get("dnmt3a_logfc_cko_vs_control")),
                    "go_bp_terms": row.get("go_bp_terms", ""),
                }
            )
    return out


def make_dotplot(test_rows: list[dict[str, object]], figdir: Path) -> None:
    figdir.mkdir(parents=True, exist_ok=True)
    plot_modules = [
        "early_neurogenesis_neuron_differentiation",
        "cortical_patterning_region_identity",
        "axon_dendrite_neuron_projection",
        "synapse_formation_maturation",
        "glial_maturation_myelination",
        "neurodevelopment_signaling",
        "cell_cycle_progenitor_proliferation",
        "final16_neurodevelopment_candidates",
        "final_synapse_window_candidates",
        "primary5_candidates",
    ]
    labels = {
        "early_neurogenesis_neuron_differentiation": "Early neurogenesis /\nneuron differentiation",
        "cortical_patterning_region_identity": "Cortical patterning /\nregion identity",
        "axon_dendrite_neuron_projection": "Axon/dendrite /\nneuron projection",
        "synapse_formation_maturation": "Synapse formation /\nmaturation",
        "glial_maturation_myelination": "Glial maturation /\nmyelination",
        "neurodevelopment_signaling": "Neurodevelopment\nsignaling",
        "cell_cycle_progenitor_proliferation": "Cell cycle /\nprogenitor proliferation",
        "final16_neurodevelopment_candidates": "Final 16 neurodevelopment\ncandidates",
        "final_synapse_window_candidates": "Final synapse-window\ncandidates",
        "primary5_candidates": "Primary 5\ncandidates",
    }
    contrasts = [contrast["contrast"] for contrast in CONTRASTS]
    row_by_key = {(row["contrast"], row["module"]): row for row in test_rows}

    fig, ax = plt.subplots(figsize=(11.5, 5.8))
    for x, contrast in enumerate(contrasts):
        for y, module in enumerate(plot_modules):
            row = row_by_key.get((contrast, module))
            if not row:
                continue
            delta = safe_float(row["delta_median_logfc"])
            fdr = safe_float(row["mannwhitney_fdr_two_sided"])
            size = 45 + min(250, max(0, -math.log10(max(fdr, 1e-12))) * 38) if math.isfinite(fdr) else 45
            ax.scatter(
                x,
                y,
                s=size,
                c=delta,
                cmap="RdBu_r",
                vmin=-0.16,
                vmax=0.16,
                edgecolor="#374151",
                linewidth=0.5,
            )
            if math.isfinite(fdr) and fdr < 0.05:
                ax.text(x, y, "*", ha="center", va="center", fontsize=10, color="black", fontweight="bold")
    ax.set_xticks(range(len(contrasts)))
    ax.set_xticklabels(["Dnmt3l OE\nP28 mFC", "Dnmt3a cKO\nP39 neurons"], fontsize=9)
    ax.set_yticks(range(len(plot_modules)))
    ax.set_yticklabels([labels[module] for module in plot_modules], fontsize=8.6)
    ax.invert_yaxis()
    ax.set_title("Rank-based developmental module shifts", fontsize=13, fontweight="bold", pad=16)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.grid(axis="both", color="#E5E7EB", linewidth=0.8)
    mappable = plt.cm.ScalarMappable(cmap="RdBu_r")
    mappable.set_clim(-0.16, 0.16)
    cbar = fig.colorbar(mappable, ax=ax, fraction=0.035, pad=0.035)
    cbar.set_label("Module median logFC - background median logFC", fontsize=8.5)
    cbar.ax.tick_params(labelsize=8)
    fig.text(0.14, 0.035, "Dot size scales with -log10(FDR). Asterisk marks FDR < 0.05.", fontsize=8.3, color="#4B5563")
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    fig.savefig(figdir / "developmental_module_shift_dotplot.png", dpi=260, bbox_inches="tight")
    fig.savefig(figdir / "developmental_module_shift_dotplot.pdf", bbox_inches="tight")
    plt.close(fig)


def write_note(path: Path, test_rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    def row_for(contrast_contains: str, module: str) -> dict[str, object] | None:
        for row in test_rows:
            if contrast_contains in str(row["contrast"]) and row["module"] == module:
                return row
        return None

    syn_d3l = row_for("Dnmt3l_OE", "synapse_formation_maturation")
    axon_d3l = row_for("Dnmt3l_OE", "axon_dendrite_neuron_projection")
    neuro_d3l = row_for("Dnmt3l_OE", "early_neurogenesis_neuron_differentiation")
    syn_d3a = row_for("Dnmt3a_cKO", "synapse_formation_maturation")
    final16_d3l = row_for("Dnmt3l_OE", "final16_neurodevelopment_candidates")

    def brief(row: dict[str, object] | None) -> str:
        if not row:
            return "not tested"
        return (
            f"delta median logFC {row['delta_median_logfc']}, "
            f"rank effect {row['rank_biserial_effect']}, FDR {row['mannwhitney_fdr_two_sided']}"
        )

    lines = [
        "# Developmental Module Shift Tests",
        "",
        "Date: 2026-06-09",
        "",
        "This analysis follows up the observation that most final neurodevelopmental candidates do not pass single-gene DE thresholds in P28/P39 perturbation data. Because these ages are closer to neuronal maturation, synapse formation, and circuit refinement than early progenitor neurogenesis, the test uses rank-based module shifts for developmental and synapse-window gene sets.",
        "",
        "## Method",
        "",
        "Gene modules were built from GSE268951 GO biological-process annotations using keyword groups for early neurogenesis/neuron differentiation, cortical patterning, axon/dendrite/neuron projection, synapse formation/maturation, glial maturation/myelination, neurodevelopmental signaling, and cell-cycle/progenitor proliferation. Final-candidate modules were also tested: all 24 candidates, the 16 neurodevelopment-related candidates, synapse-window candidates, and the primary five candidates.",
        "",
        "For each module, gene logFC values were compared with the background gene universe by a two-sided Mann-Whitney rank test. This asks whether module genes shift as a group even if individual genes do not pass hard DE cutoffs.",
        "",
        "## Key first-pass readouts",
        "",
        f"- Dnmt3l OE synapse formation/maturation module: {brief(syn_d3l)}.",
        f"- Dnmt3l OE axon/dendrite/neuron projection module: {brief(axon_d3l)}.",
        f"- Dnmt3l OE early neurogenesis/neuron differentiation module: {brief(neuro_d3l)}.",
        f"- Dnmt3a cKO synapse formation/maturation module: {brief(syn_d3a)}.",
        f"- Dnmt3l OE final 16 neurodevelopment-candidate module: {brief(final16_d3l)}.",
        "",
        "## Interpretation",
        "",
        "This module-level analysis is the right test for the P28/P39 developmental window. A lack of strong single-gene DE among final candidates does not mean neurodevelopment is unaffected. Instead, the question becomes whether synapse, projection, maturation, or signaling modules shift as a coordinated distribution.",
        "",
        "The results should be interpreted as a ranked transcriptomic module screen. They do not prove direct methylation mediation, and they remain sensitive to differences between Dnmt3l overexpression in medial frontal cortex and Dnmt3a cKO in excitatory neurons.",
        "",
        "## Generated outputs",
        "",
        "- `results/neurogenesis_mechanism/developmental_module_rank_tests.tsv`",
        "- `results/neurogenesis_mechanism/developmental_module_gene_membership.tsv`",
        "- `figures/neurogenesis_mechanism/developmental_module_shift_dotplot.png`",
        "- `figures/neurogenesis_mechanism/developmental_module_shift_dotplot.pdf`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    figdir = Path(args.figdir)
    rows = load_expression_rows(Path(args.expr), Path(args.overlap))
    final_rows = read_tsv(Path(args.final))
    modules = build_modules(rows, final_rows)
    tests = run_tests(rows, modules, args.min_module_size)
    memberships = membership_rows(modules, rows)

    test_fields = [
        "contrast",
        "source",
        "module",
        "module_gene_count",
        "background_gene_count",
        "module_mean_logfc",
        "background_mean_logfc",
        "delta_mean_logfc",
        "module_median_logfc",
        "background_median_logfc",
        "delta_median_logfc",
        "rank_biserial_effect",
        "shift_direction",
        "mannwhitney_p_two_sided",
        "mannwhitney_fdr_two_sided",
        "mannwhitney_p_up",
        "mannwhitney_p_down",
        "leading_up_genes",
        "leading_down_genes",
    ]
    membership_fields = [
        "module",
        "gene_symbol",
        "description",
        "gse268951_log2fc_oe_vs_control",
        "dnmt3a_logfc_cko_vs_control",
        "go_bp_terms",
    ]
    write_tsv(outdir / "developmental_module_rank_tests.tsv", tests, test_fields)
    write_tsv(outdir / "developmental_module_gene_membership.tsv", memberships, membership_fields)
    make_dotplot(tests, figdir)
    write_note(Path(args.note), tests)

    print(f"Wrote {outdir / 'developmental_module_rank_tests.tsv'}")
    print(f"Wrote {outdir / 'developmental_module_gene_membership.tsv'}")
    print(f"Wrote {figdir / 'developmental_module_shift_dotplot.png'}")
    print(f"Wrote {figdir / 'developmental_module_shift_dotplot.pdf'}")
    print(f"Wrote {args.note}")


if __name__ == "__main__":
    main()
