#!/usr/bin/env python3
"""First-pass Drosophila DNMT3L-OE expression/pathway analysis."""

from __future__ import annotations

import gzip
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy import stats


GSE263069_EXPR = Path("data/geo/GSE263069/supplementary/GSE263069_all.fpkm_anno.txt.gz")
GSE68861_NORM = Path("data/geo/GSE68861/supplementary/GSE68861_normalized.txt.gz")
GSE68861_SOFT = Path("data/geo/GSE68861/metadata/GPL17506_family.soft.gz")
OUTDIR = Path("results/drosophila_dnmt3l")
FIGDIR = Path("figures/drosophila_dnmt3l")
NOTE = Path("notes/drosophila_dnmt3l_analysis_first_pass.md")

CTRL_FPKM = ["C1_FPKM", "C2_FPKM", "C3_FPKM"]
OE_FPKM = ["D1_FPKM", "D2_FPKM", "D3_FPKM"]
CTRL_COUNTS = ["C1", "C2", "C3"]
OE_COUNTS = ["D1", "D2", "D3"]

NOMINAL_P = 0.05
LOG2FC_CUTOFF = 1.0


MODULE_KEYWORDS = {
    "neural development": r"neuro|nervous system|axon|dendrit|mushroom body|brain|synap|glia|photoreceptor",
    "synapse / neurotransmission": r"synap|neurotrans|cholinergic|serotonergic|glutamate|GABA|neuromuscular|axon guidance",
    "ion transport / calcium": r"calcium|voltage[- ]gated|ion transport|ion channel|potassium|sodium|calmodulin|transmembrane transport",
    "chromatin / transcription": r"chromatin|histone|DNA binding|transcription|RNA polymerase|enhancer|nucleus",
    "cell cycle / growth": r"cell cycle|mitotic|proliferation|growth|cell fate|apoptotic|differentiation",
    "stress / immune": r"stress|immune|defense|heat|oxidative|antimicrobial|response to",
    "metabolism / detox": r"metabolic|metabolism|drug|oxidation|monooxygenase|cytochrome|glutathione|hydrolase",
}

SELECTED_MARKERS = [
    ("JAK/STAT", "Stat92E"),
    ("JAK/STAT", "hop"),
    ("JAK/STAT", "upd1"),
    ("JAK/STAT", "upd2"),
    ("JAK/STAT", "upd3"),
    ("JAK/STAT", "Socs36E"),
    ("JAK/STAT", "Socs44A"),
    ("neurogenic bHLH/Notch", "ac"),
    ("neurogenic bHLH/Notch", "sc"),
    ("neurogenic bHLH/Notch", "ase"),
    ("neurogenic bHLH/Notch", "ato"),
    ("neurogenic bHLH/Notch", "da"),
    ("neurogenic bHLH/Notch", "sens"),
    ("neurogenic bHLH/Notch", "pros"),
    ("neurogenic bHLH/Notch", "elav"),
    ("neurogenic bHLH/Notch", "dpn"),
    ("neurogenic bHLH/Notch", "N"),
    ("neurogenic bHLH/Notch", "Dl"),
    ("neurogenic bHLH/Notch", "E(spl)m8-HLH"),
    ("neurogenic bHLH/Notch", "H"),
    ("neurogenic bHLH/Notch", "neur"),
    ("axon/synapse", "Abl"),
]

GSE68861_CONDITION_MAP = {
    "Drosophila_1": "Tub-3L G1 DNMT3L-OE",
    "Drosophila_2": "UAS-3L G1 control",
    "Drosophila_3": "Tub-3L G5 phenotype DNMT3L-OE",
    "Drosophila_4": "w1118 wild-type",
    "Drosophila_6": "Tub-3L G5 DNMT3L-OE",
    "Drosophila_7": "Tub-3L G4 DNMT3L-OE",
    "Drosophila_8": "Tub-3L G2 DNMT3L-OE",
    "Drosophila_9": "UAS-3L G5 control",
    "Drosophila_10": "Tub-3L G20 driver-removed",
}


def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    p = np.where(np.isfinite(p), p, 1.0)
    n = len(p)
    order = np.argsort(p)
    ranked = np.empty(n, dtype=float)
    ranked[order] = p[order] * n / np.arange(1, n + 1)
    ranked_ordered = ranked[order]
    ranked_ordered = np.minimum.accumulate(ranked_ordered[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(ranked_ordered, 0, 1)
    return out


def parse_go_terms(value: object) -> list[tuple[str, str]]:
    if pd.isna(value):
        return []
    text = str(value)
    if text == "-":
        return []
    out = []
    for item in text.split(";"):
        item = item.strip()
        match = re.match(r"(GO:\d+):(.+)", item)
        if match:
            out.append((match.group(1), match.group(2).strip()))
    return out


def read_gse263069() -> pd.DataFrame:
    df = pd.read_csv(GSE263069_EXPR, sep="\t")
    ctrl = np.log2(df[CTRL_FPKM].astype(float) + 0.1).to_numpy()
    oe = np.log2(df[OE_FPKM].astype(float) + 0.1).to_numpy()
    pvalues = []
    for idx in range(len(df)):
        with np.errstate(all="ignore"):
            res = stats.ttest_ind(oe[idx], ctrl[idx], equal_var=False, nan_policy="omit")
        pvalues.append(res.pvalue if np.isfinite(res.pvalue) else 1.0)

    df["control_mean_fpkm"] = df[CTRL_FPKM].mean(axis=1)
    df["dnmt3l_oe_mean_fpkm"] = df[OE_FPKM].mean(axis=1)
    df["control_mean_count"] = df[CTRL_COUNTS].mean(axis=1)
    df["dnmt3l_oe_mean_count"] = df[OE_COUNTS].mean(axis=1)
    df["log2fc_dnmt3l_oe_vs_control"] = np.log2(df["dnmt3l_oe_mean_fpkm"] + 0.1) - np.log2(
        df["control_mean_fpkm"] + 0.1
    )
    df["welch_pvalue_log2fpkm"] = pvalues
    df["welch_fdr_log2fpkm"] = bh_fdr(np.array(pvalues))
    df["expressed"] = df[CTRL_FPKM + OE_FPKM].sum(axis=1).gt(0)
    df["nominal_de"] = df["welch_pvalue_log2fpkm"].le(NOMINAL_P) & df[
        "log2fc_dnmt3l_oe_vs_control"
    ].abs().ge(LOG2FC_CUTOFF)
    df["direction"] = np.where(
        df["log2fc_dnmt3l_oe_vs_control"].ge(LOG2FC_CUTOFF),
        "up",
        np.where(df["log2fc_dnmt3l_oe_vs_control"].le(-LOG2FC_CUTOFF), "down", "flat"),
    )
    df["de_direction"] = np.where(df["nominal_de"], df["direction"], "not_nominal_de")
    return df


def term_enrichment(df: pd.DataFrame, column: str, label: str) -> pd.DataFrame:
    bg = df.loc[df["expressed"]].copy()
    bg_ids = set(bg["gene_id"])
    term_to_genes: dict[tuple[str, str], set[str]] = {}
    for _, row in bg.iterrows():
        for term_id, term_name in parse_go_terms(row[column]):
            term_to_genes.setdefault((term_id, term_name), set()).add(row["gene_id"])

    sets = {
        "nominal_de_all": set(bg.loc[bg["nominal_de"], "gene_id"]),
        "nominal_de_up": set(bg.loc[bg["nominal_de"] & bg["direction"].eq("up"), "gene_id"]),
        "nominal_de_down": set(bg.loc[bg["nominal_de"] & bg["direction"].eq("down"), "gene_id"]),
    }
    rows = []
    M = len(bg_ids)
    for set_name, query in sets.items():
        n = len(query)
        if n == 0:
            continue
        for (term_id, term_name), genes in term_to_genes.items():
            K = len(genes)
            k = len(query & genes)
            if K < 10 or k < 2:
                continue
            p = stats.hypergeom.sf(k - 1, M, K, n)
            rows.append(
                {
                    "annotation_namespace": label,
                    "query_set": set_name,
                    "term_id": term_id,
                    "term_name": term_name,
                    "background_genes": M,
                    "query_genes": n,
                    "term_background_genes": K,
                    "term_query_genes": k,
                    "fold_enrichment": (k / n) / (K / M),
                    "pvalue": p,
                    "example_genes": ";".join(
                        bg.loc[bg["gene_id"].isin(sorted(query & genes)), "GeneSymbol"].head(20).astype(str)
                    ),
                }
            )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["fdr"] = out.groupby(["annotation_namespace", "query_set"])["pvalue"].transform(
            lambda s: bh_fdr(s.to_numpy())
        )
        out = out.sort_values(["query_set", "fdr", "pvalue", "term_query_genes"], ascending=[True, True, True, False])
    return out


def module_enrichment(df: pd.DataFrame) -> pd.DataFrame:
    bg = df.loc[df["expressed"]].copy()
    query = bg["nominal_de"].to_numpy()
    up = bg["nominal_de"].to_numpy() & bg["direction"].eq("up").to_numpy()
    down = bg["nominal_de"].to_numpy() & bg["direction"].eq("down").to_numpy()
    searchable = (
        bg[["Description", "GO_BP", "GO_CC", "GO_MF"]].fillna("").agg(" ".join, axis=1).str.replace("%3B", ";")
    )
    rows = []
    for module, pattern in MODULE_KEYWORDS.items():
        hits = searchable.str.contains(pattern, case=False, regex=True)
        table = [[int((query & hits).sum()), int((query & ~hits).sum())], [int((~query & hits).sum()), int((~query & ~hits).sum())]]
        odds, p = stats.fisher_exact(table, alternative="greater")
        rows.append(
            {
                "module": module,
                "background_genes": len(bg),
                "module_background_genes": int(hits.sum()),
                "nominal_de_genes": int(query.sum()),
                "module_nominal_de_genes": int((query & hits).sum()),
                "module_nominal_de_up": int((up & hits).sum()),
                "module_nominal_de_down": int((down & hits).sum()),
                "fold_enrichment": (int((query & hits).sum()) / max(int(query.sum()), 1))
                / (int(hits.sum()) / max(len(bg), 1)),
                "fisher_oddsratio": odds,
                "pvalue": p,
                "example_nominal_de_genes": ";".join(
                    bg.loc[query & hits]
                    .sort_values("welch_pvalue_log2fpkm")["GeneSymbol"]
                    .head(25)
                    .astype(str)
                ),
            }
        )
    out = pd.DataFrame(rows)
    out["fdr"] = bh_fdr(out["pvalue"].to_numpy())
    return out.sort_values(["fdr", "pvalue", "module_nominal_de_genes"], ascending=[True, True, False])


def selected_marker_effects(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    lower_to_rows = {}
    for idx, symbol in df["GeneSymbol"].astype(str).str.lower().items():
        lower_to_rows.setdefault(symbol, []).append(idx)
    for group, symbol in SELECTED_MARKERS:
        for idx in lower_to_rows.get(symbol.lower(), []):
            row = df.loc[idx]
            rows.append(
                {
                    "marker_group": group,
                    "requested_symbol": symbol,
                    "gene_id": row["gene_id"],
                    "gene_symbol": row["GeneSymbol"],
                    "description": row["Description"],
                    "control_mean_fpkm": row["control_mean_fpkm"],
                    "dnmt3l_oe_mean_fpkm": row["dnmt3l_oe_mean_fpkm"],
                    "log2fc_dnmt3l_oe_vs_control": row["log2fc_dnmt3l_oe_vs_control"],
                    "welch_pvalue_log2fpkm": row["welch_pvalue_log2fpkm"],
                    "welch_fdr_log2fpkm": row["welch_fdr_log2fpkm"],
                    "nominal_de": row["nominal_de"],
                }
            )
    return pd.DataFrame(rows).sort_values(["marker_group", "requested_symbol", "gene_id"])


def read_platform_map() -> pd.DataFrame:
    if not GSE68861_SOFT.exists():
        return pd.DataFrame()
    rows = []
    in_table = False
    with gzip.open(GSE68861_SOFT, "rt", errors="replace") as handle:
        header = None
        for line in handle:
            if line.startswith("!platform_table_begin"):
                in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if in_table:
                fields = line.rstrip("\n").split("\t")
                if header is None:
                    header = fields
                else:
                    rows.append(fields)
    if not rows or header is None:
        return pd.DataFrame()
    return pd.DataFrame(rows, columns=header)


def analyze_gse68861() -> tuple[pd.DataFrame, pd.DataFrame]:
    norm = pd.read_csv(GSE68861_NORM, sep="\t", comment="#")
    norm = norm.loc[:, ~norm.columns.str.startswith("Unnamed")]
    norm["ID_REF"] = norm["ID_REF"].astype(str)
    for col in norm.columns:
        if col != "ID_REF":
            norm[col] = pd.to_numeric(norm[col], errors="coerce")

    contrasts = {
        "G1_Tub3L_vs_UAS3L": ("Drosophila_1", "Drosophila_2"),
        "G5_Tub3L_vs_UAS3L": ("Drosophila_6", "Drosophila_9"),
        "G5_phenotype_vs_G5_nonphenotype": ("Drosophila_3", "Drosophila_6"),
        "G20_driver_removed_vs_G5": ("Drosophila_10", "Drosophila_6"),
        "G5_vs_G1_Tub3L_progression": ("Drosophila_6", "Drosophila_1"),
    }
    rows = []
    delta_df = norm[["ID_REF"]].copy()
    for name, (a, b) in contrasts.items():
        delta = norm[a] - norm[b]
        delta_df[name] = delta
        rows.append(
            {
                "contrast": name,
                "numerator_condition": GSE68861_CONDITION_MAP.get(a, a),
                "denominator_condition": GSE68861_CONDITION_MAP.get(b, b),
                "n_probes": int(delta.notna().sum()),
                "median_delta_log2_signal": float(delta.median()),
                "mean_delta_log2_signal": float(delta.mean()),
                "n_abs_delta_ge_1": int(delta.abs().ge(1).sum()),
                "n_up_delta_ge_1": int(delta.ge(1).sum()),
                "n_down_delta_le_minus_1": int(delta.le(-1).sum()),
                "p95_abs_delta": float(delta.abs().quantile(0.95)),
            }
        )
    platform = read_platform_map()
    if not platform.empty:
        platform = platform.rename(columns={"ID": "ID_REF"})
        platform["ID_REF"] = platform["ID_REF"].astype(str)
        delta_df = delta_df.merge(platform[["ID_REF", "ORF", "SPOT_ID", "seqname", "start", "stop"]], on="ID_REF", how="left")
    return pd.DataFrame(rows), delta_df


def plot_volcano(df: pd.DataFrame) -> None:
    plot_df = df.loc[df["expressed"]].copy()
    colors = np.where(
        plot_df["nominal_de"] & plot_df["direction"].eq("up"),
        "#D45750",
        np.where(plot_df["nominal_de"] & plot_df["direction"].eq("down"), "#3F7CAC", "#C8C8C8"),
    )
    y = -np.log10(plot_df["welch_pvalue_log2fpkm"].clip(lower=1e-300))
    y_plot = y.clip(upper=8)
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.scatter(plot_df["log2fc_dnmt3l_oe_vs_control"], y_plot, c=colors, s=9, alpha=0.75, linewidths=0)
    ax.axvline(LOG2FC_CUTOFF, color="0.35", lw=0.8, ls="--")
    ax.axvline(-LOG2FC_CUTOFF, color="0.35", lw=0.8, ls="--")
    ax.axhline(-math.log10(NOMINAL_P), color="0.35", lw=0.8, ls="--")
    ax.set_xlabel("log2FC DNMT3L-OE / control, FPKM")
    ax.set_ylabel("-log10 nominal P, Welch test (display capped at 8)")
    ax.set_title("GSE263069 Drosophila larval brain DNMT3L-OE")
    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#D45750", markersize=6, label="nominal up"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#3F7CAC", markersize=6, label="nominal down"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#C8C8C8", markersize=6, label="other"),
    ]
    ax.legend(handles=handles, frameon=False, loc="upper right")
    out = FIGDIR / "gse263069_volcano"
    fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_go(go: pd.DataFrame) -> None:
    if go.empty:
        return
    top = (
        go.loc[go["query_set"].isin(["nominal_de_up", "nominal_de_down"])]
        .sort_values(["fdr", "pvalue", "fold_enrichment"])
        .groupby("query_set")
        .head(10)
        .copy()
    )
    if top.empty:
        return
    top["label"] = top["term_name"].str.slice(0, 64) + np.where(top["term_name"].str.len().gt(64), "...", "")
    top["plot_value"] = -np.log10(top["pvalue"].clip(lower=1e-300))
    top["signed_value"] = np.where(top["query_set"].eq("nominal_de_up"), top["plot_value"], -top["plot_value"])
    top = top.sort_values("signed_value")
    fig, ax = plt.subplots(figsize=(8.5, max(5, 0.34 * len(top))))
    colors = np.where(top["signed_value"].gt(0), "#D45750", "#3F7CAC")
    ax.barh(np.arange(len(top)), top["signed_value"], color=colors)
    ax.axvline(0, color="0.2", lw=0.8)
    ax.set_yticks(np.arange(len(top)), top["label"])
    ax.set_xlabel("signed -log10 enrichment P (right = up, left = down)")
    ax.set_title("Top GO biological-process enrichments among nominal DE genes")
    out = FIGDIR / "gse263069_go_bp_enrichment"
    fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_modules(modules: pd.DataFrame) -> None:
    if modules.empty:
        return
    m = modules.sort_values("fold_enrichment")
    fig, ax = plt.subplots(figsize=(8, 4.8))
    y = np.arange(len(m))
    sizes = 30 + m["module_nominal_de_genes"] * 7
    color = -np.log10(m["pvalue"].clip(lower=1e-300))
    sc = ax.scatter(m["fold_enrichment"], y, s=sizes, c=color, cmap="viridis", edgecolor="0.25", linewidth=0.3)
    ax.axvline(1, color="0.4", lw=0.8, ls="--")
    ax.set_yticks(y, m["module"])
    ax.set_xlabel("Nominal DE fold enrichment over expressed-gene background")
    ax.set_title("Keyword module enrichment in GSE263069 nominal DE genes")
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("-log10 Fisher P")
    out = FIGDIR / "gse263069_module_enrichment"
    fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_markers(markers: pd.DataFrame) -> None:
    if markers.empty:
        return
    m = markers.sort_values(["marker_group", "log2fc_dnmt3l_oe_vs_control"])
    fig, ax = plt.subplots(figsize=(8, max(5, 0.28 * len(m))))
    colors = np.where(m["marker_group"].eq("JAK/STAT"), "#59A14F", np.where(m["marker_group"].eq("axon/synapse"), "#8E6C8A", "#4C78A8"))
    ax.barh(np.arange(len(m)), m["log2fc_dnmt3l_oe_vs_control"], color=colors)
    ax.axvline(0, color="0.2", lw=0.8)
    ax.axvline(1, color="0.5", lw=0.7, ls="--")
    ax.axvline(-1, color="0.5", lw=0.7, ls="--")
    labels = m["gene_symbol"].astype(str) + " (" + m["marker_group"].astype(str) + ")"
    ax.set_yticks(np.arange(len(m)), labels)
    ax.set_xlabel("log2FC DNMT3L-OE / control, FPKM")
    ax.set_title("Selected fly JAK/STAT and neurogenic marker effects")
    out = FIGDIR / "gse263069_selected_marker_effects"
    fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_gse68861(summary: pd.DataFrame) -> None:
    if summary.empty:
        return
    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    x = np.arange(len(summary))
    ax.bar(x, summary["n_abs_delta_ge_1"], color="#8E6C8A")
    ax.set_xticks(x, summary["contrast"].str.replace("_", "\n"), rotation=0)
    ax.set_ylabel("Probes with |delta log2 signal| >= 1")
    ax.set_title("GSE68861 whole-larva array condition-level shifts")
    for i, row in summary.iterrows():
        ax.text(i, row["n_abs_delta_ge_1"] + max(summary["n_abs_delta_ge_1"]) * 0.02, str(int(row["n_abs_delta_ge_1"])), ha="center", fontsize=8)
    out = FIGDIR / "gse68861_probe_level_contrast_summary"
    fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def write_note(df: pd.DataFrame, modules: pd.DataFrame, markers: pd.DataFrame, gse68861_summary: pd.DataFrame) -> None:
    expressed = int(df["expressed"].sum())
    de = df.loc[df["nominal_de"]]
    up = int((de["direction"] == "up").sum())
    down = int((de["direction"] == "down").sum())
    max_fdr_hits = int((df["welch_fdr_log2fpkm"].le(0.05) & df["log2fc_dnmt3l_oe_vs_control"].abs().ge(1)).sum())
    top_modules = modules.sort_values(["pvalue", "module_nominal_de_genes"]).head(5)
    fdr_modules = modules.loc[modules["fdr"].le(0.05)]
    top_go_bp = pd.read_csv(OUTDIR / "gse263069_go_enrichment.tsv", sep="\t")
    top_go_bp = (
        top_go_bp.loc[top_go_bp["annotation_namespace"].eq("GO_BP") & top_go_bp["query_set"].eq("nominal_de_all")]
        .sort_values(["pvalue", "fdr"])
        .head(8)
    )
    marker_hits = markers.loc[markers["log2fc_dnmt3l_oe_vs_control"].abs().ge(1)].copy()
    gse68861_max = gse68861_summary.sort_values("n_abs_delta_ge_1", ascending=False).head(1)

    lines = [
        "# Drosophila DNMT3L-OE First-Pass Analysis",
        "",
        "Date: 2026-06-11",
        "",
        "## Datasets",
        "",
        "- `GSE263069`: third instar larval brain RNA-seq, 3 elavC155-Gal4 controls versus 3 human DNMT3L-OE larval brain samples.",
        "- `GSE68861`: older whole-larva Affymetrix array series; public normalized file contains condition-level averaged log2 signals, so it is treated as secondary context.",
        "",
        "## Methods",
        "",
        f"`GSE263069` was analyzed from the public processed FPKM/count annotation table. For each gene, log2FC was computed from mean FPKM as log2(mean DNMT3L-OE + 0.1) - log2(mean control + 0.1). A Welch test was run on log2(FPKM + 0.1) values across the three biological replicates per group. Because this processed-table reanalysis does not use DESeq2 raw-count dispersion modeling, the main exploratory set was defined as nominal P <= 0.05 and |log2FC| >= 1. GO enrichment used the supplied GO annotations and a hypergeometric test over expressed genes.",
        "",
        "## Results",
        "",
        f"- Expressed genes detected in `GSE263069`: {expressed:,}.",
        f"- Nominal first-pass DE genes at P <= {NOMINAL_P} and |log2FC| >= {LOG2FC_CUTOFF}: {len(de):,} ({up} up, {down} down). This closely matches the scale of the GEO summary but should not be called a replacement for the authors' DESeq2 q-value result.",
        f"- FDR-significant genes from this simple processed-FPKM Welch test: {max_fdr_hits}. Therefore FDR-level gene claims should rely on the authors' deposited summary, while this project should use the reanalysis for direction and pathway context.",
    ]
    if fdr_modules.empty:
        lines.append(
            "- No broad keyword module passed FDR <= 0.05 in this project-level enrichment test. This is important: the Drosophila result should be used as contextual pathway plausibility, not a strong independent module claim."
        )
    if not top_modules.empty:
        lines.append("- Broad keyword module counts/trends among nominal DE genes:")
        for _, row in top_modules.iterrows():
            lines.append(
                f"  - {row['module']}: {int(row['module_nominal_de_genes'])} nominal DE genes, fold enrichment {row['fold_enrichment']:.2f}, Fisher P {row['pvalue']:.2g}."
            )
    if not top_go_bp.empty:
        lines.append("- Top nominal GO biological-process enrichment trends:")
        for _, row in top_go_bp.iterrows():
            lines.append(
                f"  - {row['term_name']}: {int(row['term_query_genes'])} nominal DE genes, fold enrichment {row['fold_enrichment']:.2f}, P {row['pvalue']:.2g}, FDR {row['fdr']:.2g}."
            )
    if not marker_hits.empty:
        lines.append("- Selected fly marker genes with |log2FC| >= 1 include:")
        for _, row in marker_hits.sort_values("log2fc_dnmt3l_oe_vs_control", key=lambda s: s.abs(), ascending=False).iterrows():
            lines.append(
                f"  - `{row['gene_symbol']}` ({row['marker_group']}): log2FC {row['log2fc_dnmt3l_oe_vs_control']:.2f}, nominal P {row['welch_pvalue_log2fpkm']:.2g}."
            )
    if not gse68861_max.empty:
        row = gse68861_max.iloc[0]
        lines.append(
            f"- In `GSE68861`, the largest condition-level probe shift count was `{row['contrast']}` with {int(row['n_abs_delta_ge_1'])} probes at |delta log2 signal| >= 1."
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The Drosophila data are useful as methylation-independent pathway context. The safest conclusion is that human DNMT3L expression in a fly neural developmental system produces a transcriptomic shift of approximately the same size as reported by GEO, and the shifted genes include immune/defense, development, memory, ion-transport, glutamate-receptor, neuropeptide, and stress/metabolic annotations at nominal levels. However, this project-level processed-FPKM reanalysis does not produce FDR-robust single-gene or broad keyword-module claims. Therefore the dataset supports plausibility for a methylation-independent regulatory/stress arm, but it should not be treated as proof that the same genes or mechanisms operate in human NPCs or mouse cortex.",
            "",
            "## Generated Outputs",
            "",
            "- `results/drosophila_dnmt3l/gse263069_gene_expression_effects.tsv`",
            "- `results/drosophila_dnmt3l/gse263069_nominal_de_genes.tsv`",
            "- `results/drosophila_dnmt3l/gse263069_go_enrichment.tsv`",
            "- `results/drosophila_dnmt3l/gse263069_module_enrichment.tsv`",
            "- `results/drosophila_dnmt3l/gse263069_selected_fly_marker_effects.tsv`",
            "- `results/drosophila_dnmt3l/gse68861_probe_level_contrast_summary.tsv`",
            "- `figures/drosophila_dnmt3l/gse263069_volcano.png`",
            "- `figures/drosophila_dnmt3l/gse263069_go_bp_enrichment.png`",
            "- `figures/drosophila_dnmt3l/gse263069_module_enrichment.png`",
            "- `figures/drosophila_dnmt3l/gse263069_selected_marker_effects.png`",
            "- `figures/drosophila_dnmt3l/gse68861_probe_level_contrast_summary.png`",
            "",
        ]
    )
    NOTE.write_text("\n".join(lines))


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    FIGDIR.mkdir(parents=True, exist_ok=True)

    df = read_gse263069()
    keep_cols = [
        "gene_id",
        "GeneSymbol",
        "Description",
        "Chr",
        "Start",
        "End",
        "Strand",
        "control_mean_fpkm",
        "dnmt3l_oe_mean_fpkm",
        "control_mean_count",
        "dnmt3l_oe_mean_count",
        "log2fc_dnmt3l_oe_vs_control",
        "welch_pvalue_log2fpkm",
        "welch_fdr_log2fpkm",
        "expressed",
        "nominal_de",
        "direction",
        "de_direction",
        "GO_BP",
        "GO_CC",
        "GO_MF",
        "KO",
    ]
    df[keep_cols].sort_values("welch_pvalue_log2fpkm").to_csv(
        OUTDIR / "gse263069_gene_expression_effects.tsv", sep="\t", index=False
    )
    df.loc[df["nominal_de"], keep_cols].sort_values("welch_pvalue_log2fpkm").to_csv(
        OUTDIR / "gse263069_nominal_de_genes.tsv", sep="\t", index=False
    )

    go_bp = term_enrichment(df, "GO_BP", "GO_BP")
    go_mf = term_enrichment(df, "GO_MF", "GO_MF")
    go_cc = term_enrichment(df, "GO_CC", "GO_CC")
    go = pd.concat([go_bp, go_mf, go_cc], ignore_index=True)
    go.to_csv(OUTDIR / "gse263069_go_enrichment.tsv", sep="\t", index=False)

    modules = module_enrichment(df)
    modules.to_csv(OUTDIR / "gse263069_module_enrichment.tsv", sep="\t", index=False)

    markers = selected_marker_effects(df)
    markers.to_csv(OUTDIR / "gse263069_selected_fly_marker_effects.tsv", sep="\t", index=False)

    gse68861_summary, gse68861_delta = analyze_gse68861()
    gse68861_summary.to_csv(OUTDIR / "gse68861_probe_level_contrast_summary.tsv", sep="\t", index=False)
    gse68861_delta.to_csv(OUTDIR / "gse68861_probe_level_contrasts.tsv", sep="\t", index=False)

    summary_rows = [
        {
            "dataset": "GSE263069",
            "analysis_layer": "processed larval brain RNA-seq",
            "primary_use": "methylation-independent neural pathway comparator",
            "n_rows": len(df),
            "n_expressed": int(df["expressed"].sum()),
            "n_nominal_de": int(df["nominal_de"].sum()),
            "n_nominal_up": int((df["nominal_de"] & df["direction"].eq("up")).sum()),
            "n_nominal_down": int((df["nominal_de"] & df["direction"].eq("down")).sum()),
            "n_simple_fdr_de": int((df["welch_fdr_log2fpkm"].le(0.05) & df["log2fc_dnmt3l_oe_vs_control"].abs().ge(1)).sum()),
            "interpretation": "Use for pathway context; do not treat as mammalian methylation evidence.",
        },
        {
            "dataset": "GSE68861",
            "analysis_layer": "whole-larva array condition-level normalized signals",
            "primary_use": "lower-priority methylation-independent generational comparator",
            "n_rows": int(gse68861_summary["n_probes"].max()) if not gse68861_summary.empty else 0,
            "n_expressed": "",
            "n_nominal_de": "",
            "n_nominal_up": "",
            "n_nominal_down": "",
            "n_simple_fdr_de": "",
            "interpretation": "Use only as secondary context because public normalized file lacks per-replicate gene-level contrasts.",
        },
    ]
    pd.DataFrame(summary_rows).to_csv(OUTDIR / "drosophila_dnmt3l_first_pass_summary.tsv", sep="\t", index=False)

    plot_volcano(df)
    plot_go(go_bp)
    plot_modules(modules)
    plot_markers(markers)
    plot_gse68861(gse68861_summary)
    write_note(df, modules, markers, gse68861_summary)

    print(OUTDIR / "drosophila_dnmt3l_first_pass_summary.tsv")
    print(FIGDIR / "gse263069_volcano.png")
    print(NOTE)


if __name__ == "__main__":
    main()
