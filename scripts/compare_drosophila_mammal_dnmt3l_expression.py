#!/usr/bin/env python3
"""Compare DNMT3L-OE expression responses across fly, human, and mouse datasets."""

from __future__ import annotations

import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd
from scipy import stats


HUMAN_EXPR = Path("data/geo/GSE157080/supplementary/GSE157080_all.fpkm_anno.txt.gz")
HUMAN_AUTHOR_DE = Path("data/geo/GSE157080/supplementary/GSE157080_Table_S3_author_count_de.tsv")
MOUSE_EXPR_EFFECTS = Path("results/neurogenesis_mechanism/gse268951_rederived_gene_expression.tsv")
FLY_EXPR_EFFECTS = Path("results/drosophila_dnmt3l/gse263069_gene_expression_effects.tsv")
FLY_COUNT_EFFECTS = Path("results/drosophila_dnmt3l_count_de/gse263069_deseq2_edger_count_de_results.tsv")

OUTDIR = Path("results/cross_species_dnmt3l_expression_conservation")
FIGDIR = Path("figures/cross_species_dnmt3l_expression_conservation")
NOTE = Path("notes/cross_species_dnmt3l_expression_conservation_first_pass.md")

PSEUDOCOUNT = 0.1
RELAXED_LOG2FC = 0.58
LARGE_LOG2FC = 1.0
PVALUE_CUTOFF = 0.05

DATASET_ORDER = [
    "human_GSE157080_fetal_cortical_neuron_DNMT3L_OE",
    "fly_GSE263069_third_instar_larval_brain_DNMT3L_OE",
    "mouse_GSE268951_P28_mFC_Dnmt3l_OE",
]

DATASET_LABELS = {
    "human_GSE157080_fetal_cortical_neuron_DNMT3L_OE": "Human fetal neuron\nGSE157080",
    "fly_GSE263069_third_instar_larval_brain_DNMT3L_OE": "Fly larval brain\nGSE263069",
    "mouse_GSE268951_P28_mFC_Dnmt3l_OE": "Mouse P28 mFC\nGSE268951",
}

DATASET_META = [
    {
        "dataset_id": "human_GSE157080_fetal_cortical_neuron_DNMT3L_OE",
        "species": "human",
        "source": "GSE157080",
        "comparison": "DNMT3L-ZsGreen overexpression vs ZsGreen control",
        "stage_or_context": "cultured cortical neurons differentiated from 21 GW fetal frontal-cortex neuroprogenitors for 4 days",
        "control_n": 3,
        "dnmt3l_oe_n": 3,
        "analysis_level": "author supplementary Table S3 count-based DE statistics plus GEO processed FPKM abundance",
        "role": "primary mammalian fetal/neurogenic comparator",
    },
    {
        "dataset_id": "fly_GSE263069_third_instar_larval_brain_DNMT3L_OE",
        "species": "Drosophila",
        "source": "GSE263069",
        "comparison": "human DNMT3L overexpression vs control larval brain",
        "stage_or_context": "third-instar larval brain",
        "control_n": 3,
        "dnmt3l_oe_n": 3,
        "analysis_level": "count-based DESeq2/edgeR statistics from processed expected counts plus GEO processed FPKM abundance",
        "role": "methylation-independent developmental comparator",
    },
    {
        "dataset_id": "mouse_GSE268951_P28_mFC_Dnmt3l_OE",
        "species": "mouse",
        "source": "GSE268951",
        "comparison": "Emx1-cre;Dnmt3l overexpression vs Emx1-cre control",
        "stage_or_context": "P28 medial frontal cortex",
        "control_n": 6,
        "dnmt3l_oe_n": 6,
        "analysis_level": "public GEO processed expected-count DESeq2/edgeR reanalysis plus FPKM abundance",
        "role": "stage-mismatched postnatal maturation comparator",
    },
]

MODULES = {
    "JAK_STAT_signaling": {
        "pattern": r"\bjak\b|stat|interferon|cytokine|socs|suppressor of cytokine|unpaired|upd[123]|\bhop\b",
        "markers": {
            "human": ["STAT1", "STAT3", "JAK1", "JAK2", "SOCS1", "SOCS3", "IL6ST", "IFNGR1"],
            "mouse": ["Stat1", "Stat3", "Jak1", "Jak2", "Socs1", "Socs3", "Il6st", "Ifngr1"],
            "fly": ["Stat92E", "hop", "upd1", "upd2", "upd3", "Socs36E", "Socs44A"],
        },
    },
    "neurogenic_bHLH_Notch_TF": {
        "pattern": r"bhlh|basic helix-loop-helix|notch|proneural|neurogenin|neurod|ascl|hes|hairy|enhancer of split|achaete|scute|atonal|neural precursor",
        "markers": {
            "human": ["HES1", "HES3", "HES5", "ASCL1", "NEUROG1", "NEUROG2", "NEUROD1", "NEUROD2", "NOTCH1", "DLL1"],
            "mouse": ["Hes1", "Hes3", "Hes5", "Ascl1", "Neurog1", "Neurog2", "Neurod1", "Neurod2", "Notch1", "Dll1"],
            "fly": ["H", "E(spl)m8-HLH", "ac", "sc", "ase", "ato", "da", "dpn", "N", "Dl"],
        },
    },
    "neuronal_differentiation_axon_synapse": {
        "pattern": r"neurogenesis|generation of neurons|neuron differentiation|neuronal differentiation|neuron development|axon|dendrit|synap|neurotrans|glutamate|gaba|cholinergic|neurite|neuron projection",
        "markers": {
            "human": ["DCX", "TUBB3", "MAP2", "RBFOX3", "ELAVL4", "PROX1", "NEURL1"],
            "mouse": ["Dcx", "Tubb3", "Map2", "Rbfox3", "Elavl4", "Prox1", "Neurl1a"],
            "fly": ["elav", "pros", "neur", "Abl"],
        },
    },
    "cell_cycle_CDK_p21": {
        "pattern": r"cell cycle|mitotic|proliferation|dna replication|cdk|cyclin|cdkn|p21|growth arrest|g1/s",
        "markers": {
            "human": ["CDKN1A", "CDKN1B", "CCND1", "CDK2", "MKI67"],
            "mouse": ["Cdkn1a", "Cdkn1b", "Ccnd1", "Cdk2", "Mki67"],
            "fly": ["dap", "CycD", "Cdk1", "Cdk2", "stg"],
        },
    },
    "chromatin_transcription": {
        "pattern": r"chromatin|histone|nucleosome|transcription|rna polymerase|dna binding|methyltransferase|epigen|enhancer|polycomb|trithorax",
        "markers": {
            "human": ["DNMT3A", "DNMT3B", "DNMT3L", "EZH2", "SUZ12", "CTCF", "MEF2C", "RFX3"],
            "mouse": ["Dnmt3a", "Dnmt3b", "Dnmt3l", "Ezh2", "Suz12", "Ctcf", "Mef2c", "Rfx3"],
            "fly": ["E(z)", "Su(z)12", "Trl", "Mef2", "Rfx"],
        },
    },
    "immune_stress_defense": {
        "pattern": r"immune|defense|inflamm|interferon|cytokine|antimicrobial|stress|virus|bacter|oxidative|heat shock|response to",
        "markers": {
            "human": ["ISG15", "IFI27", "IFI44", "MX1", "IRF1", "NFKBIA"],
            "mouse": ["Isg15", "Ifi27l2a", "Ifi44", "Mx1", "Irf1", "Nfkbia"],
            "fly": ["DptA", "DptB", "Mtk", "AttA", "Hsp70Aa", "Rel"],
        },
    },
    "glial_myelination": {
        "pattern": r"glia|astrocyte|oligodendrocyte|myelin|myelination|schwann|ensheath",
        "markers": {
            "human": ["GFAP", "AQP4", "OLIG1", "OLIG2", "MBP", "PLP1", "MOG"],
            "mouse": ["Gfap", "Aqp4", "Olig1", "Olig2", "Mbp", "Plp1", "Mog"],
            "fly": ["repo", "gcm", "wrapper", "Nrg"],
        },
    },
    "metabolism_detox": {
        "pattern": r"metabolic|metabolism|oxidoreductase|cytochrome|glutathione|detox|hydrolase|monooxygenase|lipid|mitochond",
        "markers": {
            "human": ["CYP1B1", "GSTM1", "SOD2", "NQO1"],
            "mouse": ["Cyp1b1", "Gstm1", "Sod2", "Nqo1"],
            "fly": ["Cyp6a2", "Cyp6g1", "GstD1", "Sod2"],
        },
    },
    "MAPK_WNT_BMP_developmental_signaling": {
        "pattern": r"mapk|erk|wnt|bmp|fgf|fibroblast growth factor|tgf|retinoic acid|pi3k|akt|growth factor|receptor signaling",
        "markers": {
            "human": ["FGFR2", "FGF2", "WNT7B", "BMP4", "MAPK1", "AKT1"],
            "mouse": ["Fgfr2", "Fgf2", "Wnt7b", "Bmp4", "Mapk1", "Akt1"],
            "fly": ["Egfr", "rl", "wg", "dpp", "Akt1"],
        },
    },
}

MARKER_CROSSWALK = [
    ("STAT core", "STAT1", "Stat1", "Stat92E", "family-level STAT comparison; fly has one major STAT"),
    ("STAT core", "STAT3", "Stat3", "Stat92E", "family-level STAT comparison; fly has one major STAT"),
    ("JAK kinase", "JAK1", "Jak1", "hop", "family-level kinase comparison"),
    ("JAK kinase", "JAK2", "Jak2", "hop", "family-level kinase comparison"),
    ("STAT feedback", "SOCS1", "Socs1", "Socs36E", "feedback inhibitor family"),
    ("STAT feedback", "SOCS3", "Socs3", "Socs44A", "feedback inhibitor family"),
    ("HES/Notch effector", "HES3", "Hes3", "E(spl)m8-HLH", "Notch/HES-family functional comparison"),
    ("HES/Notch effector", "HES1", "Hes1", "H", "hairy/HES-family functional comparison"),
    ("Proneural bHLH", "ASCL1", "Ascl1", "ac", "achaete-scute family functional comparison"),
    ("Proneural bHLH", "NEUROG2", "Neurog2", "ato", "proneural bHLH-family comparison, not exact one-to-one orthology"),
    ("Neuronal differentiation TF", "NEUROD2", "Neurod2", "elav", "neuronal differentiation marker family/context comparison"),
    ("CDK inhibitor", "CDKN1A", "Cdkn1a", "dap", "cell-cycle inhibitor functional comparison"),
    ("Notch receptor/ligand", "NOTCH1", "Notch1", "N", "Notch receptor family"),
    ("Notch receptor/ligand", "DLL1", "Dll1", "Dl", "Delta ligand family"),
    ("Neuronal fate", "PROX1", "Prox1", "pros", "Pros/Prox neuronal fate family"),
    ("Neuronal differentiation", "ELAVL4", "Elavl4", "elav", "ELAV-family neuronal marker"),
    ("Neurogenic ubiquitin ligase", "NEURL1", "Neurl1a", "neur", "neuralized family"),
]


def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    p = np.where(np.isfinite(p), p, 1.0)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n, dtype=float)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q[order] = np.clip(ranked, 0, 1)
    return q


def safe_float(value: object, default: float = math.nan) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default
    return number if math.isfinite(number) else default


def clean_go_text(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).replace("%3B", ";")
    text = re.sub(r"GO:\d+:", " ", text)
    return text


def signed_direction(value: float, threshold: float = 0.02) -> str:
    if not math.isfinite(value) or abs(value) < threshold:
        return "near_zero"
    return "up" if value > 0 else "down"


def add_common_flags(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["pvalue"] = pd.to_numeric(df["pvalue"], errors="coerce").fillna(1.0)
    df["fdr"] = pd.to_numeric(df["fdr"], errors="coerce").fillna(1.0)
    df["log2fc"] = pd.to_numeric(df["log2fc"], errors="coerce")
    df["control_mean_fpkm"] = pd.to_numeric(df["control_mean_fpkm"], errors="coerce").fillna(0.0)
    df["dnmt3l_oe_mean_fpkm"] = pd.to_numeric(df["dnmt3l_oe_mean_fpkm"], errors="coerce").fillna(0.0)
    df["expressed"] = (df["control_mean_fpkm"] + df["dnmt3l_oe_mean_fpkm"]).gt(0)
    df["nominal_relaxed_de"] = df["pvalue"].le(PVALUE_CUTOFF) & df["log2fc"].abs().ge(RELAXED_LOG2FC)
    df["nominal_large_de"] = df["pvalue"].le(PVALUE_CUTOFF) & df["log2fc"].abs().ge(LARGE_LOG2FC)
    df["fdr_relaxed_de"] = df["fdr"].le(0.05) & df["log2fc"].abs().ge(RELAXED_LOG2FC)
    df["direction_relaxed"] = np.where(
        df["log2fc"].ge(RELAXED_LOG2FC),
        "up",
        np.where(df["log2fc"].le(-RELAXED_LOG2FC), "down", "flat"),
    )
    df["search_text"] = (
        df[["gene_symbol", "description", "go_bp_terms", "go_cc_terms", "go_mf_terms"]]
        .fillna("")
        .agg(" ".join, axis=1)
        .map(clean_go_text)
    )
    df["symbol_key"] = df["gene_symbol"].fillna("").astype(str).str.lower()
    return df


def compute_effect_from_matrix(
    df: pd.DataFrame,
    control_cols: list[str],
    oe_cols: list[str],
) -> tuple[pd.Series, pd.Series, pd.Series, pd.Series, pd.Series]:
    control = df[control_cols].astype(float)
    oe = df[oe_cols].astype(float)
    control_mean = control.mean(axis=1)
    oe_mean = oe.mean(axis=1)
    log2fc = np.log2(oe_mean + PSEUDOCOUNT) - np.log2(control_mean + PSEUDOCOUNT)
    pvalues: list[float] = []
    for idx in range(len(df)):
        ctrl_values = np.log2(control.iloc[idx].to_numpy(dtype=float) + PSEUDOCOUNT)
        oe_values = np.log2(oe.iloc[idx].to_numpy(dtype=float) + PSEUDOCOUNT)
        test = stats.ttest_ind(oe_values, ctrl_values, equal_var=False, nan_policy="omit")
        pvalues.append(float(test.pvalue) if math.isfinite(float(test.pvalue)) else 1.0)
    fdr = bh_fdr(np.asarray(pvalues, dtype=float))
    return control_mean, oe_mean, log2fc, pd.Series(pvalues, index=df.index), pd.Series(fdr, index=df.index)


def load_human_author_de() -> pd.DataFrame:
    if not HUMAN_AUTHOR_DE.exists():
        return pd.DataFrame(columns=["gene_id", "author_log2fc", "author_pvalue", "author_fdr"])
    raw = pd.read_csv(HUMAN_AUTHOR_DE, sep="\t", low_memory=False)
    out = pd.DataFrame(
        {
            "gene_id": raw["gene_id"].astype(str),
            "author_log2fc": pd.to_numeric(raw["log2FoldChange"], errors="coerce"),
            "author_pvalue": pd.to_numeric(raw["pvalue"], errors="coerce"),
            "author_fdr": pd.to_numeric(raw["qvalue"], errors="coerce"),
        }
    )
    return out.drop_duplicates("gene_id", keep="first")


def load_human() -> pd.DataFrame:
    raw = pd.read_csv(HUMAN_EXPR, sep="\t")
    control_cols = ["Z1_FPKM", "Z2_FPKM", "Z3_FPKM"]
    oe_cols = ["D1_FPKM", "D2_FPKM", "D3_FPKM"]
    control_mean, oe_mean, log2fc, pvalue, fdr = compute_effect_from_matrix(raw, control_cols, oe_cols)
    base = pd.DataFrame(
        {
            "dataset_id": DATASET_ORDER[0],
            "species": "human",
            "gene_id": raw["gene_id"].astype(str),
            "gene_symbol": raw["GeneSymbol"].astype(str),
            "description": raw["Description"].astype(str),
            "go_bp_terms": raw.get("GO_BP", "").astype(str),
            "go_cc_terms": raw.get("GO_CC", "").astype(str),
            "go_mf_terms": raw.get("GO_MF", "").astype(str),
            "control_mean_fpkm": control_mean,
            "dnmt3l_oe_mean_fpkm": oe_mean,
            "fpkm_welch_log2fc": log2fc,
            "fpkm_welch_pvalue": pvalue,
            "fpkm_welch_fdr": fdr,
        }
    )
    out = base.merge(load_human_author_de(), on="gene_id", how="left")
    out["log2fc"] = out["author_log2fc"].combine_first(out["fpkm_welch_log2fc"])
    out["pvalue"] = out["author_pvalue"].combine_first(out["fpkm_welch_pvalue"])
    out["fdr"] = out["author_fdr"].combine_first(out["fpkm_welch_fdr"])
    out["effect_source"] = np.where(
        out["author_log2fc"].notna(),
        "author_Table_S3_count_based_DE",
        "processed_FPKM_Welch_fallback",
    )
    return add_common_flags(out)


def load_fly() -> pd.DataFrame:
    if FLY_COUNT_EFFECTS.exists():
        raw = pd.read_csv(FLY_COUNT_EFFECTS, sep="\t", low_memory=False)
        out = pd.DataFrame(
            {
                "dataset_id": DATASET_ORDER[1],
                "species": "fly",
                "gene_id": raw["gene_id"],
                "gene_symbol": raw["GeneSymbol"].astype(str),
                "description": raw["Description"].astype(str),
                "go_bp_terms": raw.get("GO_BP", "").astype(str),
                "go_cc_terms": raw.get("GO_CC", "").astype(str),
                "go_mf_terms": raw.get("GO_MF", "").astype(str),
                "control_mean_fpkm": raw["control_mean_fpkm"],
                "dnmt3l_oe_mean_fpkm": raw["dnmt3l_oe_mean_fpkm"],
                "log2fc": raw["deseq2_log2fc"],
                "pvalue": raw["deseq2_pvalue"],
                "fdr": raw["deseq2_padj"],
                "effect_source": "processed_expected_count_DESeq2",
            }
        )
        return add_common_flags(out)

    raw = pd.read_csv(FLY_EXPR_EFFECTS, sep="\t")
    out = pd.DataFrame(
        {
            "dataset_id": DATASET_ORDER[1],
            "species": "fly",
            "gene_id": raw["gene_id"],
            "gene_symbol": raw["GeneSymbol"].astype(str),
            "description": raw["Description"].astype(str),
            "go_bp_terms": raw.get("GO_BP", "").astype(str),
            "go_cc_terms": raw.get("GO_CC", "").astype(str),
            "go_mf_terms": raw.get("GO_MF", "").astype(str),
            "control_mean_fpkm": raw["control_mean_fpkm"],
            "dnmt3l_oe_mean_fpkm": raw["dnmt3l_oe_mean_fpkm"],
            "log2fc": raw["log2fc_dnmt3l_oe_vs_control"],
            "pvalue": raw["welch_pvalue_log2fpkm"],
            "fdr": raw["welch_fdr_log2fpkm"],
            "effect_source": "processed_FPKM_Welch_fallback",
        }
    )
    return add_common_flags(out)


def load_mouse() -> pd.DataFrame:
    raw = pd.read_csv(MOUSE_EXPR_EFFECTS, sep="\t")
    out = pd.DataFrame(
        {
            "dataset_id": DATASET_ORDER[2],
            "species": "mouse",
            "gene_id": raw["gene_id"],
            "gene_symbol": raw["gene_symbol"].astype(str),
            "description": raw["description"].astype(str),
            "go_bp_terms": raw.get("go_bp_terms", "").astype(str),
            "go_cc_terms": "",
            "go_mf_terms": "",
            "control_mean_fpkm": raw["gse268951_control_mean_fpkm"],
            "dnmt3l_oe_mean_fpkm": raw["gse268951_dnmt3l_oe_mean_fpkm"],
            "log2fc": raw["gse268951_log2fc_oe_vs_control"],
            "pvalue": raw.get("gse268951_primary_pvalue", raw.get("gse268951_welch_pvalue")),
            "fdr": raw.get("gse268951_primary_fdr", raw.get("gse268951_bh_fdr")),
            "effect_source": raw.get("gse268951_effect_source", "validated_edgeR_best_transcript"),
        }
    )
    return add_common_flags(out)


def deduplicate_by_symbol(df: pd.DataFrame) -> pd.DataFrame:
    expressed = df.loc[df["expressed"]].copy()
    expressed["mean_total"] = expressed["control_mean_fpkm"] + expressed["dnmt3l_oe_mean_fpkm"]
    expressed = expressed.sort_values(["symbol_key", "mean_total"], ascending=[True, False])
    return expressed.drop_duplicates("symbol_key", keep="first")


def module_hits(df: pd.DataFrame, module: str, module_def: dict) -> pd.Series:
    pattern_hits = df["search_text"].str.contains(module_def["pattern"], case=False, regex=True, na=False)
    markers = [symbol.lower() for symbol in module_def["markers"].get(df["species"].iloc[0], [])]
    marker_hits = df["symbol_key"].isin(markers)
    return pattern_hits | marker_hits


def summarize_modules(all_df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for dataset_id, df in all_df.groupby("dataset_id", sort=False):
        bg = deduplicate_by_symbol(df)
        bg = bg.loc[np.isfinite(bg["log2fc"])]
        for module, module_def in MODULES.items():
            hits = module_hits(bg, module, module_def)
            module_df = bg.loc[hits]
            other_df = bg.loc[~hits]
            if len(module_df) >= 2 and len(other_df) >= 2:
                rank = stats.mannwhitneyu(module_df["log2fc"], other_df["log2fc"], alternative="two-sided")
                rank_pvalue = float(rank.pvalue)
            else:
                rank_pvalue = math.nan
            nominal_up = int((module_df["nominal_relaxed_de"] & module_df["log2fc"].gt(0)).sum())
            nominal_down = int((module_df["nominal_relaxed_de"] & module_df["log2fc"].lt(0)).sum())
            rows.append(
                {
                    "dataset_id": dataset_id,
                    "species": bg["species"].iloc[0],
                    "module": module,
                    "background_expressed_genes": len(bg),
                    "module_expressed_genes": int(hits.sum()),
                    "module_median_log2fc": float(module_df["log2fc"].median()) if len(module_df) else math.nan,
                    "module_mean_log2fc": float(module_df["log2fc"].mean()) if len(module_df) else math.nan,
                    "module_positive_fraction": float(module_df["log2fc"].gt(0).mean()) if len(module_df) else math.nan,
                    "module_nominal_relaxed_up": nominal_up,
                    "module_nominal_relaxed_down": nominal_down,
                    "module_nominal_relaxed_total": nominal_up + nominal_down,
                    "module_fdr_relaxed_total": int(module_df["fdr_relaxed_de"].sum()) if len(module_df) else 0,
                    "rank_sum_pvalue_vs_other_expressed": rank_pvalue,
                    "top_nominal_genes": ";".join(
                        module_df.loc[module_df["nominal_relaxed_de"]]
                        .assign(abs_log2fc=lambda x: x["log2fc"].abs())
                        .sort_values(["pvalue", "abs_log2fc"], ascending=[True, False])["gene_symbol"]
                        .head(25)
                        .astype(str)
                    ),
                }
            )
    out = pd.DataFrame(rows)
    out["rank_sum_fdr_within_dataset"] = out.groupby("dataset_id")["rank_sum_pvalue_vs_other_expressed"].transform(
        lambda series: bh_fdr(series.fillna(1.0).to_numpy(dtype=float))
    )
    out["module_direction"] = out["module_median_log2fc"].map(signed_direction)
    return out


def summarize_concordance(module_summary: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for module, sub in module_summary.groupby("module", sort=False):
        lookup = {row["dataset_id"]: row for _, row in sub.iterrows()}
        human = lookup.get(DATASET_ORDER[0])
        fly = lookup.get(DATASET_ORDER[1])
        mouse = lookup.get(DATASET_ORDER[2])
        def value(row: pd.Series | None, col: str) -> object:
            return row[col] if row is not None else math.nan

        human_dir = value(human, "module_direction")
        fly_dir = value(fly, "module_direction")
        mouse_dir = value(mouse, "module_direction")
        primary_pair_call = "not_testable"
        if isinstance(human_dir, str) and isinstance(fly_dir, str):
            if "near_zero" in {human_dir, fly_dir}:
                primary_pair_call = "weak_or_near_zero_direction"
            elif human_dir == fly_dir:
                primary_pair_call = "human_fly_same_direction"
            else:
                primary_pair_call = "human_fly_opposite_direction"
        rows.append(
            {
                "module": module,
                "human_median_log2fc": value(human, "module_median_log2fc"),
                "fly_median_log2fc": value(fly, "module_median_log2fc"),
                "mouse_median_log2fc": value(mouse, "module_median_log2fc"),
                "human_direction": human_dir,
                "fly_direction": fly_dir,
                "mouse_direction": mouse_dir,
                "primary_human_fly_concordance": primary_pair_call,
                "human_nominal_changed": value(human, "module_nominal_relaxed_total"),
                "fly_nominal_changed": value(fly, "module_nominal_relaxed_total"),
                "mouse_nominal_changed": value(mouse, "module_nominal_relaxed_total"),
                "interpretation": interpret_module(module, primary_pair_call, human, fly, mouse),
            }
        )
    return pd.DataFrame(rows)


def interpret_module(module: str, primary_call: str, human: pd.Series | None, fly: pd.Series | None, mouse: pd.Series | None) -> str:
    if module == "JAK_STAT_signaling":
        return (
            "Use as pathway-level support only: expression changes can reflect STAT pathway components or feedback, "
            "but STAT1/STAT3 activity is better tested by phosphorylation or motif/activity assays."
        )
    if module == "neurogenic_bHLH_Notch_TF":
        return (
            "Relevant to the STAT-paper neurogenic TF axis, but cross-species mapping is family-level and should not be "
            "treated as exact ASCL1/NEUROG2 conservation."
        )
    if module == "glial_myelination":
        return (
            "Most relevant to mouse P28 postnatal maturation; weak fly/human matching would not imply identical stages."
        )
    if primary_call == "human_fly_same_direction":
        return "Directional human-fly similarity at module level; still descriptive because orthology and stage are coarse."
    if primary_call == "human_fly_opposite_direction":
        return "Directional divergence between fetal human and larval fly expression responses."
    return "No strong directional primary-pair signal in this expression-only comparison."


def select_symbol(df: pd.DataFrame, symbol: str) -> pd.Series | None:
    if not symbol:
        return None
    matches = df.loc[df["symbol_key"].eq(symbol.lower())].copy()
    if matches.empty:
        return None
    matches["mean_total"] = matches["control_mean_fpkm"] + matches["dnmt3l_oe_mean_fpkm"]
    return matches.sort_values("mean_total", ascending=False).iloc[0]


def marker_effects(all_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    by_dataset = {dataset_id: df for dataset_id, df in all_df.groupby("dataset_id", sort=False)}
    rows: list[dict[str, object]] = []
    crosswalk_rows: list[dict[str, str]] = []
    for family, human_symbol, mouse_symbol, fly_symbol, note in MARKER_CROSSWALK:
        marker_id = f"{family}: {human_symbol}/{mouse_symbol}/{fly_symbol}"
        crosswalk_rows.append(
            {
                "marker_id": marker_id,
                "marker_family": family,
                "human_symbol": human_symbol,
                "mouse_symbol": mouse_symbol,
                "fly_symbol": fly_symbol,
                "mapping_note": note,
            }
        )
        for dataset_id, species, symbol in [
            (DATASET_ORDER[0], "human", human_symbol),
            (DATASET_ORDER[1], "fly", fly_symbol),
            (DATASET_ORDER[2], "mouse", mouse_symbol),
        ]:
            hit = select_symbol(by_dataset[dataset_id], symbol)
            rows.append(
                {
                    "marker_id": marker_id,
                    "marker_family": family,
                    "dataset_id": dataset_id,
                    "species": species,
                    "requested_symbol": symbol,
                    "matched_gene_id": "" if hit is None else hit["gene_id"],
                    "matched_gene_symbol": "" if hit is None else hit["gene_symbol"],
                    "control_mean_fpkm": math.nan if hit is None else hit["control_mean_fpkm"],
                    "dnmt3l_oe_mean_fpkm": math.nan if hit is None else hit["dnmt3l_oe_mean_fpkm"],
                    "log2fc": math.nan if hit is None else hit["log2fc"],
                    "pvalue": math.nan if hit is None else hit["pvalue"],
                    "fdr": math.nan if hit is None else hit["fdr"],
                    "effect_source": "" if hit is None else hit.get("effect_source", ""),
                    "nominal_relaxed_de": "" if hit is None else bool(hit["nominal_relaxed_de"]),
                    "fdr_relaxed_de": "" if hit is None else bool(hit["fdr_relaxed_de"]),
                    "mapping_note": note,
                }
            )

    marker_long = pd.DataFrame(rows)
    marker_summary_rows = []
    for (dataset_id, family), sub in marker_long.groupby(["dataset_id", "marker_family"], sort=False):
        present = sub.loc[np.isfinite(sub["log2fc"])]
        marker_summary_rows.append(
            {
                "dataset_id": dataset_id,
                "marker_family": family,
                "markers_tested": len(present),
                "median_log2fc": float(present["log2fc"].median()) if len(present) else math.nan,
                "mean_log2fc": float(present["log2fc"].mean()) if len(present) else math.nan,
                "positive_fraction": float(present["log2fc"].gt(0).mean()) if len(present) else math.nan,
                "nominal_relaxed_de": int(present["nominal_relaxed_de"].astype(bool).sum()) if len(present) else 0,
                "members": ";".join(present["matched_gene_symbol"].astype(str)),
            }
        )
    return pd.DataFrame(crosswalk_rows), marker_long, pd.DataFrame(marker_summary_rows)


def write_tables(
    all_df: pd.DataFrame,
    module_summary: pd.DataFrame,
    concordance: pd.DataFrame,
    marker_crosswalk: pd.DataFrame,
    marker_long: pd.DataFrame,
    marker_summary: pd.DataFrame,
) -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(DATASET_META).to_csv(OUTDIR / "cross_species_dataset_summary.tsv", sep="\t", index=False)
    all_df[
        [
            "dataset_id",
            "species",
            "gene_id",
            "gene_symbol",
            "description",
            "control_mean_fpkm",
            "dnmt3l_oe_mean_fpkm",
            "log2fc",
            "pvalue",
            "fdr",
            "effect_source",
            "nominal_relaxed_de",
            "nominal_large_de",
            "fdr_relaxed_de",
            "go_bp_terms",
        ]
    ].to_csv(OUTDIR / "cross_species_gene_expression_effects_long.tsv", sep="\t", index=False)
    module_summary.to_csv(OUTDIR / "cross_species_pathway_module_summary.tsv", sep="\t", index=False)
    concordance.to_csv(OUTDIR / "cross_species_pathway_concordance.tsv", sep="\t", index=False)
    marker_crosswalk.to_csv(OUTDIR / "cross_species_marker_family_crosswalk.tsv", sep="\t", index=False)
    marker_long.to_csv(OUTDIR / "cross_species_marker_symbol_effects.tsv", sep="\t", index=False)
    marker_summary.to_csv(OUTDIR / "cross_species_marker_family_summary.tsv", sep="\t", index=False)


def plot_heatmap(table: pd.DataFrame, value_col: str, row_col: str, out_prefix: Path, title: str, note: str = "") -> None:
    FIGDIR.mkdir(parents=True, exist_ok=True)
    pivot = table.pivot(index=row_col, columns="dataset_id", values=value_col)
    pivot = pivot.reindex(columns=DATASET_ORDER)
    values = pivot.to_numpy(dtype=float)
    finite = values[np.isfinite(values)]
    vmax = max(0.25, float(np.nanpercentile(np.abs(finite), 95))) if finite.size else 1.0
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    fig_height = max(4.5, 0.44 * len(pivot.index) + 1.7)
    fig, ax = plt.subplots(figsize=(9.4, fig_height))
    im = ax.imshow(values, cmap="RdBu_r", norm=norm, aspect="auto")
    ax.set_xticks(range(len(DATASET_ORDER)))
    ax.set_xticklabels([DATASET_LABELS[col] for col in DATASET_ORDER], fontsize=9)
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels([label.replace("_", " ") for label in pivot.index], fontsize=9)
    ax.set_title(title, fontsize=12, weight="bold", pad=14)
    for y in range(values.shape[0]):
        for x in range(values.shape[1]):
            value = values[y, x]
            if np.isfinite(value):
                ax.text(x, y, f"{value:.2f}", ha="center", va="center", fontsize=8, color="black")
            else:
                ax.text(x, y, "NA", ha="center", va="center", fontsize=8, color="0.45")
    ax.tick_params(axis="both", length=0)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label(value_col.replace("_", " "), fontsize=9)
    if note:
        fig.text(0.01, 0.01, note, ha="left", va="bottom", fontsize=8, color="0.35")
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(out_prefix.with_suffix(".png"), dpi=300)
    fig.savefig(out_prefix.with_suffix(".pdf"))
    plt.close(fig)


def plot_nominal_counts(module_summary: pd.DataFrame) -> None:
    FIGDIR.mkdir(parents=True, exist_ok=True)
    modules = list(MODULES)
    x = np.arange(len(modules))
    width = 0.24
    fig, axes = plt.subplots(1, 3, figsize=(14, 5.4), sharey=True)
    for ax, dataset_id in zip(axes, DATASET_ORDER):
        sub = module_summary.loc[module_summary["dataset_id"].eq(dataset_id)].set_index("module").reindex(modules)
        up = sub["module_nominal_relaxed_up"].fillna(0).to_numpy(dtype=float)
        down = -sub["module_nominal_relaxed_down"].fillna(0).to_numpy(dtype=float)
        ax.bar(x, up, color="#b33a3a", width=0.72, label="up")
        ax.bar(x, down, color="#3a6ea5", width=0.72, label="down")
        ax.axhline(0, color="0.25", linewidth=0.8)
        ax.set_title(DATASET_LABELS[dataset_id], fontsize=10)
        ax.set_xticks(x)
        ax.set_xticklabels([m.replace("_", "\n") for m in modules], rotation=90, fontsize=7)
        ax.tick_params(axis="y", labelsize=8)
    axes[0].set_ylabel(
        f"Number of nominally changed genes\n(P <= {PVALUE_CUTOFF}, |log2FC| >= {RELAXED_LOG2FC})",
        fontsize=9,
    )
    axes[-1].legend(frameon=False, fontsize=8, loc="upper right")
    fig.suptitle("Direction of nominal pathway-module expression changes", fontsize=12, weight="bold", y=0.98)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(FIGDIR / "cross_species_pathway_nominal_direction_counts.png", dpi=300)
    fig.savefig(FIGDIR / "cross_species_pathway_nominal_direction_counts.pdf")
    plt.close(fig)


def plot_nominal_fractions(module_summary: pd.DataFrame) -> None:
    FIGDIR.mkdir(parents=True, exist_ok=True)
    modules = list(MODULES)
    x = np.arange(len(modules))
    fig, axes = plt.subplots(1, 3, figsize=(14, 5.4), sharey=True)
    for ax, dataset_id in zip(axes, DATASET_ORDER):
        sub = module_summary.loc[module_summary["dataset_id"].eq(dataset_id)].set_index("module").reindex(modules)
        denom = sub["module_expressed_genes"].replace(0, np.nan).to_numpy(dtype=float)
        up = 100 * sub["module_nominal_relaxed_up"].fillna(0).to_numpy(dtype=float) / denom
        down = -100 * sub["module_nominal_relaxed_down"].fillna(0).to_numpy(dtype=float) / denom
        ax.bar(x, up, color="#b33a3a", width=0.72, label="up")
        ax.bar(x, down, color="#3a6ea5", width=0.72, label="down")
        ax.axhline(0, color="0.25", linewidth=0.8)
        ax.set_title(DATASET_LABELS[dataset_id], fontsize=10)
        ax.set_xticks(x)
        ax.set_xticklabels([m.replace("_", "\n") for m in modules], rotation=90, fontsize=7)
        ax.tick_params(axis="y", labelsize=8)
    axes[0].set_ylabel(
        f"Percent of expressed module genes\n(P <= {PVALUE_CUTOFF}, |log2FC| >= {RELAXED_LOG2FC})",
        fontsize=9,
    )
    axes[-1].legend(frameon=False, fontsize=8, loc="upper right")
    fig.suptitle("Normalized direction of nominal pathway-module expression changes", fontsize=12, weight="bold", y=0.98)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(FIGDIR / "cross_species_pathway_nominal_direction_fraction.png", dpi=300)
    fig.savefig(FIGDIR / "cross_species_pathway_nominal_direction_fraction.pdf")
    plt.close(fig)


def write_note(module_summary: pd.DataFrame, concordance: pd.DataFrame, marker_summary: pd.DataFrame) -> None:
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    primary_same = concordance.loc[concordance["primary_human_fly_concordance"].eq("human_fly_same_direction"), "module"].tolist()
    primary_opposite = concordance.loc[
        concordance["primary_human_fly_concordance"].eq("human_fly_opposite_direction"), "module"
    ].tolist()
    jak = marker_summary.loc[marker_summary["marker_family"].isin(["STAT core", "JAK kinase", "STAT feedback"])]
    neuro_tf = marker_summary.loc[
        marker_summary["marker_family"].isin(["HES/Notch effector", "Proneural bHLH", "Neuronal differentiation TF", "CDK inhibitor"])
    ]

    def fmt_modules(modules: list[str]) -> str:
        return ", ".join(f"`{module}`" for module in modules) if modules else "none"

    lines = [
        "# Cross-Species DNMT3L-OE Expression Conservation First Pass",
        "",
        "Date: 2026-06-11",
        "",
        "## Purpose",
        "",
        "This analysis compares the human fetal cortical-neuron DNMT3L-OE dataset (`GSE157080`) with the Drosophila third-instar larval-brain human-DNMT3L-OE dataset (`GSE263069`) to ask whether STAT/TF or other expression pathways show directional conservation. The human arm uses the published author supplementary Table S3 count-based DE statistics, with GEO FPKM retained only as abundance context. The fly arm uses the public processed expected-count DESeq2/edgeR reanalysis when available. The mouse P28 medial-frontal-cortex Dnmt3l-OE dataset (`GSE268951`) uses the public GEO processed expected-count DESeq2/edgeR reanalysis and is included as a stage-mismatched postnatal maturation comparator, not as an equivalent neurogenesis-stage dataset.",
        "",
        "## Main outputs",
        "",
        "- `results/cross_species_dnmt3l_expression_conservation/cross_species_pathway_concordance.tsv`",
        "- `results/cross_species_dnmt3l_expression_conservation/cross_species_pathway_module_summary.tsv`",
        "- `results/cross_species_dnmt3l_expression_conservation/cross_species_marker_symbol_effects.tsv`",
        "- `figures/cross_species_dnmt3l_expression_conservation/cross_species_pathway_module_median_log2fc.png`",
        "- `figures/cross_species_dnmt3l_expression_conservation/cross_species_marker_family_median_log2fc.png`",
        "- `figures/cross_species_dnmt3l_expression_conservation/cross_species_pathway_nominal_direction_fraction.png`",
        "",
        "## First-pass interpretation",
        "",
        f"- Human-fly pathway modules with the same median direction: {fmt_modules(primary_same)}.",
        f"- Human-fly pathway modules with opposite median direction: {fmt_modules(primary_opposite)}.",
        "- The comparison is useful at the pathway/family level, especially for methylation-independent developmental response logic, but it should not be described as one-to-one ortholog conservation without a curated DIOPT/FlyBase orthology pass.",
        "- STAT conservation is best interpreted as possible pathway responsiveness or feedback at the mRNA level. It does not prove conserved STAT1/STAT3 activation, because the direct STAT mechanism depends on protein phosphorylation and downstream activity.",
        "- Neurogenic bHLH/Notch-family comparisons are biologically relevant, but the fly markers are family/function analogs. They should not be written as direct ASCL1 or NEUROG2 conservation.",
        "- Mouse P28 is useful for asking whether DNMT3L perturbation also affects postnatal maturation/glial/synaptic programs, but it should remain separated from fetal human and larval fly neurodevelopmental timing.",
        "",
        "## Marker-family summary",
        "",
        "Median marker-family log2FC values are summarized below. Blank or NA values mean the requested symbol was not present in the processed expression table.",
        "",
        jak.pivot(index="marker_family", columns="dataset_id", values="median_log2fc").reindex(columns=DATASET_ORDER).to_markdown(),
        "",
        neuro_tf.pivot(index="marker_family", columns="dataset_id", values="median_log2fc").reindex(columns=DATASET_ORDER).to_markdown(),
        "",
        "## Recommended manuscript language",
        "",
        "A conservative sentence would be: \"Using the published human author count-based RNA-seq statistics and public processed expected-count reanalyses for fly and mouse, broad pathway modules show mostly same-direction positive median shifts between human fetal DNMT3L-OE neurons and Drosophila larval-brain DNMT3L-OE, whereas curated STAT/neurogenic marker families remain mixed; this supports descriptive pathway-level responsiveness, not strict ortholog conservation or proof of conserved STAT activity.\"",
        "",
    ]
    NOTE.write_text("\n".join(lines))


def main() -> None:
    human = load_human()
    fly = load_fly()
    mouse = load_mouse()
    all_df = pd.concat([human, fly, mouse], ignore_index=True)
    module_summary = summarize_modules(all_df)
    concordance = summarize_concordance(module_summary)
    marker_crosswalk, marker_long, marker_summary = marker_effects(all_df)
    write_tables(all_df, module_summary, concordance, marker_crosswalk, marker_long, marker_summary)

    plot_heatmap(
        module_summary,
        "module_median_log2fc",
        "module",
        FIGDIR / "cross_species_pathway_module_median_log2fc",
        "DNMT3L-OE pathway-module median expression shifts",
        "Primary comparison: human fetal neuron vs fly larval brain; mouse P28 mFC is stage-context.",
    )
    plot_heatmap(
        marker_summary,
        "median_log2fc",
        "marker_family",
        FIGDIR / "cross_species_marker_family_median_log2fc",
        "Curated STAT/TF marker-family median expression shifts",
        "Fly markers are family/function analogs, not strict one-to-one orthologs.",
    )
    plot_nominal_counts(module_summary)
    plot_nominal_fractions(module_summary)
    write_note(module_summary, concordance, marker_summary)

    print(f"Wrote {OUTDIR}")
    print(f"Wrote {FIGDIR}")
    print(f"Wrote {NOTE}")


if __name__ == "__main__":
    main()
