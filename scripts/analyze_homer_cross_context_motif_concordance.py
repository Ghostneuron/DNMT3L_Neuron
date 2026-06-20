#!/usr/bin/env python3
"""Integrate HOMER known-motif results across mechanism datasets."""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


INPUTS = [
    ("human_DNMT3L_OE_methylation", "results/homer_motif_enrichment/homer_known_motif_summary.tsv"),
    ("human_fetal_DS_methylation", "results/ds_homer_motif_enrichment/homer_known_motif_summary.tsv"),
    (
        "human_adult_DS_methylation",
        "results/adult_human_ds_homer_motif_enrichment_sample25k_s8/homer_known_motif_summary.tsv",
    ),
    ("mouse_Dnmt3l_OE_methylation", "results/mouse_homer_motif_enrichment/homer_known_motif_summary.tsv"),
    ("mouse_DS_model_methylation", "results/mouse_ds_homer_motif_enrichment/homer_known_motif_summary.tsv"),
    ("human_fetal_DS_scATAC", "results/gse305146_scATAC_homer_motif_enrichment/homer_known_motif_summary.tsv"),
]

FAMILY_ORDER = [
    "HOMEBOX_ANY",
    "LHX_HOMEBOX",
    "DLX_HOMEBOX",
    "NKX_HOMEBOX",
    "GSX_HOMEBOX",
    "EMX_HOMEBOX",
    "ISL_HOMEBOX",
    "POU_HOMEBOX",
    "NEURO_BHLH",
    "BHLH_ANY",
    "OLIG_BHLH",
    "SOX_HMG",
    "RFX_XBOX",
    "MEF2",
    "STAT",
    "KLF_SP_GC",
    "CTCF_BORIS",
    "ZINC_FINGER_OTHER",
    "ZEB_SNAI",
    "NUCLEAR_RECEPTOR",
    "FOX",
    "E2F",
    "AP1_CEBP",
    "SMAD",
    "TEAD",
]

GROUP_ORDER = [
    "human_DNMT3L_OE_methylation|all",
    "human_DNMT3L_OE_methylation|hyper",
    "human_DNMT3L_OE_methylation|hypo",
    "human_fetal_DS_methylation|hyper",
    "human_fetal_DS_methylation|hypo",
    "human_adult_DS_methylation|hyper",
    "human_adult_DS_methylation|hypo",
    "mouse_Dnmt3l_OE_methylation|hyper",
    "mouse_Dnmt3l_OE_methylation|hypo",
    "mouse_DS_model_methylation|hyper",
    "mouse_DS_model_methylation|hypo",
    "human_fetal_DS_scATAC|DS_open",
    "human_fetal_DS_scATAC|Control_open",
]

GROUP_LABELS = {
    "human_DNMT3L_OE_methylation|all": "Human DNMT3L-OE\nDMR all",
    "human_DNMT3L_OE_methylation|hyper": "Human DNMT3L-OE\nhyper-DMR",
    "human_DNMT3L_OE_methylation|hypo": "Human DNMT3L-OE\nhypo-DMR",
    "human_fetal_DS_methylation|hyper": "Fetal DS\nhyper-DMR",
    "human_fetal_DS_methylation|hypo": "Fetal DS\nhypo-DMR",
    "human_adult_DS_methylation|hyper": "Adult DS\nhyper-DMR",
    "human_adult_DS_methylation|hypo": "Adult DS\nhypo-DMR",
    "mouse_Dnmt3l_OE_methylation|hyper": "Mouse Dnmt3l-OE\nhyper-DMR",
    "mouse_Dnmt3l_OE_methylation|hypo": "Mouse Dnmt3l-OE\nhypo-DMR",
    "mouse_DS_model_methylation|hyper": "Mouse DS model\nhyper-DMC",
    "mouse_DS_model_methylation|hypo": "Mouse DS model\nhypo-DMC",
    "human_fetal_DS_scATAC|DS_open": "Fetal DS ATAC\nDS-open",
    "human_fetal_DS_scATAC|Control_open": "Fetal DS ATAC\nControl-open",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="results/mechanism_question")
    parser.add_argument("--figdir", default="figures/mechanism_question")
    parser.add_argument("--score-cap", type=float, default=10.0)
    parser.add_argument("--q-threshold", type=float, default=0.05)
    parser.add_argument("--suggestive-q-threshold", type=float, default=0.25)
    return parser.parse_args()


def norm_direction(direction: str) -> str:
    direction = (direction or "").strip()
    if direction in {"hyper_top", "hyper"}:
        return "hyper"
    if direction in {"hypo_top", "hypo"}:
        return "hypo"
    return direction


def parse_float(value: object) -> float:
    if value is None:
        return float("nan")
    text = str(value).strip().replace(",", "")
    if text.endswith("%"):
        text = text[:-1]
    if not text:
        return float("nan")
    try:
        return float(text)
    except ValueError:
        return float("nan")


def safe_score(q_value: float, p_value: float, cap: float) -> float:
    value = q_value
    if math.isnan(value):
        value = p_value
    if math.isnan(value):
        return 0.0
    if value <= 0:
        return cap
    return min(cap, -math.log10(value))


def motif_head(motif_name: str) -> str:
    return motif_name.split("/", 1)[0]


def motif_tf(motif_name: str) -> str:
    head = motif_head(motif_name)
    return re.sub(r"\(.*", "", head).strip()


def motif_parenthetical(motif_name: str) -> str:
    match = re.search(r"\(([^)]*)\)", motif_head(motif_name))
    return match.group(1) if match else ""


def classify_families(motif_name: str) -> list[str]:
    text = motif_head(motif_name)
    tf = motif_tf(motif_name)
    parenthetical = motif_parenthetical(motif_name)
    probe = f"{text} {tf} {parenthetical}"
    families: set[str] = set()

    if re.search(r"Homeobox|LHX|Lhx|DLX|Dlx|Nkx|NKX|Gsx|GSX|Emx|EMX|Isl|ISL|Hox|HOX|Pou|POU|Brn|BRN|Otx|OTX|Lmx|LMX|Mixl|En1|Pax|PAX|Pitx|PITX|Tlx|TLX", probe):
        families.add("HOMEBOX_ANY")
    if re.search(r"LHX|Lhx", probe):
        families.add("LHX_HOMEBOX")
    if re.search(r"DLX|Dlx", probe):
        families.add("DLX_HOMEBOX")
    if re.search(r"Nkx|NKX", probe):
        families.add("NKX_HOMEBOX")
    if re.search(r"Gsx|GSX", probe):
        families.add("GSX_HOMEBOX")
    if re.search(r"Emx|EMX", probe):
        families.add("EMX_HOMEBOX")
    if re.search(r"Isl|ISL", probe):
        families.add("ISL_HOMEBOX")
    if re.search(r"Pou|POU|Brn|BRN|OCT", probe):
        families.add("POU_HOMEBOX")
        families.add("HOMEBOX_ANY")

    if re.search(r"bHLH|HLH|Atoh|ATOH|NeuroD|NEUROD|NeuroG|NEUROG|TCF4|Olig|OLIG|Myc|MYC|E2A|Twist|HIF|Arnt|Ahr|BHLH", probe):
        families.add("BHLH_ANY")
    if re.search(r"Atoh|ATOH|NeuroD|NEUROD|NeuroG|NEUROG|TCF4", probe):
        families.add("NEURO_BHLH")
    if re.search(r"Olig|OLIG", probe):
        families.add("OLIG_BHLH")
        families.add("NEURO_BHLH")
        families.add("BHLH_ANY")

    if re.search(r"Sox|SOX|HMG", probe):
        families.add("SOX_HMG")
    if re.search(r"RFX|Rfx|X-box|X_BOX", probe):
        families.add("RFX_XBOX")
    if re.search(r"MEF2|Mef2", probe):
        families.add("MEF2")
    if re.search(r"STAT|Stat", probe):
        families.add("STAT")
    if re.search(r"KLF|Klf|Sp[0-9]?|SP[0-9]?|Maz|MAZ|ZNF148|Zfp281", probe):
        families.add("KLF_SP_GC")
    if re.search(r"CTCF|BORIS|CTCFL", probe):
        families.add("CTCF_BORIS")
    if re.search(r"ZEB|Zeb|SNAI|Snai|Slug", probe):
        families.add("ZEB_SNAI")
    if re.search(r"Zf|ZF|Znf|ZNF|ZBTB|Zbtb|ZFX|Zfx|PRDM|Zfp", probe) and "KLF_SP_GC" not in families:
        families.add("ZINC_FINGER_OTHER")
    if re.search(r"\bNR\b|Nuclear|RXR|Rxr|RAR|Rar|COUP|Nr2f|NR2F|EAR2|ROR|THR|Tlx\?\(NR\)", probe):
        families.add("NUCLEAR_RECEPTOR")
    if re.search(r"FOX|Fox", probe):
        families.add("FOX")
    if re.search(r"E2F", probe):
        families.add("E2F")
    if re.search(r"AP-1|AP1|JUN|Jun|FOS|Fos|CEBP|Cebp|MAF|Maf", probe):
        families.add("AP1_CEBP")
    if re.search(r"SMAD|Smad", probe):
        families.add("SMAD")
    if re.search(r"TEAD|Tead", probe):
        families.add("TEAD")

    return sorted(families) or ["OTHER"]


def read_known_motifs(path: Path, broad_context: str, score_cap: float) -> list[dict[str, object]]:
    if not path.exists():
        return []
    rows: list[dict[str, object]] = []
    with path.open(newline="", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            q_value = parse_float(row.get("q_value"))
            p_value = parse_float(row.get("p_value"))
            rank = int(parse_float(row.get("rank")))
            target_pct = parse_float(row.get("target_percent"))
            background_pct = parse_float(row.get("background_percent"))
            direction = norm_direction(row.get("direction", ""))
            set_id = row.get("set_id", "")
            dataset = row.get("dataset", "")
            motif_name = row.get("motif_name", "")
            score = safe_score(q_value, p_value, score_cap)
            for family in classify_families(motif_name):
                rows.append(
                    {
                        "source_file": str(path),
                        "broad_context": broad_context,
                        "analysis_group": f"{broad_context}|{direction}",
                        "set_id": set_id,
                        "dataset": dataset,
                        "genome": row.get("genome", ""),
                        "level": row.get("level", ""),
                        "direction": direction,
                        "rank": rank,
                        "motif_name": motif_name,
                        "motif_tf": motif_tf(motif_name),
                        "motif_family": family,
                        "consensus": row.get("consensus", ""),
                        "p_value": p_value,
                        "q_value": q_value,
                        "score_capped": score,
                        "target_percent": target_pct,
                        "background_percent": background_pct,
                        "target_minus_background_percent": target_pct - background_pct
                        if not math.isnan(target_pct) and not math.isnan(background_pct)
                        else float("nan"),
                    }
                )
    return rows


def summarize_groups(
    frame: pd.DataFrame,
    group_cols: list[str],
    q_threshold: float,
    suggestive_q_threshold: float,
) -> pd.DataFrame:
    summaries: list[dict[str, object]] = []
    for keys, group in frame.groupby(group_cols, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        key_values = dict(zip(group_cols, keys))
        group = group.copy()
        group["q_for_sort"] = group["q_value"].fillna(1.0)
        group = group.sort_values(
            ["score_capped", "q_for_sort", "rank", "target_minus_background_percent"],
            ascending=[False, True, True, False],
        )
        best = group.iloc[0]
        summaries.append(
            {
                **key_values,
                "n_motif_rows": len(group),
                "n_unique_sets": group["set_id"].nunique(),
                "n_significant_rows_q05": int((group["q_value"] <= q_threshold).sum()),
                "n_suggestive_rows_q25": int((group["q_value"] <= suggestive_q_threshold).sum()),
                "n_significant_sets_q05": int(group.loc[group["q_value"] <= q_threshold, "set_id"].nunique()),
                "best_score_capped": float(best["score_capped"]),
                "best_rank": int(best["rank"]),
                "best_q_value": float(best["q_value"]) if not math.isnan(float(best["q_value"])) else float("nan"),
                "best_p_value": float(best["p_value"]) if not math.isnan(float(best["p_value"])) else float("nan"),
                "best_target_minus_background_percent": float(best["target_minus_background_percent"])
                if not math.isnan(float(best["target_minus_background_percent"]))
                else float("nan"),
                "best_motif_name": best["motif_name"],
                "best_set_id": best["set_id"],
                "best_dataset": best["dataset"],
                "support_call": support_call(group, q_threshold, suggestive_q_threshold),
            }
        )
    return pd.DataFrame(summaries)


def support_call(group: pd.DataFrame, q_threshold: float, suggestive_q_threshold: float) -> str:
    n_sig_sets = group.loc[group["q_value"] <= q_threshold, "set_id"].nunique()
    n_suggestive_sets = group.loc[group["q_value"] <= suggestive_q_threshold, "set_id"].nunique()
    best_q = group["q_value"].min(skipna=True)
    if n_sig_sets >= 2:
        return "recurrent_significant"
    if n_sig_sets == 1 and best_q <= 0.01:
        return "strong_single_context"
    if n_sig_sets == 1:
        return "significant_single_context"
    if n_suggestive_sets >= 1:
        return "suggestive_only"
    return "not_supported"


def write_tsv(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, na_rep="")


def make_heatmap(group_summary: pd.DataFrame, outpath: Path) -> pd.DataFrame:
    plot = group_summary[group_summary["analysis_group"].isin(GROUP_ORDER)].copy()
    matrix = plot.pivot_table(
        index="motif_family",
        columns="analysis_group",
        values="best_score_capped",
        aggfunc="max",
        fill_value=0.0,
    )
    families = [family for family in FAMILY_ORDER if family in matrix.index]
    cols = [col for col in GROUP_ORDER if col in matrix.columns]
    matrix = matrix.reindex(families)[cols]

    fig_width = max(10, 0.8 * len(cols) + 3)
    fig_height = max(7, 0.32 * len(families) + 2)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    image = ax.imshow(matrix.values, aspect="auto", cmap="magma", vmin=0, vmax=10)
    ax.set_xticks(np.arange(len(cols)))
    ax.set_xticklabels([GROUP_LABELS.get(col, col) for col in cols], rotation=45, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(families)))
    ax.set_yticklabels(families, fontsize=8)
    ax.set_title("HOMER motif-family concordance across DNMT3L/DS mechanism contexts")
    cbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("-log10(q) capped at 10")
    fig.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=240)
    plt.close(fig)

    return matrix.reset_index()


def make_top_family_table(group_summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for analysis_group, group in group_summary.groupby("analysis_group"):
        group = group.sort_values(
            ["best_score_capped", "n_significant_sets_q05", "best_rank"],
            ascending=[False, False, True],
        )
        for rank, (_, row) in enumerate(group.head(8).iterrows(), start=1):
            rows.append(
                {
                    "analysis_group": analysis_group,
                    "rank_within_group": rank,
                    "motif_family": row["motif_family"],
                    "support_call": row["support_call"],
                    "best_score_capped": row["best_score_capped"],
                    "best_q_value": row["best_q_value"],
                    "best_motif_name": row["best_motif_name"],
                    "best_set_id": row["best_set_id"],
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    figdir = Path(args.figdir)

    rows: list[dict[str, object]] = []
    for broad_context, path in INPUTS:
        rows.extend(read_known_motifs(Path(path), broad_context, args.score_cap))
    frame = pd.DataFrame(rows)
    if frame.empty:
        raise SystemExit("No HOMER rows found")

    set_summary = summarize_groups(
        frame,
        ["broad_context", "analysis_group", "set_id", "dataset", "direction", "motif_family"],
        args.q_threshold,
        args.suggestive_q_threshold,
    )
    group_summary = summarize_groups(
        frame,
        ["broad_context", "analysis_group", "direction", "motif_family"],
        args.q_threshold,
        args.suggestive_q_threshold,
    )
    top_family = make_top_family_table(group_summary)
    matrix = make_heatmap(group_summary, figdir / "homer_cross_context_motif_family_heatmap.png")

    write_tsv(outdir / "homer_cross_context_motif_family_rows.tsv", frame)
    write_tsv(outdir / "homer_cross_context_motif_family_set_summary.tsv", set_summary)
    write_tsv(outdir / "homer_cross_context_motif_family_group_summary.tsv", group_summary)
    write_tsv(outdir / "homer_cross_context_motif_family_top_by_group.tsv", top_family)
    write_tsv(outdir / "homer_cross_context_motif_family_matrix.tsv", matrix)

    print(f"Read {len(frame)} motif-family rows from {len(INPUTS)} HOMER summary files")
    print(f"Wrote motif concordance tables to {outdir}")
    print(f"Wrote heatmap to {figdir / 'homer_cross_context_motif_family_heatmap.png'}")


if __name__ == "__main__":
    main()
