#!/usr/bin/env python3
"""Classify candidate methylation links by stricter distance windows.

The main candidate tables include promoter, gene-body, and nearest-gene
links within 100 kb. This script tests how much final-candidate support depends
on the broad 100 kb window by summarizing promoter/gene-body and nested TSS
distance windows for human CpG links and mouse DMR links.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]

FINAL_CANDIDATES = (
    ROOT
    / "results"
    / "final_candidate_evidence"
    / "final_candidate_evidence_matrix_with_regulatory_context.tsv"
)
HUMAN_LINKS = (
    ROOT
    / "results"
    / "human_methylation_convergence"
    / "human_methylation_probe_gene_links.tsv"
)
MOUSE_LINKS = (
    ROOT
    / "results"
    / "integrated_candidate_genes"
    / "mouse_triple_hit_dmr_gene_links.tsv"
)

OUT_DIR = ROOT / "results" / "candidate_link_window_sensitivity"
SUPP_DIR = ROOT / "manuscript" / "Supplementary Tables"


def yes(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.lower().eq("yes")


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def count_unique(df: pd.DataFrame, feature_col: str, mask: pd.Series) -> int:
    if df.empty:
        return 0
    return int(df.loc[mask.fillna(False), feature_col].astype(str).nunique())


def feature_window_class(row: pd.Series) -> str:
    if str(row.get("promoter_2kb", "")).lower() == "yes":
        return "promoter_2kb"
    if str(row.get("gene_body", "")).lower() == "yes":
        return "gene_body"

    distance = pd.to_numeric(pd.Series([row.get("abs_distance_to_tss")]), errors="coerce").iloc[0]
    if pd.isna(distance):
        return "distance_unavailable"
    if distance <= 10_000:
        return "tss_within_10kb"
    if distance <= 25_000:
        return "tss_within_25kb"
    if distance <= 50_000:
        return "tss_within_50kb"
    if distance <= 100_000:
        return "tss_within_100kb"
    return "outside_100kb"


def summarize_links(df: pd.DataFrame, feature_col: str) -> dict[str, int]:
    if df.empty:
        return {
            "total": 0,
            "promoter_2kb": 0,
            "gene_body": 0,
            "promoter_or_gene_body": 0,
            "tss_within_10kb": 0,
            "tss_within_25kb": 0,
            "tss_within_50kb": 0,
            "tss_within_100kb": 0,
            "nearest_100kb_only": 0,
            "distal_50_100kb_only": 0,
        }

    promoter = yes(df.get("promoter_2kb", pd.Series(index=df.index, dtype=object)))
    body = yes(df.get("gene_body", pd.Series(index=df.index, dtype=object)))
    nearest = yes(df.get("nearest_100kb", pd.Series(index=df.index, dtype=object)))
    relationship = df.get("relationship", pd.Series(index=df.index, dtype=object)).fillna("").astype(str)
    nearest = nearest | relationship.eq("nearest_100kb")
    distance = numeric(df.get("abs_distance_to_tss", pd.Series(index=df.index, dtype=object)))
    strong = promoter | body

    return {
        "total": int(df[feature_col].astype(str).nunique()),
        "promoter_2kb": count_unique(df, feature_col, promoter),
        "gene_body": count_unique(df, feature_col, body),
        "promoter_or_gene_body": count_unique(df, feature_col, strong),
        "tss_within_10kb": count_unique(df, feature_col, distance.le(10_000)),
        "tss_within_25kb": count_unique(df, feature_col, distance.le(25_000)),
        "tss_within_50kb": count_unique(df, feature_col, distance.le(50_000)),
        "tss_within_100kb": count_unique(df, feature_col, distance.le(100_000)),
        "nearest_100kb_only": count_unique(df, feature_col, nearest & ~strong),
        "distal_50_100kb_only": count_unique(
            df, feature_col, distance.gt(50_000) & distance.le(100_000) & ~strong
        ),
    }


def classify_gene(row: pd.Series) -> str:
    high_conf = (
        row["human_promoter_or_gene_body"] > 0
        or row["mouse_promoter_or_gene_body"] > 0
    )
    if high_conf:
        return "robust_promoter_or_gene_body"

    within_25 = row["human_tss_within_25kb"] > 0 or row["mouse_tss_within_25kb"] > 0
    if within_25:
        return "robust_within_25kb"

    within_50 = row["human_tss_within_50kb"] > 0 or row["mouse_tss_within_50kb"] > 0
    if within_50:
        return "proximal_50kb_only"

    within_100 = (
        row["human_tss_within_100kb"] > 0
        or row["mouse_tss_within_100kb"] > 0
        or row["human_nearest_100kb_only"] > 0
        or row["mouse_nearest_100kb_only"] > 0
    )
    if within_100:
        return "exploratory_100kb_only"

    return "no_strict_methylation_link_in_current_filters"


def subset_human_links(human: pd.DataFrame, gene_symbol: str) -> pd.DataFrame:
    rows = human.loc[human["gene_symbol"].astype(str).eq(gene_symbol)].copy()
    if "human_triple_phase3_any_ds" in rows:
        rows = rows.loc[yes(rows["human_triple_phase3_any_ds"])].copy()
    return rows


def subset_mouse_links(mouse: pd.DataFrame, mouse_symbol: str) -> pd.DataFrame:
    return mouse.loc[mouse["gene_symbol"].astype(str).eq(mouse_symbol)].copy()


def build_feature_details(
    final: pd.DataFrame, human: pd.DataFrame, mouse: pd.DataFrame
) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for _, cand in final.iterrows():
        gene = cand["gene_symbol"]
        mouse_symbol = cand["mouse_symbol"]
        tier = cand["tier"]

        human_rows = subset_human_links(human, gene)
        for _, row in human_rows.iterrows():
            records.append(
                {
                    "species": "human",
                    "feature_id": row.get("probe_id", ""),
                    "gene_symbol": gene,
                    "mouse_symbol": mouse_symbol,
                    "tier": tier,
                    "feature_direction": row.get("direction", ""),
                    "relationship": row.get("relationship", ""),
                    "promoter_2kb": row.get("promoter_2kb", ""),
                    "gene_body": row.get("gene_body", ""),
                    "nearest_100kb": row.get("nearest_100kb", ""),
                    "abs_distance_to_tss": row.get("abs_distance_to_tss", ""),
                    "distance_to_gene_body": row.get("distance_to_gene_body", ""),
                    "window_class": feature_window_class(row),
                    "dnmt3l_oe_methylation_delta": row.get("gse157078_delta_beta", ""),
                    "validation_methylation_delta": row.get("shsy_phase3_signed_delta", ""),
                    "ds_signed_methylation_delta": row.get("ds_mean_signed_delta", ""),
                    "expression_log2fc": "",
                    "expression_context": "",
                }
            )

        mouse_rows = subset_mouse_links(mouse, mouse_symbol)
        for _, row in mouse_rows.iterrows():
            records.append(
                {
                    "species": "mouse",
                    "feature_id": row.get("dmr_id", ""),
                    "gene_symbol": gene,
                    "mouse_symbol": mouse_symbol,
                    "tier": tier,
                    "feature_direction": row.get("dmr_direction", ""),
                    "relationship": row.get("relationship", ""),
                    "promoter_2kb": row.get("promoter_2kb", ""),
                    "gene_body": row.get("gene_body", ""),
                    "nearest_100kb": row.get("nearest_100kb", ""),
                    "abs_distance_to_tss": row.get("abs_distance_to_tss", ""),
                    "distance_to_gene_body": row.get("distance_to_gene_body", ""),
                    "window_class": feature_window_class(row),
                    "dnmt3l_oe_methylation_delta": row.get("dmr_mean_delta_beta", ""),
                    "validation_methylation_delta": "",
                    "ds_signed_methylation_delta": "",
                    "expression_log2fc": row.get("gse268951_log2fc_oe_vs_control", ""),
                    "expression_context": "GSE268951_Dnmt3l_OE_public_count_DE",
                }
            )

    return pd.DataFrame.from_records(records)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    SUPP_DIR.mkdir(parents=True, exist_ok=True)

    final = pd.read_csv(FINAL_CANDIDATES, sep="\t")
    human = pd.read_csv(HUMAN_LINKS, sep="\t", low_memory=False)
    mouse = pd.read_csv(MOUSE_LINKS, sep="\t", low_memory=False)

    rows: list[dict[str, object]] = []
    for _, cand in final.iterrows():
        gene = cand["gene_symbol"]
        mouse_symbol = cand["mouse_symbol"]
        human_counts = summarize_links(subset_human_links(human, gene), "probe_id")
        mouse_counts = summarize_links(subset_mouse_links(mouse, mouse_symbol), "dmr_id")
        record: dict[str, object] = {
            "gene_symbol": gene,
            "mouse_symbol": mouse_symbol,
            "tier": cand["tier"],
            "manuscript_use": cand.get("manuscript_use", ""),
        }
        record.update({f"human_{key}": value for key, value in human_counts.items()})
        record.update({f"mouse_{key}": value for key, value in mouse_counts.items()})
        rows.append(record)

    sensitivity = pd.DataFrame.from_records(rows)
    sensitivity["any_promoter_or_gene_body_link"] = (
        (sensitivity["human_promoter_or_gene_body"] > 0)
        | (sensitivity["mouse_promoter_or_gene_body"] > 0)
    ).map({True: "yes", False: "no"})
    sensitivity["any_25kb_or_stronger_link"] = (
        (sensitivity["any_promoter_or_gene_body_link"] == "yes")
        | (sensitivity["human_tss_within_25kb"] > 0)
        | (sensitivity["mouse_tss_within_25kb"] > 0)
    ).map({True: "yes", False: "no"})
    sensitivity["window_sensitivity_class"] = sensitivity.apply(classify_gene, axis=1)

    detail = build_feature_details(final, human, mouse)

    summary_records: list[dict[str, object]] = []
    tier_order = ["all_final_candidates"] + list(final["tier"].drop_duplicates())
    for tier in tier_order:
        if tier == "all_final_candidates":
            sub = sensitivity
        else:
            sub = sensitivity.loc[sensitivity["tier"].eq(tier)]
        summary_records.append(
            {
                "group": tier,
                "n_genes": int(len(sub)),
                "n_with_human_support": int((sub["human_total"] > 0).sum()),
                "n_with_mouse_support": int((sub["mouse_total"] > 0).sum()),
                "n_with_any_promoter_or_gene_body_link": int(
                    (sub["any_promoter_or_gene_body_link"] == "yes").sum()
                ),
                "n_with_any_25kb_or_stronger_link": int(
                    (sub["any_25kb_or_stronger_link"] == "yes").sum()
                ),
                "n_exploratory_100kb_only": int(
                    sub["window_sensitivity_class"].eq("exploratory_100kb_only").sum()
                ),
                "window_classes": "; ".join(
                    f"{klass}:{count}"
                    for klass, count in sub["window_sensitivity_class"].value_counts().sort_index().items()
                ),
            }
        )
    summary = pd.DataFrame.from_records(summary_records)

    sensitivity_path = OUT_DIR / "candidate_methylation_link_window_sensitivity.tsv"
    detail_path = OUT_DIR / "candidate_methylation_link_window_feature_details.tsv"
    summary_path = OUT_DIR / "candidate_methylation_link_window_summary.tsv"
    sensitivity.to_csv(sensitivity_path, sep="\t", index=False)
    detail.to_csv(detail_path, sep="\t", index=False)
    summary.to_csv(summary_path, sep="\t", index=False)

    sensitivity.to_csv(
        SUPP_DIR / "ST8h_candidate_methylation_link_window_sensitivity.tsv",
        sep="\t",
        index=False,
    )
    detail.to_csv(
        SUPP_DIR / "ST8i_candidate_methylation_link_window_feature_details.tsv",
        sep="\t",
        index=False,
    )

    print(f"Wrote {sensitivity_path}")
    print(f"Wrote {detail_path}")
    print(f"Wrote {summary_path}")
    print("Summary:")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
