#!/usr/bin/env python3
"""Bridge final candidate mechanism classes to GSE305153 fetal cortex cell types."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from scipy import stats
except Exception:  # pragma: no cover - scipy is expected but keep the script usable.
    stats = None


MECHANISM_LABELS = {
    "cross_species_dnmt3a_linked_methylation_plus_developmental_cpg_tf_context": "DNMT3A-methyl. + dev/CpG",
    "cross_species_dnmt3a_linked_methylation_plus_developmental_tf_context": "DNMT3A-methyl. + dev TF",
    "cross_species_dnmt3a_linked_methylation_without_selected_tf_instance": "DNMT3A-methyl. only",
    "mouse_dnmt3a_linked_methylation_with_tf_context_human_unresolved": "Mouse methyl. + TF; human unresolved",
}

TARGET_GENES = ["DNMT3L", "DNMT3A", "DNMT3B", "STAT1", "STAT3"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--candidate-effects",
        default="results/cross_species_candidate_expression/gse305153_candidate_celltype_effects.tsv",
    )
    parser.add_argument(
        "--target-effects",
        default="results/gse305153_cellxgene_targets/target_celltype_effect_summary.tsv",
    )
    parser.add_argument(
        "--mechanism-axis",
        default="results/mechanism_question/final_candidate_mechanism_axis_summary.tsv",
    )
    parser.add_argument("--outdir", default="results/mechanism_question")
    parser.add_argument("--figdir", default="figures/mechanism_question")
    parser.add_argument(
        "--note",
        default="notes/gse305153_candidate_mechanism_celltype_bridge.md",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def write_tsv(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n")


def safe_float_series(values: pd.Series) -> pd.Series:
    return pd.to_numeric(values, errors="coerce")


def fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}" if math.isfinite(value) else ""
    if isinstance(value, (np.floating,)):
        number = float(value)
        return f"{number:.6g}" if math.isfinite(number) else ""
    if isinstance(value, (np.integer,)):
        return str(int(value))
    return str(value)


def write_table(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})


def effect_direction(delta: float) -> str:
    if not math.isfinite(delta):
        return "unresolved"
    if delta >= 0.25:
        return "up"
    if delta <= -0.25:
        return "down"
    return "flat"


def mechanism_order(frame: pd.DataFrame) -> list[str]:
    preferred = [
        "cross_species_dnmt3a_linked_methylation_plus_developmental_cpg_tf_context",
        "cross_species_dnmt3a_linked_methylation_plus_developmental_tf_context",
        "cross_species_dnmt3a_linked_methylation_without_selected_tf_instance",
        "mouse_dnmt3a_linked_methylation_with_tf_context_human_unresolved",
    ]
    seen = set(frame["final_mechanism_axis_call"])
    return [item for item in preferred if item in seen] + sorted(seen - set(preferred))


def make_heatmap(
    matrix: pd.DataFrame,
    path_prefix: Path,
    title: str,
    cbar_label: str = "DS - control delta log2(CPM + 1)",
    vmax: float | None = None,
) -> None:
    if matrix.empty:
        return
    path_prefix.parent.mkdir(parents=True, exist_ok=True)
    values = matrix.to_numpy(dtype=float)
    if vmax is None:
        finite = values[np.isfinite(values)]
        vmax = float(np.nanpercentile(np.abs(finite), 95)) if finite.size else 1.0
        vmax = max(0.25, min(1.5, vmax))
    fig_w = max(8, 0.46 * matrix.shape[1] + 3.2)
    fig_h = max(5, 0.34 * matrix.shape[0] + 2.4)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(values, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_title(title, fontsize=13, pad=12)
    ax.set_xticks(range(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(matrix.shape[0]))
    ax.set_yticklabels(matrix.index, fontsize=9)
    ax.set_xticks(np.arange(-0.5, matrix.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, matrix.shape[0], 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.026, pad=0.02)
    cbar.set_label(cbar_label)
    fig.tight_layout()
    fig.savefig(path_prefix.with_suffix(".png"), dpi=300)
    fig.savefig(path_prefix.with_suffix(".pdf"))
    plt.close(fig)


def summarize_mechanism_by_celltype(final_effects: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    group_cols = ["label_type", "label", "final_mechanism_axis_call", "mechanism_axis_short"]
    for keys, sub in final_effects.groupby(group_cols, dropna=False, sort=False):
        label_type, label, call, short = keys
        deltas = safe_float_series(sub["delta_ds_minus_control_log2_cpm_plus_1"]).dropna()
        detect_deltas = safe_float_series(sub["delta_detected_fraction"]).dropna()
        directions = [effect_direction(float(x)) for x in deltas]
        rows.append(
            {
                "label_type": label_type,
                "label": label,
                "final_mechanism_axis_call": call,
                "mechanism_axis_short": short,
                "n_genes": int(sub["human_symbol"].nunique()),
                "mean_delta_log2_cpm": float(deltas.mean()) if len(deltas) else math.nan,
                "median_delta_log2_cpm": float(deltas.median()) if len(deltas) else math.nan,
                "mean_delta_detected_fraction": float(detect_deltas.mean()) if len(detect_deltas) else math.nan,
                "n_up_delta_ge_0_25": directions.count("up"),
                "n_down_delta_le_minus_0_25": directions.count("down"),
                "n_flat": directions.count("flat"),
                "genes_up_delta_ge_0_25": ";".join(sub.loc[deltas.index[deltas >= 0.25], "human_symbol"].tolist()),
                "genes_down_delta_le_minus_0_25": ";".join(sub.loc[deltas.index[deltas <= -0.25], "human_symbol"].tolist()),
            }
        )
    return pd.DataFrame(rows)


def summarize_targets(target_effects: pd.DataFrame) -> pd.DataFrame:
    subset = target_effects[
        (target_effects["label_type"].isin(["author_cell_type", "cell_type", "cell_class", "all_cells"]))
        & (target_effects["context"] == "all_pcw")
        & (target_effects["gene_symbol"].isin(TARGET_GENES))
    ].copy()
    subset["target_delta_log2_cpm"] = safe_float_series(subset["delta_ds_minus_control_log2_cpm"])
    subset["target_delta_detected_fraction"] = safe_float_series(subset["delta_ds_minus_control_detected_fraction"])
    keep = [
        "label_type",
        "label",
        "context",
        "gene_symbol",
        "n_ds_donors",
        "n_control_donors",
        "n_ds_cells",
        "n_control_cells",
        "target_delta_log2_cpm",
        "target_delta_detected_fraction",
    ]
    return subset[keep].rename(columns={"gene_symbol": "target_gene"})


def build_bridge_table(mech_summary: pd.DataFrame, target_summary: pd.DataFrame) -> pd.DataFrame:
    author_mech = mech_summary[mech_summary["label_type"] == "author_cell_type"].copy()
    target_wide = target_summary[target_summary["label_type"] == "author_cell_type"].pivot_table(
        index=["label_type", "label"],
        columns="target_gene",
        values="target_delta_log2_cpm",
        aggfunc="first",
    )
    target_wide.columns = [f"{gene}_delta_log2_cpm" for gene in target_wide.columns]
    target_wide = target_wide.reset_index()

    detect_wide = target_summary[target_summary["label_type"] == "author_cell_type"].pivot_table(
        index=["label_type", "label"],
        columns="target_gene",
        values="target_delta_detected_fraction",
        aggfunc="first",
    )
    detect_wide.columns = [f"{gene}_delta_detected_fraction" for gene in detect_wide.columns]
    detect_wide = detect_wide.reset_index()

    bridge = author_mech.merge(target_wide, on=["label_type", "label"], how="left").merge(
        detect_wide, on=["label_type", "label"], how="left"
    )
    return bridge


def correlations(bridge: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    targets = [f"{gene}_delta_log2_cpm" for gene in TARGET_GENES]
    for call, sub in bridge.groupby("final_mechanism_axis_call", sort=False):
        y = safe_float_series(sub["mean_delta_log2_cpm"])
        for target in targets:
            if target not in sub:
                continue
            x = safe_float_series(sub[target])
            valid = x.notna() & y.notna()
            if valid.sum() < 4:
                continue
            pearson_r = float(x[valid].corr(y[valid], method="pearson"))
            spearman_r = float(x[valid].corr(y[valid], method="spearman"))
            pearson_p = math.nan
            spearman_p = math.nan
            if stats is not None:
                pearson_p = float(stats.pearsonr(x[valid], y[valid]).pvalue)
                spearman_p = float(stats.spearmanr(x[valid], y[valid]).pvalue)
            rows.append(
                {
                    "final_mechanism_axis_call": call,
                    "mechanism_axis_short": sub["mechanism_axis_short"].iloc[0],
                    "target_gene": target.replace("_delta_log2_cpm", ""),
                    "n_author_cell_types": int(valid.sum()),
                    "pearson_r": pearson_r,
                    "pearson_p": pearson_p,
                    "spearman_r": spearman_r,
                    "spearman_p": spearman_p,
                }
            )
    return pd.DataFrame(rows)


def markdown_summary(
    path: Path,
    final_effects: pd.DataFrame,
    mech_summary: pd.DataFrame,
    bridge: pd.DataFrame,
    corr: pd.DataFrame,
) -> None:
    author_mech = mech_summary[mech_summary["label_type"] == "author_cell_type"].copy()
    author_effects = final_effects[final_effects["label_type"] == "author_cell_type"].copy()
    author_effects["delta_numeric"] = safe_float_series(author_effects["delta_ds_minus_control_log2_cpm_plus_1"])
    class_counts = final_effects.drop_duplicates(["human_symbol", "final_mechanism_axis_call"])[
        "final_mechanism_axis_call"
    ].value_counts()
    all_candidates = author_effects.groupby("label", dropna=False)["delta_numeric"].mean().sort_values(ascending=False)
    top_labels = all_candidates.head(5)
    bottom_labels = all_candidates.tail(5).sort_values()

    strongest = bridge.sort_values("mean_delta_log2_cpm", ascending=False).head(8)
    strongest_rows = [
        f"- {row.label}, {row.mechanism_axis_short}: mean delta {row.mean_delta_log2_cpm:.3g}; "
        f"up genes {int(row.n_up_delta_ge_0_25)}, down genes {int(row.n_down_delta_le_minus_0_25)}"
        for row in strongest.itertuples(index=False)
    ]
    strongest_multigene = bridge[bridge["n_genes"] >= 4].sort_values("mean_delta_log2_cpm", ascending=False).head(8)
    strongest_multigene_rows = [
        f"- {row.label}, {row.mechanism_axis_short}: mean delta {row.mean_delta_log2_cpm:.3g}; "
        f"up genes {int(row.n_up_delta_ge_0_25)}, down genes {int(row.n_down_delta_le_minus_0_25)}"
        for row in strongest_multigene.itertuples(index=False)
    ]

    corr_focus = corr[corr["target_gene"].isin(["DNMT3L", "STAT3", "STAT1"])].copy()
    corr_focus["abs_spearman"] = corr_focus["spearman_r"].abs()
    corr_focus = corr_focus.sort_values("abs_spearman", ascending=False).head(6)
    corr_rows = [
        f"- {row.target_gene} vs {row.mechanism_axis_short}: Spearman r={row.spearman_r:.3g}, p={row.spearman_p:.3g}"
        for row in corr_focus.itertuples(index=False)
    ]

    lines = [
        "# GSE305153 Candidate Mechanism Cell-Type Bridge",
        "",
        "Question: do final candidate mechanism classes show fetal-cortex cell-type expression shifts that track DNMT3L or STAT-axis expression in DS?",
        "",
        "## Inputs",
        "",
        "- Candidate expression effects: `results/cross_species_candidate_expression/gse305153_candidate_celltype_effects.tsv`.",
        "- Target expression effects: `results/gse305153_cellxgene_targets/target_celltype_effect_summary.tsv`.",
        "- Mechanism classes: `results/mechanism_question/final_candidate_mechanism_axis_summary.tsv`.",
        "",
        "## Candidate Coverage",
        "",
        f"- Final candidates tested: {final_effects['human_symbol'].nunique()} of 24.",
        f"- Author cell types tested: {final_effects.loc[final_effects['label_type'] == 'author_cell_type', 'label'].nunique()}.",
    ]
    for call, count in class_counts.items():
        plural = "gene" if int(count) == 1 else "genes"
        lines.append(f"- {MECHANISM_LABELS.get(call, call)}: {int(count)} {plural}.")

    lines.extend(
        [
            "",
            "## Main Cell-Type Pattern",
            "",
            "Across final candidates, the most DS-positive gene-weighted author-cell-type means were:",
        ]
    )
    lines.extend([f"- {label}: mean candidate delta {value:.3g}" for label, value in top_labels.items()])
    lines.extend(["", "The most DS-negative author-cell-type means were:"])
    lines.extend([f"- {label}: mean candidate delta {value:.3g}" for label, value in bottom_labels.items()])

    lines.extend(["", "## Strongest Mechanism-Class/Cell-Type Shifts", ""])
    lines.append("All mechanism-class rows, including singleton classes:")
    lines.extend(strongest_rows)
    lines.extend(["", "Restricting to multi-gene classes (`n >= 4`):"])
    lines.extend(strongest_multigene_rows)
    lines.extend(["", "## Target-Axis Correlations", ""])
    if corr_rows:
        lines.extend(corr_rows)
    else:
        lines.append("- No target-axis correlations had enough author-cell-type observations.")

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This is a cell-type-resolved expression bridge, not a methylation or chromVAR test. It asks whether candidate-gene expression shifts concentrate in the same fetal cortical cell classes where DNMT3L or STAT genes shift. Singleton mechanism classes should be read as individual-gene behavior, not module behavior. The target-axis correlations are descriptive across only 15 author cell types and do not show a strong DNMT3L/STAT tracking pattern. The analysis supports cell-context prioritization for the mechanism model, but it does not prove methylation independence. The strongest methylation inference still comes from DNMT3L/Dnmt3l methylation, DNMT3A-reference overlap, and HOMER candidate-locus evidence.",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    figdir = Path(args.figdir)

    mechanism = read_tsv(Path(args.mechanism_axis))
    mechanism["mechanism_axis_short"] = mechanism["final_mechanism_axis_call"].map(
        lambda value: MECHANISM_LABELS.get(value, value)
    )
    final_genes = set(mechanism["gene_symbol"])

    candidate_effects = read_tsv(Path(args.candidate_effects))
    final_effects = candidate_effects[candidate_effects["human_symbol"].isin(final_genes)].copy()
    final_effects = final_effects.merge(
        mechanism[
            [
                "gene_symbol",
                "tier",
                "candidate_group",
                "base_mechanism_class",
                "final_mechanism_axis_call",
                "mechanism_axis_short",
                "human_methylation_axis",
                "mouse_dnmt3l_methylation_axis",
                "mouse_dnmt3a_reference_axis",
                "developmental_tf_instance_axis",
                "cpg_ctcf_instance_axis",
                "stat_sequence_control_axis",
            ]
        ],
        left_on="human_symbol",
        right_on="gene_symbol",
        how="left",
    )
    final_effects["delta_ds_minus_control_log2_cpm_plus_1_numeric"] = safe_float_series(
        final_effects["delta_ds_minus_control_log2_cpm_plus_1"]
    )
    final_effects["delta_detected_fraction_numeric"] = safe_float_series(final_effects["delta_detected_fraction"])

    mech_summary = summarize_mechanism_by_celltype(final_effects)
    target_summary = summarize_targets(read_tsv(Path(args.target_effects)))
    bridge = build_bridge_table(mech_summary, target_summary)
    corr = correlations(bridge)

    final_effect_fields = [
        "dataset",
        "label_type",
        "label",
        "context",
        "human_symbol",
        "tier",
        "candidate_group",
        "base_mechanism_class",
        "final_mechanism_axis_call",
        "mechanism_axis_short",
        "n_ds_donors",
        "n_control_donors",
        "n_ds_cells",
        "n_control_cells",
        "mean_ds_log2_cpm_plus_1",
        "mean_control_log2_cpm_plus_1",
        "delta_ds_minus_control_log2_cpm_plus_1",
        "mean_ds_detected_fraction",
        "mean_control_detected_fraction",
        "delta_detected_fraction",
        "direction",
        "human_methylation_axis",
        "mouse_dnmt3l_methylation_axis",
        "mouse_dnmt3a_reference_axis",
        "developmental_tf_instance_axis",
        "cpg_ctcf_instance_axis",
        "stat_sequence_control_axis",
    ]
    write_tsv(outdir / "gse305153_final_candidate_celltype_mechanism_effects.tsv", final_effects[final_effect_fields])
    write_tsv(outdir / "gse305153_mechanism_axis_celltype_summary.tsv", mech_summary)
    write_tsv(outdir / "gse305153_target_axis_celltype_summary.tsv", target_summary)
    write_tsv(outdir / "gse305153_candidate_target_celltype_bridge.tsv", bridge)
    write_tsv(outdir / "gse305153_candidate_target_celltype_correlations.tsv", corr)

    author = final_effects[final_effects["label_type"] == "author_cell_type"].copy()
    gene_order = mechanism["gene_symbol"].tolist()
    labels = sorted(author["label"].unique())
    gene_matrix = author.pivot_table(
        index="human_symbol",
        columns="label",
        values="delta_ds_minus_control_log2_cpm_plus_1_numeric",
        aggfunc="first",
    ).reindex(index=gene_order, columns=labels)
    make_heatmap(
        gene_matrix,
        figdir / "gse305153_final_candidate_author_cell_type_effect_heatmap",
        "GSE305153 final candidates by fetal cortical cell type",
        vmax=1.2,
    )

    class_order = mechanism_order(mech_summary)
    class_matrix = mech_summary[mech_summary["label_type"] == "author_cell_type"].pivot_table(
        index="final_mechanism_axis_call",
        columns="label",
        values="mean_delta_log2_cpm",
        aggfunc="first",
    ).reindex(index=class_order, columns=labels)
    class_sizes = (
        final_effects.drop_duplicates(["human_symbol", "final_mechanism_axis_call"])
        .groupby("final_mechanism_axis_call")["human_symbol"]
        .nunique()
        .to_dict()
    )
    class_matrix.index = [
        f"{MECHANISM_LABELS.get(value, value)} (n={class_sizes.get(value, 0)})" for value in class_matrix.index
    ]
    make_heatmap(
        class_matrix,
        figdir / "gse305153_mechanism_axis_author_cell_type_heatmap",
        "GSE305153 mechanism classes by fetal cortical cell type",
        vmax=0.8,
    )

    markdown_summary(Path(args.note), final_effects, mech_summary, bridge, corr)

    print(f"Wrote {final_effects['human_symbol'].nunique()} final-candidate genes")
    print(outdir / "gse305153_mechanism_axis_celltype_summary.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
