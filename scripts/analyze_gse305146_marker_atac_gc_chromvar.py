#!/usr/bin/env python3
"""GC/content-matched chromVAR-style motif deviations for GSE305146 ATAC."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analyze_gse305146_marker_atac_motif_deviation import CLASS_LABELS, CLASS_ORDER, FAMILY_ORDER
from analyze_gse305146_marker_atac_motif_instance_peaksets import MOTIF_FAMILY_LABELS
from analyze_gse305146_marker_atac_peak_da import (
    bh_fdr,
    collect_donor_pseudobulk,
    pcw_adjusted_ols,
    read_background_regions,
    safe_float,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", default="data/geo/GSE305146/metadata/sample_metadata.tsv")
    parser.add_argument("--supplementary", default="data/geo/GSE305146/supplementary")
    parser.add_argument(
        "--background-table",
        default="results/gse305146_scATAC_homer_motif_enrichment/tables/GSE305146_recurrent_peak_background.tsv",
    )
    parser.add_argument(
        "--motif-hits",
        default="results/mechanism_question/gse305146_marker_atac_motif_deviation_region_hits.tsv",
    )
    parser.add_argument(
        "--composition-table",
        default=(
            "results/mechanism_question/gse305146_marker_atac_motif_deviation/"
            "annotatePeaks/gse305146_recurrent_background_motif_deviation.annotatePeaks.nmotifs.tsv"
        ),
    )
    parser.add_argument("--genome-fasta", default="/Users/jili/tools/homer/data/genomes/hg38/genome.fa")
    parser.add_argument("--outdir", default="results/mechanism_question")
    parser.add_argument("--figdir", default="figures/mechanism_question")
    parser.add_argument("--note", default="notes/gse305146_marker_atac_gc_chromvar_first_pass.md")
    parser.add_argument("--region-size", type=int, default=200)
    parser.add_argument("--bin-size", type=int, default=100_000)
    parser.add_argument("--n-background-iterations", type=int, default=50)
    parser.add_argument("--min-total-rna", type=float, default=500.0)
    parser.add_argument("--min-detected-genes", type=int, default=200)
    parser.add_argument("--min-top-score", type=float, default=0.05)
    parser.add_argument("--min-confidence-delta", type=float, default=0.05)
    parser.add_argument("--min-cells-per-donor-class", type=int, default=20)
    parser.add_argument("--min-donors-per-group", type=int, default=3)
    parser.add_argument("--seed", type=int, default=305146)
    parser.add_argument("--include-sex-chromosomes", action="store_true")
    parser.add_argument("--limit-libraries", type=int, help="Debug: process first N final-analysis libraries.")
    parser.add_argument("--limit-families", type=int, help="Debug: process first N motif families.")
    parser.add_argument("--max-background-regions", type=int, help="Debug: use first N retained recurrent peaks.")
    parser.add_argument("--plot-only", action="store_true", help="Regenerate figures and note from existing tables.")
    return parser.parse_args()


def fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return f"{value:.6g}" if math.isfinite(value) else ""
    if isinstance(value, np.floating):
        number = float(value)
        return f"{number:.6g}" if math.isfinite(number) else ""
    if isinstance(value, np.integer):
        return str(int(value))
    return str(value)


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})


def normalize_fraction(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce").astype(float)
    finite = values[np.isfinite(values)]
    if not finite.empty and float(finite.median()) > 1.5:
        values = values / 100.0
    return values


def read_homer_composition(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["peak_id", "cpg_fraction", "gc_fraction"])
    table = pd.read_csv(path, sep="\t", usecols=[0, 19, 20])
    table = table.rename(columns={table.columns[0]: "peak_id", table.columns[1]: "cpg_fraction", table.columns[2]: "gc_fraction"})
    table["peak_id"] = table["peak_id"].astype(str)
    table["cpg_fraction"] = normalize_fraction(table["cpg_fraction"])
    table["gc_fraction"] = normalize_fraction(table["gc_fraction"])
    return table[["peak_id", "cpg_fraction", "gc_fraction"]].drop_duplicates("peak_id")


def compute_gc_for_missing(regions: pd.DataFrame, fasta_path: Path) -> pd.DataFrame:
    missing = regions[regions["gc_fraction"].isna()].copy()
    if missing.empty:
        return regions
    if not fasta_path.exists():
        raise FileNotFoundError(f"Missing genome FASTA for GC calculation: {fasta_path}")

    requested: dict[str, list[tuple[int, int, int]]] = defaultdict(list)
    for row in missing.itertuples(index=True):
        requested[str(row.chrom)].append((int(row.start), int(row.end), int(row.Index)))
    for values in requested.values():
        values.sort()

    gc_values: dict[int, float] = {}
    cpg_values: dict[int, float] = {}

    current_chrom = ""
    chunks: list[str] = []

    def process_current() -> None:
        if not current_chrom or current_chrom not in requested:
            return
        sequence = "".join(chunks).upper()
        for start, end, index in requested[current_chrom]:
            subseq = sequence[start:end]
            if not subseq:
                gc_values[index] = math.nan
                cpg_values[index] = math.nan
                continue
            n_base = sum(1 for base in subseq if base in {"A", "C", "G", "T"})
            gc = subseq.count("G") + subseq.count("C")
            cpg = subseq.count("CG")
            gc_values[index] = gc / n_base if n_base else math.nan
            cpg_values[index] = cpg / max(len(subseq) - 1, 1)

    with fasta_path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                process_current()
                current_chrom = line[1:].strip().split()[0]
                chunks = []
            elif current_chrom in requested:
                chunks.append(line.strip())
        process_current()

    out = regions.copy()
    for index, value in gc_values.items():
        out.at[index, "gc_fraction"] = value
    for index, value in cpg_values.items():
        out.at[index, "cpg_fraction"] = value
    return out


def make_qbin(values: pd.Series, q: int) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce").astype(float)
    if numeric.notna().sum() == 0 or numeric.nunique(dropna=True) <= 1:
        return pd.Series(np.zeros(len(numeric), dtype=np.int16), index=values.index)
    ranked = numeric.rank(method="first")
    binned = pd.qcut(ranked, q=min(q, int(numeric.notna().sum())), labels=False, duplicates="drop")
    return binned.fillna(0).astype(np.int16)


def add_bias_columns(regions: pd.DataFrame, donor_store: dict[tuple[str, str], dict[str, object]]) -> pd.DataFrame:
    records = [record for record in donor_store.values() if float(record["total_peak_counts"]) > 0]
    if not records:
        raise RuntimeError("No donor/class ATAC records available for accessibility bias")
    counts = np.vstack([np.asarray(record["counts"], dtype=np.float64) for record in records])
    depths = np.array([float(record["total_peak_counts"]) for record in records], dtype=np.float64)
    total_counts = counts.sum(axis=0)
    total_depth = float(depths.sum())
    mean_cpm = total_counts / total_depth * 1_000_000.0 if total_depth > 0 else np.zeros(counts.shape[1])

    out = regions.copy()
    out["total_recurrent_peak_counts"] = total_counts
    out["mean_peak_cpm"] = mean_cpm
    out["mean_log2_peak_cpm"] = np.log2(mean_cpm + 1.0)
    out["peak_expected_fraction"] = total_counts / total_depth if total_depth > 0 else 0.0
    out["gc_bin"] = make_qbin(out["gc_fraction"], 20)
    out["cpg_bin"] = make_qbin(out["cpg_fraction"], 12)
    out["accessibility_bin"] = make_qbin(out["mean_log2_peak_cpm"], 20)
    out["width_bin"] = make_qbin(out["source_width"], 8)
    out["frequency_bin"] = make_qbin(out["max_frequency"], 8)
    return out


def read_motif_foregrounds(path: Path, regions: pd.DataFrame, families: list[str]) -> dict[str, np.ndarray]:
    peak_to_index = dict(zip(regions["peak_id"].astype(str), regions.index.to_numpy()))
    usecols = ["peak_id", "motif_family"]
    foregrounds: dict[str, set[int]] = {family: set() for family in families}
    for chunk in pd.read_csv(path, sep="\t", usecols=usecols, chunksize=200_000):
        chunk = chunk.dropna()
        for family, sub in chunk.groupby("motif_family", observed=True):
            if family not in foregrounds:
                continue
            foregrounds[family].update(peak_to_index[peak_id] for peak_id in sub["peak_id"].astype(str).unique() if peak_id in peak_to_index)
    return {family: np.array(sorted(indices), dtype=np.int32) for family, indices in foregrounds.items()}


def build_pool_maps(
    regions: pd.DataFrame,
    candidate_indices: np.ndarray,
    tiers: list[tuple[str, ...]],
) -> list[dict[tuple[int, ...], np.ndarray]]:
    bin_arrays = {column: regions[column].to_numpy(dtype=np.int16) for tier in tiers for column in tier}
    pools_by_tier: list[dict[tuple[int, ...], list[int]]] = []
    for tier in tiers:
        pools: dict[tuple[int, ...], list[int]] = defaultdict(list)
        for index in candidate_indices:
            key = tuple(int(bin_arrays[column][index]) for column in tier)
            pools[key].append(int(index))
        pools_by_tier.append(pools)
    return [{key: np.array(values, dtype=np.int32) for key, values in pools.items()} for pools in pools_by_tier]


def choose_background_indices(
    regions: pd.DataFrame,
    foreground_indices: np.ndarray,
    n_iterations: int,
    seed: int,
) -> tuple[np.ndarray, dict[str, object]]:
    rng = np.random.default_rng(seed)
    all_indices = np.arange(len(regions), dtype=np.int32)
    foreground_mask = np.zeros(len(regions), dtype=bool)
    foreground_mask[foreground_indices] = True
    candidate_indices = all_indices[~foreground_mask]
    if candidate_indices.size == 0:
        raise RuntimeError("No candidate background peaks remain after excluding foreground peaks")

    tiers = [
        ("gc_bin", "accessibility_bin", "cpg_bin", "width_bin", "frequency_bin"),
        ("gc_bin", "accessibility_bin", "cpg_bin", "width_bin"),
        ("gc_bin", "accessibility_bin", "cpg_bin"),
        ("gc_bin", "accessibility_bin"),
        ("gc_bin",),
        tuple(),
    ]
    pools_by_tier = build_pool_maps(regions, candidate_indices, tiers)
    bin_arrays = {column: regions[column].to_numpy(dtype=np.int16) for tier in tiers for column in tier}
    background = np.empty((n_iterations, foreground_indices.size), dtype=np.int32)
    tier_counter: Counter[str] = Counter()

    for col, peak_index in enumerate(foreground_indices):
        choices_by_tier: list[tuple[str, np.ndarray]] = []
        for tier, pools in zip(tiers, pools_by_tier):
            key = tuple(int(bin_arrays[column][peak_index]) for column in tier)
            pool = pools.get(key)
            if pool is not None and pool.size:
                choices_by_tier.append(("+".join(tier) if tier else "global_candidate", pool))
                break
        if not choices_by_tier:
            choices_by_tier.append(("global_candidate", candidate_indices))
        tier_label, pool = choices_by_tier[0]
        sampled_positions = rng.integers(0, pool.size, size=n_iterations)
        background[:, col] = pool[sampled_positions]
        tier_counter[tier_label] += n_iterations

    metric_columns = ["gc_fraction", "cpg_fraction", "mean_log2_peak_cpm", "source_width", "max_frequency"]
    fg = regions.iloc[foreground_indices]
    bg = regions.iloc[background.reshape(-1)]
    qc: dict[str, object] = {
        "n_foreground_peaks": int(foreground_indices.size),
        "n_candidate_background_peaks": int(candidate_indices.size),
        "n_background_iterations": int(n_iterations),
    }
    for column in metric_columns:
        fg_mean = float(pd.to_numeric(fg[column], errors="coerce").mean())
        bg_mean = float(pd.to_numeric(bg[column], errors="coerce").mean())
        qc[f"foreground_mean_{column}"] = fg_mean
        qc[f"background_mean_{column}"] = bg_mean
        qc[f"mean_diff_{column}"] = fg_mean - bg_mean
    total_assignments = max(sum(tier_counter.values()), 1)
    for tier, count in sorted(tier_counter.items()):
        safe_tier = tier.replace("+", "_")
        qc[f"match_fraction_{safe_tier}"] = count / total_assignments
    return background, qc


def build_count_matrix(
    donor_store: dict[tuple[str, str], dict[str, object]],
) -> tuple[list[dict[str, object]], np.ndarray, np.ndarray]:
    records = sorted(
        (record for record in donor_store.values() if float(record["total_peak_counts"]) > 0),
        key=lambda record: (
            str(record["biological_sample_id"]),
            CLASS_ORDER.index(str(record["inferred_cell_class"])) if str(record["inferred_cell_class"]) in CLASS_ORDER else 999,
        ),
    )
    counts = np.vstack([np.asarray(record["counts"], dtype=np.float64) for record in records])
    depths = np.array([float(record["total_peak_counts"]) for record in records], dtype=np.float64)
    return records, counts, depths


def calculate_chromvar_rows(
    records: list[dict[str, object]],
    counts: np.ndarray,
    depths: np.ndarray,
    regions: pd.DataFrame,
    foregrounds: dict[str, np.ndarray],
    families: list[str],
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    peak_expected = regions["peak_expected_fraction"].to_numpy(dtype=np.float64)
    donor_rows: list[dict[str, object]] = []
    qc_rows: list[dict[str, object]] = []

    for family_index, family in enumerate(families):
        foreground_indices = foregrounds.get(family, np.array([], dtype=np.int32))
        if foreground_indices.size == 0:
            continue
        print(f"[chromVAR] {family}: {foreground_indices.size} foreground peaks", flush=True)
        background_indices, qc = choose_background_indices(
            regions,
            foreground_indices,
            args.n_background_iterations,
            args.seed + family_index * 1009,
        )
        qc_rows.append(
            {
                "motif_family": family,
                "motif_family_label": MOTIF_FAMILY_LABELS.get(family, family),
                **qc,
            }
        )

        fg_observed = counts[:, foreground_indices].sum(axis=1)
        fg_expected_fraction = float(peak_expected[foreground_indices].sum())
        fg_expected = depths * fg_expected_fraction
        fg_raw_deviation = np.divide(
            fg_observed - fg_expected,
            fg_expected,
            out=np.full_like(fg_observed, np.nan, dtype=np.float64),
            where=fg_expected > 0,
        )

        bg_deviations = np.empty((counts.shape[0], args.n_background_iterations), dtype=np.float64)
        for iteration in range(args.n_background_iterations):
            bg_indices = background_indices[iteration]
            bg_observed = counts[:, bg_indices].sum(axis=1)
            bg_expected_fraction = float(peak_expected[bg_indices].sum())
            bg_expected = depths * bg_expected_fraction
            bg_deviations[:, iteration] = np.divide(
                bg_observed - bg_expected,
                bg_expected,
                out=np.full_like(bg_observed, np.nan, dtype=np.float64),
                where=bg_expected > 0,
            )
        bg_mean = np.nanmean(bg_deviations, axis=1)
        bg_sd = np.nanstd(bg_deviations, axis=1, ddof=1)
        corrected = fg_raw_deviation - bg_mean
        z = np.divide(corrected, bg_sd, out=np.full_like(corrected, np.nan), where=bg_sd > 0)

        for record, observed, expected, raw, mean_bg, sd_bg, corr, z_value in zip(
            records,
            fg_observed,
            fg_expected,
            fg_raw_deviation,
            bg_mean,
            bg_sd,
            corrected,
            z,
        ):
            donor_rows.append(
                {
                    **{key: value for key, value in record.items() if key != "counts"},
                    "motif_family": family,
                    "motif_family_label": MOTIF_FAMILY_LABELS.get(family, family),
                    "n_foreground_peaks": int(foreground_indices.size),
                    "foreground_observed_counts": float(observed),
                    "foreground_expected_counts": float(expected),
                    "foreground_raw_deviation": float(raw),
                    "background_mean_raw_deviation": float(mean_bg),
                    "background_sd_raw_deviation": float(sd_bg),
                    "chromvar_deviation": float(corr),
                    "chromvar_deviation_z": float(z_value),
                }
            )
    return pd.DataFrame(donor_rows), pd.DataFrame(qc_rows)


def effect_rows(donor: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    donor = donor[donor["inferred_cell_class"].isin(CLASS_ORDER)].copy()
    donor["dev_stage_pcw_numeric"] = pd.to_numeric(donor["dev_stage_pcw"], errors="coerce")
    for (cell_class, family), sub in donor.groupby(["inferred_cell_class", "motif_family"], dropna=False, observed=True):
        sub = sub[(sub["n_cells"].astype(int) >= args.min_cells_per_donor_class) & (sub["chromvar_deviation_z"].notna())].copy()
        ds = sub[sub["group"] == "DS"]
        control = sub[sub["group"] == "Control"]
        if (
            ds["biological_sample_id"].nunique() < args.min_donors_per_group
            or control["biological_sample_id"].nunique() < args.min_donors_per_group
        ):
            continue
        y = sub["chromvar_deviation_z"].astype(float).to_numpy()[:, None]
        groups = np.array([1 if group == "DS" else 0 for group in sub["group"]], dtype=np.int8)
        pcw = sub["dev_stage_pcw_numeric"].astype(float).to_numpy()
        beta, t_stat, pvalue = pcw_adjusted_ols(y, groups, pcw)
        row = {
            "inferred_cell_class": cell_class,
            "motif_family": family,
            "motif_family_label": sub["motif_family_label"].iloc[0],
            "n_foreground_peaks": int(sub["n_foreground_peaks"].iloc[0]),
            "n_ds_donors": int(ds["biological_sample_id"].nunique()),
            "n_control_donors": int(control["biological_sample_id"].nunique()),
            "n_ds_cells": int(ds["n_cells"].sum()),
            "n_control_cells": int(control["n_cells"].sum()),
            "mean_ds_chromvar_z": float(ds["chromvar_deviation_z"].mean()),
            "mean_control_chromvar_z": float(control["chromvar_deviation_z"].mean()),
            "delta_mean_ds_minus_control_chromvar_z": float(ds["chromvar_deviation_z"].mean() - control["chromvar_deviation_z"].mean()),
            "pcw_adjusted_group_beta_chromvar_z": float(beta[0]),
            "pcw_adjusted_group_t_chromvar_z": float(t_stat[0]),
            "pcw_adjusted_group_p_chromvar_z": float(pvalue[0]),
            "mean_ds_chromvar_deviation": float(ds["chromvar_deviation"].mean()),
            "mean_control_chromvar_deviation": float(control["chromvar_deviation"].mean()),
            "delta_mean_ds_minus_control_chromvar_deviation": float(
                ds["chromvar_deviation"].mean() - control["chromvar_deviation"].mean()
            ),
        }
        rows.append(row)
    effects = pd.DataFrame(rows)
    if effects.empty:
        return effects
    effects["fdr_bh_chromvar_z"] = bh_fdr(pd.to_numeric(effects["pcw_adjusted_group_p_chromvar_z"], errors="coerce").to_numpy())
    effects["nominal_p_lt_0_05"] = effects["pcw_adjusted_group_p_chromvar_z"] < 0.05
    effects["fdr_lt_0_05"] = effects["fdr_bh_chromvar_z"] < 0.05
    return effects


def family_summary(effects: pd.DataFrame, match_qc: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if effects.empty:
        return pd.DataFrame()
    for family, sub in effects.groupby("motif_family", dropna=False, observed=True):
        values = pd.to_numeric(sub["pcw_adjusted_group_beta_chromvar_z"], errors="coerce")
        fdr = pd.to_numeric(sub["fdr_bh_chromvar_z"], errors="coerce")
        row: dict[str, object] = {
            "motif_family": family,
            "motif_family_label": sub["motif_family_label"].iloc[0],
            "n_foreground_peaks": int(sub["n_foreground_peaks"].iloc[0]),
            "n_inferred_classes": int(values.notna().sum()),
            "mean_pcw_adjusted_beta_chromvar_z": float(values.mean()),
            "median_pcw_adjusted_beta_chromvar_z": float(values.median()),
            "n_classes_beta_positive": int((values > 0).sum()),
            "n_classes_beta_negative": int((values < 0).sum()),
            "n_classes_fdr_positive": int(((values > 0) & (fdr < 0.05)).sum()),
            "n_classes_fdr_negative": int(((values < 0) & (fdr < 0.05)).sum()),
            "max_positive_class": str(sub.loc[values.idxmax(), "inferred_cell_class"]) if len(values) else "",
            "max_pcw_adjusted_beta_chromvar_z": float(values.max()) if len(values) else math.nan,
            "max_negative_class": str(sub.loc[values.idxmin(), "inferred_cell_class"]) if len(values) else "",
            "min_pcw_adjusted_beta_chromvar_z": float(values.min()) if len(values) else math.nan,
        }
        qc_sub = match_qc[match_qc["motif_family"] == family]
        if not qc_sub.empty:
            for column in [
                "foreground_mean_gc_fraction",
                "background_mean_gc_fraction",
                "mean_diff_gc_fraction",
                "foreground_mean_cpg_fraction",
                "background_mean_cpg_fraction",
                "mean_diff_cpg_fraction",
                "foreground_mean_mean_log2_peak_cpm",
                "background_mean_mean_log2_peak_cpm",
                "mean_diff_mean_log2_peak_cpm",
            ]:
                if column in qc_sub:
                    row[column] = qc_sub[column].iloc[0]
        rows.append(row)
    summary = pd.DataFrame(rows)
    order = {family: index for index, family in enumerate(FAMILY_ORDER)}
    summary["_order"] = summary["motif_family"].map(order).fillna(999).astype(int)
    return summary.sort_values(["_order", "motif_family"]).drop(columns=["_order"])


def plot_effect_heatmap(effects: pd.DataFrame, figdir: Path) -> None:
    if effects.empty:
        return
    figdir.mkdir(parents=True, exist_ok=True)
    row_order = [family for family in FAMILY_ORDER if family in set(effects["motif_family"])]
    matrix = (
        effects.pivot_table(
            index="motif_family",
            columns="inferred_cell_class",
            values="pcw_adjusted_group_beta_chromvar_z",
            aggfunc="mean",
        )
        .reindex(index=row_order, columns=CLASS_ORDER)
    )
    labels = [MOTIF_FAMILY_LABELS.get(family, family) for family in row_order]
    values = matrix.to_numpy(dtype=float)
    finite = values[np.isfinite(values)]
    vmax = max(0.5, min(4.0, float(np.nanpercentile(np.abs(finite), 95)) if finite.size else 1.0))
    fig, ax = plt.subplots(figsize=(11.2, max(4.8, 0.42 * len(row_order) + 2.3)))
    image = ax.imshow(matrix, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(np.arange(len(CLASS_ORDER)))
    ax.set_xticklabels([CLASS_LABELS.get(name, name) for name in CLASS_ORDER], rotation=45, ha="right")
    ax.set_yticks(np.arange(len(row_order)))
    ax.set_yticklabels(labels)
    ax.set_title("GSE305146 GC/content-matched chromVAR motif deviation")
    cbar = fig.colorbar(image, ax=ax, fraction=0.028, pad=0.02)
    cbar.set_label("PCW-adjusted DS beta, chromVAR z")
    ax.set_xticks(np.arange(-0.5, len(CLASS_ORDER), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(row_order), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)
    fig.tight_layout()
    fig.savefig(figdir / "gse305146_marker_atac_gc_chromvar_beta_heatmap.png", dpi=300)
    fig.savefig(figdir / "gse305146_marker_atac_gc_chromvar_beta_heatmap.pdf")
    plt.close(fig)


def plot_fdr_heatmap(effects: pd.DataFrame, figdir: Path) -> None:
    if effects.empty:
        return
    figdir.mkdir(parents=True, exist_ok=True)
    row_order = [family for family in FAMILY_ORDER if family in set(effects["motif_family"])]
    score = effects.copy()
    score["signed_neglog10_fdr"] = np.sign(score["pcw_adjusted_group_beta_chromvar_z"].astype(float)) * (
        -np.log10(score["fdr_bh_chromvar_z"].astype(float).clip(lower=1e-300))
    ).clip(upper=10.0)
    matrix = (
        score.pivot_table(
            index="motif_family",
            columns="inferred_cell_class",
            values="signed_neglog10_fdr",
            aggfunc="mean",
        )
        .reindex(index=row_order, columns=CLASS_ORDER)
    )
    labels = [MOTIF_FAMILY_LABELS.get(family, family) for family in row_order]
    fig, ax = plt.subplots(figsize=(11.2, max(4.8, 0.42 * len(row_order) + 2.3)))
    image = ax.imshow(matrix, aspect="auto", cmap="RdBu_r", vmin=-10, vmax=10)
    ax.set_xticks(np.arange(len(CLASS_ORDER)))
    ax.set_xticklabels([CLASS_LABELS.get(name, name) for name in CLASS_ORDER], rotation=45, ha="right")
    ax.set_yticks(np.arange(len(row_order)))
    ax.set_yticklabels(labels)
    ax.set_title("GSE305146 GC/content-matched chromVAR FDR")
    cbar = fig.colorbar(image, ax=ax, fraction=0.028, pad=0.02)
    cbar.set_label("signed -log10(FDR)")
    ax.set_xticks(np.arange(-0.5, len(CLASS_ORDER), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(row_order), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)
    fig.tight_layout()
    fig.savefig(figdir / "gse305146_marker_atac_gc_chromvar_fdr_heatmap.png", dpi=300)
    fig.savefig(figdir / "gse305146_marker_atac_gc_chromvar_fdr_heatmap.pdf")
    plt.close(fig)


def write_note(path: Path, args: argparse.Namespace, summary: pd.DataFrame, effects: pd.DataFrame, match_qc: pd.DataFrame) -> None:
    top_pos = effects.sort_values("pcw_adjusted_group_beta_chromvar_z", ascending=False).head(12) if not effects.empty else pd.DataFrame()
    top_neg = effects.sort_values("pcw_adjusted_group_beta_chromvar_z", ascending=True).head(12) if not effects.empty else pd.DataFrame()
    fdr_hits = effects[effects["fdr_lt_0_05"]].copy() if not effects.empty else pd.DataFrame()
    lines = [
        "# GSE305146 GC/Content-Matched chromVAR-Style Motif Deviation",
        "",
        "Question: do HOMER motif-family accessibility shifts survive a chromVAR-style sequence-content matched background?",
        "",
        "## Methods Snapshot",
        "",
        "- Unit of analysis: biological-sample pseudobulk within marker-inferred fetal cortical cell classes.",
        "- Peak universe: recurrent 200 bp centered GSE305146 ATAC peaks; autosomes only unless the script is run with `--include-sex-chromosomes`.",
        "- Motif foregrounds: recurrent peaks containing each HOMER motif family from the prior motif scan.",
        f"- Backgrounds: {args.n_background_iterations} random peak sets per motif family, matched by hg38 GC fraction, CpG fraction, average accessibility, source peak width, and recurrence frequency bins.",
        "- Deviation: foreground raw accessibility deviation from expected recurrent-peak accessibility, centered and scaled by the matched background deviation distribution, giving a chromVAR-style z score.",
        "- Differential test: chromVAR z ~ DS status + post-conception week, fit across donors within each inferred class.",
        "",
        "## Motif-Family Summary",
        "",
    ]
    if summary.empty:
        lines.append("- No motif-family effects were available.")
    else:
        for _, row in summary.iterrows():
            lines.append(
                f"- {row['motif_family_label']}: mean PCW-adjusted beta "
                f"{float(row['mean_pcw_adjusted_beta_chromvar_z']):.3g}; "
                f"{int(row['n_classes_beta_positive'])}/{int(row['n_inferred_classes'])} classes positive; "
                f"FDR-positive {int(row['n_classes_fdr_positive'])}, FDR-negative {int(row['n_classes_fdr_negative'])}; "
                f"foreground peaks {int(row['n_foreground_peaks'])}."
            )
    lines.extend(["", "## FDR-Significant Class/Motif Effects", ""])
    if fdr_hits.empty:
        lines.append("- No class/motif effect passed FDR < 0.05 in this first GC/content-matched chromVAR-style pass.")
    else:
        fdr_hits = fdr_hits.assign(abs_beta=fdr_hits["pcw_adjusted_group_beta_chromvar_z"].abs()).sort_values(
            ["fdr_bh_chromvar_z", "abs_beta"], ascending=[True, False]
        )
        for _, row in fdr_hits.head(20).iterrows():
            lines.append(
                f"- {row['motif_family_label']} in {row['inferred_cell_class']}: "
                f"beta {float(row['pcw_adjusted_group_beta_chromvar_z']):.3g}, "
                f"FDR {float(row['fdr_bh_chromvar_z']):.3g}; "
                f"DS donors {int(row['n_ds_donors'])}, control donors {int(row['n_control_donors'])}."
            )
    lines.extend(["", "## Strongest Positive Effects", ""])
    for _, row in top_pos.iterrows():
        lines.append(
            f"- {row['motif_family_label']} in {row['inferred_cell_class']}: "
            f"beta {float(row['pcw_adjusted_group_beta_chromvar_z']):.3g}, "
            f"FDR {float(row['fdr_bh_chromvar_z']):.3g}."
        )
    lines.extend(["", "## Strongest Negative Effects", ""])
    for _, row in top_neg.iterrows():
        lines.append(
            f"- {row['motif_family_label']} in {row['inferred_cell_class']}: "
            f"beta {float(row['pcw_adjusted_group_beta_chromvar_z']):.3g}, "
            f"FDR {float(row['fdr_bh_chromvar_z']):.3g}."
        )
    if not match_qc.empty:
        lines.extend(["", "## Background Matching QC", ""])
        for _, row in match_qc.iterrows():
            lines.append(
                f"- {row['motif_family_label']}: mean GC diff {float(row['mean_diff_gc_fraction']):.4g}, "
                f"mean CpG diff {float(row['mean_diff_cpg_fraction']):.4g}, "
                f"mean accessibility diff {float(row['mean_diff_mean_log2_peak_cpm']):.4g}."
            )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This is the strongest motif-activity control so far because the background is matched by actual hg38 sequence composition and average accessibility, not only by recurrent-peak structure. It still remains a first-pass analysis because cell classes are marker-inferred and not author-curated clusters.",
            "",
            "Effects that remain after this correction are more compatible with TF/regulatory activity shifts. Effects that shrink or fail FDR should be described as broader accessibility or developmental peak-context signals rather than direct motif activity.",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def run_plot_only(args: argparse.Namespace) -> int:
    outdir = Path(args.outdir)
    figdir = Path(args.figdir)
    effects = pd.read_csv(outdir / "gse305146_marker_atac_gc_chromvar_effects.tsv", sep="\t")
    summary = pd.read_csv(outdir / "gse305146_marker_atac_gc_chromvar_summary.tsv", sep="\t")
    match_qc = pd.read_csv(outdir / "gse305146_marker_atac_gc_chromvar_background_qc.tsv", sep="\t")
    plot_effect_heatmap(effects, figdir)
    plot_fdr_heatmap(effects, figdir)
    write_note(Path(args.note), args, summary, effects, match_qc)
    print("Regenerated GC/content-matched chromVAR figures and note")
    return 0


def main() -> int:
    args = parse_args()
    if args.plot_only:
        return run_plot_only(args)

    outdir = Path(args.outdir)
    figdir = Path(args.figdir)
    outdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    families = [family for family in FAMILY_ORDER if family != "BACKGROUND"]
    if args.limit_families:
        families = families[: args.limit_families]

    regions = read_background_regions(args).reset_index(drop=True)
    composition = read_homer_composition(Path(args.composition_table))
    regions = regions.merge(composition, on="peak_id", how="left")
    regions = compute_gc_for_missing(regions, Path(args.genome_fasta))

    donor_store, mapping_qc, library_class_qc = collect_donor_pseudobulk(args, regions)
    records, counts, depths = build_count_matrix(donor_store)
    regions = add_bias_columns(regions, donor_store)
    regions.to_csv(
        outdir / "gse305146_marker_atac_gc_chromvar_region_bias.tsv.gz",
        sep="\t",
        index=False,
        lineterminator="\n",
        compression="gzip",
    )
    mapping_qc.to_csv(
        outdir / "gse305146_marker_atac_gc_chromvar_mapping_qc.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    library_class_qc.to_csv(
        outdir / "gse305146_marker_atac_gc_chromvar_library_class_qc.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )

    foregrounds = read_motif_foregrounds(Path(args.motif_hits), regions, families)
    donor, match_qc = calculate_chromvar_rows(records, counts, depths, regions, foregrounds, families, args)
    donor_fields = [
        "library_prefix",
        "geo_sample_accession",
        "biological_sample_id",
        "sample_name",
        "group",
        "genotype",
        "dev_stage_pcw",
        "sex",
        "platform_refs",
        "inferred_cell_class",
        "n_cells",
        "total_peak_counts",
        "motif_family",
        "motif_family_label",
        "n_foreground_peaks",
        "foreground_observed_counts",
        "foreground_expected_counts",
        "foreground_raw_deviation",
        "background_mean_raw_deviation",
        "background_sd_raw_deviation",
        "chromvar_deviation",
        "chromvar_deviation_z",
    ]
    write_tsv(
        outdir / "gse305146_marker_atac_gc_chromvar_donor.tsv",
        donor.to_dict(orient="records"),
        donor_fields,
    )
    match_qc.to_csv(
        outdir / "gse305146_marker_atac_gc_chromvar_background_qc.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    effects = effect_rows(donor, args)
    effects.to_csv(
        outdir / "gse305146_marker_atac_gc_chromvar_effects.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    summary = family_summary(effects, match_qc)
    summary.to_csv(
        outdir / "gse305146_marker_atac_gc_chromvar_summary.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    plot_effect_heatmap(effects, figdir)
    plot_fdr_heatmap(effects, figdir)
    write_note(Path(args.note), args, summary, effects, match_qc)

    print("Wrote GC/content-matched chromVAR-style outputs")
    print(outdir / "gse305146_marker_atac_gc_chromvar_effects.tsv")
    print(outdir / "gse305146_marker_atac_gc_chromvar_summary.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
