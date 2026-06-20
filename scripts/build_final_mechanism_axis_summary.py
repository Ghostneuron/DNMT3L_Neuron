#!/usr/bin/env python3
"""Build a manuscript-facing mechanism-axis summary for final candidates."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap


FINAL_CANDIDATES = Path("results/manuscript_summary/final_candidate_manuscript_summary.tsv")
CANDIDATE_LOCUS = Path("results/mechanism_question/candidate_locus_mechanism_integration.tsv")
MOTIF_INSTANCES = Path(
    "results/mechanism_question/candidate_locus_homer_instance_scan/candidate_locus_homer_motif_instance_gene_summary.tsv"
)
GSE305146_SUPP = Path("data/geo/GSE305146/metadata/supplementary_files.tsv")
GSE305146_LOCAL = Path("data/geo/GSE305146/supplementary")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="results/mechanism_question")
    parser.add_argument("--figdir", default="figures/mechanism_question")
    return parser.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t")


def fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        return f"{value:.6g}"
    if pd.isna(value):
        return ""
    return str(value)


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fieldnames})


def to_int(value: object) -> int:
    if value is None or pd.isna(value) or value == "":
        return 0
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def yesno(flag: bool) -> str:
    return "yes" if flag else "no"


def direction_class(text: object) -> str:
    text = str(text).lower()
    if "opposite_possible_repression" in text or "opposite_possible_derepression" in text:
        return "repression_or_derepression_like"
    if "same_direction" in text:
        return "directionally_complex_same_direction"
    return "unresolved_or_mixed"


def mechanism_axis_call(row: dict[str, object]) -> str:
    human = row["human_methylation_axis"] == "yes"
    mouse = row["mouse_dnmt3l_methylation_axis"] == "yes"
    dnmt3a = row["mouse_dnmt3a_reference_axis"] == "yes"
    dev = row["developmental_tf_instance_axis"] == "yes"
    cpg = row["cpg_ctcf_instance_axis"] == "yes"

    if human and mouse and dnmt3a and dev and cpg:
        return "cross_species_dnmt3a_linked_methylation_plus_developmental_cpg_tf_context"
    if human and mouse and dnmt3a and dev:
        return "cross_species_dnmt3a_linked_methylation_plus_developmental_tf_context"
    if human and mouse and dnmt3a and cpg:
        return "cross_species_dnmt3a_linked_methylation_plus_cpg_ctcf_context"
    if human and mouse and dnmt3a:
        return "cross_species_dnmt3a_linked_methylation_without_selected_tf_instance"
    if mouse and dnmt3a and (dev or cpg):
        return "mouse_dnmt3a_linked_methylation_with_tf_context_human_unresolved"
    if mouse and dnmt3a:
        return "mouse_dnmt3a_linked_methylation_human_unresolved"
    if human or mouse:
        return "methylation_associated_dnmt3a_dependence_unresolved"
    if dev or cpg:
        return "tf_context_without_current_methylation_anchor"
    return "mechanism_unresolved"


def concise_interpretation(row: dict[str, object]) -> str:
    call = row["final_mechanism_axis_call"]
    stat = row["stat_axis_interpretation"]
    if call.startswith("cross_species_dnmt3a_linked") and "developmental_cpg" in call:
        base = "strongest combined support: DNMT3A-linked methylation plus local developmental and CpG/CTCF motif context"
    elif call.startswith("cross_species_dnmt3a_linked") and "developmental" in call:
        base = "DNMT3A-linked methylation with local developmental TF sequence context"
    elif call.startswith("cross_species_dnmt3a_linked") and "cpg" in call:
        base = "DNMT3A-linked methylation with local CpG/CTCF sequence context"
    elif call.startswith("cross_species_dnmt3a_linked"):
        base = "DNMT3A-linked methylation support, but selected local motif instances are weak or absent"
    elif call.startswith("mouse_dnmt3a_linked"):
        base = "mouse DNMT3A-linked methylation support; human methylation support remains unresolved"
    else:
        base = "mechanism remains unresolved with current public data"
    if stat == "stat_sequence_only_not_enriched":
        return f"{base}; STAT motif presence is sequence-only and not enriched"
    return base


def audit_gse305146_inputs() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    expected = read_tsv(GSE305146_SUPP) if GSE305146_SUPP.exists() else pd.DataFrame(columns=["url"])
    expected_names = expected["url"].astype(str).map(lambda url: url.rsplit("/", 1)[-1]).tolist() if "url" in expected else []
    local_names = {path.name for path in GSE305146_LOCAL.glob("*")} if GSE305146_LOCAL.exists() else set()
    for suffix in ["barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"]:
        expected_suffix = sorted(name for name in expected_names if name.endswith(suffix))
        local_suffix = sorted(name for name in local_names if name.endswith(suffix))
        rows.append(
            {
                "input_type": suffix,
                "expected_from_geo": len(expected_suffix),
                "present_locally": len(local_suffix),
                "missing": max(0, len(expected_suffix) - len(local_suffix)),
                "status": "ready" if len(expected_suffix) and len(expected_suffix) == len(local_suffix) else "incomplete",
            }
        )
    return rows


def make_heatmap(rows: list[dict[str, object]], figdir: Path) -> None:
    figdir.mkdir(parents=True, exist_ok=True)
    columns = [
        ("Human\nmethyl.", "human_methylation_axis"),
        ("Mouse\nmethyl.", "mouse_dnmt3l_methylation_axis"),
        ("DNMT3A\nref.", "mouse_dnmt3a_reference_axis"),
        ("Dev. TF\ninstance", "developmental_tf_instance_axis"),
        ("CpG/CTCF\ninstance", "cpg_ctcf_instance_axis"),
        ("DS ATAC\nlocal", "fetal_ds_atac_local_axis"),
        ("STAT ctrl\nonly", "stat_sequence_control_axis"),
        ("Repress.\nlike", "repression_like_axis"),
        ("Human\nunresolved", "human_methylation_unresolved_axis"),
    ]
    matrix = np.array([[1 if str(row[column]) == "yes" else 0 for _, column in columns] for row in rows], dtype=float)
    fig_height = max(6.0, 0.34 * len(rows) + 1.8)
    fig, ax = plt.subplots(figsize=(9.7, fig_height))
    cmap = ListedColormap(["#f4f4f2", "#305f72"])
    im = ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(range(len(columns)))
    ax.set_xticklabels([label for label, _ in columns], rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([str(row["gene_symbol"]) for row in rows], fontsize=9)
    ax.set_title("Final candidate mechanism axes", fontsize=12, pad=12)
    ax.set_xticks(np.arange(-0.5, len(columns), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(rows), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.0)
    ax.tick_params(which="minor", bottom=False, left=False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02, ticks=[0.25, 0.75])
    cbar.ax.set_yticklabels(["absent", "present"])
    fig.tight_layout()
    fig.savefig(figdir / "final_candidate_mechanism_axis_heatmap.png", dpi=300)
    fig.savefig(figdir / "final_candidate_mechanism_axis_heatmap.pdf")
    plt.close(fig)


def write_note(path: Path, rows: list[dict[str, object]], audit_rows: list[dict[str, object]]) -> None:
    df = pd.DataFrame(rows)
    counts = df["final_mechanism_axis_call"].value_counts()
    n = len(rows)
    human_mouse = int(((df["human_methylation_axis"] == "yes") & (df["mouse_dnmt3l_methylation_axis"] == "yes")).sum())
    dev = int((df["developmental_tf_instance_axis"] == "yes").sum())
    cpg = int((df["cpg_ctcf_instance_axis"] == "yes").sum())
    stat = int((df["stat_sequence_control_axis"] == "yes").sum())
    matrix_audit = {row["input_type"]: row for row in audit_rows}
    lines = [
        "# Final candidate mechanism-axis summary",
        "",
        "This table merges candidate methylation evidence, mouse DNMT3A-reference overlap, local fetal DS ATAC context, and candidate-linked HOMER motif-instance evidence.",
        "",
        "## Main outputs",
        "",
        "- Candidate summary: `results/mechanism_question/final_candidate_mechanism_axis_summary.tsv`",
        "- Class counts: `results/mechanism_question/final_candidate_mechanism_axis_class_counts.tsv`",
        "- GSE305146 input audit: `results/mechanism_question/gse305146_cell_resolved_input_audit.tsv`",
        "- Heatmap: `figures/mechanism_question/final_candidate_mechanism_axis_heatmap.png`",
        "",
        "## High-level result",
        "",
        f"- Human plus mouse methylation support: {human_mouse}/{n} candidates.",
        f"- Developmental TF motif instances: {dev}/{n} candidates.",
        f"- CpG/CTCF/GC motif instances: {cpg}/{n} candidates.",
        f"- STAT1/STAT3 control motif instances: {stat}/{n} candidates, but these are non-enriched sequence motifs.",
        "",
        "## Mechanism-axis class counts",
        "",
    ]
    for call, count in counts.items():
        genes = ";".join(df.loc[df["final_mechanism_axis_call"] == call, "gene_symbol"].astype(str).tolist())
        lines.append(f"- {call}: {count} ({genes})")
    matrices_ready = (
        str(matrix_audit.get("matrix.mtx.gz", {}).get("present_locally", "0"))
        == str(matrix_audit.get("matrix.mtx.gz", {}).get("expected_from_geo", ""))
        and str(matrix_audit.get("matrix.mtx.gz", {}).get("expected_from_geo", "0")) != "0"
    )
    if matrices_ready:
        audit_interpretation = (
            "Interpretation: GSE305146 matrix recovery is now complete. True cell-resolved chromVAR/"
            "differential-accessibility analysis is no longer blocked by Matrix Market files, but still requires "
            "a cell/cell-type label strategy, such as barcode annotation recovery, de novo clustering, or label transfer."
        )
    else:
        audit_interpretation = (
            "Interpretation: true GSE305146 cell-resolved chromVAR/differential-accessibility analysis still "
            "requires missing matrix files. The current final mechanism summary therefore uses pseudobulk ATAC, "
            "recurrent peak sets, HOMER motif enrichment, and candidate-linked motif-instance scans."
        )

    lines.extend(
        [
            "",
            "## GSE305146 Cell-Resolved Input Audit",
            "",
            f"- Barcodes present: {matrix_audit.get('barcodes.tsv.gz', {}).get('present_locally', 0)} of {matrix_audit.get('barcodes.tsv.gz', {}).get('expected_from_geo', 0)}.",
            f"- Features present: {matrix_audit.get('features.tsv.gz', {}).get('present_locally', 0)} of {matrix_audit.get('features.tsv.gz', {}).get('expected_from_geo', 0)}.",
            f"- Matrices present: {matrix_audit.get('matrix.mtx.gz', {}).get('present_locally', 0)} of {matrix_audit.get('matrix.mtx.gz', {}).get('expected_from_geo', 0)}.",
            "",
            audit_interpretation,
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines))


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    figdir = Path(args.figdir)
    outdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    final = read_tsv(FINAL_CANDIDATES)
    locus = read_tsv(CANDIDATE_LOCUS).set_index("gene_symbol")
    motif = read_tsv(MOTIF_INSTANCES).set_index("gene_symbol")
    rows: list[dict[str, object]] = []

    for base in final.to_dict(orient="records"):
        gene = str(base["gene_symbol"])
        loc = locus.loc[gene].to_dict() if gene in locus.index else {}
        mot = motif.loc[gene].to_dict() if gene in motif.index else {}

        human_methyl = (
            to_int(loc.get("human_dnmt3l_oe_effective_probe_or_region_count", 0))
            + to_int(loc.get("human_fetal_ds_effective_probe_count", 0))
            + to_int(loc.get("human_adult_ds_probe_count", 0))
            > 0
        )
        mouse_methyl = to_int(loc.get("mouse_dnmt3l_dmr_count", 0)) > 0
        dnmt3a_ref = to_int(loc.get("mouse_dnmt3a_p39_cko_hypo_overlap_dmr_count", 0)) > 0
        fetal_atac = (
            to_int(loc.get("gse305146_ds_open_peak_count_100kb", 0))
            + to_int(loc.get("gse305146_control_open_peak_count_100kb", 0))
            > 0
        )
        dev_instance = to_int(mot.get("n_regions_with_developmental_tf_motif", 0)) > 0
        cpg_instance = to_int(mot.get("n_regions_with_cpg_ctcf_gc_motif", 0)) > 0
        stat_instance = to_int(mot.get("n_regions_with_stat_motif", 0)) > 0
        direction = direction_class(base.get("directionality", ""))

        row = {
            "gene_symbol": gene,
            "tier": base.get("tier", ""),
            "candidate_group": base.get("candidate_group", ""),
            "paper_role": base.get("paper_role", ""),
            "base_mechanism_class": base.get("mechanism_class", ""),
            "human_methylation_axis": yesno(human_methyl),
            "human_dnmt3l_oe_support_count": loc.get("human_dnmt3l_oe_effective_probe_or_region_count", 0),
            "human_fetal_ds_methylation_count": loc.get("human_fetal_ds_effective_probe_count", 0),
            "human_adult_ds_methylation_count": loc.get("human_adult_ds_probe_count", 0),
            "mouse_dnmt3l_methylation_axis": yesno(mouse_methyl),
            "mouse_dnmt3l_dmr_count": loc.get("mouse_dnmt3l_dmr_count", 0),
            "mouse_dnmt3a_reference_axis": yesno(dnmt3a_ref),
            "mouse_dnmt3a_reference_dmr_count": loc.get("mouse_dnmt3a_p39_cko_hypo_overlap_dmr_count", 0),
            "mouse_h3k27me3_overlap_count": loc.get("mouse_h3k27me3_p39_cko_up_overlap_dmr_count", 0),
            "fetal_ds_atac_local_axis": yesno(fetal_atac),
            "gse305146_ds_open_peak_count_100kb": loc.get("gse305146_ds_open_peak_count_100kb", 0),
            "gse305146_control_open_peak_count_100kb": loc.get("gse305146_control_open_peak_count_100kb", 0),
            "developmental_tf_instance_axis": yesno(dev_instance),
            "developmental_tf_instance_region_count": mot.get("n_regions_with_developmental_tf_motif", 0),
            "developmental_tf_instance_families": mot.get("developmental_tf_motif_families", ""),
            "cpg_ctcf_instance_axis": yesno(cpg_instance),
            "cpg_ctcf_instance_region_count": mot.get("n_regions_with_cpg_ctcf_gc_motif", 0),
            "cpg_ctcf_instance_families": mot.get("cpg_ctcf_gc_motif_families", ""),
            "stat_sequence_control_axis": yesno(stat_instance),
            "stat_control_instance_region_count": mot.get("n_regions_with_stat_motif", 0),
            "stat_axis_interpretation": "stat_sequence_only_not_enriched" if stat_instance else "not_supported_by_motif_instance_scan",
            "repression_like_axis": yesno(direction == "repression_or_derepression_like"),
            "directionality_axis": direction,
            "human_methylation_unresolved_axis": yesno(mouse_methyl and not human_methyl),
            "candidate_locus_mechanism_call": loc.get("candidate_locus_mechanism_call", ""),
            "short_interpretation": base.get("short_interpretation", ""),
        }
        row["final_mechanism_axis_call"] = mechanism_axis_call(row)
        row["final_mechanism_axis_interpretation"] = concise_interpretation(row)
        rows.append(row)

    fields = [
        "gene_symbol",
        "tier",
        "candidate_group",
        "paper_role",
        "base_mechanism_class",
        "final_mechanism_axis_call",
        "final_mechanism_axis_interpretation",
        "human_methylation_axis",
        "human_dnmt3l_oe_support_count",
        "human_fetal_ds_methylation_count",
        "human_adult_ds_methylation_count",
        "mouse_dnmt3l_methylation_axis",
        "mouse_dnmt3l_dmr_count",
        "mouse_dnmt3a_reference_axis",
        "mouse_dnmt3a_reference_dmr_count",
        "mouse_h3k27me3_overlap_count",
        "fetal_ds_atac_local_axis",
        "gse305146_ds_open_peak_count_100kb",
        "gse305146_control_open_peak_count_100kb",
        "developmental_tf_instance_axis",
        "developmental_tf_instance_region_count",
        "developmental_tf_instance_families",
        "cpg_ctcf_instance_axis",
        "cpg_ctcf_instance_region_count",
        "cpg_ctcf_instance_families",
        "stat_sequence_control_axis",
        "stat_control_instance_region_count",
        "stat_axis_interpretation",
        "repression_like_axis",
        "directionality_axis",
        "human_methylation_unresolved_axis",
        "candidate_locus_mechanism_call",
        "short_interpretation",
    ]
    write_tsv(outdir / "final_candidate_mechanism_axis_summary.tsv", rows, fields)

    df = pd.DataFrame(rows)
    class_rows = []
    for call, group in df.groupby("final_mechanism_axis_call", sort=False):
        class_rows.append(
            {
                "final_mechanism_axis_call": call,
                "n_candidates": int(len(group)),
                "genes": ";".join(group["gene_symbol"].astype(str).tolist()),
                "n_human_methylation": int((group["human_methylation_axis"] == "yes").sum()),
                "n_mouse_methylation": int((group["mouse_dnmt3l_methylation_axis"] == "yes").sum()),
                "n_dnmt3a_reference": int((group["mouse_dnmt3a_reference_axis"] == "yes").sum()),
                "n_developmental_tf_instance": int((group["developmental_tf_instance_axis"] == "yes").sum()),
                "n_cpg_ctcf_instance": int((group["cpg_ctcf_instance_axis"] == "yes").sum()),
                "n_stat_sequence_control": int((group["stat_sequence_control_axis"] == "yes").sum()),
            }
        )
    write_tsv(
        outdir / "final_candidate_mechanism_axis_class_counts.tsv",
        class_rows,
        [
            "final_mechanism_axis_call",
            "n_candidates",
            "genes",
            "n_human_methylation",
            "n_mouse_methylation",
            "n_dnmt3a_reference",
            "n_developmental_tf_instance",
            "n_cpg_ctcf_instance",
            "n_stat_sequence_control",
        ],
    )

    audit_rows = audit_gse305146_inputs()
    write_tsv(
        outdir / "gse305146_cell_resolved_input_audit.tsv",
        audit_rows,
        ["input_type", "expected_from_geo", "present_locally", "missing", "status"],
    )
    make_heatmap(rows, figdir)
    write_note(outdir / "final_candidate_mechanism_axis_summary.md", rows, audit_rows)
    print(f"Wrote {len(rows)} candidate mechanism-axis rows")
    print(outdir / "final_candidate_mechanism_axis_summary.tsv")


if __name__ == "__main__":
    main()
