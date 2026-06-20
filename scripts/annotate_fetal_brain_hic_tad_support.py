#!/usr/bin/env python3
"""Annotate candidate CpG-gene links with fetal brain Hi-C TAD support.

This is intentionally conservative. Won et al. GSE77565 provides processed
CP/GZ TADs in hg19. The analysis asks whether DNMT3L-sensitive CpGs and the
candidate gene TSS fall in the same fetal cortex TAD. Same-TAD support is not
treated as a direct enhancer-promoter loop.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from bisect import bisect_right
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


DEFAULT_CANDIDATES = Path(
    "manuscript/Supplementary Tables/"
    "ST6b_final_candidate_evidence_matrix_with_regulatory_context.tsv"
)
DEFAULT_LINKS = Path("results/human_methylation_convergence/human_methylation_probe_gene_links.tsv")
DEFAULT_DMC = Path("results/human_candidate_dmrs/GSE157078_candidate_dmcs.tsv")
DEFAULT_GTF = Path("data/reference/gencode/Homo_sapiens.GRCh37.75.gtf.gz")
DEFAULT_CP_TAD = Path("data/external/won2016_developing_human_brain_hic/GSE77565_CP_TAD.bed.gz")
DEFAULT_GZ_TAD = Path("data/external/won2016_developing_human_brain_hic/GSE77565_GZ_TAD.bed.gz")
DEFAULT_OUTDIR = Path("results/fetal_brain_3d_contact_support")
DEFAULT_ST = Path(
    "manuscript/Supplementary Tables/"
    "ST6d_final_candidate_fetal_brain_hic_tad_support.tsv"
)
DEFAULT_FIGDIR = Path("figures/fetal_brain_3d_contact_support")
SYMBOL_ALIASES_GRCH37 = {
    "PLPP3": {"PLPP3", "PPAP2B"},
}


@dataclass(frozen=True)
class Gene:
    symbol: str
    chrom: str
    start: int
    end: int
    strand: str

    @property
    def tss(self) -> int:
        if self.strand == "-":
            return self.end
        return self.start


@dataclass(frozen=True)
class Tad:
    chrom: str
    start: int
    end: int
    tad_id: str


def parse_attrs(attr: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for match in re.finditer(r'(\S+) "([^"]+)"', attr):
        out[match.group(1)] = match.group(2)
    return out


def norm_chrom(chrom: str) -> str:
    chrom = str(chrom)
    if chrom.startswith("chr"):
        return chrom
    return f"chr{chrom}"


def read_genes(gtf: Path, target_symbols: set[str]) -> dict[str, Gene]:
    genes: dict[str, Gene] = {}
    lookup: dict[str, str] = {}
    for symbol in target_symbols:
        for alias in SYMBOL_ALIASES_GRCH37.get(symbol, {symbol}):
            lookup[alias] = symbol
    with gzip.open(gtf, "rt") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = parse_attrs(fields[8])
            gtf_symbol = attrs.get("gene_name", "")
            symbol = lookup.get(gtf_symbol)
            if symbol is None:
                continue
            chrom = norm_chrom(fields[0])
            # GTF is 1-based inclusive; use 0-based half-open for BED-like logic.
            start = int(fields[3]) - 1
            end = int(fields[4])
            candidate = Gene(symbol=symbol, chrom=chrom, start=start, end=end, strand=fields[6])
            current = genes.get(symbol)
            if current is None or (candidate.end - candidate.start) > (current.end - current.start):
                genes[symbol] = candidate
    missing = sorted(target_symbols - set(genes))
    if missing:
        raise SystemExit(f"Missing candidate genes in GRCh37 GTF: {', '.join(missing)}")
    return genes


def read_tads(path: Path, label: str) -> dict[str, list[Tad]]:
    by_chrom: dict[str, list[Tad]] = defaultdict(list)
    with gzip.open(path, "rt") as handle:
        for i, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, *_ = line.rstrip("\n").split("\t")
            tad = Tad(norm_chrom(chrom), int(start), int(end), f"{label}_TAD_{i:04d}")
            by_chrom[tad.chrom].append(tad)
    for chrom in by_chrom:
        by_chrom[chrom].sort(key=lambda x: x.start)
    return dict(by_chrom)


def tad_lookup(tads_by_chrom: dict[str, list[Tad]], chrom: str, pos: int) -> Tad | None:
    tads = tads_by_chrom.get(norm_chrom(chrom), [])
    if not tads:
        return None
    starts = [t.start for t in tads]
    idx = bisect_right(starts, pos) - 1
    if idx >= 0:
        tad = tads[idx]
        if tad.start <= pos < tad.end:
            return tad
    return None


def read_probe_hg19(dmc: Path) -> dict[str, dict[str, object]]:
    probes: dict[str, dict[str, object]] = {}
    with dmc.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            start = int(row["start"])
            end = int(row["end"])
            probes[row["probe_id"]] = {
                "chrom": norm_chrom(row["chrom"]),
                "start": start,
                "end": end,
                "midpoint": (start + end) // 2,
                "direction": row.get("direction", ""),
                "delta_beta": row.get("delta_beta", ""),
            }
    return probes


def yesno(value: bool) -> str:
    return "yes" if value else "no"


def unique_join(values) -> str:
    cleaned = []
    seen = set()
    for value in values:
        if value is None:
            continue
        value = str(value)
        if value == "" or value in seen:
            continue
        seen.add(value)
        cleaned.append(value)
    return ",".join(cleaned)


def make_outputs(args: argparse.Namespace) -> None:
    candidates = pd.read_csv(args.candidates, sep="\t")
    candidate_symbols = list(candidates["gene_symbol"])
    target_symbols = set(candidate_symbols)
    genes = read_genes(args.gtf, target_symbols)
    probe_hg19 = read_probe_hg19(args.dmc)
    cp_tads = read_tads(args.cp_tad, "CP")
    gz_tads = read_tads(args.gz_tad, "GZ")

    links = pd.read_csv(args.links, sep="\t", low_memory=False)
    links = links[links["gene_symbol"].isin(target_symbols)].copy()
    links = links[links["probe_id"].isin(probe_hg19)].copy()

    detailed_rows: list[dict[str, object]] = []
    for _, row in links.iterrows():
        symbol = row["gene_symbol"]
        gene = genes[symbol]
        probe = probe_hg19[row["probe_id"]]
        if probe["chrom"] != gene.chrom:
            same_cp = same_gz = False
            cp_probe_tad = cp_gene_tad = gz_probe_tad = gz_gene_tad = None
        else:
            cp_probe_tad = tad_lookup(cp_tads, gene.chrom, int(probe["midpoint"]))
            cp_gene_tad = tad_lookup(cp_tads, gene.chrom, gene.tss)
            gz_probe_tad = tad_lookup(gz_tads, gene.chrom, int(probe["midpoint"]))
            gz_gene_tad = tad_lookup(gz_tads, gene.chrom, gene.tss)
            same_cp = cp_probe_tad is not None and cp_gene_tad is not None and cp_probe_tad.tad_id == cp_gene_tad.tad_id
            same_gz = gz_probe_tad is not None and gz_gene_tad is not None and gz_probe_tad.tad_id == gz_gene_tad.tad_id

        direct = str(row.get("relationship", "")) in {"promoter_2kb", "gene_body"}
        contexts = []
        if same_cp:
            contexts.append("CP")
        if same_gz:
            contexts.append("GZ")
        distance = int(probe["midpoint"]) - gene.tss if probe["chrom"] == gene.chrom else ""
        detailed_rows.append(
            {
                "gene_symbol": symbol,
                "probe_id": row["probe_id"],
                "probe_chrom_hg19": probe["chrom"],
                "probe_start_hg19": probe["start"],
                "probe_end_hg19": probe["end"],
                "probe_direction": probe["direction"],
                "gse157078_delta_beta": probe["delta_beta"],
                "gene_chrom_hg19": gene.chrom,
                "gene_start_hg19": gene.start,
                "gene_end_hg19": gene.end,
                "gene_strand": gene.strand,
                "gene_tss_hg19": gene.tss,
                "signed_distance_to_tss_hg19": distance,
                "abs_distance_to_tss_hg19": abs(distance) if distance != "" else "",
                "original_hg38_link_relationship": row.get("relationship", ""),
                "direct_promoter_or_gene_body_link": yesno(direct),
                "cp_probe_tad": cp_probe_tad.tad_id if cp_probe_tad else "",
                "cp_gene_tad": cp_gene_tad.tad_id if cp_gene_tad else "",
                "cp_same_tad": yesno(same_cp),
                "gz_probe_tad": gz_probe_tad.tad_id if gz_probe_tad else "",
                "gz_gene_tad": gz_gene_tad.tad_id if gz_gene_tad else "",
                "gz_same_tad": yesno(same_gz),
                "same_fetal_brain_tad_contexts": ",".join(contexts),
                "same_fetal_brain_tad_any": yesno(bool(contexts)),
            }
        )

    detailed = pd.DataFrame(detailed_rows)
    if detailed.empty:
        raise SystemExit("No candidate probe-gene links were available for fetal Hi-C TAD annotation.")

    rows = []
    for _, cand in candidates.iterrows():
        symbol = cand["gene_symbol"]
        gene_rows = detailed[detailed["gene_symbol"].eq(symbol)]
        strict_ids = {
            x.strip()
            for x in str(cand.get("human_strict_fetal_positive_probe_ids", "")).split(",")
            if x.strip() and x.strip().lower() != "nan"
        }
        strict_rows = gene_rows[gene_rows["probe_id"].isin(strict_ids)] if strict_ids else gene_rows.iloc[0:0]
        same_any = gene_rows["same_fetal_brain_tad_any"].eq("yes")
        same_cp = gene_rows["cp_same_tad"].eq("yes")
        same_gz = gene_rows["gz_same_tad"].eq("yes")
        direct_any = gene_rows["direct_promoter_or_gene_body_link"].eq("yes")
        strict_same = strict_rows["same_fetal_brain_tad_any"].eq("yes") if not strict_rows.empty else pd.Series(dtype=bool)
        rows.append(
            {
                "tier": cand["tier"],
                "gene_symbol": symbol,
                "evidence_class": cand["evidence_class"],
                "n_human_dnmt3l_sensitive_probe_gene_links": int(gene_rows["probe_id"].nunique()),
                "n_direct_promoter_or_gene_body_links": int(gene_rows.loc[direct_any, "probe_id"].nunique()),
                "n_same_fetal_brain_tad_links_any": int(gene_rows.loc[same_any, "probe_id"].nunique()),
                "n_same_cp_tad_links": int(gene_rows.loc[same_cp, "probe_id"].nunique()),
                "n_same_gz_tad_links": int(gene_rows.loc[same_gz, "probe_id"].nunique()),
                "n_strict_fetal_positive_probe_links": len(strict_ids),
                "n_strict_same_fetal_brain_tad_links_any": int(strict_rows.loc[strict_same, "probe_id"].nunique()) if not strict_rows.empty else 0,
                "fraction_same_fetal_brain_tad_links": (
                    round(float(gene_rows.loc[same_any, "probe_id"].nunique()) / float(gene_rows["probe_id"].nunique()), 4)
                    if gene_rows["probe_id"].nunique()
                    else ""
                ),
                "same_tad_contexts": unique_join(gene_rows.loc[same_any, "same_fetal_brain_tad_contexts"]),
                "same_tad_probe_ids": unique_join(gene_rows.loc[same_any, "probe_id"]),
                "strict_same_tad_probe_ids": unique_join(strict_rows.loc[strict_same, "probe_id"]) if not strict_rows.empty else "",
                "interpretation": (
                    "same_fetal_brain_TAD_support"
                    if gene_rows.loc[same_any, "probe_id"].nunique()
                    else ("direct_or_local_link_only" if gene_rows.loc[direct_any, "probe_id"].nunique() else "no_human_methylation_link")
                ),
            }
        )
    summary = pd.DataFrame(rows)

    args.outdir.mkdir(parents=True, exist_ok=True)
    args.figdir.mkdir(parents=True, exist_ok=True)
    args.supplementary_table.parent.mkdir(parents=True, exist_ok=True)
    detailed_path = args.outdir / "final_candidate_human_cpg_gene_fetal_brain_hic_tad_links.tsv"
    summary_path = args.outdir / "final_candidate_human_fetal_brain_hic_tad_support_summary.tsv"
    detailed.to_csv(detailed_path, sep="\t", index=False)
    summary.to_csv(summary_path, sep="\t", index=False)
    summary.to_csv(args.supplementary_table, sep="\t", index=False)

    # Preserve the final-candidate order used in the other manuscript panels.
    ordered = summary.copy()
    plot_df = ordered.set_index("gene_symbol")[
        [
            "n_direct_promoter_or_gene_body_links",
            "n_same_cp_tad_links",
            "n_same_gz_tad_links",
            "n_strict_same_fetal_brain_tad_links_any",
        ]
    ].astype(float)
    plt.figure(figsize=(7.8, max(5.2, 0.34 * len(plot_df))))
    ax = sns.heatmap(
        plot_df,
        cmap="YlGnBu",
        linewidths=0.4,
        linecolor="white",
        cbar_kws={"label": "Probe-gene links"},
        annot=True,
        fmt=".0f",
        annot_kws={"fontsize": 8, "color": "black"},
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Fetal Human Brain Hi-C TAD Support For Candidate CpG-Gene Links", fontsize=11)
    ax.set_xticklabels(
        [
            "Direct promoter/\ngene body",
            "Same CP\nTAD",
            "Same GZ\nTAD",
            "Strict fetal\nsame TAD",
        ],
        rotation=0,
        fontsize=8,
    )
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=8, fontstyle="italic")
    plt.tight_layout()
    fig_png = args.figdir / "final_candidate_fetal_brain_hic_tad_support_heatmap.png"
    fig_pdf = args.figdir / "final_candidate_fetal_brain_hic_tad_support_heatmap.pdf"
    plt.savefig(fig_png, dpi=300)
    plt.savefig(fig_pdf)
    plt.close()

    plt.figure(figsize=(7.2, max(5.2, 0.34 * len(plot_df))))
    ax = sns.heatmap(
        plot_df,
        cmap="YlGnBu",
        linewidths=0.4,
        linecolor="white",
        cbar_kws={"label": "Probe-gene links"},
        annot=True,
        fmt=".0f",
        annot_kws={"fontsize": 8, "color": "black"},
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Fetal CP/GZ Hi-C Same-TAD Support", fontsize=11)
    ax.set_xticklabels(
        [
            "Direct promoter/\ngene body",
            "Same CP\nTAD",
            "Same GZ\nTAD",
            "Strict fetal\nsame TAD",
        ],
        rotation=0,
        fontsize=8,
    )
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=8, fontstyle="italic")
    ax.text(-0.12, 1.04, "c", transform=ax.transAxes, fontsize=18, fontweight="bold", va="bottom", ha="left")
    plt.tight_layout()
    ed_png = args.figdir / "ED_Fig6c_final_candidate_fetal_brain_hic_tad_support_heatmap.png"
    ed_pdf = args.figdir / "ED_Fig6c_final_candidate_fetal_brain_hic_tad_support_heatmap.pdf"
    plt.savefig(ed_png, dpi=300)
    plt.savefig(ed_pdf)
    plt.close()

    manifest = args.outdir / "fetal_brain_hic_tad_support_manifest.tsv"
    with manifest.open("w") as handle:
        handle.write("field\tvalue\n")
        handle.write("source_dataset\tGSE77565\n")
        handle.write("source_publication_pmid\t27760116\n")
        handle.write("source_publication_title\tChromosome conformation elucidates regulatory relationships in developing human brain\n")
        handle.write("genome_build\tGRCh37/hg19\n")
        handle.write("contact_feature\tTopologically associating domains from fetal cortical plate and germinal zone Hi-C\n")
        handle.write("interpretation_limit\tSame-TAD support indicates domain-level regulatory plausibility, not a direct enhancer-promoter loop.\n")
        handle.write(f"detailed_table\t{detailed_path}\n")
        handle.write(f"summary_table\t{summary_path}\n")
        handle.write(f"supplementary_table\t{args.supplementary_table}\n")
        handle.write(f"figure_png\t{fig_png}\n")
        handle.write(f"ed_fig6c_panel_png\t{ed_png}\n")

    print(f"Wrote {detailed_path}")
    print(f"Wrote {summary_path}")
    print(f"Wrote {args.supplementary_table}")
    print(f"Wrote {fig_png}")
    print(f"Wrote {ed_png}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", type=Path, default=DEFAULT_CANDIDATES)
    parser.add_argument("--links", type=Path, default=DEFAULT_LINKS)
    parser.add_argument("--dmc", type=Path, default=DEFAULT_DMC)
    parser.add_argument("--gtf", type=Path, default=DEFAULT_GTF)
    parser.add_argument("--cp-tad", type=Path, default=DEFAULT_CP_TAD)
    parser.add_argument("--gz-tad", type=Path, default=DEFAULT_GZ_TAD)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--figdir", type=Path, default=DEFAULT_FIGDIR)
    parser.add_argument("--supplementary-table", type=Path, default=DEFAULT_ST)
    args = parser.parse_args()
    make_outputs(args)


if __name__ == "__main__":
    main()
