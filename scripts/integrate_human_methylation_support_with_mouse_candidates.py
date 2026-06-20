#!/usr/bin/env python3
"""Build human methylation support and compare with mouse candidate genes."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


DEFAULT_LIFTOVER = "results/shsy5y_gse168276/gse157078_candidate_dmcs_hg38_liftover.tsv"
DEFAULT_SHSY = "results/shsy5y_gse168276/gse168276_gse157078_target_phase_effects.tsv"
DEFAULT_DS = "results/ds_methylation_dnmt3l_concordance/ds_dnmt3l_candidate_probe_values.tsv"
DEFAULT_GENCODE = "data/reference/gencode/gencode.v39.annotation.gtf.gz"
DEFAULT_CROSS_SPECIES = "results/cross_species_candidate_expression/cross_species_candidate_priority.tsv"
DEFAULT_OUTDIR = "results/human_methylation_convergence"


@dataclass(frozen=True)
class Probe:
    probe_id: str
    direction: str
    delta_beta: float
    chrom: str
    start: int
    end: int
    row: dict[str, str]


@dataclass(frozen=True)
class Gene:
    gene_id: str
    gene_symbol: str
    gene_type: str
    chrom: str
    start: int
    end: int
    strand: str

    @property
    def tss(self) -> int:
        return self.end if self.strand == "-" else self.start


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--liftover", default=DEFAULT_LIFTOVER)
    parser.add_argument("--shsy-effects", default=DEFAULT_SHSY)
    parser.add_argument("--ds-values", default=DEFAULT_DS)
    parser.add_argument("--gencode-gtf", default=DEFAULT_GENCODE)
    parser.add_argument("--cross-species", default=DEFAULT_CROSS_SPECIES)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--promoter-window", type=int, default=2000)
    parser.add_argument("--nearest-window", type=int, default=100000)
    return parser.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="", errors="replace")
    return path.open(newline="", errors="replace")


def clean_key(value: str) -> str:
    return value.strip().replace("\r", "")


def safe_float(value: object) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return math.nan
    return number if math.isfinite(number) else math.nan


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


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})


def dict_reader(path: Path):
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        reader.fieldnames = [clean_key(field) for field in (reader.fieldnames or [])]
        for row in reader:
            yield {clean_key(key): value for key, value in row.items()}


def phase_key(phase: str) -> str:
    return phase.lower().replace(" ", "_")


def is_yes(value: object) -> bool:
    return str(value).strip().lower() == "yes"


def read_probes(path: Path) -> dict[str, Probe]:
    probes: dict[str, Probe] = {}
    for row in dict_reader(path):
        probe_id = row["probe_id"]
        chrom = row.get("hg38_chrom", "")
        if not probe_id or not chrom:
            continue
        probes[probe_id] = Probe(
            probe_id=probe_id,
            direction=row.get("direction", ""),
            delta_beta=safe_float(row.get("gse157078_delta_beta")),
            chrom=chrom,
            start=int(float(row["hg38_start"])),
            end=int(float(row["hg38_end"])),
            row=row,
        )
    return probes


def read_shsy(path: Path) -> dict[str, dict[str, dict[str, object]]]:
    by_probe: dict[str, dict[str, dict[str, object]]] = defaultdict(dict)
    for row in dict_reader(path):
        probe_id = row["probe_id"]
        key = phase_key(row["phase"])
        by_probe[probe_id][key] = {
            "delta": safe_float(row.get("shsy5y_delta_beta")),
            "signed_delta": safe_float(row.get("signed_shsy5y_delta_by_gse157078_direction")),
            "concordant": is_yes(row.get("concordant")),
            "n_gfp": row.get("n_gfp", ""),
            "n_dnmt3l_gfp": row.get("n_dnmt3l_gfp", ""),
        }
    return by_probe


def context_kind(context: str) -> str:
    text = context.lower()
    if "fetal" in text or "fetus" in text:
        return "fetal"
    if "adult" in text:
        return "adult"
    return "mixed"


def read_ds(path: Path) -> dict[str, dict[str, object]]:
    values: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in dict_reader(path):
        probe_id = row["probe_id"]
        context = row["context"]
        values[probe_id].append(
            {
                "dataset": row["dataset"],
                "context": context,
                "kind": context_kind(context),
                "delta": safe_float(row.get("ds_delta_beta")),
                "signed_delta": safe_float(row.get("signed_ds_delta_by_dnmt3l_direction")),
                "concordant": is_yes(row.get("concordant")),
                "control_n": row.get("control_n_detected", ""),
                "ds_n": row.get("ds_n_detected", ""),
            }
        )

    summary: dict[str, dict[str, object]] = {}
    for probe_id, rows in values.items():
        concordant = [row for row in rows if row["concordant"]]
        fetal = [row for row in rows if row["kind"] == "fetal"]
        fetal_concordant = [row for row in fetal if row["concordant"]]
        adult = [row for row in rows if row["kind"] == "adult"]
        adult_concordant = [row for row in adult if row["concordant"]]
        signed = [safe_float(row["signed_delta"]) for row in rows]
        fetal_signed = [safe_float(row["signed_delta"]) for row in fetal]
        adult_signed = [safe_float(row["signed_delta"]) for row in adult]
        summary[probe_id] = {
            "contexts_tested": len(rows),
            "concordant_contexts": len(concordant),
            "concordant_fraction": len(concordant) / len(rows) if rows else math.nan,
            "fetal_contexts_tested": len(fetal),
            "fetal_concordant_contexts": len(fetal_concordant),
            "fetal_concordant_fraction": len(fetal_concordant) / len(fetal) if fetal else math.nan,
            "adult_contexts_tested": len(adult),
            "adult_concordant_contexts": len(adult_concordant),
            "adult_concordant_fraction": len(adult_concordant) / len(adult) if adult else math.nan,
            "mean_signed_delta": mean(signed),
            "fetal_mean_signed_delta": mean(fetal_signed),
            "adult_mean_signed_delta": mean(adult_signed),
            "max_signed_delta": max(signed) if signed else math.nan,
            "min_signed_delta": min(signed) if signed else math.nan,
            "concordant_context_list": ";".join(f"{row['dataset']}:{row['context']}" for row in concordant),
            "fetal_concordant_context_list": ";".join(
                f"{row['dataset']}:{row['context']}" for row in fetal_concordant
            ),
        }
    return summary


def parse_gtf_attributes(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for match in re.finditer(r'([A-Za-z0-9_]+) "([^"]+)"', text):
        attrs[match.group(1)] = match.group(2)
    return attrs


def strip_version(value: str) -> str:
    return value.split(".", 1)[0]


def read_genes(path: Path) -> dict[str, list[Gene]]:
    genes_by_chrom: dict[str, list[Gene]] = defaultdict(list)
    with open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = parse_gtf_attributes(parts[8])
            gene_id = strip_version(attrs.get("gene_id", ""))
            symbol = attrs.get("gene_name", "")
            gene_type = attrs.get("gene_type", "")
            if not gene_id or not symbol:
                continue
            genes_by_chrom[parts[0]].append(
                Gene(
                    gene_id=gene_id,
                    gene_symbol=symbol.upper(),
                    gene_type=gene_type,
                    chrom=parts[0],
                    start=int(parts[3]) - 1,
                    end=int(parts[4]),
                    strand=parts[6],
                )
            )
    for chrom in genes_by_chrom:
        genes_by_chrom[chrom].sort(key=lambda gene: (gene.start, gene.end, gene.gene_id))
    return genes_by_chrom


def overlap(start_a: int, end_a: int, start_b: int, end_b: int) -> bool:
    return min(end_a, end_b) > max(start_a, start_b)


def signed_distance_to_tss(probe: Probe, gene: Gene) -> int:
    midpoint = (probe.start + probe.end) // 2
    return midpoint - gene.tss


def body_distance(probe: Probe, gene: Gene) -> int:
    if overlap(probe.start, probe.end, gene.start, gene.end):
        return 0
    if probe.end <= gene.start:
        return gene.start - probe.end
    return probe.start - gene.end


def relation_for(probe: Probe, gene: Gene, promoter_window: int) -> dict[str, object]:
    promoter_start = gene.tss - promoter_window
    promoter_end = gene.tss + promoter_window
    promoter = overlap(probe.start, probe.end, promoter_start, promoter_end)
    body = overlap(probe.start, probe.end, gene.start, gene.end)
    if promoter:
        relation = "promoter_2kb"
    elif body:
        relation = "gene_body"
    else:
        relation = "nearest_100kb"
    signed_tss = signed_distance_to_tss(probe, gene)
    return {
        "relationship": relation,
        "promoter_2kb": "yes" if promoter else "no",
        "gene_body": "yes" if body else "no",
        "nearest_100kb": "yes" if relation == "nearest_100kb" else "no",
        "signed_distance_to_tss": signed_tss,
        "abs_distance_to_tss": abs(signed_tss),
        "distance_to_gene_body": body_distance(probe, gene),
    }


def probe_support_rows(
    probes: dict[str, Probe],
    shsy: dict[str, dict[str, dict[str, object]]],
    ds: dict[str, dict[str, object]],
) -> dict[str, dict[str, object]]:
    rows: dict[str, dict[str, object]] = {}
    for probe_id, probe in probes.items():
        shsy_probe = shsy.get(probe_id, {})
        ds_probe = ds.get(probe_id, {})
        shsy_concordant_n = sum(1 for phase in shsy_probe.values() if phase.get("concordant"))
        phase3 = shsy_probe.get("phase_3", {})
        row = {
            "probe_id": probe_id,
            "direction": probe.direction,
            "gse157078_delta_beta": probe.delta_beta,
            "hg38_chrom": probe.chrom,
            "hg38_start": probe.start,
            "hg38_end": probe.end,
            "ucsc_refgene_name": probe.row.get("ucsc_refgene_name", ""),
            "ucsc_refgene_group": probe.row.get("ucsc_refgene_group", ""),
            "gencode_name_array_annotation": probe.row.get("gencode_name", ""),
            "gencode_group_array_annotation": probe.row.get("gencode_group", ""),
            "shsy_growth_delta_beta": shsy_probe.get("growth", {}).get("delta"),
            "shsy_growth_signed_delta": shsy_probe.get("growth", {}).get("signed_delta"),
            "shsy_growth_concordant": "yes" if shsy_probe.get("growth", {}).get("concordant") else "no",
            "shsy_phase1_delta_beta": shsy_probe.get("phase_1", {}).get("delta"),
            "shsy_phase1_signed_delta": shsy_probe.get("phase_1", {}).get("signed_delta"),
            "shsy_phase1_concordant": "yes" if shsy_probe.get("phase_1", {}).get("concordant") else "no",
            "shsy_phase3_delta_beta": phase3.get("delta"),
            "shsy_phase3_signed_delta": phase3.get("signed_delta"),
            "shsy_phase3_concordant": "yes" if phase3.get("concordant") else "no",
            "shsy_concordant_phase_count": shsy_concordant_n,
            "shsy_any_phase_concordant": "yes" if shsy_concordant_n > 0 else "no",
            "ds_contexts_tested": ds_probe.get("contexts_tested", 0),
            "ds_concordant_contexts": ds_probe.get("concordant_contexts", 0),
            "ds_concordant_fraction": ds_probe.get("concordant_fraction"),
            "ds_fetal_contexts_tested": ds_probe.get("fetal_contexts_tested", 0),
            "ds_fetal_concordant_contexts": ds_probe.get("fetal_concordant_contexts", 0),
            "ds_fetal_concordant_fraction": ds_probe.get("fetal_concordant_fraction"),
            "ds_adult_contexts_tested": ds_probe.get("adult_contexts_tested", 0),
            "ds_adult_concordant_contexts": ds_probe.get("adult_concordant_contexts", 0),
            "ds_adult_concordant_fraction": ds_probe.get("adult_concordant_fraction"),
            "ds_mean_signed_delta": ds_probe.get("mean_signed_delta"),
            "ds_fetal_mean_signed_delta": ds_probe.get("fetal_mean_signed_delta"),
            "ds_adult_mean_signed_delta": ds_probe.get("adult_mean_signed_delta"),
            "ds_max_signed_delta": ds_probe.get("max_signed_delta"),
            "ds_min_signed_delta": ds_probe.get("min_signed_delta"),
            "ds_concordant_context_list": ds_probe.get("concordant_context_list", ""),
            "ds_fetal_concordant_context_list": ds_probe.get("fetal_concordant_context_list", ""),
        }
        row["human_triple_phase3_any_ds"] = (
            "yes" if phase3.get("concordant") and int(row["ds_concordant_contexts"]) >= 1 else "no"
        )
        row["human_triple_phase3_fetal_ds"] = (
            "yes" if phase3.get("concordant") and int(row["ds_fetal_concordant_contexts"]) >= 1 else "no"
        )
        row["human_triple_phase3_two_ds_contexts"] = (
            "yes" if phase3.get("concordant") and int(row["ds_concordant_contexts"]) >= 2 else "no"
        )
        row["human_triple_phase3_ds_majority"] = (
            "yes"
            if phase3.get("concordant") and safe_float(row["ds_concordant_fraction"]) >= 0.5
            else "no"
        )
        row["human_triple_phase3_fetal_ds_majority"] = (
            "yes"
            if phase3.get("concordant") and safe_float(row["ds_fetal_concordant_fraction"]) >= 0.5
            else "no"
        )
        row["human_triple_phase3_fetal_mean_positive"] = (
            "yes"
            if phase3.get("concordant") and safe_float(row["ds_fetal_mean_signed_delta"]) > 0
            else "no"
        )
        rows[probe_id] = row
    return rows


def link_probes_to_genes(
    probes: dict[str, Probe],
    probe_rows: dict[str, dict[str, object]],
    genes_by_chrom: dict[str, list[Gene]],
    promoter_window: int,
    nearest_window: int,
) -> list[dict[str, object]]:
    links: dict[tuple[str, str], dict[str, object]] = {}
    for probe in probes.values():
        chrom_genes = genes_by_chrom.get(probe.chrom, [])
        nearest_gene = None
        nearest_distance = math.inf
        for gene in chrom_genes:
            distance = abs(signed_distance_to_tss(probe, gene))
            if distance < nearest_distance:
                nearest_gene = gene
                nearest_distance = distance
            promoter_start = gene.tss - promoter_window
            promoter_end = gene.tss + promoter_window
            if not overlap(probe.start, probe.end, promoter_start, promoter_end) and not overlap(
                probe.start, probe.end, gene.start, gene.end
            ):
                continue
            links[(probe.probe_id, gene.gene_id)] = build_link(probe, gene, probe_rows[probe.probe_id], promoter_window)

        if nearest_gene is not None and nearest_distance <= nearest_window:
            key = (probe.probe_id, nearest_gene.gene_id)
            if key not in links:
                links[key] = build_link(probe, nearest_gene, probe_rows[probe.probe_id], promoter_window)
            else:
                links[key]["nearest_100kb"] = "yes"

    return sorted(
        links.values(),
        key=lambda row: (row["hg38_chrom"], int(row["probe_start"]), int(row["probe_end"]), row["gene_symbol"]),
    )


def build_link(probe: Probe, gene: Gene, support: dict[str, object], promoter_window: int) -> dict[str, object]:
    relation = relation_for(probe, gene, promoter_window)
    return {
        "probe_id": probe.probe_id,
        "hg38_chrom": probe.chrom,
        "probe_start": probe.start,
        "probe_end": probe.end,
        "gene_id": gene.gene_id,
        "gene_symbol": gene.gene_symbol,
        "gene_type": gene.gene_type,
        "gene_chrom": gene.chrom,
        "gene_start": gene.start,
        "gene_end": gene.end,
        "gene_strand": gene.strand,
        "gene_tss": gene.tss,
        **relation,
        **support,
    }


def read_cross_species(path: Path) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    if not path.exists():
        return rows
    for row in dict_reader(path):
        symbol = row.get("human_symbol", "")
        if symbol:
            rows[symbol] = row
    return rows


def aggregate_genes(links: list[dict[str, object]], cross: dict[str, dict[str, str]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in links:
        grouped[(str(row["gene_id"]), str(row["gene_symbol"]))].append(row)

    gene_rows: list[dict[str, object]] = []
    for (gene_id, symbol), rows in grouped.items():
        by_probe = {str(row["probe_id"]): row for row in rows}
        triple_any = [row for row in by_probe.values() if row.get("human_triple_phase3_any_ds") == "yes"]
        triple_fetal = [row for row in by_probe.values() if row.get("human_triple_phase3_fetal_ds") == "yes"]
        triple_fetal_majority = [
            row for row in by_probe.values() if row.get("human_triple_phase3_fetal_ds_majority") == "yes"
        ]
        triple_fetal_mean_positive = [
            row for row in by_probe.values() if row.get("human_triple_phase3_fetal_mean_positive") == "yes"
        ]
        phase3 = [row for row in by_probe.values() if row.get("shsy_phase3_concordant") == "yes"]
        ds_concordant = [row for row in by_probe.values() if int(row.get("ds_concordant_contexts") or 0) > 0]
        promoter_triple = {
            row["probe_id"]
            for row in rows
            if row.get("promoter_2kb") == "yes" and row.get("human_triple_phase3_any_ds") == "yes"
        }
        body_triple = {
            row["probe_id"]
            for row in rows
            if row.get("gene_body") == "yes" and row.get("human_triple_phase3_any_ds") == "yes"
        }
        nearest_triple = {
            row["probe_id"]
            for row in rows
            if row.get("nearest_100kb") == "yes" and row.get("human_triple_phase3_any_ds") == "yes"
        }
        cross_row = cross.get(symbol, {})
        methylation_score = (
            min(10, 3 * len(triple_fetal_mean_positive))
            + min(5, 2 * len(triple_fetal_majority))
            + min(3, len(triple_any))
            + min(3, len(promoter_triple))
            + min(2, len(body_triple))
        )
        cross_score = safe_float(cross_row.get("cross_species_priority_score"))
        combined_score = methylation_score + (cross_score if math.isfinite(cross_score) else 0)
        gene_rows.append(
            {
                "gene_id": gene_id,
                "gene_symbol": symbol,
                "gene_type": rows[0].get("gene_type", ""),
                "chrom": rows[0].get("gene_chrom", ""),
                "start": rows[0].get("gene_start", ""),
                "end": rows[0].get("gene_end", ""),
                "strand": rows[0].get("gene_strand", ""),
                "n_linked_gse157078_probes": len(by_probe),
                "n_shsy_phase3_concordant_probes": len(phase3),
                "n_ds_concordant_probes": len(ds_concordant),
                "n_human_triple_phase3_any_ds_probes": len(triple_any),
                "n_human_triple_phase3_fetal_ds_probes": len(triple_fetal),
                "n_human_triple_phase3_fetal_ds_majority_probes": len(triple_fetal_majority),
                "n_human_triple_phase3_fetal_mean_positive_probes": len(triple_fetal_mean_positive),
                "n_promoter_triple_any_ds_probes": len(promoter_triple),
                "n_gene_body_triple_any_ds_probes": len(body_triple),
                "n_nearest_triple_any_ds_probes": len(nearest_triple),
                "mean_gse157078_delta_beta": mean([safe_float(row.get("gse157078_delta_beta")) for row in by_probe.values()]),
                "mean_shsy_phase3_signed_delta": mean(
                    [safe_float(row.get("shsy_phase3_signed_delta")) for row in by_probe.values()]
                ),
                "mean_ds_signed_delta": mean([safe_float(row.get("ds_mean_signed_delta")) for row in by_probe.values()]),
                "triple_probe_ids": ",".join(sorted(str(row["probe_id"]) for row in triple_any)[:30])
                + ("..." if len(triple_any) > 30 else ""),
                "human_methylation_support_score": methylation_score,
                "mouse_candidate": "yes" if cross_row else "no",
                "mouse_rank": cross_row.get("mouse_rank", ""),
                "mouse_gene_symbol": cross_row.get("mouse_gene_symbol", ""),
                "mouse_evidence_score": safe_float(cross_row.get("mouse_evidence_score")),
                "mouse_gse268951_direction": cross_row.get("mouse_gse268951_direction", ""),
                "mouse_dp10_direction": cross_row.get("mouse_dp10_direction", ""),
                "human_expression_support_level": cross_row.get("human_support_level", ""),
                "human_expression_support_score": safe_float(cross_row.get("human_support_score")),
                "cross_species_expression_priority_score": cross_score,
                "combined_methylation_expression_score": combined_score,
                "description": cross_row.get("description", ""),
            }
        )

    gene_rows.sort(
        key=lambda row: (
            -safe_float(row["combined_methylation_expression_score"]),
            -safe_float(row["n_human_triple_phase3_fetal_mean_positive_probes"]),
            -safe_float(row["n_human_triple_phase3_fetal_ds_majority_probes"]),
            -safe_float(row["n_human_triple_phase3_fetal_ds_probes"]),
            -safe_float(row["n_human_triple_phase3_any_ds_probes"]),
            row["gene_symbol"],
        )
    )
    return gene_rows


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)

    probes = read_probes(Path(args.liftover))
    shsy = read_shsy(Path(args.shsy_effects))
    ds = read_ds(Path(args.ds_values))
    genes_by_chrom = read_genes(Path(args.gencode_gtf))
    cross = read_cross_species(Path(args.cross_species))

    probe_rows_by_id = probe_support_rows(probes, shsy, ds)
    probe_rows = [probe_rows_by_id[probe_id] for probe_id in sorted(probe_rows_by_id)]
    links = link_probes_to_genes(
        probes,
        probe_rows_by_id,
        genes_by_chrom,
        args.promoter_window,
        args.nearest_window,
    )
    gene_rows = aggregate_genes(links, cross)
    convergent_rows = [
        row
        for row in gene_rows
        if row.get("mouse_candidate") == "yes" and int(row.get("n_human_triple_phase3_any_ds_probes") or 0) > 0
    ]

    probe_fields = [
        "probe_id",
        "direction",
        "gse157078_delta_beta",
        "hg38_chrom",
        "hg38_start",
        "hg38_end",
        "shsy_growth_delta_beta",
        "shsy_growth_signed_delta",
        "shsy_growth_concordant",
        "shsy_phase1_delta_beta",
        "shsy_phase1_signed_delta",
        "shsy_phase1_concordant",
        "shsy_phase3_delta_beta",
        "shsy_phase3_signed_delta",
        "shsy_phase3_concordant",
        "shsy_concordant_phase_count",
        "shsy_any_phase_concordant",
        "ds_contexts_tested",
        "ds_concordant_contexts",
        "ds_concordant_fraction",
        "ds_fetal_contexts_tested",
        "ds_fetal_concordant_contexts",
        "ds_fetal_concordant_fraction",
        "ds_adult_contexts_tested",
        "ds_adult_concordant_contexts",
        "ds_adult_concordant_fraction",
        "ds_mean_signed_delta",
        "ds_fetal_mean_signed_delta",
        "ds_adult_mean_signed_delta",
        "human_triple_phase3_any_ds",
        "human_triple_phase3_fetal_ds",
        "human_triple_phase3_two_ds_contexts",
        "human_triple_phase3_ds_majority",
        "human_triple_phase3_fetal_ds_majority",
        "human_triple_phase3_fetal_mean_positive",
        "ds_concordant_context_list",
        "ds_fetal_concordant_context_list",
        "ucsc_refgene_name",
        "ucsc_refgene_group",
        "gencode_name_array_annotation",
        "gencode_group_array_annotation",
    ]
    support_link_fields = [
        field for field in probe_fields if field not in {"probe_id", "hg38_chrom", "hg38_start", "hg38_end"}
    ]
    link_fields = [
        "probe_id",
        "hg38_chrom",
        "probe_start",
        "probe_end",
        "gene_id",
        "gene_symbol",
        "gene_type",
        "gene_chrom",
        "gene_start",
        "gene_end",
        "gene_strand",
        "gene_tss",
        "relationship",
        "promoter_2kb",
        "gene_body",
        "nearest_100kb",
        "signed_distance_to_tss",
        "abs_distance_to_tss",
        "distance_to_gene_body",
        *support_link_fields,
    ]
    gene_fields = [
        "gene_id",
        "gene_symbol",
        "gene_type",
        "chrom",
        "start",
        "end",
        "strand",
        "n_linked_gse157078_probes",
        "n_shsy_phase3_concordant_probes",
        "n_ds_concordant_probes",
        "n_human_triple_phase3_any_ds_probes",
        "n_human_triple_phase3_fetal_ds_probes",
        "n_human_triple_phase3_fetal_ds_majority_probes",
        "n_human_triple_phase3_fetal_mean_positive_probes",
        "n_promoter_triple_any_ds_probes",
        "n_gene_body_triple_any_ds_probes",
        "n_nearest_triple_any_ds_probes",
        "mean_gse157078_delta_beta",
        "mean_shsy_phase3_signed_delta",
        "mean_ds_signed_delta",
        "triple_probe_ids",
        "human_methylation_support_score",
        "mouse_candidate",
        "mouse_rank",
        "mouse_gene_symbol",
        "mouse_evidence_score",
        "mouse_gse268951_direction",
        "mouse_dp10_direction",
        "human_expression_support_level",
        "human_expression_support_score",
        "cross_species_expression_priority_score",
        "combined_methylation_expression_score",
        "description",
    ]

    write_tsv(outdir / "human_methylation_probe_support.tsv", probe_rows, probe_fields)
    write_tsv(outdir / "human_methylation_probe_gene_links.tsv", links, link_fields)
    write_tsv(outdir / "human_methylation_gene_support.tsv", gene_rows, gene_fields)
    write_tsv(outdir / "human_mouse_convergent_candidate_genes.tsv", convergent_rows, gene_fields)

    phase3_tested = [row for row in probe_rows if row.get("shsy_phase3_concordant") in {"yes", "no"}]
    summary_rows = [
        {"metric": "gse157078_candidate_probes_lifted_to_hg38", "value": len(probes)},
        {"metric": "shsy_phase3_tested_probes", "value": len(phase3_tested)},
        {
            "metric": "shsy_phase3_concordant_probes",
            "value": sum(1 for row in phase3_tested if row["shsy_phase3_concordant"] == "yes"),
        },
        {
            "metric": "ds_tested_probes",
            "value": sum(1 for row in probe_rows if int(row.get("ds_contexts_tested") or 0) > 0),
        },
        {
            "metric": "ds_any_context_concordant_probes",
            "value": sum(1 for row in probe_rows if int(row.get("ds_concordant_contexts") or 0) > 0),
        },
        {
            "metric": "ds_fetal_context_concordant_probes",
            "value": sum(1 for row in probe_rows if int(row.get("ds_fetal_concordant_contexts") or 0) > 0),
        },
        {
            "metric": "human_triple_phase3_any_ds_probes",
            "value": sum(1 for row in probe_rows if row["human_triple_phase3_any_ds"] == "yes"),
        },
        {
            "metric": "human_triple_phase3_fetal_ds_probes",
            "value": sum(1 for row in probe_rows if row["human_triple_phase3_fetal_ds"] == "yes"),
        },
        {
            "metric": "human_triple_phase3_fetal_ds_majority_probes",
            "value": sum(1 for row in probe_rows if row["human_triple_phase3_fetal_ds_majority"] == "yes"),
        },
        {
            "metric": "human_triple_phase3_fetal_mean_positive_probes",
            "value": sum(1 for row in probe_rows if row["human_triple_phase3_fetal_mean_positive"] == "yes"),
        },
        {"metric": "probe_gene_links", "value": len(links)},
        {"metric": "linked_genes", "value": len(gene_rows)},
        {
            "metric": "genes_with_human_triple_phase3_any_ds_probe",
            "value": sum(1 for row in gene_rows if int(row["n_human_triple_phase3_any_ds_probes"]) > 0),
        },
        {
            "metric": "genes_with_human_triple_phase3_fetal_ds_probe",
            "value": sum(1 for row in gene_rows if int(row["n_human_triple_phase3_fetal_ds_probes"]) > 0),
        },
        {
            "metric": "genes_with_human_triple_phase3_fetal_mean_positive_probe",
            "value": sum(1 for row in gene_rows if int(row["n_human_triple_phase3_fetal_mean_positive_probes"]) > 0),
        },
        {"metric": "mouse_candidate_genes_with_human_triple_any_ds", "value": len(convergent_rows)},
        {
            "metric": "mouse_candidate_genes_with_human_triple_fetal_ds",
            "value": sum(1 for row in convergent_rows if int(row["n_human_triple_phase3_fetal_ds_probes"]) > 0),
        },
        {
            "metric": "mouse_candidate_genes_with_human_triple_fetal_mean_positive",
            "value": sum(
                1 for row in convergent_rows if int(row["n_human_triple_phase3_fetal_mean_positive_probes"]) > 0
            ),
        },
    ]
    write_tsv(outdir / "human_methylation_convergence_summary.tsv", summary_rows, ["metric", "value"])


if __name__ == "__main__":
    main()
