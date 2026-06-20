#!/usr/bin/env python3
"""Validate GSE157078 DNMT3L-sensitive CpGs in GSE168276 SH-SY5Y WGBS."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from pyliftover import LiftOver


DEFAULT_DNMT3L_DMCS = "results/human_candidate_dmrs/GSE157078_candidate_dmcs.tsv"
DEFAULT_METADATA = "data/geo/GSE168276/metadata/sample_metadata.tsv"
DEFAULT_SUPPLEMENTARY_DIR = "data/geo/GSE168276/supplementary"
DEFAULT_CHAIN = "data/reference/liftover/hg19ToHg38.over.chain.gz"
DEFAULT_OUTDIR = "results/shsy5y_gse168276"
PHASE_ORDER = {"Growth": 0, "Phase 1": 1, "Phase 3": 2}
GENOTYPE_ORDER = {"GFP": 0, "DNMT3L+GFP": 1}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dnmt3l-dmcs", default=DEFAULT_DNMT3L_DMCS)
    parser.add_argument("--metadata", default=DEFAULT_METADATA)
    parser.add_argument("--supplementary-dir", default=DEFAULT_SUPPLEMENTARY_DIR)
    parser.add_argument("--chain", default=DEFAULT_CHAIN)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--min-coverage", type=int, default=5)
    parser.add_argument("--min-group-samples", type=int, default=2)
    parser.add_argument("--match-window", type=int, default=0, help="Allow +/- this many bp around lifted CpG start.")
    parser.add_argument("--workers", type=int, default=4)
    return parser.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="", errors="replace")
    return path.open(newline="", errors="replace")


def safe_float(value: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return number if math.isfinite(number) else float("nan")


def fmt(value: float) -> str:
    return f"{value:.6g}" if math.isfinite(value) else ""


def mean(values: list[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    return sum(finite) / len(finite) if finite else float("nan")


def parse_characteristics(value: str) -> dict[str, str]:
    fields = {}
    for token in value.split(";"):
        if "=" not in token:
            continue
        key, raw = token.split("=", 1)
        fields[key.strip().lower()] = raw.strip()
    return fields


def chrom_sort_key(chrom: str) -> tuple[int, int | str]:
    name = chrom.removeprefix("chr")
    if name.isdigit():
        return (0, int(name))
    special = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    return (0, special[name]) if name in special else (1, name)


def read_samples(metadata_path: Path, supplementary_dir: Path) -> list[dict[str, str]]:
    samples = []
    with metadata_path.open(newline="", errors="replace") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            chars = parse_characteristics(row["characteristics"])
            filename = Path(row["supplementary_data"]).name
            samples.append(
                {
                    "accession": row["accession"],
                    "sample": row["title"],
                    "genotype": chars.get("genotype", ""),
                    "phase": chars.get("differentiation phase", ""),
                    "replicate": chars.get("cell line replicate", ""),
                    "path": str(supplementary_dir / filename),
                }
            )
    return samples


def lift_candidate_dmcs(path: Path, chain_path: Path) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    liftover = LiftOver(str(chain_path))
    lifted = []
    unmapped = []
    with open_text(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            chrom = row["chrom"]
            start0 = int(row["start"])
            converted = liftover.convert_coordinate(chrom, start0)
            converted = [item for item in converted if item[0].startswith("chr") and "_" not in item[0]]
            if len(converted) != 1:
                unmapped.append({**row, "lift_status": f"ambiguous_or_unmapped:{len(converted)}"})
                continue
            hg38_chrom, hg38_start0, hg38_strand, _ = converted[0]
            hg38_start0 = int(hg38_start0)
            lifted.append(
                {
                    "probe_id": row["probe_id"],
                    "direction": row["direction"],
                    "gse157078_delta_beta": row["delta_beta"],
                    "hg19_chrom": row["chrom"],
                    "hg19_start": row["start"],
                    "hg19_end": row["end"],
                    "hg38_chrom": hg38_chrom,
                    "hg38_start": str(hg38_start0),
                    "hg38_end": str(hg38_start0 + 1),
                    "hg38_cpg_start_1based": str(hg38_start0 + 1),
                    "liftover_strand": hg38_strand,
                    "ucsc_refgene_name": row.get("ucsc_refgene_name", ""),
                    "ucsc_refgene_group": row.get("ucsc_refgene_group", ""),
                    "relation_to_cpg_island": row.get("relation_to_cpg_island", ""),
                    "regulatory_feature_group": row.get("regulatory_feature_group", ""),
                    "gencode_name": row.get("gencode_name", ""),
                    "gencode_group": row.get("gencode_group", ""),
                }
            )
    lifted.sort(key=lambda row: (chrom_sort_key(row["hg38_chrom"]), int(row["hg38_start"]), row["probe_id"]))
    return lifted, unmapped


def target_lookup(lifted: list[dict[str, str]], match_window: int) -> dict[tuple[str, int], list[tuple[str, int]]]:
    lookup: dict[tuple[str, int], list[tuple[str, int]]] = {}
    for row in lifted:
        chrom = row["hg38_chrom"]
        pos = int(row["hg38_cpg_start_1based"])
        for offset in range(-match_window, match_window + 1):
            lookup.setdefault((chrom, pos + offset), []).append((row["probe_id"], offset))
    return lookup


def normalized_cpg_start_1based(position_1based: int, strand: str) -> int:
    # Bismark cytosine reports list both cytosines in a CpG pair.  Collapse the
    # reverse-strand cytosine back to the plus-strand CpG start.
    return position_1based - 1 if strand == "-" else position_1based


def summarize_sample_at_targets(
    sample: dict[str, str],
    lookup: dict[tuple[str, int], list[tuple[str, int]]],
    min_coverage: int,
) -> list[dict[str, str]]:
    counts: dict[str, dict[str, int]] = {}
    offsets: dict[str, set[int]] = {}
    with gzip.open(sample["path"], "rt", errors="replace") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7 or parts[5] != "CG":
                continue
            chrom = parts[0]
            pos = int(parts[1])
            cpg_start = normalized_cpg_start_1based(pos, parts[2])
            matches = lookup.get((chrom, cpg_start))
            if not matches:
                continue
            methylated = int(parts[3])
            unmethylated = int(parts[4])
            for probe_id, offset in matches:
                rec = counts.setdefault(probe_id, {"methylated": 0, "unmethylated": 0})
                rec["methylated"] += methylated
                rec["unmethylated"] += unmethylated
                offsets.setdefault(probe_id, set()).add(offset)

    rows = []
    for probe_id, rec in counts.items():
        coverage = rec["methylated"] + rec["unmethylated"]
        if coverage < min_coverage:
            continue
        rows.append(
            {
                "accession": sample["accession"],
                "sample": sample["sample"],
                "genotype": sample["genotype"],
                "phase": sample["phase"],
                "replicate": sample["replicate"],
                "probe_id": probe_id,
                "methylated_reads": str(rec["methylated"]),
                "unmethylated_reads": str(rec["unmethylated"]),
                "coverage": str(coverage),
                "beta": fmt(rec["methylated"] / coverage),
                "match_offsets": ";".join(str(value) for value in sorted(offsets.get(probe_id, []))),
            }
        )
    rows.sort(key=lambda row: row["probe_id"])
    return rows


def phase_probe_effects(
    sample_rows: list[dict[str, str]],
    lifted_by_probe: dict[str, dict[str, str]],
    min_group_samples: int,
) -> list[dict[str, str]]:
    grouped: dict[tuple[str, str], dict[str, list[dict[str, str]]]] = {}
    for row in sample_rows:
        grouped.setdefault((row["phase"], row["probe_id"]), {}).setdefault(row["genotype"], []).append(row)

    effects = []
    for (phase, probe_id), by_genotype in sorted(grouped.items(), key=lambda item: (PHASE_ORDER.get(item[0][0], 99), item[0][1])):
        gfp = by_genotype.get("GFP", [])
        dnmt3l = by_genotype.get("DNMT3L+GFP", [])
        if len(gfp) < min_group_samples or len(dnmt3l) < min_group_samples:
            continue
        gfp_betas = [safe_float(row["beta"]) for row in gfp]
        dnmt3l_betas = [safe_float(row["beta"]) for row in dnmt3l]
        gfp_mean = mean(gfp_betas)
        dnmt3l_mean = mean(dnmt3l_betas)
        delta = dnmt3l_mean - gfp_mean
        meta = lifted_by_probe[probe_id]
        expected_sign = 1 if meta["direction"] == "hyper" else -1
        signed_delta = delta * expected_sign
        effects.append(
            {
                "phase": phase,
                "probe_id": probe_id,
                "direction": meta["direction"],
                "gse157078_delta_beta": meta["gse157078_delta_beta"],
                "hg38_chrom": meta["hg38_chrom"],
                "hg38_start": meta["hg38_start"],
                "hg38_end": meta["hg38_end"],
                "hg38_cpg_start_1based": meta["hg38_cpg_start_1based"],
                "n_gfp": str(len(gfp)),
                "n_dnmt3l_gfp": str(len(dnmt3l)),
                "gfp_mean_beta": fmt(gfp_mean),
                "dnmt3l_gfp_mean_beta": fmt(dnmt3l_mean),
                "shsy5y_delta_beta": fmt(delta),
                "signed_shsy5y_delta_by_gse157078_direction": fmt(signed_delta),
                "concordant": "yes" if signed_delta > 0 else "no",
                "ucsc_refgene_name": meta["ucsc_refgene_name"],
                "ucsc_refgene_group": meta["ucsc_refgene_group"],
                "relation_to_cpg_island": meta["relation_to_cpg_island"],
                "regulatory_feature_group": meta["regulatory_feature_group"],
                "gencode_name": meta["gencode_name"],
                "gencode_group": meta["gencode_group"],
            }
        )
    return effects


def summarize_effects(effects: list[dict[str, str]], lifted_count: int) -> list[dict[str, str]]:
    rows = []
    for phase in sorted({row["phase"] for row in effects}, key=lambda value: PHASE_ORDER.get(value, 99)):
        phase_rows = [row for row in effects if row["phase"] == phase]
        signed = [safe_float(row["signed_shsy5y_delta_by_gse157078_direction"]) for row in phase_rows]
        deltas = [safe_float(row["shsy5y_delta_beta"]) for row in phase_rows]
        gse_deltas = [safe_float(row["gse157078_delta_beta"]) for row in phase_rows]
        hyper = [row for row in phase_rows if row["direction"] == "hyper"]
        hypo = [row for row in phase_rows if row["direction"] == "hypo"]
        concordant = [value for value in signed if value > 0]
        discordant = [value for value in signed if value < 0]
        rows.append(
            {
                "phase": phase,
                "lifted_candidate_probes": str(lifted_count),
                "tested_probes": str(len(phase_rows)),
                "tested_fraction_of_lifted": fmt(len(phase_rows) / lifted_count if lifted_count else float("nan")),
                "n_hyper_signature_probes": str(len(hyper)),
                "n_hypo_signature_probes": str(len(hypo)),
                "mean_gse157078_delta_beta": fmt(mean(gse_deltas)),
                "mean_shsy5y_delta_beta": fmt(mean(deltas)),
                "mean_signed_shsy5y_delta_by_gse157078_direction": fmt(mean(signed)),
                "concordant_n": str(len(concordant)),
                "discordant_n": str(len(discordant)),
                "concordant_fraction": fmt(len(concordant) / len(phase_rows) if phase_rows else float("nan")),
                "same_direction_ge_0.02": str(sum(1 for value in signed if value >= 0.02)),
                "same_direction_ge_0.05": str(sum(1 for value in signed if value >= 0.05)),
                "opposite_direction_ge_0.02": str(sum(1 for value in signed if value <= -0.02)),
                "opposite_direction_ge_0.05": str(sum(1 for value in signed if value <= -0.05)),
                "mean_abs_shsy5y_delta_beta": fmt(mean([abs(value) for value in deltas])),
            }
        )
    return rows


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    lifted, unmapped = lift_candidate_dmcs(Path(args.dnmt3l_dmcs), Path(args.chain))
    lifted_by_probe = {row["probe_id"]: row for row in lifted}
    lookup = target_lookup(lifted, args.match_window)
    samples = read_samples(Path(args.metadata), Path(args.supplementary_dir))
    missing = [sample["path"] for sample in samples if not Path(sample["path"]).exists()]
    if missing:
        raise FileNotFoundError(f"Missing {len(missing)} GSE168276 CpG reports: {missing[:3]}")

    sample_rows = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(summarize_sample_at_targets, sample, lookup, args.min_coverage): sample for sample in samples}
        for future in as_completed(futures):
            sample = futures[future]
            rows = future.result()
            print(f"Summarized {sample['sample']}: {len(rows)} covered target probes", flush=True)
            sample_rows.extend(rows)
    sample_rows.sort(
        key=lambda row: (
            PHASE_ORDER.get(row["phase"], 99),
            GENOTYPE_ORDER.get(row["genotype"], 99),
            int(row["replicate"]),
            row["probe_id"],
        )
    )

    effects = phase_probe_effects(sample_rows, lifted_by_probe, args.min_group_samples)
    summary = summarize_effects(effects, len(lifted))

    lifted_fields = [
        "probe_id",
        "direction",
        "gse157078_delta_beta",
        "hg19_chrom",
        "hg19_start",
        "hg19_end",
        "hg38_chrom",
        "hg38_start",
        "hg38_end",
        "hg38_cpg_start_1based",
        "liftover_strand",
        "ucsc_refgene_name",
        "ucsc_refgene_group",
        "relation_to_cpg_island",
        "regulatory_feature_group",
        "gencode_name",
        "gencode_group",
    ]
    sample_fields = [
        "accession",
        "sample",
        "genotype",
        "phase",
        "replicate",
        "probe_id",
        "methylated_reads",
        "unmethylated_reads",
        "coverage",
        "beta",
        "match_offsets",
    ]
    effect_fields = [
        "phase",
        "probe_id",
        "direction",
        "gse157078_delta_beta",
        "hg38_chrom",
        "hg38_start",
        "hg38_end",
        "hg38_cpg_start_1based",
        "n_gfp",
        "n_dnmt3l_gfp",
        "gfp_mean_beta",
        "dnmt3l_gfp_mean_beta",
        "shsy5y_delta_beta",
        "signed_shsy5y_delta_by_gse157078_direction",
        "concordant",
        "ucsc_refgene_name",
        "ucsc_refgene_group",
        "relation_to_cpg_island",
        "regulatory_feature_group",
        "gencode_name",
        "gencode_group",
    ]
    summary_fields = [
        "phase",
        "lifted_candidate_probes",
        "tested_probes",
        "tested_fraction_of_lifted",
        "n_hyper_signature_probes",
        "n_hypo_signature_probes",
        "mean_gse157078_delta_beta",
        "mean_shsy5y_delta_beta",
        "mean_signed_shsy5y_delta_by_gse157078_direction",
        "concordant_n",
        "discordant_n",
        "concordant_fraction",
        "same_direction_ge_0.02",
        "same_direction_ge_0.05",
        "opposite_direction_ge_0.02",
        "opposite_direction_ge_0.05",
        "mean_abs_shsy5y_delta_beta",
    ]
    write_tsv(outdir / "gse157078_candidate_dmcs_hg38_liftover.tsv", lifted, lifted_fields)
    write_tsv(outdir / "gse157078_candidate_dmcs_hg38_liftover_unmapped.tsv", unmapped, list(unmapped[0].keys()) if unmapped else ["probe_id"])
    write_tsv(outdir / "gse168276_gse157078_target_sample_values.tsv", sample_rows, sample_fields)
    write_tsv(outdir / "gse168276_gse157078_target_phase_effects.tsv", effects, effect_fields)
    write_tsv(outdir / "gse168276_gse157078_target_concordance_summary.tsv", summary, summary_fields)
    print(f"Lifted {len(lifted)} candidate probes; unmapped/ambiguous {len(unmapped)}")
    print(f"Wrote {len(sample_rows)} target sample rows")
    print(f"Wrote {len(effects)} phase effect rows")
    print(f"Wrote {len(summary)} summary rows")


if __name__ == "__main__":
    main()
