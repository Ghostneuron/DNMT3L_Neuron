#!/usr/bin/env python3
"""Call first-pass human EPIC candidate DMCs/DMRs for GSE157078."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from pathlib import Path


SAMPLES = [
    {"sample": "Z1", "condition": "control", "meth": "SAMPLE 1 Methylated signal", "unmeth": "SAMPLE 1 Unmethylated Signal", "pval": "SAMPLE 1 Detection Pval"},
    {"sample": "Z2", "condition": "control", "meth": "SAMPLE 2 Methylated signal", "unmeth": "SAMPLE 2 Unmethylated Signal", "pval": "SAMPLE 2 Detection Pval"},
    {"sample": "Z3", "condition": "control", "meth": "SAMPLE 3 Methylated signal", "unmeth": "SAMPLE 3 Unmethylated Signal", "pval": "SAMPLE 3 Detection Pval"},
    {"sample": "D1", "condition": "DNMT3L_overexpression", "meth": "SAMPLE 4 Methylated signal", "unmeth": "SAMPLE 4 Unmethylated Signal", "pval": "SAMPLE 4 Detection Pval"},
    {"sample": "D2", "condition": "DNMT3L_overexpression", "meth": "SAMPLE 5 Methylated signal", "unmeth": "SAMPLE 5 Unmethylated Signal", "pval": "SAMPLE 5 Detection Pval"},
    {"sample": "D3", "condition": "DNMT3L_overexpression", "meth": "SAMPLE 6 Methylated signal", "unmeth": "SAMPLE 6 Unmethylated Signal", "pval": "SAMPLE 6 Detection Pval"},
]

CONTROL = ["Z1", "Z2", "Z3"]
CASE = ["D1", "D2", "D3"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--signals", default="data/geo/GSE157078/supplementary/GSE157078_signals.txt.gz")
    parser.add_argument(
        "--manifest",
        default="data/geo/GPL21145/GPL21145_MethylationEPIC_15073387_v-1-0.csv.gz",
        help="GEO GPL21145 EPIC annotation CSV.",
    )
    parser.add_argument("--outdir", default="results/human_candidate_dmrs")
    parser.add_argument("--min-delta", type=float, default=0.05)
    parser.add_argument("--max-gap", type=int, default=1000)
    parser.add_argument("--min-probes", type=int, default=2)
    parser.add_argument("--detection-p", type=float, default=0.01)
    parser.add_argument("--beta-offset", type=float, default=100.0)
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


def mean(values: list[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    return sum(finite) / len(finite) if finite else float("nan")


def fmt(value: float) -> str:
    return f"{value:.6g}" if math.isfinite(value) else ""


def chrom_name(value: str) -> str:
    value = value or ""
    value = value.strip()
    if not value:
        return ""
    return value if value.startswith("chr") else f"chr{value}"


def chrom_sort_key(chrom: str) -> tuple[int, int | str]:
    name = chrom.removeprefix("chr")
    if name.isdigit():
        return (0, int(name))
    special = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    return (0, special[name]) if name in special else (1, name)


def read_epic_annotation(path: Path) -> dict[str, dict[str, str]]:
    annotations: dict[str, dict[str, str]] = {}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("IlmnID,"):
                header = next(csv.reader([line]))
                break
        else:
            raise RuntimeError(f"Could not find IlmnID header in {path}")

        reader = csv.DictReader(handle, fieldnames=header)
        for row in reader:
            probe = row.get("Name") or row.get("IlmnID")
            chrom = chrom_name(row.get("CHR", ""))
            pos = (row.get("MAPINFO") or "").strip()
            if not probe or not chrom or not pos.isdigit():
                continue
            annotations[probe] = {
                "chrom": chrom,
                "position": pos,
                "ucsc_refgene_name": row.get("UCSC_RefGene_Name", ""),
                "ucsc_refgene_group": row.get("UCSC_RefGene_Group", ""),
                "cpg_island": row.get("UCSC_CpG_Islands_Name", ""),
                "relation_to_cpg_island": row.get("Relation_to_UCSC_CpG_Island", ""),
                "regulatory_feature_name": row.get("Regulatory_Feature_Name", ""),
                "regulatory_feature_group": row.get("Regulatory_Feature_Group", ""),
                "gencode_name": row.get("GencodeBasicV12_NAME", ""),
                "gencode_group": row.get("GencodeBasicV12_Group", ""),
                "snp_id": row.get("SNP_ID", ""),
                "snp_distance": row.get("SNP_DISTANCE", ""),
                "snp_maf": row.get("SNP_MinorAlleleFrequency", ""),
            }
    return annotations


def beta_values(row: dict[str, str], detection_p: float, beta_offset: float) -> dict[str, float] | None:
    values = {}
    for sample in SAMPLES:
        pvalue = safe_float(row.get(sample["pval"], ""))
        if not math.isfinite(pvalue) or pvalue > detection_p:
            return None
        methylated = safe_float(row.get(sample["meth"], ""))
        unmethylated = safe_float(row.get(sample["unmeth"], ""))
        denominator = methylated + unmethylated + beta_offset
        if not math.isfinite(denominator) or denominator <= 0:
            return None
        values[sample["sample"]] = methylated / denominator
    return values


def direction(delta: float, min_delta: float) -> str:
    if delta >= min_delta:
        return "hyper"
    if delta <= -min_delta:
        return "hypo"
    return "none"


def read_candidate_dmcs(
    signal_path: Path,
    annotations: dict[str, dict[str, str]],
    min_delta: float,
    detection_p: float,
    beta_offset: float,
) -> tuple[list[dict[str, object]], dict[str, int]]:
    candidates: list[dict[str, object]] = []
    counts = {"seen": 0, "passed_detection": 0, "annotated": 0}
    with open_text(signal_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            counts["seen"] += 1
            probe = row["ID_REF"]
            values = beta_values(row, detection_p=detection_p, beta_offset=beta_offset)
            if values is None:
                continue
            counts["passed_detection"] += 1
            anno = annotations.get(probe)
            if not anno:
                continue
            counts["annotated"] += 1
            control_mean = mean([values[sample] for sample in CONTROL])
            case_mean = mean([values[sample] for sample in CASE])
            delta = case_mean - control_mean
            call = direction(delta, min_delta)
            if call == "none":
                continue
            pos = int(anno["position"])
            candidates.append(
                {
                    "chrom": anno["chrom"],
                    "position": pos,
                    "start": pos - 1,
                    "end": pos,
                    "probe_id": probe,
                    "direction": call,
                    "mean_control_beta": control_mean,
                    "mean_case_beta": case_mean,
                    "delta_beta": delta,
                    **anno,
                }
            )
    candidates.sort(key=lambda row: (chrom_sort_key(str(row["chrom"])), int(row["position"])))
    return candidates, counts


def dmc_for_output(row: dict[str, object]) -> dict[str, str]:
    return {
        "chrom": str(row["chrom"]),
        "start": str(row["start"]),
        "end": str(row["end"]),
        "probe_id": str(row["probe_id"]),
        "direction": str(row["direction"]),
        "mean_control_beta": fmt(float(row["mean_control_beta"])),
        "mean_case_beta": fmt(float(row["mean_case_beta"])),
        "delta_beta": fmt(float(row["delta_beta"])),
        "ucsc_refgene_name": str(row.get("ucsc_refgene_name", "")),
        "ucsc_refgene_group": str(row.get("ucsc_refgene_group", "")),
        "cpg_island": str(row.get("cpg_island", "")),
        "relation_to_cpg_island": str(row.get("relation_to_cpg_island", "")),
        "regulatory_feature_name": str(row.get("regulatory_feature_name", "")),
        "regulatory_feature_group": str(row.get("regulatory_feature_group", "")),
        "gencode_name": str(row.get("gencode_name", "")),
        "gencode_group": str(row.get("gencode_group", "")),
        "snp_id": str(row.get("snp_id", "")),
        "snp_distance": str(row.get("snp_distance", "")),
        "snp_maf": str(row.get("snp_maf", "")),
    }


def should_group(prev: dict[str, object], row: dict[str, object], max_gap: int) -> bool:
    return (
        prev["chrom"] == row["chrom"]
        and prev["direction"] == row["direction"]
        and int(row["position"]) - int(prev["position"]) <= max_gap
    )


def unique_join(values: list[str]) -> str:
    seen = []
    for value in values:
        for token in value.split(";"):
            token = token.strip()
            if token and token not in seen:
                seen.append(token)
    return ";".join(seen)


def summarize_region(rows: list[dict[str, object]], region_id: int) -> dict[str, str]:
    chrom = str(rows[0]["chrom"])
    start = min(int(row["start"]) for row in rows)
    end = max(int(row["end"]) for row in rows)
    deltas = [float(row["delta_beta"]) for row in rows]
    control = [float(row["mean_control_beta"]) for row in rows]
    case = [float(row["mean_case_beta"]) for row in rows]
    return {
        "dmr_id": f"GSE157078_DMR_{region_id:06d}",
        "chrom": chrom,
        "start": str(start),
        "end": str(end),
        "direction": str(rows[0]["direction"]),
        "n_probes": str(len(rows)),
        "width_bp": str(end - start),
        "mean_control_beta": fmt(mean(control)),
        "mean_case_beta": fmt(mean(case)),
        "mean_delta_beta": fmt(mean(deltas)),
        "max_abs_delta_beta": fmt(max(abs(delta) for delta in deltas)),
        "probe_ids": ";".join(str(row["probe_id"]) for row in rows),
        "ucsc_refgene_name": unique_join([str(row.get("ucsc_refgene_name", "")) for row in rows]),
        "ucsc_refgene_group": unique_join([str(row.get("ucsc_refgene_group", "")) for row in rows]),
        "relation_to_cpg_island": unique_join([str(row.get("relation_to_cpg_island", "")) for row in rows]),
        "regulatory_feature_group": unique_join([str(row.get("regulatory_feature_group", "")) for row in rows]),
        "gencode_name": unique_join([str(row.get("gencode_name", "")) for row in rows]),
        "gencode_group": unique_join([str(row.get("gencode_group", "")) for row in rows]),
    }


def call_dmrs(candidates: list[dict[str, object]], max_gap: int, min_probes: int) -> list[dict[str, str]]:
    regions = []
    current: list[dict[str, object]] = []
    region_id = 1
    for row in candidates:
        if not current:
            current = [row]
            continue
        if should_group(current[-1], row, max_gap):
            current.append(row)
            continue
        if len(current) >= min_probes:
            regions.append(summarize_region(current, region_id))
            region_id += 1
        current = [row]
    if len(current) >= min_probes:
        regions.append(summarize_region(current, region_id))
    return regions


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_bed(path: Path, regions: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        for region in regions:
            score = min(1000, int(abs(float(region["mean_delta_beta"])) * 1000))
            handle.write(
                "\t".join(
                    [
                        region["chrom"],
                        region["start"],
                        region["end"],
                        region["dmr_id"],
                        str(score),
                        ".",
                        region["direction"],
                        region["n_probes"],
                        region["mean_delta_beta"],
                    ]
                )
                + "\n"
            )


def write_summary(
    path: Path,
    candidates: list[dict[str, object]],
    regions: list[dict[str, str]],
    counts: dict[str, int],
    args: argparse.Namespace,
) -> None:
    rows = [
        {
            "dataset": "GSE157078",
            "genome_build": "hg19_GRCh37_from_GPL21145",
            "min_delta": fmt(args.min_delta),
            "max_gap": str(args.max_gap),
            "min_probes": str(args.min_probes),
            "detection_p": fmt(args.detection_p),
            "signals_seen": str(counts["seen"]),
            "signals_passed_detection": str(counts["passed_detection"]),
            "signals_annotated": str(counts["annotated"]),
            "candidate_probes": str(len(candidates)),
            "hyper_candidate_probes": str(sum(1 for row in candidates if row["direction"] == "hyper")),
            "hypo_candidate_probes": str(sum(1 for row in candidates if row["direction"] == "hypo")),
            "candidate_dmrs": str(len(regions)),
            "hyper_candidate_dmrs": str(sum(1 for row in regions if row["direction"] == "hyper")),
            "hypo_candidate_dmrs": str(sum(1 for row in regions if row["direction"] == "hypo")),
        }
    ]
    write_tsv(
        path,
        rows,
        [
            "dataset",
            "genome_build",
            "min_delta",
            "max_gap",
            "min_probes",
            "detection_p",
            "signals_seen",
            "signals_passed_detection",
            "signals_annotated",
            "candidate_probes",
            "hyper_candidate_probes",
            "hypo_candidate_probes",
            "candidate_dmrs",
            "hyper_candidate_dmrs",
            "hypo_candidate_dmrs",
        ],
    )


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    annotations = read_epic_annotation(Path(args.manifest))
    candidates, counts = read_candidate_dmcs(
        Path(args.signals),
        annotations,
        min_delta=args.min_delta,
        detection_p=args.detection_p,
        beta_offset=args.beta_offset,
    )
    regions = call_dmrs(candidates, max_gap=args.max_gap, min_probes=args.min_probes)

    dmc_fields = [
        "chrom",
        "start",
        "end",
        "probe_id",
        "direction",
        "mean_control_beta",
        "mean_case_beta",
        "delta_beta",
        "ucsc_refgene_name",
        "ucsc_refgene_group",
        "cpg_island",
        "relation_to_cpg_island",
        "regulatory_feature_name",
        "regulatory_feature_group",
        "gencode_name",
        "gencode_group",
        "snp_id",
        "snp_distance",
        "snp_maf",
    ]
    dmr_fields = [
        "dmr_id",
        "chrom",
        "start",
        "end",
        "direction",
        "n_probes",
        "width_bp",
        "mean_control_beta",
        "mean_case_beta",
        "mean_delta_beta",
        "max_abs_delta_beta",
        "probe_ids",
        "ucsc_refgene_name",
        "ucsc_refgene_group",
        "relation_to_cpg_island",
        "regulatory_feature_group",
        "gencode_name",
        "gencode_group",
    ]
    write_tsv(outdir / "GSE157078_candidate_dmcs.tsv", [dmc_for_output(row) for row in candidates], dmc_fields)
    write_tsv(outdir / "GSE157078_candidate_dmrs.tsv", regions, dmr_fields)
    write_bed(outdir / "GSE157078_candidate_dmrs_hg19.bed", regions)
    write_summary(outdir / "GSE157078_candidate_dmr_summary.tsv", candidates, regions, counts, args)
    print(
        f"Wrote {len(candidates)} candidate probes and {len(regions)} candidate DMRs to {outdir}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
