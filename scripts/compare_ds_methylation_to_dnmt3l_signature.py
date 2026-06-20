#!/usr/bin/env python3
"""Compare DS brain methylation datasets with the DNMT3L-OE methylation signature.

This is a candidate-probe first pass.  The direct DNMT3L-overexpression anchor
is GSE157078, called as EPIC candidate DMCs in
results/human_candidate_dmrs/GSE157078_candidate_dmcs.tsv.  DS datasets are
450K-era arrays, so probe IDs provide the cleanest conservative intersection.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
from collections import defaultdict
from pathlib import Path


DEFAULT_DNMT3L_DMCS = "results/human_candidate_dmrs/GSE157078_candidate_dmcs.tsv"
DEFAULT_OUTDIR = "results/ds_methylation_dnmt3l_concordance"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dnmt3l-dmcs", default=DEFAULT_DNMT3L_DMCS)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--min-group-samples", type=int, default=2)
    parser.add_argument("--detection-p", type=float, default=0.01)
    parser.add_argument("--gse74486-detection-p", type=float, default=0.05)
    parser.add_argument("--beta-offset", type=float, default=100.0)
    parser.add_argument("--min-absolute-ds-delta", type=float, default=0.02)
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


def median(values: list[float]) -> float:
    finite = sorted(value for value in values if math.isfinite(value))
    if not finite:
        return float("nan")
    midpoint = len(finite) // 2
    if len(finite) % 2:
        return finite[midpoint]
    return (finite[midpoint - 1] + finite[midpoint]) / 2


def slug(value: str) -> str:
    value = value.lower().strip()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return value.strip("_") or "unknown"


def parse_characteristics(value: str) -> dict[str, str]:
    fields = {}
    for token in value.split(";"):
        if "=" not in token:
            continue
        key, raw = token.split("=", 1)
        fields[key.strip().lower()] = raw.strip()
    return fields


def ds_group(value: str) -> str:
    normalized = value.lower().strip()
    if "down" in normalized or "trisomy" in normalized or normalized == "ds":
        return "DS"
    if "normal" in normalized or "control" in normalized or normalized == "n":
        return "Control"
    return ""


def read_dnmt3l_signature(path: Path) -> dict[str, dict[str, str]]:
    signature = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            probe = row["probe_id"]
            signature[probe] = {
                "probe_id": probe,
                "dnmt3l_direction": row["direction"],
                "dnmt3l_delta_beta": row["delta_beta"],
                "chrom": row.get("chrom", ""),
                "start": row.get("start", ""),
                "end": row.get("end", ""),
                "ucsc_refgene_name": row.get("ucsc_refgene_name", ""),
                "ucsc_refgene_group": row.get("ucsc_refgene_group", ""),
                "relation_to_cpg_island": row.get("relation_to_cpg_island", ""),
                "regulatory_feature_group": row.get("regulatory_feature_group", ""),
                "gencode_name": row.get("gencode_name", ""),
                "gencode_group": row.get("gencode_group", ""),
            }
    return signature


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def gse73747_samples(path: Path) -> list[dict[str, str]]:
    samples = []
    with path.open(newline="", errors="replace") as handle:
        for index, row in enumerate(csv.DictReader(handle, delimiter="\t"), start=1):
            chars = parse_characteristics(row["characteristics"])
            group = ds_group(chars.get("disease state", ""))
            tissue = chars.get("tissue", row["source_name"])
            stage = "fetal" if "fetal" in tissue.lower() else "adult"
            contexts = ["GSE73747_all_brain"]
            if stage == "fetal":
                contexts.append("GSE73747_fetal_all")
            if "frontal" in tissue.lower():
                contexts.append(f"GSE73747_{stage}_frontal")
            if "temporal" in tissue.lower():
                contexts.append(f"GSE73747_{stage}_temporal")
            samples.append(
                {
                    "dataset": "GSE73747",
                    "sample_key": str(index),
                    "matrix_id": f"Sample {index}",
                    "accession": row["accession"],
                    "group": group,
                    "context_labels": ";".join(contexts),
                    "tissue": tissue,
                    "stage": stage,
                    "cell_type": "bulk brain",
                    "age": chars.get("gestational.week") or chars.get("adult.age", ""),
                    "sex": chars.get("gender", ""),
                }
            )
    return samples


def gse63347_samples(path: Path) -> list[dict[str, str]]:
    samples = []
    with path.open(newline="", errors="replace") as handle:
        for index, row in enumerate(csv.DictReader(handle, delimiter="\t"), start=1):
            chars = parse_characteristics(row["characteristics"])
            tissue = chars.get("tissue", row["source_name"])
            group = ds_group(chars.get("downsyndrome.status", ""))
            tissue_context = slug(tissue.replace("Brain", "").strip())
            contexts = ["GSE63347_adult_all_brain", f"GSE63347_adult_{tissue_context}"]
            samples.append(
                {
                    "dataset": "GSE63347",
                    "sample_key": str(index),
                    "matrix_id": f"Sample{index}",
                    "accession": row["accession"],
                    "group": group,
                    "context_labels": ";".join(contexts),
                    "tissue": tissue,
                    "stage": "adult",
                    "cell_type": "bulk brain",
                    "age": chars.get("age", ""),
                    "sex": chars.get("sex", ""),
                }
            )
    return samples


def gse74486_xml_samples(path: Path) -> list[dict[str, str]]:
    text = path.read_text(errors="replace")
    blocks = re.findall(r"<Sample iid=\".*?</Sample>", text, flags=re.S)
    samples = []
    for block in blocks:
        accession = extract_xml_text(block, "Accession")
        title = extract_xml_text(block, "Title")
        source = extract_xml_text(block, "Source")
        description = extract_xml_text(block, "Description")
        matrix_ids = re.findall(r"\b(?:ds|tb)\d+\b", description)
        matrix_id = matrix_ids[-1] if matrix_ids else ""
        if not matrix_id:
            continue
        chars = {
            match.group(1).strip().lower(): clean_xml_text(match.group(2))
            for match in re.finditer(
                r"<Characteristics tag=\"([^\"]+)\">\s*(.*?)\s*</Characteristics>",
                block,
                flags=re.S,
            )
        }
        tissue = chars.get("tissue", source)
        cell_type = chars.get("cell type", "")
        stage = chars.get("developmental stage", "")
        group = ds_group(chars.get("disease state", ""))
        if "blood" in tissue.lower():
            continue
        contexts = ["GSE74486_brain_all"]
        if stage:
            contexts.append(f"GSE74486_{slug(stage)}_all_brain")
        if tissue and cell_type:
            contexts.append(f"GSE74486_{slug(stage)}_{slug(tissue)}_{slug(cell_type)}")
        elif tissue:
            contexts.append(f"GSE74486_{slug(stage)}_{slug(tissue)}")
        samples.append(
            {
                "dataset": "GSE74486",
                "sample_key": matrix_id,
                "matrix_id": matrix_id,
                "accession": accession,
                "group": group,
                "context_labels": ";".join(contexts),
                "tissue": tissue,
                "stage": stage,
                "cell_type": cell_type,
                "age": "",
                "sex": chars.get("gender", ""),
                "title": title,
            }
        )
    return samples


def extract_xml_text(block: str, tag: str) -> str:
    match = re.search(rf"<{tag}(?:\s[^>]*)?>\s*(.*?)\s*</{tag}>", block, flags=re.S)
    return clean_xml_text(match.group(1)) if match else ""


def clean_xml_text(value: str) -> str:
    return re.sub(r"\s+", " ", value.replace("&gt;", ">").replace("&lt;", "<").replace("&amp;", "&")).strip()


def eligible_contexts(samples: list[dict[str, str]], min_group_samples: int) -> dict[str, dict[str, list[dict[str, str]]]]:
    contexts: dict[str, dict[str, list[dict[str, str]]]] = defaultdict(lambda: {"Control": [], "DS": []})
    for sample in samples:
        group = sample["group"]
        if group not in {"Control", "DS"}:
            continue
        for context in sample["context_labels"].split(";"):
            contexts[context][group].append(sample)
    return {
        context: groups
        for context, groups in contexts.items()
        if len(groups["Control"]) >= min_group_samples and len(groups["DS"]) >= min_group_samples
    }


def sample_context_rows(contexts: dict[str, dict[str, list[dict[str, str]]]]) -> list[dict[str, str]]:
    rows = []
    for context, groups in sorted(contexts.items()):
        examples = groups["Control"][:2] + groups["DS"][:2]
        rows.append(
            {
                "context": context,
                "n_control": str(len(groups["Control"])),
                "n_ds": str(len(groups["DS"])),
                "control_accessions": ";".join(sample["accession"] for sample in groups["Control"]),
                "ds_accessions": ";".join(sample["accession"] for sample in groups["DS"]),
                "example_tissues": ";".join(sorted({sample["tissue"] for sample in examples})),
                "example_cell_types": ";".join(sorted({sample["cell_type"] for sample in examples})),
            }
        )
    return rows


def mean_beta_for_signal_row(
    row: list[str],
    specs: list[dict[str, int | str]],
    detection_p: float,
    beta_offset: float,
) -> tuple[float, int]:
    values = []
    for spec in specs:
        methylated = safe_float(row[int(spec["meth_index"])])
        unmethylated = safe_float(row[int(spec["unmeth_index"])])
        pvalue = safe_float(row[int(spec["pval_index"])])
        if not math.isfinite(pvalue) or pvalue > detection_p:
            continue
        denominator = methylated + unmethylated + beta_offset
        if not math.isfinite(denominator) or denominator <= 0:
            continue
        values.append(methylated / denominator)
    return mean(values), len(values)


def mean_beta_for_processed_row(row: list[str], specs: list[dict[str, int | str]], detection_p: float) -> tuple[float, int]:
    values = []
    for spec in specs:
        beta = safe_float(row[int(spec["beta_index"])])
        pvalue = safe_float(row[int(spec["pval_index"])])
        if not math.isfinite(beta) or not math.isfinite(pvalue) or pvalue > detection_p:
            continue
        values.append(beta)
    return mean(values), len(values)


def signal_matrix_sample_specs(header: list[str], samples: list[dict[str, str]], compact: bool) -> dict[str, dict[str, int | str]]:
    specs = {}
    for sample in samples:
        number = int(sample["sample_key"])
        if compact:
            meth_name = f"Sample{number}Methylated"
            unmeth_name = f"Sample{number}Unmethylated"
            pval_name = f"Sample{number}DectionPvalue"
        else:
            meth_name = f"Sample {number} Methylated Signal"
            unmeth_name = f"Sample {number} Unmethylated Signal"
            pval_name = "Detection Pval"
        try:
            meth_index = header.index(meth_name)
            unmeth_index = header.index(unmeth_name)
        except ValueError:
            continue
        if compact:
            try:
                pval_index = header.index(pval_name)
            except ValueError:
                continue
        else:
            pval_index = unmeth_index + 1
            if pval_index >= len(header) or header[pval_index] != pval_name:
                continue
        specs[sample["sample_key"]] = {
            "meth_index": meth_index,
            "unmeth_index": unmeth_index,
            "pval_index": pval_index,
        }
    return specs


def processed_matrix_sample_specs(header: list[str], samples: list[dict[str, str]]) -> dict[str, dict[str, int | str]]:
    specs = {}
    header_index = {name: index for index, name in enumerate(header)}
    for sample in samples:
        matrix_id = sample["matrix_id"]
        pval_name = f"{matrix_id}_detectionpval"
        if matrix_id not in header_index or pval_name not in header_index:
            continue
        specs[sample["sample_key"]] = {
            "beta_index": header_index[matrix_id],
            "pval_index": header_index[pval_name],
        }
    return specs


def compare_signal_dataset(
    dataset: str,
    matrix_path: Path,
    samples: list[dict[str, str]],
    signature: dict[str, dict[str, str]],
    min_group_samples: int,
    detection_p: float,
    beta_offset: float,
    compact_headers: bool,
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    contexts = eligible_contexts(samples, min_group_samples)
    detail_rows = []
    with open_text(matrix_path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        sample_specs = signal_matrix_sample_specs(header, samples, compact=compact_headers)
        context_specs = build_context_specs(contexts, sample_specs)
        for row in reader:
            probe = row[0]
            sig = signature.get(probe)
            if not sig:
                continue
            for context, groups in context_specs.items():
                control_mean, control_n = mean_beta_for_signal_row(row, groups["Control"], detection_p, beta_offset)
                ds_mean, ds_n = mean_beta_for_signal_row(row, groups["DS"], detection_p, beta_offset)
                if control_n < min_group_samples or ds_n < min_group_samples:
                    continue
                detail_rows.append(build_detail_row(dataset, context, sig, control_mean, ds_mean, control_n, ds_n))
    summary_rows = summarize_details(detail_rows)
    return summary_rows, detail_rows, sample_context_rows(contexts)


def compare_processed_dataset(
    dataset: str,
    matrix_path: Path,
    samples: list[dict[str, str]],
    signature: dict[str, dict[str, str]],
    min_group_samples: int,
    detection_p: float,
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    contexts = eligible_contexts(samples, min_group_samples)
    detail_rows = []
    with open_text(matrix_path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        sample_specs = processed_matrix_sample_specs(header, samples)
        context_specs = build_context_specs(contexts, sample_specs)
        for row in reader:
            probe = row[0]
            sig = signature.get(probe)
            if not sig:
                continue
            for context, groups in context_specs.items():
                control_mean, control_n = mean_beta_for_processed_row(row, groups["Control"], detection_p)
                ds_mean, ds_n = mean_beta_for_processed_row(row, groups["DS"], detection_p)
                if control_n < min_group_samples or ds_n < min_group_samples:
                    continue
                detail_rows.append(build_detail_row(dataset, context, sig, control_mean, ds_mean, control_n, ds_n))
    summary_rows = summarize_details(detail_rows)
    return summary_rows, detail_rows, sample_context_rows(contexts)


def build_context_specs(
    contexts: dict[str, dict[str, list[dict[str, str]]]],
    sample_specs: dict[str, dict[str, int | str]],
) -> dict[str, dict[str, list[dict[str, int | str]]]]:
    context_specs: dict[str, dict[str, list[dict[str, int | str]]]] = {}
    for context, groups in contexts.items():
        context_specs[context] = {"Control": [], "DS": []}
        for group_name, group_samples in groups.items():
            for sample in group_samples:
                spec = sample_specs.get(sample["sample_key"])
                if spec:
                    context_specs[context][group_name].append(spec)
    return context_specs


def build_detail_row(
    dataset: str,
    context: str,
    sig: dict[str, str],
    control_mean: float,
    ds_mean: float,
    control_n: int,
    ds_n: int,
) -> dict[str, str]:
    ds_delta = ds_mean - control_mean
    dnmt3l_delta = safe_float(sig["dnmt3l_delta_beta"])
    expected_sign = 1 if sig["dnmt3l_direction"] == "hyper" else -1
    signed_ds_delta = ds_delta * expected_sign
    return {
        "dataset": dataset,
        "context": context,
        "probe_id": sig["probe_id"],
        "chrom": sig["chrom"],
        "start": sig["start"],
        "end": sig["end"],
        "dnmt3l_direction": sig["dnmt3l_direction"],
        "dnmt3l_delta_beta": fmt(dnmt3l_delta),
        "ds_control_mean_beta": fmt(control_mean),
        "ds_case_mean_beta": fmt(ds_mean),
        "ds_delta_beta": fmt(ds_delta),
        "signed_ds_delta_by_dnmt3l_direction": fmt(signed_ds_delta),
        "concordant": "yes" if signed_ds_delta > 0 else "no",
        "control_n_detected": str(control_n),
        "ds_n_detected": str(ds_n),
        "ucsc_refgene_name": sig["ucsc_refgene_name"],
        "ucsc_refgene_group": sig["ucsc_refgene_group"],
        "relation_to_cpg_island": sig["relation_to_cpg_island"],
        "regulatory_feature_group": sig["regulatory_feature_group"],
        "gencode_name": sig["gencode_name"],
        "gencode_group": sig["gencode_group"],
    }


def summarize_details(detail_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    by_context: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in detail_rows:
        by_context[(row["dataset"], row["context"])].append(row)

    summaries = []
    for (dataset, context), rows in sorted(by_context.items()):
        dnmt3l_deltas = [safe_float(row["dnmt3l_delta_beta"]) for row in rows]
        ds_deltas = [safe_float(row["ds_delta_beta"]) for row in rows]
        signed = [safe_float(row["signed_ds_delta_by_dnmt3l_direction"]) for row in rows]
        hyper_rows = [row for row in rows if row["dnmt3l_direction"] == "hyper"]
        hypo_rows = [row for row in rows if row["dnmt3l_direction"] == "hypo"]
        concordant = [value for value in signed if value > 0]
        discordant = [value for value in signed if value < 0]
        summaries.append(
            {
                "dataset": dataset,
                "context": context,
                "n_overlap_probes": str(len(rows)),
                "n_dnmt3l_hyper_probes": str(len(hyper_rows)),
                "n_dnmt3l_hypo_probes": str(len(hypo_rows)),
                "mean_dnmt3l_delta_beta": fmt(mean(dnmt3l_deltas)),
                "mean_ds_delta_beta": fmt(mean(ds_deltas)),
                "median_ds_delta_beta": fmt(median(ds_deltas)),
                "mean_signed_ds_delta_by_dnmt3l_direction": fmt(mean(signed)),
                "median_signed_ds_delta_by_dnmt3l_direction": fmt(median(signed)),
                "concordant_n": str(len(concordant)),
                "discordant_n": str(len(discordant)),
                "concordant_fraction": fmt(len(concordant) / len(rows)) if rows else "",
                "same_direction_ge_0.02": str(sum(1 for value in signed if value >= 0.02)),
                "same_direction_ge_0.05": str(sum(1 for value in signed if value >= 0.05)),
                "opposite_direction_ge_0.02": str(sum(1 for value in signed if value <= -0.02)),
                "opposite_direction_ge_0.05": str(sum(1 for value in signed if value <= -0.05)),
                "mean_abs_ds_delta_beta": fmt(mean([abs(value) for value in ds_deltas])),
            }
        )
    return summaries


def add_counts_to_summary(
    summary_rows: list[dict[str, str]],
    context_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    counts = {row["context"]: row for row in context_rows}
    out = []
    for row in summary_rows:
        count_row = counts.get(row["context"], {})
        merged = dict(row)
        merged["n_control"] = count_row.get("n_control", "")
        merged["n_ds"] = count_row.get("n_ds", "")
        out.append(merged)
    return out


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    signature = read_dnmt3l_signature(Path(args.dnmt3l_dmcs))

    all_summary = []
    all_detail = []
    all_contexts = []

    dataset_specs = [
        (
            "GSE73747",
            compare_signal_dataset,
            {
                "matrix_path": Path("data/geo/GSE73747/supplementary/GSE73747_unmethylated_methylated.txt.gz"),
                "samples": gse73747_samples(Path("data/geo/GSE73747/metadata/sample_metadata.tsv")),
                "detection_p": args.detection_p,
                "beta_offset": args.beta_offset,
                "compact_headers": False,
            },
        ),
        (
            "GSE63347",
            compare_signal_dataset,
            {
                "matrix_path": Path("data/geo/GSE63347/supplementary/GSE63347_datSignal.csv.gz"),
                "samples": gse63347_samples(Path("data/geo/GSE63347/metadata/sample_metadata.tsv")),
                "detection_p": args.detection_p,
                "beta_offset": args.beta_offset,
                "compact_headers": True,
            },
        ),
        (
            "GSE74486",
            compare_processed_dataset,
            {
                "matrix_path": Path(
                    "data/geo/GSE74486/supplementary/"
                    "GSE74486_DS_processed_norm_batch_corrected_filtered_matrix_2.txt.gz"
                ),
                "samples": gse74486_xml_samples(Path("data/geo/GSE74486/metadata/GSE74486_family.xml")),
                "detection_p": args.gse74486_detection_p,
            },
        ),
    ]

    for dataset, fn, kwargs in dataset_specs:
        if not kwargs["matrix_path"].exists():
            continue
        summary_rows, detail_rows, context_rows = fn(
            dataset=dataset,
            signature=signature,
            min_group_samples=args.min_group_samples,
            **kwargs,
        )
        all_summary.extend(add_counts_to_summary(summary_rows, context_rows))
        all_detail.extend(detail_rows)
        for row in context_rows:
            row = dict(row)
            row["dataset"] = dataset
            all_contexts.append(row)

    summary_fields = [
        "dataset",
        "context",
        "n_control",
        "n_ds",
        "n_overlap_probes",
        "n_dnmt3l_hyper_probes",
        "n_dnmt3l_hypo_probes",
        "mean_dnmt3l_delta_beta",
        "mean_ds_delta_beta",
        "median_ds_delta_beta",
        "mean_signed_ds_delta_by_dnmt3l_direction",
        "median_signed_ds_delta_by_dnmt3l_direction",
        "concordant_n",
        "discordant_n",
        "concordant_fraction",
        "same_direction_ge_0.02",
        "same_direction_ge_0.05",
        "opposite_direction_ge_0.02",
        "opposite_direction_ge_0.05",
        "mean_abs_ds_delta_beta",
    ]
    detail_fields = [
        "dataset",
        "context",
        "probe_id",
        "chrom",
        "start",
        "end",
        "dnmt3l_direction",
        "dnmt3l_delta_beta",
        "ds_control_mean_beta",
        "ds_case_mean_beta",
        "ds_delta_beta",
        "signed_ds_delta_by_dnmt3l_direction",
        "concordant",
        "control_n_detected",
        "ds_n_detected",
        "ucsc_refgene_name",
        "ucsc_refgene_group",
        "relation_to_cpg_island",
        "regulatory_feature_group",
        "gencode_name",
        "gencode_group",
    ]
    context_fields = [
        "dataset",
        "context",
        "n_control",
        "n_ds",
        "control_accessions",
        "ds_accessions",
        "example_tissues",
        "example_cell_types",
    ]

    write_tsv(outdir / "ds_dnmt3l_candidate_probe_concordance_summary.tsv", all_summary, summary_fields)
    write_tsv(outdir / "ds_dnmt3l_candidate_probe_values.tsv", all_detail, detail_fields)
    write_tsv(outdir / "ds_methylation_dataset_contexts.tsv", all_contexts, context_fields)

    print(f"Wrote {len(all_summary)} context summaries")
    print(f"Wrote {len(all_detail)} candidate-probe rows")
    print(f"Wrote {len(all_contexts)} eligible dataset contexts")


if __name__ == "__main__":
    main()
