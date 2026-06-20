#!/usr/bin/env python3
"""Compare DS/control target-gene expression across downloaded GEO datasets.

The broad cross-dataset view reports within-dataset, within-context case-control
effects instead of pooling values across platforms. Fetal bulk GSE301886 and
GSE305139 use precomputed edgeR count-based logFC values when available.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
from collections import defaultdict
from pathlib import Path


TARGET_ORDER = ["DNMT3L", "DNMT3A", "DNMT3B", "STAT1", "STAT3"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-genes", default="config/target_genes.tsv")
    parser.add_argument("--outdir", default="results/ds_expression_comparison")
    parser.add_argument(
        "--human-ds-bulk-count-de",
        default="results/human_ds_bulk_count_de/human_ds_bulk_count_de_all_contexts.tsv",
    )
    return parser.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="", errors="replace")
    return path.open(newline="", errors="replace")


def strip_version(value: str) -> str:
    return re.sub(r"\.\d+$", "", value.strip().strip('"'))


def parse_float(value: str) -> float:
    try:
        number = float(value)
    except ValueError:
        return math.nan
    return number if math.isfinite(number) else math.nan


def fmt(value: float) -> str:
    if value is None or not math.isfinite(value):
        return ""
    return f"{value:.6g}"


def mean(values: list[float]) -> float:
    clean = [value for value in values if math.isfinite(value)]
    return sum(clean) / len(clean) if clean else math.nan


def read_targets(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    symbol_by_ensembl: dict[str, str] = {}
    symbol_aliases: dict[str, str] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["species"] != "human":
                continue
            symbol = row["gene_symbol"]
            symbol_aliases[symbol.upper()] = symbol
            for alias in re.split(r"[|;,]\s*", row.get("aliases", "")):
                if alias:
                    symbol_aliases[alias.upper()] = symbol
            ensembl = strip_version(row.get("ensembl_hint", ""))
            if ensembl:
                symbol_by_ensembl[ensembl] = symbol
                symbol_aliases[ensembl.upper()] = symbol
    return symbol_by_ensembl, symbol_aliases


def parse_characteristics(value: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for part in value.split(";"):
        part = part.strip()
        if not part or "=" not in part:
            continue
        key, val = part.split("=", 1)
        result[key.strip().lower()] = val.strip()
    return result


def load_geo_metadata(accession: str) -> dict[str, dict[str, str]]:
    path = Path("data/geo") / accession / "metadata" / "sample_metadata.tsv"
    rows: dict[str, dict[str, str]] = {}
    if not path.exists():
        return rows
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            chars = parse_characteristics(row.get("characteristics", ""))
            merged = {**row, **chars}
            rows[row["accession"]] = merged
    return rows


def metadata_by_field(accession: str, field: str) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for row in load_geo_metadata(accession).values():
        key = row.get(field, "")
        if key:
            result[key] = row
    return result


def canonical_group(*values: str) -> str:
    text = " ".join(value or "" for value in values).lower()
    if any(token in text for token in ["down syndrome", "trisomy 21", "trisomic", "t21", " dsd", "dsd"]):
        return "DS"
    if re.search(r"\bds\b", text) and not re.search(r"\bds2u\b", text):
        return "DS"
    if any(token in text for token in ["control", "euploid", "disomic", "non-demented", "normal", "ctr"]):
        return "Control"
    if re.search(r"\bn_\d|\bn$", text):
        return "Control"
    return "Other"


def add_row(
    rows: list[dict[str, str]],
    dataset: str,
    context: str,
    sample_label: str,
    group: str,
    gene_symbol: str,
    value: float,
    value_type: str,
    source_file: str,
    library_size: float = math.nan,
    count: float = math.nan,
    notes: str = "",
) -> None:
    normalized = value
    log2_value = math.log2(value + 1.0) if math.isfinite(value) and value >= 0 else math.nan
    rows.append(
        {
            "dataset": dataset,
            "context": context,
            "sample_label": sample_label,
            "group": group,
            "gene_symbol": gene_symbol,
            "value_type": value_type,
            "value": fmt(value),
            "log2_value": fmt(log2_value),
            "library_size": fmt(library_size),
            "count": fmt(count),
            "source_file": source_file,
            "notes": notes,
        }
    )


def target_symbol(gene_id: str, symbol_by_ensembl: dict[str, str], symbol_aliases: dict[str, str]) -> str | None:
    gene = strip_version(gene_id)
    if gene in symbol_by_ensembl:
        return symbol_by_ensembl[gene]
    return symbol_aliases.get(gene.upper())


def detect_delimiter(path: Path) -> str:
    suffix = path.with_suffix("").suffix if path.suffix == ".gz" else path.suffix
    return "," if suffix == ".csv" else "\t"


def read_matrix_targets(
    path: Path,
    dataset: str,
    context: str,
    sample_info: dict[str, dict[str, str]],
    symbol_by_ensembl: dict[str, str],
    symbol_aliases: dict[str, str],
    gene_col: str | None = None,
    sample_start: int = 1,
    value_type: str = "raw_count_cpm",
    first_header_is_samples: bool = False,
    group_from_sample: bool = False,
    group_from_info=None,
    context_from_sample=None,
) -> list[dict[str, str]]:
    delimiter = detect_delimiter(path)
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        header = next(reader)
        if first_header_is_samples:
            sample_cols = header
            first_col_name = "gene_id"
            sample_start = 1
        else:
            first_col_name = header[0]
            sample_cols = header[sample_start:]

        target_rows: dict[str, list[str]] = {}
        library_sizes = {sample: 0.0 for sample in sample_cols}
        for row in reader:
            if not row:
                continue
            gene_id = row[0]
            symbol = target_symbol(gene_id, symbol_by_ensembl, symbol_aliases)
            values = row[sample_start : sample_start + len(sample_cols)] if not first_header_is_samples else row[1:]
            for sample, value in zip(sample_cols, values):
                number = parse_float(value)
                if math.isfinite(number):
                    library_sizes[sample] += number
            if symbol:
                target_rows[symbol] = row

    rows: list[dict[str, str]] = []
    for symbol, row in target_rows.items():
        values = row[sample_start : sample_start + len(sample_cols)] if not first_header_is_samples else row[1:]
        for sample, value in zip(sample_cols, values):
            raw = parse_float(value)
            if not math.isfinite(raw):
                continue
            info = sample_info.get(sample, {})
            if group_from_info:
                group = group_from_info(sample, info)
            elif group_from_sample:
                group = canonical_group(sample)
            else:
                group = canonical_group(
                    info.get("group", ""),
                    info.get("genotype", ""),
                    info.get("condition", ""),
                    info.get("diagnosis", ""),
                    info.get("genotype/variation", ""),
                    sample,
                )
            row_context = context_from_sample(sample, info, context) if context_from_sample else context
            if value_type == "raw_count_cpm":
                lib = library_sizes.get(sample, math.nan)
                normalized = raw / lib * 1_000_000 if lib and math.isfinite(lib) else math.nan
                count = raw
            else:
                normalized = raw
                count = math.nan
                lib = math.nan
            add_row(rows, dataset, row_context, sample, group, symbol, normalized, value_type, str(path), lib, count)
    return rows


def read_existing_ds_bulk() -> list[dict[str, str]]:
    path = Path("results/ds_fetal_expression/ds_target_sample_values.tsv")
    rows: list[dict[str, str]] = []
    if not path.exists():
        return rows
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("included_for_summary", "TRUE").upper() == "FALSE":
                continue
            context = row.get("region") or "fetal brain"
            add_row(
                rows,
                row["dataset"],
                context,
                row.get("sample_label", row.get("geo_sample_accession", "")),
                row["group"],
                row["gene_symbol"],
                parse_float(row["cpm"]),
                "raw_count_cpm",
                str(path),
                parse_float(row.get("library_size", "")),
                parse_float(row.get("count", "")),
            )
    return rows


def read_gse59630_array() -> list[dict[str, str]]:
    path = Path("results/ds_fetal_expression/GSE59630_target_sample_values.tsv")
    rows: list[dict[str, str]] = []
    if not path.exists():
        return rows
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            context = f"{row.get('life_stage', '')}|{row.get('analysis_region', row.get('region', ''))}"
            add_row(
                rows,
                "GSE59630",
                context,
                row.get("source_name", row.get("geo_sample_accession", "")),
                row["group"],
                row["gene_symbol"],
                parse_float(row["expression_value"]),
                "array_expression_native",
                str(path),
            )
    return rows


def read_gse305146_multiome() -> list[dict[str, str]]:
    path = Path("results/ds_fetal_expression/GSE305146_target_library_values.tsv")
    rows: list[dict[str, str]] = []
    if not path.exists():
        return rows
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            contexts = [
                "fetal_single_nucleus_multiome_all_pcw",
                f"fetal_single_nucleus_multiome|PCW{row.get('dev_stage_pcw', '')}",
            ]
            for context in contexts:
                add_row(
                    rows,
                    "GSE305146",
                    context,
                    row.get("library_prefix", row.get("sample_name", "")),
                    row["group"],
                    row["gene_symbol"],
                    parse_float(row["target_cpm"]),
                    "raw_count_cpm",
                    str(path),
                    parse_float(row.get("total_gene_expression_umi", "")),
                    parse_float(row.get("target_umi", "")),
                    f"detected_fraction={row.get('detected_fraction', '')}",
                )
    return rows


def read_star_sample_files(
    dataset: str,
    pattern: str,
    context_from_meta,
    symbol_by_ensembl: dict[str, str],
    symbol_aliases: dict[str, str],
) -> list[dict[str, str]]:
    metadata = load_geo_metadata(dataset)
    rows: list[dict[str, str]] = []
    for path in sorted(Path("data/geo").glob(pattern)):
        accession_match = re.search(r"(GSM\d+)", path.name)
        accession = accession_match.group(1) if accession_match else ""
        info = metadata.get(accession, {})
        counts: dict[str, float] = {}
        library_size = 0.0
        with open_text(path) as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                gene_id = parts[0]
                value = parse_float(parts[1])
                if not math.isfinite(value):
                    continue
                if not gene_id.startswith("N_"):
                    library_size += value
                symbol = target_symbol(gene_id, symbol_by_ensembl, symbol_aliases)
                if symbol:
                    counts[symbol] = value
        context = context_from_meta(info, path)
        group = canonical_group(info.get("genotype", ""), info.get("condition", ""), info.get("strain", ""), info.get("title", ""))
        for symbol, count in counts.items():
            cpm = count / library_size * 1_000_000 if library_size else math.nan
            add_row(rows, dataset, context, accession or path.stem, group, symbol, cpm, "raw_count_cpm", str(path), library_size, count)
    return rows


def read_htseq_sample_files(
    dataset: str,
    pattern: str,
    context_from_meta,
    symbol_by_ensembl: dict[str, str],
    symbol_aliases: dict[str, str],
) -> list[dict[str, str]]:
    metadata = load_geo_metadata(dataset)
    rows: list[dict[str, str]] = []
    for path in sorted(Path("data/geo").glob(pattern)):
        accession_match = re.search(r"(GSM\d+)", path.name)
        accession = accession_match.group(1) if accession_match else ""
        info = metadata.get(accession, {})
        counts: dict[str, float] = {}
        library_size = 0.0
        with open_text(path) as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                gene_id, value_text = parts[0], parts[1]
                value = parse_float(value_text)
                if not math.isfinite(value):
                    continue
                if not gene_id.startswith("__"):
                    library_size += value
                symbol = target_symbol(gene_id, symbol_by_ensembl, symbol_aliases)
                if symbol:
                    counts[symbol] = value
        context = context_from_meta(info, path)
        group = canonical_group(info.get("genotype", ""), info.get("condition", ""), info.get("strain", ""), info.get("title", ""))
        for symbol, count in counts.items():
            cpm = count / library_size * 1_000_000 if library_size else math.nan
            add_row(rows, dataset, context, accession or path.stem, group, symbol, cpm, "raw_count_cpm", str(path), library_size, count)
    return rows


def sample_info_from_metadata(accession: str, field: str = "sample id") -> dict[str, dict[str, str]]:
    return metadata_by_field(accession, field)


def sample_info_from_titles(accession: str, sample_id_pattern: str) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for row in load_geo_metadata(accession).values():
        match = re.search(sample_id_pattern, row.get("title", ""))
        if match:
            result[match.group(1)] = row
    return result


def compact_label(value: str, fallback: str = "unknown") -> str:
    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", value.strip()).strip("_")
    return cleaned or fallback


def context_with_treatment(sample: str, info: dict[str, str], base: str) -> str:
    treatment = compact_label(info.get("treatment", ""), "untagged_treatment")
    return f"{base}|{treatment}"


def context_with_cell_type_and_protocol(sample: str, info: dict[str, str], base: str) -> str:
    cell_type = compact_label(info.get("cell type", ""), "unknown_cell_type")
    protocol = compact_label(info.get("protocol", ""), "unknown_protocol")
    return f"{cell_type}|protocol_{protocol}"


def context_gse144857(sample: str, info: dict[str, str], base: str) -> str:
    lower = sample.lower()
    if "neuron" in lower:
        return "iPSC_derived_neuron"
    if "stemcell" in lower:
        return "iPSC_stemcell"
    return base


def group_gse222365_expanded(sample: str, info: dict[str, str]) -> str:
    if re.search(r"_par[A-Za-z]*_", sample):
        return "DS"
    if re.search(r"_dis[A-Za-z]*_", sample):
        return "Control"
    return canonical_group(info.get("genotype", ""), sample)


def read_gse52249(symbol_by_ensembl: dict[str, str], symbol_aliases: dict[str, str]) -> list[dict[str, str]]:
    path = Path("data/geo/GSE52249/supplementary/GSE52249_RPKM_expression_profile.txt.gz")
    rows: list[dict[str, str]] = []
    if not path.exists():
        return rows
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            symbol = target_symbol(row.get("gene_name", ""), symbol_by_ensembl, symbol_aliases)
            if not symbol:
                continue
            for key, value_text in row.items():
                if not key.startswith("rpkm_"):
                    continue
                group = "DS" if "_DS_" in key else "Control"
                add_row(rows, "GSE52249", "twin_iPSC", key, group, symbol, parse_float(value_text), "rpkm", str(path))
    return rows


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def de_context_key(dataset: str, context: str) -> tuple[str, str]:
    aliases = {
        ("GSE305139", "fetal brain"): ("GSE305139", "fetal_brain_final"),
    }
    return aliases.get((dataset, context), (dataset, context))


def read_count_de_target_effects(path: Path, target_genes: set[str]) -> dict[tuple[str, str, str], dict[str, str]]:
    effects: dict[tuple[str, str, str], dict[str, str]] = {}
    if not path.exists():
        return effects
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            gene = row.get("gene_symbol", "")
            if gene not in target_genes:
                continue
            key = (row.get("dataset", ""), row.get("context", ""), gene)
            current = effects.get(key)
            if current is None or parse_float(row.get("mean_ds_cpm", "")) + parse_float(row.get("mean_control_cpm", "")) > parse_float(current.get("mean_ds_cpm", "")) + parse_float(current.get("mean_control_cpm", "")):
                effects[key] = row
    return effects


def summarize_effects(
    rows: list[dict[str, str]],
    count_de: dict[tuple[str, str, str], dict[str, str]] | None = None,
) -> list[dict[str, str]]:
    grouped: dict[tuple[str, str, str], dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    value_types: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    for row in rows:
        if row["group"] not in {"DS", "Control"}:
            continue
        key = (row["dataset"], row["context"], row["gene_symbol"])
        value = parse_float(row["log2_value"])
        if math.isfinite(value):
            grouped[key][row["group"]].append(value)
            value_types[key].add(row["value_type"])

    out: list[dict[str, str]] = []
    for key, values in sorted(grouped.items()):
        dataset, context, gene = key
        ds = values.get("DS", [])
        control = values.get("Control", [])
        if not ds or not control:
            continue
        ds_mean = mean(ds)
        control_mean = mean(control)
        descriptive_delta = ds_mean - control_mean
        de_key = (*de_context_key(dataset, context), gene)
        de_row = (count_de or {}).get(de_key, {})
        de_effect = parse_float(de_row.get("effect_log2fc_for_claim", ""))
        effect_for_claim = de_effect if math.isfinite(de_effect) else descriptive_delta
        out.append(
            {
                "dataset": dataset,
                "context": context,
                "gene_symbol": gene,
                "n_ds": str(len(ds)),
                "n_control": str(len(control)),
                "mean_ds_log2": fmt(ds_mean),
                "mean_control_log2": fmt(control_mean),
                "delta_ds_minus_control_log2": fmt(descriptive_delta),
                "effect_log2fc_for_claim": fmt(effect_for_claim),
                "effect_source": de_row.get("effect_source", "mean_log2_value_descriptive"),
                "edgeR_logFC": fmt(parse_float(de_row.get("logFC", ""))),
                "edgeR_PValue": fmt(parse_float(de_row.get("PValue", ""))),
                "edgeR_FDR": fmt(parse_float(de_row.get("FDR", ""))),
                "tested_by_edgeR": de_row.get("tested_by_edgeR", ""),
                "value_type": "|".join(sorted(value_types[key])),
            }
        )
    return out


def plot_effects(effects: list[dict[str, str]], figure_dir: Path) -> list[Path]:
    if not effects:
        return []
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return []

    figure_dir.mkdir(parents=True, exist_ok=True)
    contexts = sorted({f"{row['dataset']}|{row['context']}" for row in effects})
    genes = TARGET_ORDER
    values: dict[tuple[str, str], float] = {}
    for row in effects:
        context = f"{row['dataset']}|{row['context']}"
        effect = parse_float(row.get("effect_log2fc_for_claim", ""))
        values[(context, row["gene_symbol"])] = effect if math.isfinite(effect) else parse_float(row["delta_ds_minus_control_log2"])

    matrix = []
    finite_values = []
    for context in contexts:
        row_values = []
        for gene in genes:
            value = values.get((context, gene), math.nan)
            row_values.append(value)
            if math.isfinite(value):
                finite_values.append(value)
        matrix.append(row_values)
    if not finite_values:
        return []

    vmax = max(0.25, min(3.0, max(abs(value) for value in finite_values)))
    height = max(4.5, min(18.0, 0.32 * len(contexts) + 1.8))
    fig, ax = plt.subplots(figsize=(8.5, height))
    image = ax.imshow(matrix, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=35, ha="right")
    ax.set_yticks(range(len(contexts)))
    ax.set_yticklabels(contexts, fontsize=6)
    ax.set_title("DS-control expression difference")
    ax.set_xlabel("Target gene")
    ax.set_ylabel("Dataset and context")
    for y, context in enumerate(contexts):
        for x, gene in enumerate(genes):
            value = values.get((context, gene), math.nan)
            if math.isfinite(value):
                ax.text(x, y, f"{value:+.2f}", ha="center", va="center", fontsize=5)
    colorbar = fig.colorbar(image, ax=ax, shrink=0.75)
    colorbar.set_label("Mean log2(value+1), DS - control")
    fig.tight_layout()
    heatmap_png = figure_dir / "ds_target_gene_effect_heatmap.png"
    heatmap_pdf = figure_dir / "ds_target_gene_effect_heatmap.pdf"
    fig.savefig(heatmap_png, dpi=220)
    fig.savefig(heatmap_pdf)
    plt.close(fig)

    dnmt3l = [row for row in effects if row["gene_symbol"] == "DNMT3L"]
    outputs = [heatmap_png, heatmap_pdf]
    if dnmt3l:
        labels = [f"{row['dataset']}|{row['context']}" for row in dnmt3l]
        deltas = [
            parse_float(row.get("effect_log2fc_for_claim", ""))
            if math.isfinite(parse_float(row.get("effect_log2fc_for_claim", "")))
            else parse_float(row["delta_ds_minus_control_log2"])
            for row in dnmt3l
        ]
        valid = [(label, delta) for label, delta in zip(labels, deltas) if math.isfinite(delta)]
        if valid:
            valid.sort(key=lambda item: item[1])
            labels = [item[0] for item in valid]
            deltas = [item[1] for item in valid]
            fig_height = max(3.5, min(12.0, 0.36 * len(labels) + 1.5))
            fig, ax = plt.subplots(figsize=(8.5, fig_height))
            colors = ["#2f6f9f" if value >= 0 else "#b65d3b" for value in deltas]
            ax.barh(range(len(labels)), deltas, color=colors)
            ax.axvline(0, color="#333333", linewidth=0.8)
            ax.set_yticks(range(len(labels)))
            ax.set_yticklabels(labels, fontsize=7)
            ax.set_xlabel("Mean log2(value+1), DS - control")
            ax.set_title("DNMT3L DS-control expression difference")
            fig.tight_layout()
            dnmt3l_png = figure_dir / "dnmt3l_effects_by_dataset_context.png"
            dnmt3l_pdf = figure_dir / "dnmt3l_effects_by_dataset_context.pdf"
            fig.savefig(dnmt3l_png, dpi=220)
            fig.savefig(dnmt3l_pdf)
            plt.close(fig)
            outputs.extend([dnmt3l_png, dnmt3l_pdf])
    return outputs


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    symbol_by_ensembl, symbol_aliases = read_targets(Path(args.target_genes))

    rows: list[dict[str, str]] = []
    rows.extend(read_existing_ds_bulk())
    rows.extend(read_gse59630_array())
    rows.extend(read_gse305146_multiome())

    rows.extend(
        read_star_sample_files(
            "GSE247990",
            "GSE247990/supplementary/*ReadsPerGene.out.tab.gz",
            lambda info, _path: f"{info.get('time', '')}|{info.get('source_name', '')}",
            symbol_by_ensembl,
            symbol_aliases,
        )
    )
    rows.extend(
        read_htseq_sample_files(
            "GSE185192",
            "GSE185192/supplementary/*HTSeq_RawCounts.txt.gz",
            lambda info, _path: f"{info.get('cell type', '')}|{info.get('treatment', '')}",
            symbol_by_ensembl,
            symbol_aliases,
        )
    )

    rows.extend(
        read_matrix_targets(
            Path("data/geo/GSE305132/supplementary/GSE305132_counts_comb_bulk_in_vitro_ASO_exp.tsv.gz"),
            "GSE305132",
            "iPSC_NPC_ASO_primary",
            sample_info_from_metadata("GSE305132", "sample id"),
            symbol_by_ensembl,
            symbol_aliases,
            sample_start=2,
            context_from_sample=context_with_treatment,
        )
    )
    rows.extend(
        read_matrix_targets(
            Path("data/geo/GSE305132/supplementary/GSE305132_counts_comb_bulk_in_vitro_ASO_exp_E22.tsv.gz"),
            "GSE305132",
            "iPSC_NPC_ASO_E22",
            sample_info_from_metadata("GSE305132", "sample id"),
            symbol_by_ensembl,
            symbol_aliases,
            sample_start=2,
            context_from_sample=context_with_treatment,
        )
    )
    for path in [
        Path("data/geo/GSE305135/supplementary/GSE305135_Bulk_Singleron_in_vitro_basic_char_counts_comb.csv.gz"),
        Path("data/geo/GSE305135/supplementary/GSE305135_Bulk_Singleron_in_vitro_basic_char_v06_add_samples_counts_comb.csv.gz"),
    ]:
        rows.extend(
            read_matrix_targets(
                path,
                "GSE305135",
                "iPSC_neural_cells",
                sample_info_from_metadata("GSE305135", "sample id"),
                symbol_by_ensembl,
                symbol_aliases,
                sample_start=1,
                context_from_sample=context_with_cell_type_and_protocol,
            )
        )

    rows.extend(
        read_matrix_targets(
            Path("data/geo/GSE221178/supplementary/GSE221178_raw_counts_human_cortex_matrix_coding.txt.gz"),
            "GSE221178",
            "adult_cortex_totalRNA",
            {},
            symbol_by_ensembl,
            symbol_aliases,
            sample_start=1,
            group_from_sample=True,
            context_from_sample=context_gse144857,
        )
    )
    rows.extend(
        read_matrix_targets(
            Path("data/geo/GSE221178/supplementary/GSE221178_raw_counts_human_hippo_matrix_coding.txt.gz"),
            "GSE221178",
            "adult_hippocampus_totalRNA",
            {},
            symbol_by_ensembl,
            symbol_aliases,
            sample_start=1,
            group_from_sample=True,
        )
    )
    rows.extend(
        read_matrix_targets(
            Path("data/geo/GSE251939/supplementary/GSE251939_count_matrix_l3.txt.gz"),
            "GSE251939",
            "adult_frontal_cortex_pyramidal_neurons_layer3",
            sample_info_from_titles("GSE251939", r"(h[0-9A-Z]+)_Layer3"),
            symbol_by_ensembl,
            symbol_aliases,
            first_header_is_samples=True,
        )
    )
    rows.extend(
        read_matrix_targets(
            Path("data/geo/GSE251939/supplementary/GSE251939_count_matrix_l5.txt.gz"),
            "GSE251939",
            "adult_frontal_cortex_pyramidal_neurons_layer5",
            sample_info_from_titles("GSE251939", r"(h[0-9A-Z]+)_Layer5"),
            symbol_by_ensembl,
            symbol_aliases,
            first_header_is_samples=True,
        )
    )
    rows.extend(
        read_matrix_targets(
            Path("data/geo/GSE144857/supplementary/GSE144857_counts_raw-2019-12-06.csv.gz"),
            "GSE144857",
            "iPSC_neuron_and_stemcell",
            {},
            symbol_by_ensembl,
            symbol_aliases,
            sample_start=1,
            group_from_sample=True,
        )
    )
    rows.extend(
        read_matrix_targets(
            Path("data/geo/GSE222365/supplementary/GSE222365_counts_pilot_study.txt.gz"),
            "GSE222365",
            "pilot_cerebral_organoid",
            sample_info_from_titles("GSE222365", r"(organoid_\d+)"),
            symbol_by_ensembl,
            symbol_aliases,
            sample_start=1,
        )
    )
    rows.extend(
        read_matrix_targets(
            Path("data/geo/GSE222365/supplementary/GSE222365_counts_expanded_study.txt.gz"),
            "GSE222365",
            "expanded_cerebral_organoid",
            {},
            symbol_by_ensembl,
            symbol_aliases,
            sample_start=1,
            group_from_info=group_gse222365_expanded,
        )
    )
    rows.extend(
        read_matrix_targets(
            Path("data/geo/GSE318560/supplementary/GSE318560_raw_counts_matrix.txt.gz"),
            "GSE318560",
            "adult_hippocampus_dementia_series",
            sample_info_from_titles("GSE318560", r"(HC_\d+)"),
            symbol_by_ensembl,
            symbol_aliases,
            first_header_is_samples=True,
        )
    )
    rows.extend(read_gse52249(symbol_by_ensembl, symbol_aliases))

    sample_fields = [
        "dataset",
        "context",
        "sample_label",
        "group",
        "gene_symbol",
        "value_type",
        "value",
        "log2_value",
        "library_size",
        "count",
        "source_file",
        "notes",
    ]
    write_tsv(outdir / "ds_target_gene_sample_values.tsv", rows, sample_fields)

    count_de = read_count_de_target_effects(Path(args.human_ds_bulk_count_de), set(TARGET_ORDER))
    effects = summarize_effects(rows, count_de)
    effect_fields = [
        "dataset",
        "context",
        "gene_symbol",
        "n_ds",
        "n_control",
        "mean_ds_log2",
        "mean_control_log2",
        "delta_ds_minus_control_log2",
        "effect_log2fc_for_claim",
        "effect_source",
        "edgeR_logFC",
        "edgeR_PValue",
        "edgeR_FDR",
        "tested_by_edgeR",
        "value_type",
    ]
    write_tsv(outdir / "ds_target_gene_effect_summary.tsv", effects, effect_fields)
    figures = plot_effects(effects, Path("figures/ds_expression_comparison"))
    figure_note = f" and {len(figures)} figure files" if figures else ""
    print(f"Wrote {len(rows)} sample-level target rows, {len(effects)} DS-control effects{figure_note} to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
