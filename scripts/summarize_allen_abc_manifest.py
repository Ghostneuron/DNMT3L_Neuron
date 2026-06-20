#!/usr/bin/env python3
"""Summarize Allen Brain Cell Atlas manifest files for this project.

The ABC matrices are often multi-GB H5AD files. This script inventories the
release manifest, highlights project-relevant datasets, and checks whether the
target genes are present in small gene metadata tables that have been downloaded.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Iterable


RELEVANT_DATASETS = {
    "WHB-10Xv3": {
        "label": "Whole human brain 10Xv3 snRNA-seq",
        "species_scope": "human",
        "project_use": "adult human cell-type and region reference",
    },
    "WHB-taxonomy": {
        "label": "Whole human brain taxonomy",
        "species_scope": "human",
        "project_use": "adult human cell-type annotation reference",
    },
    "WMB-10X": {
        "label": "Whole mouse brain 10X metadata",
        "species_scope": "mouse",
        "project_use": "adult mouse metadata and cell-type reference",
    },
    "WMB-10Xv2": {
        "label": "Whole mouse brain 10Xv2 scRNA-seq",
        "species_scope": "mouse",
        "project_use": "adult mouse region-partitioned expression matrices",
    },
    "WMB-10Xv3": {
        "label": "Whole mouse brain 10Xv3 scRNA-seq",
        "species_scope": "mouse",
        "project_use": "adult mouse region-partitioned expression matrices",
    },
    "WMB-taxonomy": {
        "label": "Whole mouse brain taxonomy",
        "species_scope": "mouse",
        "project_use": "adult mouse cell-type annotation reference",
    },
    "Developing-Mouse-Vis-Cortex-10X": {
        "label": "Developing mouse visual cortex 10X scRNA-seq",
        "species_scope": "mouse",
        "project_use": "developmental mouse cortical cell-type reference",
    },
    "Developing-Mouse-Vis-Cortex-taxonomy": {
        "label": "Developing mouse visual cortex taxonomy",
        "species_scope": "mouse",
        "project_use": "developmental mouse cortical annotation reference",
    },
    "Zeng-Aging-Mouse-10Xv3": {
        "label": "Aging mouse whole-brain 10Xv3 scRNA-seq",
        "species_scope": "mouse",
        "project_use": "adult/aging mouse cell-type reference",
    },
    "HMBA-10xMultiome-BG": {
        "label": "Human-mammalian basal ganglia 10x multiome",
        "species_scope": "human/mammalian",
        "project_use": "cross-species basal ganglia expression and chromatin context",
    },
    "HMBA-10xMultiome-BG-Aligned": {
        "label": "Human-mammalian basal ganglia aligned 10x multiome",
        "species_scope": "human/mammalian",
        "project_use": "aligned cross-species basal ganglia target-gene comparison",
    },
}

GENE_TABLES = [
    {
        "atlas": "WHB-10Xv3",
        "species": "human",
        "path": "WHB_10Xv3_gene.csv",
    },
    {
        "atlas": "WMB-10X",
        "species": "mouse",
        "path": "WMB_10X_gene.csv",
    },
    {
        "atlas": "Developing-Mouse-Vis-Cortex-10X",
        "species": "mouse",
        "path": "Developing_Mouse_Vis_Cortex_10X_gene.csv",
    },
    {
        "atlas": "HMBA-10xMultiome-BG-Aligned",
        "species": "human",
        "path": "HMBA_10xMultiome_BG_Aligned_gene.csv",
    },
]

INDEX_FIELDS = [
    "dataset",
    "dataset_label",
    "species_scope",
    "project_use",
    "manifest_path",
    "section",
    "file_role",
    "version",
    "relative_path",
    "url",
    "size_bytes",
    "size_gb",
    "download_strategy",
]

PRESENCE_FIELDS = [
    "atlas",
    "atlas_label",
    "species",
    "target_gene_symbol",
    "present",
    "matched_gene_symbol",
    "gene_identifier",
    "matched_by",
    "gene_table",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        default="data/allen/abc_atlas/manifest_20260415.json",
        help="Downloaded ABC manifest JSON.",
    )
    parser.add_argument(
        "--target-genes",
        default="config/target_genes.tsv",
        help="Project target gene table.",
    )
    parser.add_argument(
        "--gene-metadata-dir",
        default="data/allen/abc_atlas/gene_metadata",
        help="Directory containing selected ABC gene.csv files.",
    )
    parser.add_argument(
        "--outdir",
        default="results/allen_abc",
        help="Output directory for manifest summaries.",
    )
    return parser.parse_args()


def iter_file_entries(
    dataset: str,
    node: Any,
    path_parts: list[str],
) -> Iterable[dict[str, str]]:
    if isinstance(node, dict):
        if "url" in node and "relative_path" in node:
            yield {
                "dataset": dataset,
                "manifest_path": "/".join(path_parts),
                "version": str(node.get("version", "")),
                "relative_path": str(node.get("relative_path", "")),
                "url": str(node.get("url", "")),
                "size_bytes": str(node.get("size", "")),
            }
            return
        for key, value in node.items():
            yield from iter_file_entries(dataset, value, path_parts + [str(key)])
    elif isinstance(node, list):
        for idx, value in enumerate(node):
            yield from iter_file_entries(dataset, value, path_parts + [str(idx)])


def file_role(relative_path: str, manifest_path: str) -> str:
    suffix = Path(relative_path).suffix.lower()
    if suffix == ".h5ad":
        return "expression_matrix_h5ad"
    if suffix == ".h5":
        return "precomputed_stats_h5"
    if relative_path.endswith("gene.csv"):
        return "gene_metadata_csv"
    if relative_path.endswith("cell_metadata.csv"):
        return "cell_metadata_csv"
    if "taxonomy" in relative_path or "cluster" in relative_path:
        return "taxonomy_metadata_csv"
    if suffix == ".csv":
        return "metadata_csv"
    return manifest_path.split("/")[-1] if manifest_path else "file"


def download_strategy(role: str, size_bytes: int) -> str:
    if role == "gene_metadata_csv":
        return "download_now_small_gene_table"
    if role in {"taxonomy_metadata_csv", "metadata_csv"} and size_bytes <= 50_000_000:
        return "download_as_needed_small_metadata"
    if role == "cell_metadata_csv" and size_bytes <= 150_000_000:
        return "download_as_needed_medium_cell_metadata"
    if role == "expression_matrix_h5ad" and size_bytes <= 2_000_000_000:
        return "candidate_target_gene_matrix_download"
    if role == "expression_matrix_h5ad":
        return "defer_large_matrix_until_needed"
    if size_bytes >= 1_000_000_000:
        return "defer_large_file_until_needed"
    return "download_as_needed"


def load_targets(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_gene_table(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def summarize_manifest(manifest_path: Path) -> list[dict[str, str]]:
    manifest = json.loads(manifest_path.read_text())
    file_listing = manifest.get("file_listing", {})
    rows: list[dict[str, str]] = []
    for dataset, node in file_listing.items():
        dataset_meta = RELEVANT_DATASETS.get(dataset)
        for entry in iter_file_entries(dataset, node, []):
            relative_path = entry["relative_path"]
            role = file_role(relative_path, entry["manifest_path"])
            try:
                size_bytes = int(entry["size_bytes"])
            except ValueError:
                size_bytes = 0
            rows.append(
                {
                    "dataset": dataset,
                    "dataset_label": dataset_meta["label"] if dataset_meta else "",
                    "species_scope": dataset_meta["species_scope"] if dataset_meta else "",
                    "project_use": dataset_meta["project_use"] if dataset_meta else "",
                    "manifest_path": entry["manifest_path"],
                    "section": entry["manifest_path"].split("/")[0] if entry["manifest_path"] else "",
                    "file_role": role,
                    "version": entry["version"],
                    "relative_path": relative_path,
                    "url": entry["url"],
                    "size_bytes": str(size_bytes),
                    "size_gb": f"{size_bytes / 1_000_000_000:.3f}",
                    "download_strategy": download_strategy(role, size_bytes),
                }
            )
    return rows


def summarize_target_presence(
    target_rows: list[dict[str, str]],
    gene_metadata_dir: Path,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for table_info in GENE_TABLES:
        atlas = table_info["atlas"]
        species = table_info["species"]
        table_path = gene_metadata_dir / table_info["path"]
        genes = load_gene_table(table_path)
        by_symbol = {
            row.get("gene_symbol", "").upper(): row
            for row in genes
            if row.get("gene_symbol")
        }
        by_identifier = {
            row.get("gene_identifier", ""): row
            for row in genes
            if row.get("gene_identifier")
        }
        for target in target_rows:
            if target["species"].lower() != species:
                continue
            symbol = target["gene_symbol"]
            ensembl_hint = target.get("ensembl_hint", "")
            match = by_symbol.get(symbol.upper())
            matched_by = "symbol" if match else ""
            if not match and ensembl_hint:
                match = by_identifier.get(ensembl_hint)
                matched_by = "ensembl_hint" if match else ""
            rows.append(
                {
                    "atlas": atlas,
                    "atlas_label": RELEVANT_DATASETS[atlas]["label"],
                    "species": species,
                    "target_gene_symbol": symbol,
                    "present": "yes" if match else "no",
                    "matched_gene_symbol": match.get("gene_symbol", "") if match else "",
                    "gene_identifier": match.get("gene_identifier", "") if match else "",
                    "matched_by": matched_by,
                    "gene_table": str(table_path),
                }
            )
    return rows


def write_tsv(rows: list[dict[str, str]], path: Path, fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    manifest_path = Path(args.manifest)
    outdir = Path(args.outdir)
    target_rows = load_targets(Path(args.target_genes))

    index_rows = summarize_manifest(manifest_path)
    relevant_rows = [
        row for row in index_rows if row["dataset"] in RELEVANT_DATASETS
    ]
    presence_rows = summarize_target_presence(target_rows, Path(args.gene_metadata_dir))

    write_tsv(
        index_rows,
        Path("data/allen/abc_atlas/manifest_20260415_file_index.tsv"),
        INDEX_FIELDS,
    )
    write_tsv(
        relevant_rows,
        outdir / "allen_abc_relevant_file_index.tsv",
        INDEX_FIELDS,
    )
    write_tsv(
        presence_rows,
        outdir / "allen_abc_target_gene_presence.tsv",
        PRESENCE_FIELDS,
    )

    h5ad_rows = [
        row
        for row in relevant_rows
        if row["file_role"] == "expression_matrix_h5ad"
    ]
    print(f"Indexed {len(index_rows)} ABC files.")
    print(f"Project-relevant files: {len(relevant_rows)}.")
    print(f"Project-relevant H5AD matrices: {len(h5ad_rows)}.")
    for row in h5ad_rows[:12]:
        print(
            f"- {row['dataset']}: {Path(row['relative_path']).name} "
            f"({row['size_gb']} GB; {row['download_strategy']})"
        )
    missing = [
        row
        for row in presence_rows
        if row["present"] != "yes"
    ]
    print(
        "Target-gene presence checks: "
        f"{len(presence_rows) - len(missing)}/{len(presence_rows)} present."
    )
    if missing:
        print("Missing target genes:")
        for row in missing:
            print(f"- {row['atlas']} {row['target_gene_symbol']}")


if __name__ == "__main__":
    main()
