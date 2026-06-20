#!/usr/bin/env python3
"""Download Allen Institute ISH metadata and quantified target-gene expression.

This intentionally avoids downloading full-resolution image stacks. It keeps
experiment metadata, section image IDs/URLs, and StructureUnionize rows
containing expression density, intensity, and energy by anatomical structure.
"""

from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


API_BASE = "https://api.brain-map.org"
QUERY_URL = f"{API_BASE}/api/v2/data/query.json"

PRODUCTS = {
    1: {
        "name": "Mouse Brain adult ISH",
        "species": "mouse",
        "portal": "https://mouse.brain-map.org",
    },
    3: {
        "name": "Developing Mouse Brain ISH",
        "species": "mouse",
        "portal": "https://developingmouse.brain-map.org",
    },
    22: {
        "name": "Developing Human Brain ISH",
        "species": "human",
        "portal": "https://brainspan.org",
    },
    9: {
        "name": "Human Brain ISH Cortex Study",
        "species": "human",
        "portal": "https://human.brain-map.org",
    },
    10: {
        "name": "Human Brain ISH SubCortex Study",
        "species": "human",
        "portal": "https://human.brain-map.org",
    },
}


AVAILABILITY_FIELDS = [
    "product_id",
    "product_name",
    "species",
    "gene_symbol",
    "total_experiments",
    "delegate_experiments",
]

DATASET_FIELDS = [
    "product_id",
    "product_name",
    "species",
    "gene_symbol",
    "gene_id",
    "gene_name",
    "entrez_id",
    "section_data_set_id",
    "delegate",
    "expression",
    "failed",
    "plane_of_section",
    "reference_space_id",
    "section_thickness",
    "storage_directory",
    "specimen_id",
    "specimen_name",
    "donor_id",
    "donor_name",
    "sex",
    "strain",
    "age_name",
    "age_days",
    "age_embryonic",
    "age_tags",
    "section_images_api_url",
    "grid_data_download_url",
    "experiment_url",
]

IMAGE_FIELDS = [
    "product_id",
    "product_name",
    "species",
    "gene_symbol",
    "section_data_set_id",
    "section_image_id",
    "section_number",
    "image_type",
    "resolution",
    "width",
    "height",
    "image_width",
    "image_height",
    "x",
    "y",
    "image_download_url_downsample6",
    "expression_mask_url_downsample6",
]

UNIONIZE_FIELDS = [
    "product_id",
    "product_name",
    "species",
    "gene_symbol",
    "section_data_set_id",
    "delegate",
    "age_name",
    "age_days",
    "age_embryonic",
    "plane_of_section",
    "structure_id",
    "structure_acronym",
    "structure_name",
    "structure_graph_id",
    "structure_graph_order",
    "structure_depth",
    "expression_density",
    "expression_energy",
    "sum_expressing_pixel_intensity",
    "sum_expressing_pixels",
    "sum_pixel_intensity",
    "sum_pixels",
    "voxel_energy_mean",
    "voxel_energy_cv",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--target-genes",
        default="config/target_genes.tsv",
        help="Target gene table with human and mouse symbols.",
    )
    parser.add_argument(
        "--outdir",
        default="data/allen/ish",
        help="Output directory for compact Allen ISH files.",
    )
    parser.add_argument(
        "--products",
        default="1,3,22,9,10",
        help="Comma-separated Allen product IDs to query.",
    )
    parser.add_argument("--page-size", type=int, default=500)
    parser.add_argument("--timeout", type=float, default=60)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--sleep", type=float, default=0.05)
    parser.add_argument(
        "--delegate-only-unionizes",
        action="store_true",
        help="Only fetch StructureUnionize rows for delegate experiments.",
    )
    return parser.parse_args()


def fetch_url(url: str, timeout: float, retries: int) -> bytes:
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            request = Request(url, headers={"User-Agent": "DNMT3L_Neuron/0.1"})
            with urlopen(request, timeout=timeout) as response:
                return response.read()
        except (HTTPError, URLError, TimeoutError) as exc:
            last_error = exc
            if attempt == retries:
                break
            time.sleep(min(2.0 * attempt, 10.0))
    raise RuntimeError(f"failed to fetch {url}: {last_error}")


def fetch_json(criteria: str, timeout: float, retries: int) -> dict[str, Any]:
    url = f"{QUERY_URL}?{urlencode({'criteria': criteria})}"
    payload = fetch_url(url, timeout=timeout, retries=retries)
    data = json.loads(payload.decode("utf-8"))
    if not data.get("success"):
        raise RuntimeError(f"Allen API query failed: {data.get('msg')}")
    return data


def fetch_paged(
    base_criteria: str,
    page_size: int,
    timeout: float,
    retries: int,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    start = 0
    total: int | None = None
    while total is None or start < total:
        criteria = (
            f"{base_criteria},"
            f"rma::options[start_row$eq{start}][num_rows$eq{page_size}]"
        )
        data = fetch_json(criteria, timeout=timeout, retries=retries)
        total = int(data["total_rows"])
        page = data.get("msg", [])
        rows.extend(page)
        if not page:
            break
        start += len(page)
    return rows


def count_experiments(
    product_id: int,
    gene_symbol: str,
    delegate_only: bool,
    timeout: float,
    retries: int,
) -> int:
    delegate_filter = "[delegate$eqtrue]" if delegate_only else ""
    criteria = (
        "model::SectionDataSet,"
        f"rma::criteria,[failed$eqfalse]{delegate_filter},"
        f"products[id$eq{product_id}],"
        f"genes[acronym$eq'{gene_symbol}'],"
        "rma::options[num_rows$eq0]"
    )
    return int(fetch_json(criteria, timeout=timeout, retries=retries)["total_rows"])


def load_targets(path: Path) -> dict[str, list[str]]:
    targets: dict[str, list[str]] = {"human": [], "mouse": []}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            species = row["species"].lower()
            if species in targets:
                targets[species].append(row["gene_symbol"])
    return targets


def flatten_dataset(product: dict[str, str], product_id: int, row: dict[str, Any]) -> dict[str, str]:
    gene = (row.get("genes") or [{}])[0]
    specimen = row.get("specimen") or {}
    donor = specimen.get("donor") or {}
    age = donor.get("age") or {}
    plane = row.get("plane_of_section") or {}
    dataset_id = str(row.get("id", ""))
    return {
        "product_id": str(product_id),
        "product_name": product["name"],
        "species": product["species"],
        "gene_symbol": str(gene.get("acronym", "")),
        "gene_id": str(gene.get("id", "")),
        "gene_name": str(gene.get("name", "")),
        "entrez_id": str(gene.get("entrez_id", "")),
        "section_data_set_id": dataset_id,
        "delegate": str(row.get("delegate", "")),
        "expression": str(row.get("expression", "")),
        "failed": str(row.get("failed", "")),
        "plane_of_section": str(plane.get("name", "")),
        "reference_space_id": str(row.get("reference_space_id", "")),
        "section_thickness": str(row.get("section_thickness", "")),
        "storage_directory": str(row.get("storage_directory", "")),
        "specimen_id": str(row.get("specimen_id", "")),
        "specimen_name": str(specimen.get("name", "")),
        "donor_id": str(specimen.get("donor_id", "")),
        "donor_name": str(donor.get("name", "")),
        "sex": str(donor.get("sex", "")),
        "strain": str(donor.get("strain", "")),
        "age_name": str(age.get("name", "")),
        "age_days": str(age.get("days", "")),
        "age_embryonic": str(age.get("embryonic", "")),
        "age_tags": str(age.get("tags", "")),
        "section_images_api_url": (
            f"{API_BASE}/api/v2/data/query.json?criteria="
            f"model::SectionImage,rma::criteria,data_set[id$eq{dataset_id}]"
        ),
        "grid_data_download_url": f"{API_BASE}/grid_data/download/{dataset_id}",
        "experiment_url": f"{product['portal']}/experiment/show/{dataset_id}",
    }


def flatten_image(dataset: dict[str, str], row: dict[str, Any]) -> dict[str, str]:
    image_id = str(row.get("id", ""))
    return {
        "product_id": dataset["product_id"],
        "product_name": dataset["product_name"],
        "species": dataset["species"],
        "gene_symbol": dataset["gene_symbol"],
        "section_data_set_id": dataset["section_data_set_id"],
        "section_image_id": image_id,
        "section_number": str(row.get("section_number", "")),
        "image_type": str(row.get("image_type", "")),
        "resolution": str(row.get("resolution", "")),
        "width": str(row.get("width", "")),
        "height": str(row.get("height", "")),
        "image_width": str(row.get("image_width", "")),
        "image_height": str(row.get("image_height", "")),
        "x": str(row.get("x", "")),
        "y": str(row.get("y", "")),
        "image_download_url_downsample6": (
            f"{API_BASE}/api/v2/image_download/{image_id}?downsample=6"
        ),
        "expression_mask_url_downsample6": (
            f"{API_BASE}/api/v2/image_download/{image_id}?downsample=6&view=expression"
        ),
    }


def flatten_unionize(dataset: dict[str, str], row: dict[str, Any]) -> dict[str, str]:
    structure = row.get("structure") or {}
    return {
        "product_id": dataset["product_id"],
        "product_name": dataset["product_name"],
        "species": dataset["species"],
        "gene_symbol": dataset["gene_symbol"],
        "section_data_set_id": dataset["section_data_set_id"],
        "delegate": dataset["delegate"],
        "age_name": dataset["age_name"],
        "age_days": dataset["age_days"],
        "age_embryonic": dataset["age_embryonic"],
        "plane_of_section": dataset["plane_of_section"],
        "structure_id": str(row.get("structure_id", "")),
        "structure_acronym": str(structure.get("acronym", "")),
        "structure_name": str(structure.get("name", "")),
        "structure_graph_id": str(structure.get("graph_id", "")),
        "structure_graph_order": str(structure.get("graph_order", "")),
        "structure_depth": str(structure.get("depth", "")),
        "expression_density": str(row.get("expression_density", "")),
        "expression_energy": str(row.get("expression_energy", "")),
        "sum_expressing_pixel_intensity": str(
            row.get("sum_expressing_pixel_intensity", "")
        ),
        "sum_expressing_pixels": str(row.get("sum_expressing_pixels", "")),
        "sum_pixel_intensity": str(row.get("sum_pixel_intensity", "")),
        "sum_pixels": str(row.get("sum_pixels", "")),
        "voxel_energy_mean": str(row.get("voxel_energy_mean", "")),
        "voxel_energy_cv": str(row.get("voxel_energy_cv", "")),
    }


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    targets = load_targets(Path(args.target_genes))
    product_ids = [int(pid) for pid in args.products.split(",") if pid.strip()]

    availability_rows: list[dict[str, str]] = []
    dataset_rows: list[dict[str, str]] = []
    image_rows: list[dict[str, str]] = []
    unionize_rows: list[dict[str, str]] = []

    for product_id in product_ids:
        product = PRODUCTS[product_id]
        for gene_symbol in targets[product["species"]]:
            total = count_experiments(
                product_id,
                gene_symbol,
                delegate_only=False,
                timeout=args.timeout,
                retries=args.retries,
            )
            delegate_total = count_experiments(
                product_id,
                gene_symbol,
                delegate_only=True,
                timeout=args.timeout,
                retries=args.retries,
            )
            availability_rows.append(
                {
                    "product_id": str(product_id),
                    "product_name": product["name"],
                    "species": product["species"],
                    "gene_symbol": gene_symbol,
                    "total_experiments": str(total),
                    "delegate_experiments": str(delegate_total),
                }
            )
            if total == 0:
                continue

            criteria = (
                "model::SectionDataSet,"
                "rma::criteria,[failed$eqfalse],"
                f"products[id$eq{product_id}],"
                f"genes[acronym$eq'{gene_symbol}'],"
                "rma::include,genes,specimen(donor(age)),plane_of_section"
            )
            raw_datasets = fetch_paged(
                criteria,
                page_size=args.page_size,
                timeout=args.timeout,
                retries=args.retries,
            )
            for raw_dataset in raw_datasets:
                dataset = flatten_dataset(product, product_id, raw_dataset)
                dataset_rows.append(dataset)

                image_criteria = (
                    "model::SectionImage,"
                    f"rma::criteria,data_set[id$eq{dataset['section_data_set_id']}]"
                )
                raw_images = fetch_paged(
                    image_criteria,
                    page_size=args.page_size,
                    timeout=args.timeout,
                    retries=args.retries,
                )
                image_rows.extend(flatten_image(dataset, image) for image in raw_images)

                if (
                    args.delegate_only_unionizes
                    and dataset["delegate"].lower() != "true"
                ):
                    continue
                unionize_criteria = (
                    "model::StructureUnionize,"
                    f"rma::criteria,section_data_set[id$eq{dataset['section_data_set_id']}],"
                    "rma::include,structure"
                )
                raw_unionizes = fetch_paged(
                    unionize_criteria,
                    page_size=args.page_size,
                    timeout=args.timeout,
                    retries=args.retries,
                )
                unionize_rows.extend(
                    flatten_unionize(dataset, unionize) for unionize in raw_unionizes
                )
                if args.sleep:
                    time.sleep(args.sleep)

    write_tsv(outdir / "allen_ish_target_availability.tsv", AVAILABILITY_FIELDS, availability_rows)
    write_tsv(outdir / "allen_ish_section_datasets.tsv", DATASET_FIELDS, dataset_rows)
    write_tsv(outdir / "allen_ish_section_images.tsv", IMAGE_FIELDS, image_rows)
    write_tsv(outdir / "allen_ish_structure_unionizes.tsv", UNIONIZE_FIELDS, unionize_rows)

    print(f"Wrote {len(availability_rows)} availability rows")
    print(f"Wrote {len(dataset_rows)} section dataset rows")
    print(f"Wrote {len(image_rows)} section image rows")
    print(f"Wrote {len(unionize_rows)} structure unionize rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
