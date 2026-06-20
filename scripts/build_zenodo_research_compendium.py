#!/usr/bin/env python3
"""Build the manuscript-facing Zenodo research compendium."""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATE_TAG = "20260620"
DEST = ROOT / f"Zenodo_Research_Compendium_{DATE_TAG}"
ZIP_PATH = ROOT / "release" / f"DNMT3L_Neuron_research_compendium_{DATE_TAG}.zip"
CODE_MANIFEST = ROOT / "docs" / "code_file_manifest.tsv"
README_TEMPLATE = ROOT / "docs" / "ZENODO_RESEARCH_COMPENDIUM_README.md"

RESULT_EXTENSIONS = {".tsv", ".csv", ".txt", ".md", ".log", ".bed"}
MAX_DERIVED_RESULT_BYTES = 20 * 1024 * 1024
EXCLUDED_NAMES = {".DS_Store"}
EXCLUDED_RESULT_PARTS = {"homer_out", "backgrounds", "__pycache__"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def copy_file(source: Path, archive_rel: Path, rows: list[dict[str, object]]) -> None:
    target = DEST / archive_rel
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    rows.append(
        {
            "archive_path": archive_rel.as_posix(),
            "source_path": source.relative_to(ROOT).as_posix(),
            "bytes": source.stat().st_size,
            "sha256": sha256(source),
        }
    )


def iter_code_release_paths() -> list[Path]:
    paths: list[Path] = []
    with CODE_MANIFEST.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            paths.append(ROOT / row["path"])
    paths.append(CODE_MANIFEST)
    return sorted(set(paths), key=lambda path: path.relative_to(ROOT).as_posix())


def iter_manuscript_files() -> list[Path]:
    names = [
        "Manuscript_DNMT3L_Neuron_MolNeurobio_submission_Final.docx",
        "Supplementary Information.docx",
        "graphical_abstract_DNMT3L_Neuron_v2.png",
        "graphical_abstract_DNMT3L_Neuron_v2.pdf",
    ]
    return [ROOT / "manuscript" / name for name in names if (ROOT / "manuscript" / name).is_file()]


def iter_tree_files(base: Path) -> list[Path]:
    if not base.exists():
        return []
    return sorted(
        [path for path in base.rglob("*") if path.is_file() and path.name not in EXCLUDED_NAMES],
        key=lambda path: path.relative_to(ROOT).as_posix(),
    )


def iter_derived_result_files() -> list[Path]:
    base = ROOT / "results"
    if not base.exists():
        return []
    paths: list[Path] = []
    for path in sorted(base.rglob("*"), key=lambda item: item.as_posix()):
        if not path.is_file() or path.name in EXCLUDED_NAMES:
            continue
        rel_parts = set(path.relative_to(ROOT).parts)
        if rel_parts & EXCLUDED_RESULT_PARTS:
            continue
        if path.suffix.lower() not in RESULT_EXTENSIONS:
            continue
        if path.stat().st_size > MAX_DERIVED_RESULT_BYTES:
            continue
        paths.append(path)
    return paths


def write_manifest(rows: list[dict[str, object]]) -> None:
    manifest = DEST / "zenodo_compendium_manifest.tsv"
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["archive_path", "source_path", "bytes", "sha256"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda row: str(row["archive_path"])))


def build_zip() -> None:
    ZIP_PATH.parent.mkdir(exist_ok=True)
    tmp = ZIP_PATH.with_suffix(".zip.tmp")
    with zipfile.ZipFile(tmp, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zip_handle:
        for path in sorted(DEST.rglob("*"), key=lambda item: item.as_posix()):
            if path.is_file():
                zip_handle.write(path, arcname=f"{DEST.name}/{path.relative_to(DEST).as_posix()}")
    tmp.replace(ZIP_PATH)


def build_compendium(overwrite: bool) -> None:
    if DEST.exists():
        if not overwrite:
            raise SystemExit(f"{DEST} already exists; rerun with --overwrite to rebuild it.")
        shutil.rmtree(DEST)
    DEST.mkdir(parents=True)

    rows: list[dict[str, object]] = []
    copy_file(README_TEMPLATE, Path("README.md"), rows)

    for source in iter_code_release_paths():
        copy_file(source, Path("code") / source.relative_to(ROOT), rows)

    for source in iter_manuscript_files():
        copy_file(source, Path("manuscript") / source.name, rows)

    for source in iter_tree_files(ROOT / "manuscript" / "Final figures"):
        copy_file(source, Path("final_figures") / source.relative_to(ROOT / "manuscript" / "Final figures"), rows)

    for source in iter_tree_files(ROOT / "manuscript" / "Supplementary Tables"):
        copy_file(
            source,
            Path("supplementary_tables") / source.relative_to(ROOT / "manuscript" / "Supplementary Tables"),
            rows,
        )

    for source in iter_tree_files(ROOT / "figures"):
        copy_file(source, Path("source_figure_panels") / source.relative_to(ROOT / "figures"), rows)

    for source in iter_derived_result_files():
        copy_file(source, Path("derived_results") / source.relative_to(ROOT / "results"), rows)

    write_manifest(rows)
    build_zip()
    print(f"Wrote {DEST}")
    print(f"Wrote {ZIP_PATH}")
    print(f"Archived {len(rows)} files plus zenodo_compendium_manifest.tsv")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--overwrite", action="store_true", help="rebuild an existing compendium folder")
    args = parser.parse_args()
    build_compendium(overwrite=args.overwrite)


if __name__ == "__main__":
    main()
