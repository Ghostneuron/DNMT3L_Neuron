#!/usr/bin/env python3
"""Build GitHub/Zenodo code-release manifests."""

from __future__ import annotations

import ast
import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"

INCLUDE_DIRS = [
    "scripts",
    "config",
    "docs",
    ".github",
]

INCLUDE_FILES = [
    "README.md",
    "LICENSE",
    "CITATION.cff",
    ".zenodo.json",
    "CODE_AVAILABILITY.md",
    "requirements.txt",
    "environment.yml",
    ".gitignore",
]

EXCLUDED_SUFFIXES = {".pyc", ".DS_Store"}
EXCLUDED_PARTS = {"__pycache__"}
GENERATED_MANIFESTS = {
    "docs/code_file_manifest.tsv",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def iter_release_files() -> list[Path]:
    paths: list[Path] = []
    for rel in INCLUDE_FILES:
        path = ROOT / rel
        if path.is_file():
            paths.append(path)
    for rel_dir in INCLUDE_DIRS:
        base = ROOT / rel_dir
        if not base.exists():
            continue
        for path in sorted(base.rglob("*")):
            if not path.is_file():
                continue
            parts = set(path.relative_to(ROOT).parts)
            if parts & EXCLUDED_PARTS:
                continue
            if path.name in EXCLUDED_SUFFIXES or path.suffix in EXCLUDED_SUFFIXES:
                continue
            if path.relative_to(ROOT).as_posix() in GENERATED_MANIFESTS:
                continue
            paths.append(path)
    return sorted(set(paths), key=lambda p: p.relative_to(ROOT).as_posix())


def script_description(path: Path) -> str:
    try:
        tree = ast.parse(path.read_text(errors="replace"))
    except SyntaxError:
        return ""
    doc = ast.get_docstring(tree)
    if not doc:
        return ""
    return " ".join(doc.strip().splitlines()).strip()


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    script_rows = []
    for path in sorted((ROOT / "scripts").glob("*")):
        if path.suffix not in {".py", ".R"}:
            continue
        rel = path.relative_to(ROOT).as_posix()
        script_rows.append(
            {
                "script": rel,
                "language": "R" if path.suffix == ".R" else "Python",
                "description": script_description(path) if path.suffix == ".py" else "",
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    write_tsv(DOCS / "script_manifest.tsv", script_rows, ["script", "language", "description", "bytes", "sha256"])

    files = iter_release_files()
    file_rows = []
    for path in files:
        rel = path.relative_to(ROOT).as_posix()
        file_rows.append(
            {
                "path": rel,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    write_tsv(DOCS / "code_file_manifest.tsv", file_rows, ["path", "bytes", "sha256"])
    print(f"Wrote {DOCS / 'code_file_manifest.tsv'}")
    print(f"Wrote {DOCS / 'script_manifest.tsv'}")


if __name__ == "__main__":
    main()
