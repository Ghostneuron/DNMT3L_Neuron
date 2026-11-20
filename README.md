# DNMT3L Neuron

Analysis code and reproducibility materials for **Context-dependent concordance of DNMT3L-sensitive methylation in Down syndrome neurodevelopment**.

Jie Lu, Independent researcher, San Francisco, California, USA.

## R1 reproducibility materials

Version 1.1.0 adds the analyses for the revised results.

- Human methylation concordance compares the observed sign agreement with a stratified expectation within the archived DNMT3L-sensitive CpG set.
- The mouse segmental-trisomy comparison evaluates the 3,089 CpGs available in both comparison tracks.
- Candidate scores use gene-level GSE268951 expression and are evaluated across alternative weights and evidence components.

The public summary tables are in [`docs/r1/`](docs/r1/). See the [R1 reproducibility guide](docs/R1_REPRODUCIBILITY.md) for commands, input requirements and interpretation, and the [release notes](docs/RELEASE_NOTES_20261120.md) for changes. The original Python and R pipeline remains available below.

## Original pipeline and repository contents

- `scripts/`: Python and R analysis scripts.
- `config/`: dataset, target-gene, and regulatory-source manifests.
- `docs/`: release notes, reproducibility guide, file manifests, and data-source notes.
- `requirements.txt`: Python package requirements for most scripts.
- `environment.yml`: Conda-style environment summary.
- `.zenodo.json`: draft Zenodo metadata for a software/code release.
- `CITATION.cff`: GitHub citation metadata.

Large downloaded data, controlled-access data, manuscript drafts, and full result tables are stored outside the code repository. The R1 summary tables are in [`docs/r1/`](docs/r1/). They can be regenerated or obtained from the sources listed in `config/datasets.tsv` and `docs/DATA_SOURCES.md`.

## Quick Start

Create the analysis environment:

```bash
conda env create -f environment.yml
conda activate dnmt3l-neuron
```

or install the Python package set into an existing environment:

```bash
python -m pip install -r requirements.txt
```

Install the R/Bioconductor dependencies separately:

```r
install.packages(c("BiocManager"))
BiocManager::install(c("DESeq2", "edgeR"))
```

HOMER must be installed separately for motif analyses. Set `HOMER_HOME` to the HOMER installation root and add `$HOMER_HOME/bin` to `PATH`. For GC-matched analyses, set `HG38_FASTA` or pass `--genome-fasta`.

## Data Layout

Downloaded data are expected under:

```text
data/geo/<ACCESSION>/
data/allen/
data/cellxgene/
data/reference/
data/external/
data/zenodo/
```

These directories are ignored by git. The scripts assume this layout unless command-line arguments override paths.

## Main Reanalysis Blocks

Download and inventory GEO metadata:

```bash
python scripts/download_geo_metadata.py --manifest config/datasets.tsv
python scripts/inventory_geo_supplementary_files.py \
  --accession-table results/ds_ncbi_geo/ncbi_geo_high_priority_brain_neural_candidates.tsv \
  --out results/ds_ncbi_geo/high_priority_supplementary_inventory.tsv
```

Run count-based expression reanalyses:

```bash
Rscript scripts/analyze_gse268951_count_de.R
Rscript scripts/analyze_drosophila_dnmt3l_count_de.R
Rscript scripts/analyze_human_ds_bulk_count_de.R
```

Run methylation and candidate-link analyses:

```bash
python scripts/call_human_epic_candidate_dmrs.py
python scripts/call_mouse_candidate_dmrs.py
python scripts/integrate_human_methylation_support_with_mouse_candidates.py
python scripts/build_final_candidate_evidence_matrix.py
```

Run motif and chromatin-context analyses:

```bash
python scripts/prepare_homer_motif_inputs.py
python scripts/run_homer_motif_enrichment.py
python scripts/analyze_homer_cross_context_motif_concordance.py
python scripts/analyze_gse305146_marker_atac_peak_da.py
python scripts/analyze_gse305146_marker_atac_gc_chromvar.py
```

Build final evidence summaries:

```bash
python scripts/build_manuscript_candidate_summary.py
python scripts/build_final_mechanism_axis_summary.py
```

See `docs/REPRODUCIBILITY.md` for a fuller pipeline map.

## Release Manifests

Generate code-release manifests after changes:

```bash
python scripts/build_code_release_manifest.py
```

This writes:

- `docs/code_file_manifest.tsv`
- `docs/script_manifest.tsv`

## Important Limits

This repository provides reproducible code and accession-based provenance, not a frozen copy of every raw dataset. Some inputs are public but very large; some postnatal Down syndrome PFC resources require controlled access through dbGaP. Public processed-count reanalyses are used where raw FASTQ re-quantification was not feasible.

## License

The code is released under the MIT License. Data remain governed by the terms of the original repositories and accessions.
