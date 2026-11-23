# DNMT3L Neuron

Analysis code and reproducibility materials for **Dye bias and donor variation limit the evidence for DNMT3L-driven methylation changes in fetal Down syndrome brain**.

Jie Lu, Independent researcher, San Francisco, California, USA.

## Brief Report reproducibility materials

Zenodo v1.3.0 contains the donor-level fetal methylation analysis, matched CpG benchmark, paired comparison with the human DNMT3A W330R fibroblast signature, and the array and donor sensitivity addendum.

Download the [version 1.3.0 research compendium on Zenodo](https://doi.org/10.5281/zenodo.24431348) for the frozen signature, analysis code, native result tables and figures.

| Fetal cohort | Stratified DNMT3L OR (95% donor bootstrap CI) | Donor permutation P | Balanced DNMT3L and DNMT3A ORs | Paired ratio (95% donor bootstrap CI) |
| --- | ---: | ---: | ---: | ---: |
| GSE74486 fetal cerebrum | 1.63 (0.66–3.16) | 0.539 | 1.63 and 0.76 | 2.15 (0.81–3.91) |
| GSE73747 fetal cortex | 1.62 (0.96–2.24) | 0.236 | 1.57 and 0.75 | 2.10 (1.16–2.88) |
| GSE73747 deposited beta sensitivity | 1.93 (0.97–2.48) | — | 1.87 and 0.77 | 2.42 (1.16–3.24) |

Both fetal cohorts exceeded all 2,000 matched CpG panels with the observed donors held fixed. Both primary ORs were inside their donor-permutation null ranges. The fibroblast comparator differs from the DNMT3L model in lineage and perturbation.

The [dated post hoc sensitivity analyses](docs/brief_report/sensitivity/) include the revised Figure 2 source, donor permutation tests, and the external Spiers age-slope adjustment. The adjusted ORs were 1.85 (95% donor bootstrap CI 0.92–2.64, permutation P = 0.226) and 1.87 (1.09–2.53, P = 0.123) in GSE74486 and GSE73747, respectively.

See the [Brief Report reproducibility guide](docs/brief_report/README.md) and [full estimates](docs/brief_report/estimates.tsv) for the frozen signature, inputs, and original commands.

The [array and donor sensitivity addendum](docs/brief_report/array_sensitivity/)
contains the Table S6/S8 code and results, the Qin-rule reanalysis, the 21
November 2026 exploratory power simulation, and the Table S7G global
mean-difference analysis. Across 865,743 quality-filtered CpGs, the global
overexpression-minus-control mean beta difference was +0.011882 for raw beta,
+0.011719 after background-only Noob, and +0.001035 after ssNoob dye correction.

## Archived R1 reproducibility materials

Version 1.1.0 documents the preceding analysis, “Context-dependent concordance of DNMT3L-sensitive methylation in Down syndrome neurodevelopment.”

Download the [version 1.1.0 research compendium on Zenodo](https://doi.org/10.5281/zenodo.20770052) for the processed inputs, full supporting tables, expected outputs and figures.

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
