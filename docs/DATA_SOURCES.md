# Data Sources

This project is accession-driven. The full working dataset manifest is `config/datasets.tsv`; the manuscript-facing dataset table is `manuscript/Supplementary Tables/ST1_datasets.tsv` in the local submission package.

## Main Public Data Classes

- Human DNMT3L overexpression methylation and RNA-seq: `GSE157081`, `GSE157078`, `GSE157080`.
- Human SH-SY5Y DNMT3L overexpression WGBS validation: `GSE168276`.
- Mouse Dnmt3l overexpression methylation and RNA-seq: `GSE271949`, `GSE268951`.
- Mouse Dp10 Down syndrome model expression: `GSE214273`.
- Human fetal Down syndrome expression and multiome context: `GSE301886`, `GSE305139`, `GSE305146`, `GSE305153`.
- Human adult Down syndrome methylation context: `GSE63347`, `GSE74486`.
- Mouse Down syndrome methylation context: `GSE74505`.
- Drosophila human DNMT3L overexpression context: `GSE263069`; array comparator `GSE68861`.
- Dnmt3a neuron reference comparison: `GSE141587`.
- Human fetal single-cell target-gene context: `GSE104276`.
- Allen Institute references: BrainSpan, Allen in situ hybridization, and Allen Brain Cell Atlas resources.
- Postnatal Down syndrome PFC contextual resource: Zenodo record `17527842` and controlled-access dbGaP `phs004098.v1.p1`.

## Local Data Layout

The scripts expect downloaded inputs under:

```text
data/geo/<ACCESSION>/
data/allen/
data/cellxgene/
data/reference/
data/external/
data/zenodo/
data/dbgap/
```

These directories are intentionally ignored by git. The code release records accessions and analysis scripts, not the full 24 GB local data cache.

## Controlled-Access Data

The postnatal Down syndrome PFC matrices are not redistributed here. The public Shiny/Zenodo resources were used only as contextual or visual audit support unless controlled-access matrices are available to the user under the relevant dbGaP terms.
