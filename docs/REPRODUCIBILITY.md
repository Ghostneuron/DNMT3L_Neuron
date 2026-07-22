# Reproducibility Guide

This guide summarizes the code organization and the main commands needed to recreate the paper-facing analyses. Paths assume execution from the repository root.

## 1. Environment

Create the Python/R environment:

```bash
conda env create -f environment.yml
conda activate dnmt3l-neuron
```

If using an existing environment:

```bash
python -m pip install -r requirements.txt
```

R scripts require Bioconductor `DESeq2` and `edgeR`.

HOMER motif analyses require an external HOMER installation. Set `HOMER_HOME` to the installation root, add `$HOMER_HOME/bin` to `PATH`, and set `HG38_FASTA` (or pass `--genome-fasta`) for GC-matched analyses.

## 2. Data Inventory And Download

The dataset manifest is:

```text
config/datasets.tsv
```

Download and inventory public GEO resources:

```bash
python scripts/download_geo_metadata.py --manifest config/datasets.tsv
python scripts/collect_ncbi_ds_geo_candidates.py --outdir results/ds_ncbi_geo
python scripts/inventory_geo_supplementary_files.py \
  --accession-table results/ds_ncbi_geo/ncbi_geo_high_priority_brain_neural_candidates.tsv \
  --out results/ds_ncbi_geo/high_priority_supplementary_inventory.tsv
```

Some large GEO supplementary files were downloaded manually or with accession-specific scripts. Controlled-access dbGaP matrices are not included.

## 3. Target-Gene And Expression Analyses

Normal and Down syndrome target-gene summaries:

```bash
python scripts/summarize_brainspan_target_expression.py
python scripts/summarize_gse104276_single_cell_targets.py
python scripts/compare_ds_gene_expression_targets.py
python scripts/analyze_gse305153_cellxgene_targets.py
```

Count-based expression reanalyses:

```bash
Rscript scripts/analyze_human_ds_bulk_count_de.R
Rscript scripts/analyze_gse268951_count_de.R
Rscript scripts/analyze_drosophila_dnmt3l_count_de.R
```

## 4. Methylation And Candidate Discovery

Human DNMT3L-OE methylation:

```bash
python scripts/call_human_epic_candidate_dmrs.py
python scripts/validate_gse157078_signature_in_gse168276_shsy5y.py
```

Mouse Dnmt3l-OE methylation:

```bash
python scripts/call_mouse_candidate_dmrs.py
python scripts/rank_mouse_methylation_expression_candidate_genes.py
```

Human/mouse integration:

```bash
python scripts/integrate_human_methylation_support_with_mouse_candidates.py
python scripts/integrate_mouse_candidates_with_human_ds_expression.py
python scripts/build_final_candidate_evidence_matrix.py
python scripts/build_manuscript_candidate_summary.py
```

## 5. HOMER And Chromatin Analyses

Prepare and run HOMER inputs:

```bash
python scripts/prepare_homer_motif_inputs.py
python scripts/prepare_ds_homer_motif_inputs.py
python scripts/prepare_mouse_homer_first_pass.py
python scripts/prepare_mouse_ds_homer_inputs.py
python scripts/prepare_gse305146_scATAC_homer_inputs.py
python scripts/run_homer_motif_enrichment.py
python scripts/summarize_homer_motif_results.py
python scripts/analyze_homer_cross_context_motif_concordance.py
```

Fetal Down syndrome multiome/chromatin analyses:

```bash
python scripts/analyze_gse305146_marker_atac_peak_da.py
python scripts/analyze_gse305146_marker_atac_motif_deviation.py
python scripts/analyze_gse305146_marker_atac_matched_motif_deviation.py --plot-only
python scripts/analyze_gse305146_motif_instance_specificity.py
python scripts/analyze_gse305146_marker_atac_gc_chromvar.py
python scripts/analyze_gse305146_peak_da_ranked_motif_shift.py
```

## 6. Developmental Module And Mechanism Analyses

```bash
python scripts/test_developmental_module_shifts.py
python scripts/analyze_glial_module_methylation_overlap.py
python scripts/analyze_neurogenesis_and_dnmt3a_independence.py
python scripts/compare_dnmt3l_modules_to_ds_brain.py
python scripts/analyze_human_ds_neurogenesis_timing_modules.py
python scripts/build_final_mechanism_axis_summary.py
```

## 7. Release Manifests

```bash
python scripts/build_code_release_manifest.py
```

The code-release manifest files are written to `docs/`.

## 8. Interpretation Boundaries

- Public processed expected-count reanalyses are not raw FASTQ re-quantification.
- Methylation-to-gene links based on proximity or same-TAD support are regulatory plausibility annotations, not causal enhancer-promoter validation.
- STAT pathway activity is not inferred from motif enrichment alone; phosphorylation and perturbation experiments in prior studies are separate evidence.
