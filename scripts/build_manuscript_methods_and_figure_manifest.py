#!/usr/bin/env python3
"""Write manuscript methods text and a figure/table manifest."""

from __future__ import annotations

import csv
from pathlib import Path


METHODS_NOTE = Path("notes/manuscript_methods_candidate_mechanism.md")
FIGURE_PLAN_NOTE = Path("notes/manuscript_figure_table_plan.md")
MANIFEST = Path("results/manuscript_summary/figure_table_manifest.tsv")


METHODS_TEXT = """# Manuscript Methods Draft: DNMT3L-Linked Candidate Mechanism Analysis

Date: 2026-06-09

This draft covers the second analysis section of the project: integration of DNMT3L/Dnmt3l perturbation methylation data, Down syndrome-associated methylation and expression data, mouse model expression data, and fetal brain chromatin annotations to prioritize candidate downstream neurodevelopmental genes.

## Dataset Selection And Analysis Scope

The mechanism analysis used public datasets that could connect at least one of four evidence layers: DNMT3L/Dnmt3l-sensitive methylation, Down syndrome-associated methylation, Down syndrome or trisomy-model expression, and fetal brain regulatory chromatin. Datasets were tracked in `config/datasets.tsv` with accession, species, assay, system, design, planned use, and processing status.

The primary human DNMT3L perturbation methylation anchor was `GSE157078`, an EPIC methylation array study in cultured human fetal cortical neurons with DNMT3L overexpression versus control. `GSE168276`, a DNMT3L-overexpression WGBS dataset in SH-SY5Y neuronal differentiation, was used as an orthogonal human neuronal-model validation layer. Human Down syndrome methylation support was evaluated from fetal brain-relevant and DS methylation datasets processed in the human methylation convergence workflow. Human DS expression support was taken from the fetal bulk, fetal cortex single-cell/multiome, and CELLxGENE fetal cortex analyses described in the expression methods.

The primary mouse mechanism anchors were `GSE271949`, a P28 medial frontal cortex Dnmt3l conditional overexpression methylation dataset in the Emx1 lineage, and `GSE268951`, the matched Dnmt3l overexpression RNA-seq/noncoding RNA profiling dataset. `GSE214273` was used as a Dp(10)2Yey adult hippocampus expression context. The mouse methylation analysis also used DNMT3A-dependent neuronal methylation and postnatal methylation-gain reference tracks from the DNMT3A neuron-dependency workflow.

Regulatory annotation used human fetal telencephalon ATAC/OCR, putative regulatory element, H3K27ac, H3K4me1, and H3K27me3 tracks from `GSE149268`. Mouse regulatory annotation used embryonic or developing brain ATAC and histone-mark tracks from `GSE172703`, `GSE201400`, `GSE201487`, and `GSE201488`.

## Mouse Dnmt3l-Sensitive Methylation And Candidate Ranking

Mouse Dnmt3l-sensitive DMRs were called from the `GSE271949` methylation data and summarized by genomic coordinate, methylation direction, mean methylation difference, and overlap with reference methylation or chromatin features. A high-confidence mouse methylation mechanism signal was defined as a triple-hit DMR when a Dnmt3l overexpression DMR overlapped neuronal DNMT3A-required methylation, normal postnatal methylation gain, and a region with altered H3K27me3 after Dnmt3a loss. This definition was designed to prioritize loci where Dnmt3l-sensitive methylation is embedded in a known developmental neuronal methylation program.

DMRs were linked to genes by proximity and genomic relationship. Links were classified as promoter-proximal when the DMR was within 2 kb of a transcription start site, gene-body linked when the DMR overlapped the gene body, or nearest-neighbor linked when the closest gene was within 100 kb. These links were treated as candidate regulatory links, not proof of direct target-gene regulation.

Mouse methylation-linked genes were ranked by the number and quality of linked DMRs, methylation direction, relationship class, Dnmt3l overexpression RNA-seq direction in `GSE268951`, and Dp10 expression direction in `GSE214273`. Human-mouse comparison used a first-pass ortholog mapping by gene symbol and upper-case symbol matching, with manual checks for the final candidate set.

## Human DNMT3L And Down Syndrome Methylation Support

Human DNMT3L-sensitive CpGs from `GSE157078` were mapped to hg38 and linked to GENCODE v39 protein-coding genes using the same relationship classes used for mouse DMRs: promoter 2 kb, gene body, and nearest gene within 100 kb. Probe-level methylation direction was summarized as mean delta beta for DNMT3L overexpression versus control.

To identify human methylation convergence, DNMT3L-sensitive probes were intersected with DS-associated methylation direction. A strict human fetal-positive support probe required a DNMT3L-sensitive methylation change and a concordant positive signed DS methylation signal in fetal brain-relevant data. SH-SY5Y `GSE168276` Phase 3 WGBS validation was used as an additional neuronal-model concordance layer, but was not treated as fetal brain tissue evidence.

Human methylation support scores combined the number of linked DNMT3L-sensitive probes, strict fetal-positive probes, concordance with SH-SY5Y validation where available, and DS methylation direction. These scores were used for prioritization and visualization rather than as formal probability estimates.

## Expression Integration

Human expression evidence came from the DS/control expression workflow, including fetal bulk RNA-seq, GSE59630 array data, GSE305146 multiome target-gene summaries, and the GSE305153/CELLxGENE fetal cortex single-nucleus analysis. For each candidate gene, the integration used the strongest available within-dataset DS/control effect, with context labels preserved for tissue, brain region, or cell type. Single-cell expression was summarized at donor or sample pseudobulk level when donor metadata were available, and cell-level observations were not treated as independent biological replicates.

Mouse expression evidence used `GSE268951` Dnmt3l overexpression RNA-seq and the `GSE214273` Dp(10)2Yey hippocampus pseudobulk single-cell analysis. For `GSE268951`, primary differential-expression effects were reanalyzed from the public GEO supplementary table `GSE268951_all.fpkm_mRNA_anno.txt.gz`, which contains per-sample processed `expected_count` columns in addition to FPKM. Expected-count columns were collapsed to gene level and analyzed with DESeq2 and edgeR using sex and condition in the design where applicable; FPKM values were retained only for abundance and sex-context summaries. The local `CON-VS-OE_DE_anno.xls` table was not used for manuscript evidence because it is not listed as a GEO public supplementary file. Because the Dp10 dataset contains only two model and two control samples and lacks cell-type labels, it was interpreted descriptively.

Expression effects were compared only within each dataset or context. Absolute expression levels were not pooled across platforms.

## Regulatory Annotation

Human DNMT3L-sensitive CpG-gene links and strict human support probes were annotated for overlap with fetal telencephalon regulatory features from `GSE149268`: open chromatin/OCR, putative regulatory elements, ATAC peaks, H3K27ac, H3K4me1, and H3K27me3. Mouse triple-hit DMR-gene links were annotated for overlap with mouse embryonic/developing brain ATAC, H3K27ac, H3K4me3, and H3K27me3 tracks from `GSE172703`, `GSE201400`, `GSE201487`, and `GSE201488`.

All regulatory overlaps used BED-style interval intersection at the CpG/probe or DMR coordinate level. A gene was considered regulatory-supported when at least one linked methylation feature overlapped a relevant fetal or embryonic brain regulatory annotation. Regulatory annotation was used to ask whether candidate methylation features occur in plausible regulatory chromatin, not to assign enhancer target genes definitively.

## Candidate Tiering

Final candidate genes were grouped into four manuscript tiers.

Tier 1 primary methylation-convergent candidates required balanced cross-species support: mouse triple-hit DMR evidence, human DNMT3L-sensitive methylation, DS methylation concordance, expression evidence in human or mouse data, and regulatory annotation support. This tier contains `PLPP3`, `FOXP1`, `FGFR2`, `ZFHX3`, and `SIAH3`.

Tier 2 secondary balanced candidates had methylation and expression support but lower overall priority, less complete cross-species convergence, or less direct directionality than the primary set. This tier contains `FAM222A`, `COL26A1`, `KIRREL3`, `CAMTA1`, and `LMX1B`.

Tier 3 methylation-heavy developmental candidates had strong methylation and regulatory evidence, often with H3K27me3/Polycomb context, but weaker or more context-restricted expression evidence. This tier contains `OTX1`, `OTX2`, `LHX2`, `DRAXIN`, `NFIA`, `PRDM16`, and `NR2F1`.

Tier 4 expression/pathway candidates were expression- or pathway-relevant genes with incomplete strict human methylation support. This tier contains `FZD7`, `TGIF2`, `PRDM8`, `CRABP2`, `BCL11B`, `RXRA`, and `CYP26A1`. These genes were not labeled methylation-independent; rather, their methylation dependence remains unresolved in the current data.

## Enrichment Testing

Two complementary enrichment approaches were used. First, gene-level hypergeometric tests compared candidate-gene regulatory support against gene-level backgrounds derived from all human DNMT3L-linked probe genes, strict human triple-support probe genes, or mouse triple-hit DMR-linked genes. These tests ask whether candidate genes are enriched for at least one regulatory annotation compared with a broad gene background.

Second, matched link-level permutation tests were used as the main regulatory enrichment test. Human CpG-gene links were matched by gene relationship class and DNMT3L methylation direction. Mouse DMR-gene links were matched by relationship class and DMR methylation direction. For each candidate set and regulatory assay, 10,000 random matched candidate sets were sampled with seed `20260609`. The observed number of regulatory-overlap records was compared with the matched null distribution. The reported fold over null was calculated as observed successes divided by mean null successes, and the empirical greater-tail p value was calculated with a plus-one correction.

The matched permutation tests are stronger than the broad hypergeometric tests because they preserve the relationship and methylation-direction structure of the candidate records. They still do not fully match for all possible confounders such as CpG density, CpG island state, gene length, chromosomal neighborhood, and mappability.

## Methylation-Expression Direction Classification

For each candidate gene, methylation direction was summarized separately for human DNMT3L-sensitive probes and mouse Dnmt3l-sensitive DMRs. Human expression direction was taken from the strongest available DS/control effect among all-cell, bulk, and cell-type-specific contexts. Mouse expression direction was taken from Dnmt3l overexpression RNA-seq and Dp10 hippocampus pseudobulk summaries.

Methylation-expression direction was classified as repression-like when increased methylation was paired with decreased expression, or decreased methylation was paired with increased expression, in a relevant context. Same-direction changes such as hypermethylation with increased expression were classified as directionally complex because they are not consistent with a simple promoter-silencing model. Directionally complex candidates may reflect enhancer methylation, gene-body methylation, indirect effects, developmental timing shifts, cell-state composition changes, or non-causal correlation.

FOXP1 showed the clearest primary-candidate direction compatible with methylation-associated repression. PLPP3, FGFR2, ZFHX3, and SIAH3 showed strong regulatory convergence but directionally complex methylation-expression relationships in at least one major context.

## Reproducibility Map

| Analysis layer | Main script | Key output |
| --- | --- | --- |
| Mouse Dnmt3l DMR calling | `scripts/call_mouse_candidate_dmrs.py` | `results/mouse_dnmt3l_dmrs/` |
| DNMT3A neuronal reference overlap | `scripts/overlap_mouse_dnmt3l_dmrs_with_dnmt3a_neuron_reference.py` | `results/dnmt3a_neuron_dependency_overlap/` |
| Mouse methylation-expression ranking | `scripts/rank_mouse_methylation_expression_candidate_genes.py` | `results/integrated_candidate_genes/` |
| Human DNMT3L methylation signature | `scripts/call_human_epic_candidate_dmrs.py` | `results/human_dnmt3l_epic/` |
| SH-SY5Y validation | `scripts/validate_gse157078_signature_in_gse168276_shsy5y.py` | `results/shsy5y_gse168276/` |
| Human methylation convergence | `scripts/integrate_human_methylation_support_with_mouse_candidates.py` | `results/human_methylation_convergence/` |
| Cross-species expression integration | `scripts/integrate_mouse_candidates_with_human_ds_expression.py` | `results/cross_species_candidate_expression/` |
| Final evidence matrix | `scripts/build_final_candidate_evidence_matrix.py` | `results/final_candidate_evidence/final_candidate_evidence_matrix.tsv` |
| Regulatory context annotation | `scripts/add_regulatory_context_to_final_candidate_matrix.py` | `results/final_candidate_evidence/final_candidate_evidence_matrix_with_regulatory_context.tsv` |
| Regulatory enrichment and classification | `scripts/test_regulatory_enrichment_and_classify_mechanisms.py` | `results/mechanism_classification/` |
| Matched regulatory permutations | `scripts/run_matched_regulatory_permutation_tests.py` | `results/mechanism_classification/matched_regulatory_permutation_tests.tsv` |
| Manuscript candidate summary | `scripts/build_manuscript_candidate_summary.py` | `results/manuscript_summary/final_candidate_manuscript_summary.tsv` |
| Expanded in-silico mechanism tests | `scripts/expand_final_candidate_mechanism_tests.py` | `results/manuscript_summary/final_candidate_in_silico_mechanism_tests.tsv` |
| Neurogenesis and DNMT3A-independence screen | `scripts/analyze_neurogenesis_and_dnmt3a_independence.py` | `results/neurogenesis_mechanism/` |

## Limitations

This analysis is designed for prioritization and mechanism generation, not definitive causal proof. Most DS methylation and expression layers are not paired in the same samples, so methylation-mediated expression change cannot be inferred directly. DMR-to-gene and CpG-to-gene links are proximity-based and may miss distal enhancer targets or assign methylation features to the wrong gene. The SH-SY5Y validation dataset is a neuronal model, not fetal brain tissue. The mouse Dnmt3l overexpression model is constitutive within the Emx1 lineage and may not capture endogenous, temporally restricted Dnmt3l regulation. The matched permutation tests control relationship class and methylation direction but not all genomic covariates. Finally, available mouse regulatory overlaps support individual candidate loci but do not show enrichment above matched mouse DMR background.

## Manuscript-Level Interpretation

The integrated analysis supports a working model in which altered DNMT3L dosage is linked to methylation changes in fetal-brain regulatory chromatin near neurodevelopmental genes. The strongest current evidence is regulatory convergence, especially in human fetal brain chromatin, rather than direct proof of methylation-mediated transcriptional repression. FOXP1 is the clearest repression-like primary candidate, while PLPP3, FGFR2, ZFHX3, and SIAH3 are strong but mechanistically more complex regulatory candidates.
"""


FIGURE_PLAN_TEXT = """# Manuscript Figure And Table Plan

Date: 2026-06-09

This plan maps each major manuscript claim to the current local figures, result tables, and remaining figure-building tasks.

## Proposed Main Figures

### Figure 1: DNMT3L/DNMT3A/DNMT3B/STAT1/STAT3 Expression Landscape

Primary files:

- `figures/target_gene_summary/human_brainspan_temporal_rpkm.png`
- `figures/target_gene_summary/human_brainspan_region_heatmaps.png`
- `figures/target_gene_summary/human_gse104276_fetal_pfc_week_trends.png`
- `figures/target_gene_summary/mouse_allen_ish_temporal_energy.png`
- `figures/ds_expression_comparison/ds_target_gene_effect_heatmap.png`
- `figures/gse305153_cellxgene_targets/target_gene_author_cell_type_effect_heatmap.png`

Main claim: DNMT3L/Dnmt3l is not broadly expressed in normal brain, but sparse human fetal cortex/DS-associated cell contexts can show low-level DNMT3L signal. DNMT3A/DNMT3B and STAT1/STAT3 provide the broader developmental and DS-associated expression context.

Caveat: For DNMT3L, most evidence is sparse detection or context-specific signal, not robust bulk upregulation.

### Figure 2: Integrated Final Candidate Evidence Matrix

Primary files:

- `figures/final_candidate_evidence/final_candidate_evidence_regulatory_heatmap.png`
- `figures/final_candidate_evidence/final_candidate_regulatory_support_barplot.png`
- `results/final_candidate_evidence/final_candidate_evidence_matrix_with_regulatory_context.tsv`

Main claim: PLPP3, FOXP1, FGFR2, ZFHX3, and SIAH3 form the primary methylation-convergent candidate set, with additional secondary, methylation-heavy, and expression/pathway tiers.

Caveat: A high evidence score means convergent public-data support, not causal validation.

### Figure 3: Primary Candidate Locus Panels

Primary files:

- `figures/primary_candidate_locus_panels/PLPP3_locus_evidence.png`
- `figures/primary_candidate_locus_panels/FOXP1_locus_evidence.png`
- `figures/primary_candidate_locus_panels/FGFR2_locus_evidence.png`
- `figures/primary_candidate_locus_panels/ZFHX3_locus_evidence.png`
- `figures/primary_candidate_locus_panels/SIAH3_locus_evidence.png`
- `results/final_candidate_evidence/primary_candidate_locus_panel_summary.tsv`

Main claim: Primary candidates have locus-level methylation evidence in human and mouse, with regulatory chromatin support at or near linked methylation features.

Caveat: Locus panels show proximity and chromatin context; they do not prove the methylation feature controls that gene.

### Figure 4: Regulatory Enrichment And Mechanism Classes

Primary files:

- `figures/mechanism_classification/matched_regulatory_permutation_fold_over_null.png`
- `figures/mechanism_classification/mechanism_class_counts.png`
- `results/mechanism_classification/matched_regulatory_permutation_tests.tsv`
- `results/mechanism_classification/methylation_expression_mechanism_table.tsv`

Main claim: Candidate methylation features are enriched for human fetal brain regulatory chromatin after matched link-level permutation, and candidate genes separate into repression-like, directionally complex, Polycomb/developmental, and unresolved expression/pathway groups.

Caveat: Human regulatory enrichment is strong, but mouse DMR regulatory overlap is not enriched above matched mouse DMR background.

### Figure 5: Working Model Schematic

Primary files:

- `figures/manuscript_model/dnmt3l_candidate_mechanism_model.png`
- `figures/manuscript_model/dnmt3l_candidate_mechanism_model.pdf`

Main claim: DNMT3L dosage may perturb fetal-brain regulatory methylation programs near neurodevelopmental genes, with FOXP1 as a repression-like example and PLPP3/FGFR2/ZFHX3/SIAH3 as regulatory-convergent but directionally complex examples.

Caveat: The model should be framed as a hypothesis for validation, not as a demonstrated causal pathway.

## Proposed Main Tables

### Table 1: Dataset And Analysis Module Summary

Primary file:

- `config/datasets.tsv`

Include accession/resource, species, system, assay, design, analysis role, status, and reason for exclusion where applicable.

### Table 2: Final Candidate Manuscript Summary

Primary file:

- `results/manuscript_summary/final_candidate_manuscript_summary.tsv`

Use this as the concise main candidate table. Include gene, tier, candidate group, human methylation, human expression, mouse methylation, mouse expression, regulatory context, and short interpretation.

### Table 3: Methylation-Expression Mechanism Classification

Primary file:

- `results/mechanism_classification/methylation_expression_mechanism_table.tsv`

Use this as a detailed supplementary table or a compact main table if the mechanism section is central.

### Table 4: Regulatory Enrichment Tests

Primary files:

- `results/mechanism_classification/matched_regulatory_permutation_tests.tsv`
- `results/mechanism_classification/regulatory_enrichment_tests.tsv`

Use matched permutation results as the main test. Keep gene-level hypergeometric tests as supporting or supplementary evidence.

### Table 5: Locus-Level Evidence Summary

Primary file:

- `results/final_candidate_evidence/primary_candidate_locus_panel_summary.tsv`

Use as a supplementary table that links each primary gene to its human probes, mouse DMRs, methylation direction, expression direction, and regulatory marks.

### Table 6: In-Silico Mechanism Test Roadmap

Primary files:

- `results/manuscript_summary/final_candidate_in_silico_mechanism_tests.tsv`
- `results/manuscript_summary/nonprimary_candidate_in_silico_mechanism_tests.tsv`
- `figures/manuscript_summary/nonprimary_candidate_in_silico_mechanism_test_map.png`

Use this as a supplementary mechanism roadmap. It keeps the final candidate set unchanged but adds explicit computational hypotheses, methylation/chromatin/expression/pathway tests, expected supporting patterns, counter-evidence, and caveats.

### Table 7: Neurogenesis And DNMT3A-Independence Screen

Primary files:

- `results/neurogenesis_mechanism/final_candidate_neurogenesis_annotation.tsv`
- `results/neurogenesis_mechanism/dnmt3l_vs_dnmt3a_transcriptomic_overlap.tsv`
- `results/neurogenesis_mechanism/dnmt3l_not_dnmt3a_neurogenesis_candidates.tsv`
- `results/neurogenesis_mechanism/developmental_module_rank_tests.tsv`
- `figures/neurogenesis_mechanism/final_candidate_neurogenesis_annotation_heatmap.png`
- `figures/neurogenesis_mechanism/dnmt3l_dnmt3a_transcriptomic_overlap_counts.png`
- `figures/neurogenesis_mechanism/developmental_module_shift_dotplot.png`

Use this as a functional follow-up table/figure set. It asks which final candidates are connected to neurogenesis annotations and which Dnmt3l-responsive neurogenesis genes are not altered in the Dnmt3a cKO reference.

## Recommended Figure Order For The Paper

1. Start with expression landscape of DNMT3L/DNMT3A/DNMT3B/STAT1/STAT3 in normal and DS brain.
2. Introduce the methylation/expression integration and the final candidate evidence matrix.
3. Show primary candidate locus panels for concrete examples.
4. Quantify regulatory enrichment and classify methylation-expression mechanisms.
5. End with a working model and validation roadmap.

## Immediate Next Build Task

The next useful build task is panel assembly: combine the existing expression, evidence-matrix, locus, enrichment, and model figures into journal-style multi-panel figures with consistent lettering and export dimensions.
"""


MANIFEST_ROWS = [
    {
        "item": "Figure 1",
        "type": "main_figure",
        "title": "DNMT3L target-panel expression landscape",
        "primary_file": "figures/target_gene_summary/human_brainspan_temporal_rpkm.png",
        "supporting_files": "figures/target_gene_summary/human_brainspan_region_heatmaps.png;figures/target_gene_summary/human_gse104276_fetal_pfc_week_trends.png;figures/target_gene_summary/mouse_allen_ish_temporal_energy.png;figures/ds_expression_comparison/ds_target_gene_effect_heatmap.png;figures/gse305153_cellxgene_targets/target_gene_author_cell_type_effect_heatmap.png",
        "main_claim": "DNMT3L/Dnmt3l is not broadly expressed in normal brain, with sparse context-specific DS/fetal cortex signal.",
        "caveat": "DNMT3L evidence is sparse and context-specific, not robust bulk upregulation.",
        "status": "existing_files_ready_needs_panel_assembly",
    },
    {
        "item": "Figure 2",
        "type": "main_figure",
        "title": "Integrated candidate evidence matrix",
        "primary_file": "figures/final_candidate_evidence/final_candidate_evidence_regulatory_heatmap.png",
        "supporting_files": "figures/final_candidate_evidence/final_candidate_regulatory_support_barplot.png;results/final_candidate_evidence/final_candidate_evidence_matrix_with_regulatory_context.tsv",
        "main_claim": "PLPP3, FOXP1, FGFR2, ZFHX3, and SIAH3 are primary methylation-convergent candidates.",
        "caveat": "Integrated support prioritizes candidates but does not prove causal regulation.",
        "status": "ready",
    },
    {
        "item": "Figure 3",
        "type": "main_figure",
        "title": "Primary candidate locus evidence",
        "primary_file": "figures/primary_candidate_locus_panels/PLPP3_locus_evidence.png",
        "supporting_files": "figures/primary_candidate_locus_panels/FOXP1_locus_evidence.png;figures/primary_candidate_locus_panels/FGFR2_locus_evidence.png;figures/primary_candidate_locus_panels/ZFHX3_locus_evidence.png;figures/primary_candidate_locus_panels/SIAH3_locus_evidence.png;results/final_candidate_evidence/primary_candidate_locus_panel_summary.tsv",
        "main_claim": "Primary candidates show locus-level human and mouse methylation evidence with fetal/developmental regulatory chromatin support.",
        "caveat": "Locus proximity does not prove enhancer-target or methylation-mediated regulation.",
        "status": "ready",
    },
    {
        "item": "Figure 4",
        "type": "main_figure",
        "title": "Regulatory enrichment and mechanism classes",
        "primary_file": "figures/mechanism_classification/matched_regulatory_permutation_fold_over_null.png",
        "supporting_files": "figures/mechanism_classification/mechanism_class_counts.png;results/mechanism_classification/matched_regulatory_permutation_tests.tsv;results/mechanism_classification/methylation_expression_mechanism_table.tsv",
        "main_claim": "Candidate methylation features are enriched for human fetal brain regulatory chromatin and separate into distinct mechanism classes.",
        "caveat": "Mouse regulatory overlap supports loci individually but is not enriched above matched mouse DMR background.",
        "status": "ready",
    },
    {
        "item": "Figure 5",
        "type": "main_figure",
        "title": "Working model of DNMT3L-linked regulatory convergence",
        "primary_file": "figures/manuscript_model/dnmt3l_candidate_mechanism_model.png",
        "supporting_files": "notes/manuscript_results_candidate_mechanism.md;notes/manuscript_methods_candidate_mechanism.md",
        "main_claim": "Altered DNMT3L dosage may perturb fetal regulatory methylation near neurodevelopmental genes.",
        "caveat": "Working model remains a hypothesis requiring paired methylation-expression or perturbation validation.",
        "status": "ready",
    },
    {
        "item": "Table 1",
        "type": "main_table",
        "title": "Dataset and analysis module summary",
        "primary_file": "config/datasets.tsv",
        "supporting_files": "notes/dataset_readiness.md;notes/atac_dataset_readiness.md;notes/chip_dataset_readiness.md",
        "main_claim": "The study integrates public human, mouse, methylation, expression, and chromatin datasets.",
        "caveat": "Some audited datasets were excluded from quantitative inference because of missing or ambiguous processed data.",
        "status": "ready_needs_manuscript_formatting",
    },
    {
        "item": "Table 2",
        "type": "main_table",
        "title": "Final candidate manuscript summary",
        "primary_file": "results/manuscript_summary/final_candidate_manuscript_summary.tsv",
        "supporting_files": "results/final_candidate_evidence/final_candidate_evidence_matrix_with_regulatory_context.tsv",
        "main_claim": "Final candidates can be grouped by primary, secondary, methylation-heavy, and expression/pathway categories.",
        "caveat": "Expression/pathway candidates are methylation-unresolved rather than proven methylation-independent.",
        "status": "ready",
    },
    {
        "item": "Table 3",
        "type": "main_or_supp_table",
        "title": "Methylation-expression mechanism classification",
        "primary_file": "results/mechanism_classification/methylation_expression_mechanism_table.tsv",
        "supporting_files": "notes/mechanism_enrichment_and_direction_first_pass.md",
        "main_claim": "FOXP1 is repression-like, while several other primary genes are directionally complex regulatory candidates.",
        "caveat": "Mechanism classes are based on cross-dataset directionality, not paired mediation.",
        "status": "ready",
    },
    {
        "item": "Table 4",
        "type": "main_or_supp_table",
        "title": "Regulatory enrichment tests",
        "primary_file": "results/mechanism_classification/matched_regulatory_permutation_tests.tsv",
        "supporting_files": "results/mechanism_classification/regulatory_enrichment_tests.tsv;results/mechanism_classification/regulatory_background_summary.tsv",
        "main_claim": "Human candidate methylation features are enriched for fetal brain regulatory annotations after matched permutation.",
        "caveat": "Matched tests do not control every genomic covariate.",
        "status": "ready",
    },
    {
        "item": "Table 5",
        "type": "supp_table",
        "title": "Primary candidate locus evidence summary",
        "primary_file": "results/final_candidate_evidence/primary_candidate_locus_panel_summary.tsv",
        "supporting_files": "figures/primary_candidate_locus_panels/",
        "main_claim": "Each primary candidate has explicit linked methylation features and regulatory annotations.",
        "caveat": "Coordinate-linked evidence should be followed by locus-specific validation.",
        "status": "ready",
    },
    {
        "item": "Table 6",
        "type": "supp_table",
        "title": "In-silico mechanism test roadmap",
        "primary_file": "results/manuscript_summary/final_candidate_in_silico_mechanism_tests.tsv",
        "supporting_files": "results/manuscript_summary/nonprimary_candidate_in_silico_mechanism_tests.tsv;figures/manuscript_summary/nonprimary_candidate_in_silico_mechanism_test_map.png;notes/expanded_candidate_in_silico_mechanism_tests.md",
        "main_claim": "The 19 non-primary candidates can be assigned to testable in-silico mechanisms without changing the final candidate set.",
        "caveat": "These are computational hypotheses and validation priorities, not causal proof.",
        "status": "ready",
    },
    {
        "item": "Table 7",
        "type": "supp_table_and_figure",
        "title": "Neurogenesis and DNMT3A-independence screen",
        "primary_file": "results/neurogenesis_mechanism/final_candidate_neurogenesis_annotation.tsv",
        "supporting_files": "results/neurogenesis_mechanism/dnmt3l_vs_dnmt3a_transcriptomic_overlap.tsv;results/neurogenesis_mechanism/dnmt3l_not_dnmt3a_neurogenesis_candidates.tsv;results/neurogenesis_mechanism/developmental_module_rank_tests.tsv;figures/neurogenesis_mechanism/final_candidate_neurogenesis_annotation_heatmap.png;figures/neurogenesis_mechanism/dnmt3l_dnmt3a_transcriptomic_overlap_counts.png;figures/neurogenesis_mechanism/developmental_module_shift_dotplot.png;notes/neurogenesis_and_dnmt3a_independence_first_pass.md;notes/developmental_module_shift_tests.md",
        "main_claim": "Many final candidates are neurogenesis/neurodevelopment-context genes, and Dnmt3l OE shifts synapse-window developmental modules even when most final candidates do not pass single-gene DE thresholds.",
        "caveat": "Dnmt3l OE and Dnmt3a cKO are different perturbation systems, so DNMT3A independence and methylation independence remain candidate interpretations.",
        "status": "ready",
    },
]


MANIFEST_FIELDS = [
    "item",
    "type",
    "title",
    "primary_file",
    "supporting_files",
    "main_claim",
    "caveat",
    "status",
]


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def write_manifest(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(MANIFEST_ROWS)


def main() -> None:
    write_text(METHODS_NOTE, METHODS_TEXT)
    write_text(FIGURE_PLAN_NOTE, FIGURE_PLAN_TEXT)
    write_manifest(MANIFEST)
    print(f"Wrote {METHODS_NOTE}")
    print(f"Wrote {FIGURE_PLAN_NOTE}")
    print(f"Wrote {MANIFEST}")


if __name__ == "__main__":
    main()
