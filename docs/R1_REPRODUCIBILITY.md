# Reproduce the R1 Analyses

Version 1.1.0 accompanies *Context-dependent concordance of DNMT3L-sensitive methylation in Down syndrome neurodevelopment*. Python 3.10 or newer is sufficient for the analytical workflow. The figure script also needs ReportLab 4.4.9 to match the tested environment.

## Inputs

Extract the research compendium and use its `inputs/` directory. It holds the public-data-derived probe comparisons, candidate evidence, and gene-level effect columns needed by the R1 script. See [the dataset manifest](../config/r1_datasets.tsv) for study roles. The published source studies and accession terms govern source data.

From the compendium root:

```bash
python code/scripts/reanalyze_r1.py inputs --outdir reproduced
python -m pip install 'reportlab==4.4.9'
python code/scripts/draw_r1_figures.py reproduced/results.json --outdir reproduced/figures
```

From a GitHub checkout, use `scripts/` instead of `code/scripts/` and point the analysis at the extracted compendium's `inputs/` directory.

The analysis creates `human_concordance.tsv`, `mouse_shared_cpg_summary.tsv`, `mouse_shared_cpg_values.tsv`, `rank_sensitivity.tsv`, `all_discovery_ranks.tsv`, `curated_candidate_score_audit.tsv`, and `results.json`. The four compact summaries can be compared with [`docs/r1/`](r1/). The complete compendium includes all seven results for byte-level comparison.

## Read the outputs

- The human analysis compares observed sign agreement with the expectation given the marginal directions within strata of chromosome, CpG-island annotation and control methylation. It uses the archived perturbation-sensitive probes and has no donor-level uncertainty estimate.
- The mouse comparison uses exact genomic coordinates shared by both differential methylation tracks. It describes one sample per genotype.
- The candidate analysis recalculates an exploratory discovery score using gene-level GSE268951 results. The score's weights reflect evidence-prioritization choices. The alternative rankings quantify their effect.
- The ATAC sample manifest in the compendium records the archived fetal analysis cohort and its library counts. Source provenance and availability are listed by accession.

The expanded data-generation pipeline is described in [REPRODUCIBILITY.md](REPRODUCIBILITY.md).
