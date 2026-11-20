# Brief Report donor sensitivity analyses

The original and follow-up plans are both dated 20 November 2026. The original cohort results were known before the first plan, and the original sensitivity results were known before the addendum. Both plans are post hoc. The analysis reuses the frozen donor reconstruction, probe strata, and effects from the published v1.2.0 compendium.

## Reproduction

1. Use the [v1.2.0 archive](https://doi.org/10.5281/zenodo.24406697) or the [primary analysis at commit 10e9270](https://github.com/Ghostneuron/DNMT3L_Neuron/tree/10e927060fa2bff3140cdfce0e6664c26fe49041/docs/brief_report). Run its input preparation and keep the cache. Install the dependencies from the primary analysis `requirements.txt`.
2. Download the [GSE58885 normalized series matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE58nnn/GSE58885/matrix/GSE58885_series_matrix.txt.gz). The analyzed gzip SHA-256 was `85c4c99814a1c738690ae1f23f8248889f8666161f1e1e8f93b1ac499388a45b`.
3. From this directory, set the path variables to your own inputs and desired output locations and run the commands below. The code ran with Python 3.13.5, NumPy 2.3.4, pandas 2.3.3 and SciPy 1.16.2.

```sh
python run_sensitivity.py "$v12_compendium" "$prepared_cache" "$spiers_matrix" "$output_dir"
python run_followup.py "$v12_compendium" "$prepared_cache" "$spiers_matrix" "$output_dir"
python draw_donor_null.py "$output_dir" "$figure_output_dir"
```

The two analyses use 9,999 Freedman–Lane donor residual permutations with age, sex, and available brain-region terms, permuting donors within sex and region-availability blocks. The Spiers adjustment uses a fixed external reference estimated over 70–112 days post-conception. The adjustment is applied to every permuted DS coefficient vector before recalculating the OR. The central 95% null ranges use the 2.5th and 97.5th percentiles of the permuted ORs. The two-sided critical boundaries use the absolute log OR and the plus-one Monte Carlo P value. The original plan also reports a blocked diagnosis-label test and results for two alternate Spiers windows.

| Cohort | Primary OR (95% donor bootstrap CI) | Primary permutation P | OR boundaries for P < 0.05 | Spiers-adjusted OR (95% donor bootstrap CI) | Adjusted permutation P |
| --- | ---: | ---: | ---: | ---: | ---: |
| GSE74486 | 1.63 (0.66–3.16) | 0.539 | 0.359–2.787 | 1.85 (0.92–2.64) | 0.226 |
| GSE73747 | 1.62 (0.96–2.24) | 0.236 | 0.486–2.058 | 1.87 (1.09–2.53) | 0.123 |

The rejection rule uses ORs strictly outside the boundaries before rounding. The tests depend on exchangeability under the donor blocking and nuisance models. A check of the GSE74486 [GEO metadata](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74486), [archive inventory](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74486/suppl/filelist.txt), pre-ComBat matrix header, and [published sample table](https://doi.org/10.1186/s13059-015-0827-6) found sample identifiers and platform manifests but no mapping to physical array chips. A chip-adjusted donor permutation could not be fit from these files.

## Outputs

- `donor_permutation.tsv` and `spiers_maturation.tsv` give the original sensitivity results, with full window and resampling details.
- `followup_permutation.tsv` gives the primary and Spiers-adjusted tests, null quantiles, two-sided critical boundaries, and tail counts.
- `crude_and_stratified.tsv` pairs the pooled direction-count OR with the Mantel–Haenszel OR for each fetal comparison.
- `donor_null_draws.tsv.gz` contains the OR for each permutation and is used to regenerate the revised Figure 2. The remaining figures use the v1.2.0 results.
