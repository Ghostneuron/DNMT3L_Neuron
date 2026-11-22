# GSE157078 array-position sensitivity checks

The dated plan specifies the three post hoc checks. The main fetal results use the fixed discovery signature and donor workflow. Table S6 of the supplementary PDF reports the new checks. The limma and complete-separation results were carried forward from the preceding signature analysis.

## Reproduction

Use Python 3.13.5 and the pinned requirements from the [Brief Report source](https://github.com/Ghostneuron/DNMT3L_Neuron/tree/4768daae4bb0c5f647602d5a80e3c24c7c93a600/docs/brief_report), or the combined v1.3.0 ZIP, and retain the prepared fetal cache. The two scripts reuse the archived `run_signature_sensitivity.py` from the R5 analysis as the `signature_module` argument. R 4.5.2 reads the Hecker package data.

- [GSE157078 raw signal matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157078/suppl/GSE157078_signals.txt.gz) SHA-256 `1f2f98443b96891557adbc4fa4bee454677c6581edfb6e880dac94c9cbc8a0c5`.
- [GPL21145 EPIC annotation](https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL21nnn/GPL21145/suppl/GPL21145_MethylationEPIC_15073387_v-1-0.csv.gz) SHA-256 `bcb11dc3de47b509864238ca2f2d24a505684624b4a18e362597756da6f3bea1`.
- [GSE250556 raw signal matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE250nnn/GSE250556/suppl/GSE250556_GEO_submission_matrix_signal_intensities.csv.gz) SHA-256 `b04fe6673cf70c8bfc33847390ce92c33d49fa9109c18fc646ea13a1848bb3f4`.
- [GSE250556 SeSAMe matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE250nnn/GSE250556/suppl/GSE250556_GEO_submission_matrix_processed.csv.gz) SHA-256 `9ee58d5b2e54f117bfccef8b183418aa68136b85c589f71db43cf2883ca0967f`.
- [GSE250556 sample SOFT](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250556&targ=gsm&form=text&view=brief) gives subject, pooling, and chip placement. [Butler et al.](https://doi.org/10.1371/journal.pone.0326337) identify one excluded QC outlier.
- [Hecker et al. source](https://github.com/julianhecker/SeffCovar/tree/3ba58ea7a5a4728c76af1cb859d441fd7dc8bd01/data) provides the `S_high` set and five slide-effect categories defining its scored probe universe.

```sh
python run_chip_checks.py "$discovery_signals" "$epic_annotation" "$signature_module" "$compendium" "$fetal_cache" "$out"
python run_reference_checks.py "$discovery_signals" "$epic_annotation" "$signature_module" "$compendium" "$technical_soft" "$technical_raw" "$technical_sesame" "$hecker_data" "$out"
```

The external chamber contrast used the same row blocks as GSE157078 with fixed pooled-donor and slide effects. The 23 analyzed aliquots came from four blood donors and four chips. Probes needed complete valid values for the respective external matrix. I repeated the top-decile comparison on the same complete probe universe for raw and SeSAMe beta to examine how callability influenced the comparison. The Hecker benchmark estimates between-slide susceptibility. It does not model chamber effects within slides. Chromosome bootstraps used 10,000 draws. The external contrasts provide evidence of positional susceptibility. They do not quantify the contribution of position to the observed neural-culture contrast.

## Outputs

- `omission_signature_counts.tsv` and `omission_fetal_results.tsv` give each rebuilt signature and its fetal ORs, donor intervals, paired ratios, and permutation P values.
- `within_group_row_trends.tsv` and `row_beta_contrasts.tsv` give sign-oriented shifts, monotonicity, and the descriptive reference from the 36 within-group orders.
- `position_susceptibility_enrichment.tsv` and `external_chamber_direction.tsv` give the matched enrichment and direction comparisons with independent EPIC data.
