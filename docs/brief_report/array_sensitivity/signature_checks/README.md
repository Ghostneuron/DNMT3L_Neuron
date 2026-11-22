# Post hoc DNMT3L signature sensitivity

The archived GSE157078 signature was reconstructed before applying the new filters. The archived discovery signs and the two primary fetal ORs were reproduced. The dated plan records what was specified before running the checks.

## Reproduction

Use the pinned Python environment and prepared fetal cache from the [Brief Report source at commit 4768daa](https://github.com/Ghostneuron/DNMT3L_Neuron/tree/4768daae4bb0c5f647602d5a80e3c24c7c93a600/docs/brief_report) or the extracted v1.3.0 compendium. This run also used R 4.5.2 with limma 3.64.3.

- [GSE157078 signals](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157078/suppl/GSE157078_signals.txt.gz), SHA-256 `1f2f98443b96891557adbc4fa4bee454677c6581edfb6e880dac94c9cbc8a0c5`.
- [GPL21145 EPIC annotation](https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL21nnn/GPL21145/suppl/GPL21145_MethylationEPIC_15073387_v-1-0.csv.gz), SHA-256 `bcb11dc3de47b509864238ca2f2d24a505684624b4a18e362597756da6f3bea1`.
- [GSE157078 sample SOFT](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157078&targ=gsm&form=text&view=brief) and [GEO IDAT inventory](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157078/suppl/filelist.txt) provide sample labels and array positions.

```sh
python run_signature_sensitivity.py "$signals" "$epic_annotation" "$soft" "$geo_inventory" "$compendium" "$cache" "$output"
```

The two-sided donor test uses the absolute log OR with zero ORs counted as extreme. For the sparse BH subset, four GSE74486 donor bootstrap draws had no defined OR. Percentile summaries use the 2,496 defined draws. One GSE73747 BH permutation had no defined OR, and its permutation P value counts that draw in the tail. These choices preserve the distinction between finite zero odds ratios and undefined ratios.

## Outputs

- `array_positions.tsv` matches the six array labels to GEO accessions, the chip ID, and the slots.
- `label_split_counts.tsv` contains all ten partitions. Sample 4 is always in the group labelled putative overexpression, which anchors sign for each unlabeled partition.
- `signature_counts.tsv` gives discovery-wide counts. `discovery_filter_statistics.tsv.gz` gives the limma statistics and separation indicator for the original 15,602 CpGs. Benjamini–Hochberg adjustment precedes intersection with the original signature.
- `donor_signature_sensitivity.tsv` gives each fetal OR, paired bootstrap ratio to the primary OR, and blocked donor permutation results. These uncertainties condition on the signature selected from the discovery arrays.
