# Brief Report reproducibility package

This analysis tests whether the archived fetal neural DNMT3L overexpression signature aligns with the direction of DS methylation differences. The endpoint is the Mantel–Haenszel odds ratio (OR) for a positive fitted DS effect at DNMT3L hypermethylated CpGs compared with DNMT3L hypomethylated CpGs, conditional on probe strata. The bootstrap intervals propagate uncertainty from the DS donors, conditional on the fixed external signatures and CpG sets.

The [dated post hoc sensitivity directory](sensitivity/) gives the donor permutation and Spiers reference analyses, crude ORs, and the source for the revised Figure 2. The commands below reproduce the archived v1.2.0 primary analysis and original figures.

## Inputs

Save the public downloads below in a local data directory under the indicated names. Extract `DS Suppl Tables GB4.xlsx` from the Mendioroz supplement as `Mendioroz_2015_Tables.xlsx`.

| Local name | Public source |
| --- | --- |
| `GSE74486_precombat.txt.gz` | [GSE74486 normalized pre-ComBat matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74486/suppl/GSE74486_DS_processed_norm_no_batch_correction_no_filtering_matrix_1.txt.gz) |
| `GSE74486_brief.soft` | [GSE74486 sample SOFT](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74486&targ=gsm&form=text&view=brief) |
| `Mendioroz_2015_Tables.xlsx` | [Mendioroz Additional file 2](https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs13059-015-0827-6/MediaObjects/13059_2015_827_MOESM2_ESM.zip) |
| `GSE73747_signal.txt.gz` | [GSE73747 methylated and unmethylated signals](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE73nnn/GSE73747/suppl/GSE73747_unmethylated_methylated.txt.gz) |
| `GSE73747_series_matrix.txt.gz` | [GSE73747 sample metadata and deposited beta values](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE73nnn/GSE73747/matrix/GSE73747_series_matrix.txt.gz) |
| `GSE157078_signals.txt.gz` | [GSE157078 EPIC assay signals](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157078/suppl/GSE157078_signals.txt.gz) |
| `GSE120558_EPIC_matrix.txt.gz` | [GSE120558 human DNMT3A W330R methylation matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120558/matrix/GSE120558-GPL21145_series_matrix.txt.gz) |
| `GPL13534_annotation.soft.gz` | [GPL13534 manufacturer annotation and SOFT record](https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13534/soft/GPL13534_family.soft.gz) |
| `HM450_Zhou_manifest.tsv.gz` | [Zhou HM450 hg19 annotation and mask](https://zhouserver.research.chop.edu/InfiniumAnnotation/HM450/HM450.hg19.manifest.tsv.gz) |

The frozen DNMT3L signature input, `GSE157078_candidate_dmcs.tsv`, contains the probe IDs and beta differences from the project’s GSE157078 discovery analysis. The absolute beta-change cutoff is 0.05 [Qin et al.](https://doi.org/10.1016/j.pneurobio.2021.102028).

From this directory, install the dependencies in `requirements.txt`. The analysis ran with Python 3.13.5. Set the input, cache, and output directory variables to paths of your choice.

```sh
python -m pip install -r requirements.txt
python prepare_inputs.py "$inputs_dir" "$cache_dir"
python analyze.py "$cache_dir" GSE157078_candidate_dmcs.tsv "$output_dir"
python prepare_inputs.py "$inputs_dir" "$cache_dir" --series-only
python analyze.py "$cache_dir" GSE157078_candidate_dmcs.tsv "$series_output_dir" --cortex-series-only
python draw_figures.py . "$figure_output_dir"
```

The bundled tables combine the main and deposited-beta sensitivity outputs. Compare the outputs by context and statistic. The figures use the bundled tables. The analysis and matching seed is 20261120.

## Analysis details

- GSE74486 uses normalized pre-ComBat beta values, with sex and gestational age in the per-CpG regression. The one reported age range, 12–13 weeks, is represented by its midpoint. The supplement and GEO sample metadata support the 8 DS and 6 control donor identities. Array batch was unavailable in the sample files used here. [Mendioroz et al.](https://doi.org/10.1186/s13059-015-0827-6)
- GSE73747 uses `M / (M + U + 100)` and fits sex, gestational age and brain region. Within-region fits omit the region term. The 59 fetal arrays included 17 pairs with 100% genotype agreement across 64–65 jointly called array SNP probes. The next-highest agreement was 67.7%. The 17 pairs comprised 16 frontal–temporal matches and one same-region repeat. The analysis averages the repeat and resamples the 42 donor clusters within diagnosis and region pattern. The deposited beta-value analysis retains the same detection filter. [El Hajj et al.](https://doi.org/10.1080/15592294.2016.1192736)
- The matched CpG panels preserve hyper and hypo counts within chromosome, CpG island context, design type and the label-blind cohort mean beta bin. The 2,000-panel upper tail describes the observed donors held fixed. The primary package reports donor bootstrap intervals. Donor permutation results are available in the [sensitivity directory](sensitivity/).
- The DNMT3A comparison uses autosomal sites shared by the three assays. The two W330R donors must agree in direction, with mean absolute beta change of at least 0.05 compared with the three healthy fibroblast donors. The comparison balances genomic strata, sign and perturbation magnitude. The ratio of ORs and interval are calculated within the same 2,500 donor draws. The fibroblast comparator differs in cell type and perturbation from the neural DNMT3L model. The estimates do not account for uncertainty in the external perturbation datasets. [Heyn et al.](https://doi.org/10.1038/s41588-018-0274-x)

## Results

| Cohort | DNMT3L direction OR (95% CI) | OR ratio versus balanced DNMT3A (95% CI) |
| --- | ---: | ---: |
| GSE74486 fetal cerebrum, 8 DS and 6 controls | 1.63 (0.66–3.16) | 2.15 (0.81–3.91) |
| GSE73747 fetal cortex, 16 DS and 26 controls | 1.62 (0.96–2.24) | 2.10 (1.16–2.88) |
| GSE73747 deposited beta sensitivity | 1.93 (0.97–2.48) | 2.42 (1.16–3.24) |

In GSE74486 and GSE73747, the observed DNMT3L OR exceeded all 2,000 random matched-CpG panels (Monte Carlo upper-tail probability 1/2,001 in each cohort, conditional on the sampled donors). The donor bootstrap intervals for the primary DNMT3L OR include 1, and the corresponding donor permutation P values are 0.539 and 0.236. The balanced DNMT3L and DNMT3A ORs were 1.63 and 0.76 for GSE74486, 1.57 and 0.75 for GSE73747, and 1.87 and 0.77 for the deposited-beta sensitivity. These comparisons leave neural-context specificity unresolved.

The full estimates table also provides frontal and temporal results, thresholds on the magnitude of the DS effect, and leave-one-donor-out ranges. The regional results share donors. The temporal regression was rank deficient in 7 of 2,500 bootstrap draws. The corresponding intervals use 2,493 valid draws. Differences in tissue dissection, gestational age and cell composition limit direct comparisons between the two cohorts.

## Version and licenses

Version 1.2.0 accompanies *Directional concordance of DNMT3L-associated methylation in fetal Down syndrome brain*. The previous v1.1.0 compendium documents the earlier study and is independently citable at [10.5281/zenodo.20770052](https://doi.org/10.5281/zenodo.20770052). See [licenses](LICENSES.md) for reuse terms.
