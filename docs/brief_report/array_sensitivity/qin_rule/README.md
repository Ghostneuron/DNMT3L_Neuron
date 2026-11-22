# Qin et al. 2021 rule and preprocessing sensitivity

Qin et al. used nominal pooled Student t-test P ≤ 0.05 and an absolute mean
beta difference of at least 0.10. The deposited GSE157078 series beta
values reproduce the values for their 347 published Table S2 probes. Their
full legacy probe mask, exact gene-union construction and historical
HOMER background are unavailable. The figure and motif tables describe
reconstructions and sensitivity checks.

Place these public inputs in a work directory `QWORK`. The publisher
workbooks are accessed at their original URLs.

- `GSE157078_series_matrix.txt.gz` from [GSE157078 GEO](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157078/matrix/).
- [Table S2](https://ars.els-cdn.com/content/image/1-s2.0-S0301008221000423-mmc4.xlsx)
  as `mmc4.xlsx`, [Table S3](https://ars.els-cdn.com/content/image/1-s2.0-S0301008221000423-mmc5.xlsx)
  as `mmc5.xlsx`, and [Table S5](https://ars.els-cdn.com/content/image/1-s2.0-S0301008221000423-mmc7.xlsx)
  as `mmc7.xlsx`.
- [Zhou hg19 EPIC manifest](https://zhouserver.research.chop.edu/InfiniumAnnotation/EPIC/EPIC.hg19.manifest.tsv.gz)
  as `zhou_epic.tsv.gz`.
- For the motif runs, [UCSC hg19 2bit](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit)
  as `hg19.2bit` and the [Zhou hg19 450K manifest](https://zhouserver.research.chop.edu/InfiniumAnnotation/HM450/HM450.hg19.manifest.tsv.gz)
  in `FETAL_INPUTS/HM450_Zhou_manifest.tsv.gz`.

Use the adjacent raw IDAT instructions to generate `IDAT_BETA`. Run:

```sh
Rscript preprocess.R "$QWORK" "$IDATS"
Rscript preprocess_background.R "$QWORK" "$IDATS"
python probe_tests.py "$QWORK" "$GEO_SIGNALS" "$IDAT_BETA"
python region_tests.py "$QWORK"
python panel_checks.py "$QWORK"
python correlation_scope.py "$QWORK"
python bumps_in.py "$QWORK"
Rscript bumps.R "$QWORK"
python bumps_summary.py "$QWORK"
python motif_inputs.py "$QWORK" "$FETAL_INPUTS"
python extra_motif_targets.py "$QWORK"
```

HOMER 5.1 runs from `$QWORK/motif` with target FASTA windows named by
the set key. The fixed seed samples 60,000 windows from the 450K/EPIC union.
Use the SNP/mapping filter documented in `probe_tests.py` and the source
chromosome/coordinates for published sets. For a published Figure 1A or
ssNoob+SWAN set run:

```sh
findMotifs.pl "$set.fa" fasta "$set" -fastaBg background.fa -mset vertebrates -b -nogo -nofacts -p 2 -S 15
```

For the smaller-de-novo sensitivity use the first 12,000 background
FASTA records, `-noknown -S 8`, and a distinct output directory. Use
`-nomotif` for library-only comparisons. The pretargeted STAT3 and related
checks run `targeted_homer.py "$QWORK"`, then, before the HOMER run
cleans its intermediates, use

```sh
homer2 known -s "$set/seq.tsv" -g "$set/group.adj" -m targeted.motifs -stat binomial -p 2 -mout "$set/targeted_results.motifs"
```

Reproduce compact motif outputs with `summarize_homer.py` and
`summarize_targeted.py`. Cases without a reconstructed 5% BH result are
not evidence of independent TF occupancy. See [Qin et al.](https://doi.org/10.1016/j.pneurobio.2021.102028).
