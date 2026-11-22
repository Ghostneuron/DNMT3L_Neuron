# Post hoc DNMT3L array and donor sensitivity

These dated analyses accompany the fetal Down syndrome Brief Report. The
[previous Zenodo record](https://doi.org/10.5281/zenodo.24406697) contains
the frozen primary analysis. Use the adjacent Brief Report directory
(the archive root in the combined Zenodo ZIP) as `COMPENDIUM`.

| Directory | Analyses |
| --- | --- |
| [signature_checks](signature_checks/) | Label splits, complete separation, limma and fetal OR sensitivity |
| [chip_position](chip_position/) | Leave-one-array-out, within-group trends and public technical-reference comparisons |
| [raw_idat](raw_idat/) | Sample mapping, control channel ratios, raw IDAT QC, Noob and ssNoob |
| [qin_rule](qin_rule/) | Pooled Student t-test, effect-size cutoff, published figure and motif sensitivity |
| [results](results/) | Effect-size, context, observed donor OR and pseudo-donor power tables |

The [21 November addendum](Power_Addendum_2026-11-21.md) describes
the newest exploratory analyses. Earlier post hoc plans are dated 20 November
2026. Prepared beta matrices are regenerated from GEO using `prepare_inputs.py`
in the main compendium. IDAT files, publisher supplements and reference genome
are obtained from the original public sources.

```sh
python run_idat_bins.py "$IDAT_BETA" "$GEO_SIGNALS" "$EPIC_MANIFEST" "$OUT"
python run_power.py "$COMPENDIUM" "$FETAL_CACHE" "$OUT/ssnoob_signature.tsv" "$OUT" --cohort GSE74486
python run_power.py "$COMPENDIUM" "$FETAL_CACHE" "$OUT/ssnoob_signature.tsv" "$OUT" --cohort GSE73747
```

The beta denominator offset is 100 and the Noob intensity offset is 15.
Bins reconcile exactly with the original support totals. Crude ORs pool the
probe counts after covariate adjustment of the fetal differences. The
Mantel–Haenszel ORs use the frozen CpG strata. Observed intervals use the
2,500 original donor bootstrap draws. The simulated sample sizes each use
50 pseudo-cohorts, 9,999 total Freedman–Lane draws and 2,500 total donor
bootstrap draws. The simulated power models a location shift on log OR.
Repeated source donors can underrepresent population heterogeneity.

Code uses the repository MIT license. The author-generated result tables
use CC BY 4.0. Original data and annotations retain their own terms.
