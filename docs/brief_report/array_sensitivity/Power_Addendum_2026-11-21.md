# Exploratory donor resampling and array sensitivity

Recorded on 21 November 2026 before computing these tables. Earlier donor
concordance, chip QC and Qin-rule results were known. This addendum does not
alter the frozen analysis.

- Use the original 865,743-probe GEO detection support for all IDAT summaries.
  Tabulate gain/loss at 0.05–<0.10, 0.10–<0.15 and ≥0.15 and by EPIC probe
  chemistry and island, shore, shelf or open sea.
- Use the prepared per-CpG GEO donor-region beta matrices and frozen
  covariates and signature support. For N = 28, 56, 112 and 224, generate 50
  pseudo-cohorts at each N by sampling donor IDs with replacement within
  diagnosis × sex × region-pattern cells. Use largest-remainder allocation
  with each original nonempty cell represented. Keep all regions from each
  selected donor together and treat copies independently in inner draws.
- Pair the frozen and ssNoob signatures on the same panels, permuted residuals
  and bootstrap draws. Split 9,999 Freedman–Lane permutations and 2,500
  diagnosis × region-pattern donor-bootstrap draws across the 50 panels in
  each scenario. The permutation blocks are sex × region pattern. Seeds
  20261120 (bootstrap), 2026112002 (permutations) and 2026112101 (outer panels)
  are fixed. Keep baseline-beta bins from the source cohort.
- Estimate two-sided, 0.05 power by shifting the within-panel centered
  empirical permutation distribution by log(1.6), log(2.0) or log(2.4) and
  applying that panel's empirical |log OR| cutoff. For the Monte Carlo power
  interval, use the standard error across the 50 pseudo-cohort powers and a
  t critical value. Report pooled permutation SD, mean within-panel bootstrap
  SE on log OR and the smallest tested N at 80%. Power models a shift of the
  statistic with unchanged noise shape. The empirical donors do not model
  additional between-donor heterogeneity in a new study.
