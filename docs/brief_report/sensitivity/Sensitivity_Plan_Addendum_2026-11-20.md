# Post hoc sensitivity plan addendum, 20 November 2026

This addendum follows review of the results from the dated sensitivity plan. It was recorded before computing the tests specified below. The existing estimates, cohort definitions, probe sets, and Spiers age window are fixed.

- Repeat the original 9,999 Freedman–Lane residual permutations within the recorded donor blocks, using the original seed. For the age-adjusted test, apply the same fixed external Spiers projection to each permuted DS coefficient vector, then evaluate the stratified direction odds ratio on the same age-callable CpGs. Use the absolute log odds ratio and the plus-one P value from the original donor test.
- Report the exact Monte Carlo rejection boundary for the primary two-sided P < 0.05 test, in addition to the central 95% permutation range. Plot the two primary donor permutation distributions and the observed odds ratios.
- Compute crude odds ratios from the pooled 2 × 2 direction counts and place them alongside their stratified estimates. Check the public GSE74486 metadata for sample-to-array chip identifiers, retaining the original model if a mapping is unavailable.
