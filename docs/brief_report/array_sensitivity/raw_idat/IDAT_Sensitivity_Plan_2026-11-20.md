# Post hoc GSE157078 raw IDAT sensitivity

Recorded 20 November 2026 before reading the raw IDAT results. The earlier effect-size signature, fetal association and chip-position checks were already known.

Verify the six GSM accessions, archived Z/D identities and IDAT Sentrix headers. Use raw control metrics, negative-control backgrounds, median M and U, and detection failure fractions at P > 0.01 for each array. Plot the metrics by Sentrix row. Treat the row pattern as descriptive with one array per row and fully segregated treatments.

Run single-sample Noob with dye-bias correction, retaining the archived beta offset and 0.05 group-mean difference cutoff. Also run the same Noob background correction with dye correction disabled to distinguish the incremental change due to dye correction. Compare archived raw, Noob and dye-corrected Noob using identical callable support for each contrast, applying the original all-six GEO quality filter and independent raw-IDAT detection filter as separate analyses. Rebuild signatures with each discovery array omitted on the same support. Compare sign, gain share, overlap, and the fetal donor ORs using the fixed donor models, 2,500 paired donor bootstrap draws and original 9,999 donor residual permutations.

Preprocessing changes are evidence about sensitivity to measurement choices. The design cannot identify the effect of treatment separately from row or order.
