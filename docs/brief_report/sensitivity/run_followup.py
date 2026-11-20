"""Evaluate the 2026-11-20 post hoc addendum against the frozen v1.2.0 inputs."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import run_sensitivity as sensitivity


COHORTS = ("GSE74486", "GSE73747")
WINDOW = "overlap_70_112"
TAIL_TOL = 1e-10  # Same tail comparison as the original donor test.


def one_test(original, name, probes, beta, weights, observed_weights, slopes=None):
    if slopes is not None:
        slopes = slopes.reindex(probes.index).dropna()
        probes = probes.loc[slopes.index]
        centered = slopes - slopes.groupby(probes.stratum, sort=False).transform("mean")
        predictor = centered.to_numpy()
        projection = predictor @ beta.loc[probes.index].to_numpy() / (predictor @ predictor)
    signed, panel, n_strata = sensitivity.signature_panel(probes)
    y = np.ascontiguousarray(beta.loc[signed.index], dtype="float64")
    if slopes is not None:
        y -= centered.loc[signed.index].to_numpy()[:, None] * projection
    observed = sensitivity.or_from_weights(original, y, observed_weights, panel, n_strata)[0]
    draws = sensitivity.or_from_weights(original, y, weights, panel, n_strata)
    assert len(draws) == sensitivity.N_PERM and np.all(np.isfinite(draws) & (draws > 0))
    score = np.abs(np.log(draws))
    tail = int(np.count_nonzero(score >= abs(np.log(observed)) - TAIL_TOL))
    # P = (1 + tail) / 10000. Strictly less than .05 requires at most 498 draws in the tail.
    rank = len(draws) - 499
    critical_log_or = np.partition(score, rank)[rank] + TAIL_TOL
    assert np.count_nonzero(score >= np.nextafter(critical_log_or - TAIL_TOL, np.inf)) <= 498
    low, median, high = sensitivity.interval(draws)
    return (
        {"contrast": name, "or": observed, "draws": len(draws), "tail": tail,
         "p_two_sided": (1 + tail) / (1 + len(draws)), "null_low": low,
         "null_median": median, "null_high": high, "critical_low": np.exp(-critical_log_or),
         "critical_high": np.exp(critical_log_or), "signed_probes": len(signed)},
        draws,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("frozen", type=Path)
    parser.add_argument("cache", type=Path)
    parser.add_argument("spiers", type=Path)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    original = sensitivity.original_code(args.frozen / "analyze.py")
    signature = pd.read_csv(args.frozen / "GSE157078_candidate_dmcs.tsv", sep="\t",
                            index_col="probe_id").delta_beta
    annotation = pd.read_pickle(args.cache / "eligible_annotations.pkl")
    estimates = pd.read_csv(args.frozen / "estimates.tsv", sep="\t")
    counts = pd.read_csv(args.frozen / "direction_counts.tsv", sep="\t").query('panel == "DNMT3L"')
    crude = counts.assign(
        crude_or=counts.hyper_hyper * counts.hypo_hypo / (counts.hyper_hypo * counts.hypo_hyper)
    )
    paired = crude.merge(estimates.query('statistic == "DNMT3L_OR"')[["context", "estimate"]],
                         on="context", validate="one_to_one").rename(columns={"estimate": "stratified_or"})
    paired.drop(columns="panel").to_csv(args.out / "crude_and_stratified.tsv", sep="\t",
                                       index=False, float_format="%.9g")
    slopes, _ = sensitivity.spiers_slopes(args.spiers)
    summaries = []
    distributions = []
    for cohort in COHORTS:
        metadata, beta = original.load_cohort(args.cache, cohort, None, False)
        probes = sensitivity.probe_frame(original, beta, annotation, signature)
        design, _ = original.designs(metadata)
        blocks = sensitivity.donor_blocks(metadata)
        weights = sensitivity.residual_weights(design, np.delete(design, 1, axis=1), blocks)
        observed = np.linalg.pinv(design)[1, :, None]
        for name, reference in (("primary", None), ("spiers_adjusted", slopes[WINDOW])):
            summary, scores = one_test(original, name, probes, beta, weights, observed, reference)
            if name == "primary":
                expected = paired.loc[paired.context.eq(cohort), "stratified_or"].iloc[0]
                assert np.isclose(summary["or"], expected, rtol=2e-5)
            summaries.append({"cohort": cohort, **summary})
            distributions.append(pd.DataFrame({"cohort": cohort, "contrast": name,
                                               "draw": np.arange(1, len(scores) + 1), "or": scores}))
    pd.DataFrame(summaries).to_csv(args.out / "followup_permutation.tsv", sep="\t", index=False,
                                  float_format="%.9g")
    pd.concat(distributions, ignore_index=True).to_csv(
        args.out / "donor_null_draws.tsv.gz", sep="\t", index=False, float_format="%.17g",
        compression={"method": "gzip", "mtime": 0},
    )


if __name__ == "__main__":
    main()
