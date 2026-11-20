"""Estimate sign concordance and specificity with paired donor resampling."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, vstack


N_BOOT = 2500
N_NULL = 2000
SEED = 20261120


def load_cohort(cache, cohort, region, use_series):
    meta = pd.read_csv(cache / f"{cohort}_sample_manifest.tsv", sep="\t", index_col=0)
    filename = f"{cohort}_series_beta.pkl" if use_series else f"{cohort}_beta.pkl"
    beta = pd.read_pickle(cache / filename)
    if region:
        meta = meta.loc[meta.region.eq(region)]
    observations = []
    values = {}
    for (donor, region), arrays in meta.groupby(["donor", "region"], sort=False):
        key = f"{donor}/{region}"
        values[key] = beta[arrays.index].mean(axis=1)
        record = arrays.iloc[0].to_dict()
        record.update({"region": region, "donor": donor, "observation": key})
        observations.append(record)
    beta = pd.DataFrame(values).dropna()
    beta = beta.loc[beta.ge(0).all(axis=1) & beta.le(1).all(axis=1)]
    return pd.DataFrame(observations).set_index("observation"), beta


def designs(meta):
    terms = {"intercept": np.ones(len(meta)), "DS": meta.dx.to_numpy(dtype=float)}
    if meta.sex.nunique() > 1:
        terms["male"] = meta.sex.eq("male").to_numpy(dtype=float)
    if meta.region.nunique() > 1:
        terms["temporal"] = meta.region.eq("temporal").to_numpy(dtype=float)
    if "age" in meta.columns and meta.age.notna().all():
        terms["age"] = (meta.age - meta.age.mean()).to_numpy()
    return np.column_stack(list(terms.values())), ", ".join(k for k in terms if k not in ("intercept", "DS"))


def fit_weights(meta, rng, draws):
    groups = meta.groupby("donor", sort=False)
    donors = [(key, block.dx.iloc[0], "+".join(sorted(block.region))) for key, block in groups]
    by_label = pd.DataFrame(donors, columns=["donor", "dx", "region_pattern"]).set_index("donor")
    design, terms = designs(meta)
    n = len(donors)
    member = by_label.index.get_indexer(meta.donor)
    unit = np.zeros((n, 1 + draws + n))
    unit[:, 0] = 1
    unit[:, 1 + draws:] = 1 - np.eye(n)
    for key, block in by_label.groupby(["dx", "region_pattern"]):
        index = by_label.index.get_indexer(block.index)
        probabilities = np.full(len(index), 1 / len(index))
        unit[index, 1:1+draws] = rng.multinomial(len(index), probabilities, size=draws).T
    weights = np.empty((len(meta), unit.shape[1]))
    for index, vector in enumerate(unit.T):
        weighted = design.T * vector[member]
        gram = weighted @ design
        if np.linalg.matrix_rank(gram) != design.shape[1]:
            if index == 0:
                raise ValueError("Rank deficient observed design")
            weights[:, index] = np.nan
        else:
            weights[:, index] = np.linalg.solve(gram, weighted)[1]
    return weights, by_label, terms


def stratum_codes(probes):
    key = pd.MultiIndex.from_frame(probes[["CpG_chrm", "island", "designType", "mean_bin"]])
    return pd.factorize(key, sort=True)[0]


def make_panels(probes, n_strata):
    """Build sparse signed weights, including both jointly balanced signatures."""
    n = len(probes)
    strata = probes.stratum.to_numpy()
    l = probes.l_delta.to_numpy()
    a = probes.a_delta.to_numpy()
    names = []
    blocks = []
    summary = {}
    matched_sets = []

    def add(name, h, d):
        h_idx, d_idx = np.flatnonzero(h), np.flatnonzero(d)
        signed = csr_matrix(
            (np.r_[h[h_idx], d[d_idx]],
             (np.r_[strata[h_idx] * 2, strata[d_idx] * 2 + 1], np.r_[h_idx, d_idx])),
            shape=(n_strata * 2, n),
        )
        names.append(name)
        blocks.append(signed)
        summary[name] = (float(h.sum()), float(d.sum()))

    selected_l = np.isfinite(l)
    add("DNMT3L", (l > 0).astype(float), (l < 0).astype(float))
    null_weights = np.zeros((4, n))
    for s, positions in probes.groupby("stratum", sort=True).indices.items():
        target = positions[selected_l[positions]]
        candidates = positions[~selected_l[positions]]
        hyper, hypo = np.count_nonzero(l[target] > 0), np.count_nonzero(l[target] < 0)
        if len(target) == 0 or len(candidates) < len(target):
            continue
        null_weights[0, target] = l[target] > 0
        null_weights[1, target] = l[target] < 0
        null_weights[2, candidates] = hyper / len(candidates)
        null_weights[3, candidates] = hypo / len(candidates)
        matched_sets.append((s, candidates, hyper, hypo))
    add("DNMT3L_matchable", *null_weights[:2])
    add("matched_CpG_expected", *null_weights[2:])

    # Equivalent to drawing equal numbers without replacement within each
    # stratum, direction and perturbation-magnitude bin, then averaging counts.
    contrast_weights = np.zeros((4, n))
    bin_l = np.digitize(np.abs(l), [0.10, 0.15])
    bin_a = np.digitize(np.abs(a), [0.10, 0.15])
    shared = probes.a_callable.to_numpy()
    for positions in probes.groupby("stratum", sort=True).indices.values():
        positions = positions[shared[positions]]
        for sign_idx, sign in enumerate((1, -1)):
            for strength in range(3):
                a_set = positions[(sign * a[positions] > 0) & (bin_a[positions] == strength)]
                l_set = positions[(sign * l[positions] > 0) & (bin_l[positions] == strength)]
                size = min(len(a_set), len(l_set))
                if size:
                    contrast_weights[sign_idx, l_set] = size / len(l_set)
                    contrast_weights[2 + sign_idx, a_set] = size / len(a_set)
    add("DNMT3L_balanced", *contrast_weights[:2])
    add("DNMT3A_balanced", *contrast_weights[2:])
    return names, vstack(blocks, format="csr"), summary, matched_sets


def mh_or(h_plus, h_minus, l_plus, l_minus):
    total = h_plus + h_minus + l_plus + l_minus
    divisor = np.where(total > 0, total, 1)
    numerator = ((h_plus * l_minus) / divisor).sum(axis=-2)
    denominator = ((h_minus * l_plus) / divisor).sum(axis=-2)
    return np.divide(numerator, denominator, out=np.full_like(numerator, np.nan), where=denominator > 0)


def count_effects(beta, weights, panel, n_strata, n_panels):
    estimates = np.full((n_panels, weights.shape[1]), np.nan)
    sensitivities = np.full((2, weights.shape[1]), np.nan)
    observed_counts = None
    observed_delta = None
    for begin in range(0, weights.shape[1], 48):
        end = min(weights.shape[1], begin + 48)
        delta = beta @ weights[:, begin:end].astype("float32")
        active = np.isfinite(weights[:, begin:end]).all(axis=0)
        # Float32 products can leave tiny residuals at exact ties.
        tol = 1e-7
        pos = (panel @ (delta > tol)).reshape(n_panels, n_strata, 2, -1)
        neg = (panel @ (delta < -tol)).reshape(n_panels, n_strata, 2, -1)
        estimates[:, begin:end] = np.where(active, mh_or(pos[:, :, 0], neg[:, :, 0], pos[:, :, 1], neg[:, :, 1]), np.nan)
        if begin == 0:
            observed_counts = np.stack((pos[:, :, 0, 0], neg[:, :, 0, 0], pos[:, :, 1, 0], neg[:, :, 1, 0]), axis=-1)
            observed_delta = delta[:, 0].copy()
        for i, cutoff in enumerate((0.01, 0.02)):
            sig = panel[:2 * n_strata]
            sp = (sig @ (delta > cutoff)).reshape(n_strata, 2, -1)
            sn = (sig @ (delta < -cutoff)).reshape(n_strata, 2, -1)
            sensitivities[i, begin:end] = np.where(active, mh_or(sp[:, 0], sn[:, 0], sp[:, 1], sn[:, 1]), np.nan)
    return estimates, sensitivities, observed_counts, observed_delta


def matched_null(observed_delta, matched_sets, rng, draws):
    counts = np.zeros((4, len(matched_sets), draws))
    for j, (_, candidates, hyper, hypo) in enumerate(matched_sets):
        calls = observed_delta[candidates]
        for i in range(draws):
            sample = rng.choice(calls, hyper + hypo, replace=False)
            h, l = sample[:hyper], sample[hyper:]
            counts[:, j, i] = np.count_nonzero(h > 1e-7), np.count_nonzero(h < -1e-7), np.count_nonzero(l > 1e-7), np.count_nonzero(l < -1e-7)
    return mh_or(*counts)


def result_row(label, statistic, estimate, draws, info):
    boot = estimate[1:draws+1]
    valid = np.isfinite(boot)
    low, high = np.quantile(boot[valid], [0.025, 0.975]) if valid.any() else (np.nan, np.nan)
    jack = estimate[draws+1:]
    return {
        "context": label, "statistic": statistic, "estimate": estimate[0],
        "ci_low": low, "ci_high": high, "bootstrap_valid": int(valid.sum()),
        "jackknife_min": np.nanmin(jack) if np.isfinite(jack).any() else np.nan,
        "jackknife_max": np.nanmax(jack) if np.isfinite(jack).any() else np.nan,
        **info,
    }


def run(cache, annotation, dnmt3l, dnmt3a, cohort, region, null_draws, draws, use_series):
    label = f"{cohort}_{region}" if region else cohort
    if use_series:
        label += "_deposited_beta"
    meta, beta = load_cohort(cache, cohort, region, use_series)
    probes = annotation.loc[beta.index].copy()
    probes["mean_bin"] = np.digitize(beta.mean(axis=1), [0.20, 0.40, 0.60, 0.80])
    probes["stratum"] = stratum_codes(probes)
    probes["l_delta"] = dnmt3l.reindex(probes.index)
    probes["a_callable"] = probes.index.isin(dnmt3a.index)
    common_sign = (np.sign(dnmt3a.delta_a) == np.sign(dnmt3a.delta_b)) & dnmt3a.delta.abs().ge(0.05)
    probes["a_delta"] = dnmt3a.delta.where(common_sign).reindex(probes.index)
    n_strata = probes.stratum.nunique()
    names, panel, sizes, matched_sets = make_panels(probes, n_strata)
    rng = np.random.default_rng(SEED)
    weights, donors, terms = fit_weights(meta, rng, draws)
    estimates, sensitivities, counts, observed_delta = count_effects(
        np.ascontiguousarray(beta, dtype="float32"), weights, panel, n_strata, len(names)
    )
    print(label, "donors", donors.dx.value_counts().to_dict(), "CpGs", len(probes), "terms", terms, flush=True)
    h_plus, h_minus, l_plus, l_minus = counts[0].sum(axis=0)
    n = h_plus + h_minus + l_plus + l_minus
    agreement = (h_plus + l_minus) / n
    expected = ((h_plus + h_minus) * (h_plus + l_plus) + (l_plus+l_minus) * (h_minus+l_minus)) / n**2
    common = {"ds_donors": int(donors.dx.sum()), "control_donors": int(len(donors)-donors.dx.sum()),
              "region_observations": len(meta), "eligible_cpgs": len(probes), "adjusted_covariates": terms}
    rows = []
    cpg_counts = []
    for i, name in enumerate(names):
        rows.append(result_row(label, name + "_OR", estimates[i], draws, {
            **common, "n_hyper": sizes[name][0], "n_hypo": sizes[name][1],
        }))
        cpg_counts.append({
            "context": label, "panel": name,
            "hyper_hyper": counts[i, :, 0].sum(), "hyper_hypo": counts[i, :, 1].sum(),
            "hypo_hyper": counts[i, :, 2].sum(), "hypo_hypo": counts[i, :, 3].sum(),
        })
    for i, cutoff in enumerate((0.01, 0.02)):
        rows.append(result_row(label, f"DNMT3L_OR_absDSge{cutoff:.2f}", sensitivities[i], draws, common))
    rows.append(result_row(label, "DNMT3L_vs_matched_CpGs_ratio_of_ORs", estimates[1] / estimates[2], draws, common))
    rows.append(result_row(label, "DNMT3L_vs_DNMT3A_ratio_of_ORs", estimates[3] / estimates[4], draws, common))
    summary = {"context": label, "same_direction": agreement, "global_sign_expected": expected, "ties": sum(sizes["DNMT3L"]) - n,
               "joint_DNMT3A_background_cpgs": int(probes.a_callable.sum()), "DNMT3A_signature_eligible_cpgs": int(probes.a_delta.notna().sum())}
    if null_draws:
        null = matched_null(observed_delta, matched_sets, rng, null_draws)
        summary.update({
            "matched_CpG_null_draws": null_draws,
            "matched_CpG_null_OR_2.5%": np.nanquantile(null, .025),
            "matched_CpG_null_OR_50%": np.nanquantile(null, .5),
            "matched_CpG_null_OR_97.5%": np.nanquantile(null, .975),
            "matched_CpG_upper_tail": (1 + np.count_nonzero(null >= estimates[1, 0])) / (1 + np.isfinite(null).sum()),
        })
    return rows, cpg_counts, summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("cache", type=Path)
    parser.add_argument("signature", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--bootstrap", type=int, default=N_BOOT)
    parser.add_argument("--matched-nulls", type=int, default=N_NULL)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--cortex-series-only", action="store_true")
    modes.add_argument("--cerebrum-only", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    dnmt3l = pd.read_csv(args.signature, sep="\t", index_col="probe_id")["delta_beta"]
    assert dnmt3l.index.is_unique and dnmt3l.abs().ge(.0499999).all()
    dnmt3a = pd.read_pickle(args.cache / "dnmt3a_effects.pkl")
    annotation = pd.read_pickle(args.cache / "eligible_annotations.pkl")
    rows, counts, summaries = [], [], []
    contexts = [("GSE74486", None), ("GSE73747", None), ("GSE73747", "frontal"), ("GSE73747", "temporal")]
    if args.cortex_series_only:
        contexts = [("GSE73747", None)]
    if args.cerebrum_only:
        contexts = [("GSE74486", None)]
    for cohort, region in contexts:
        result, observed, summary = run(
            args.cache, annotation, dnmt3l, dnmt3a, cohort, region, args.matched_nulls, args.bootstrap, args.cortex_series_only
        )
        rows.extend(result)
        counts.extend(observed)
        summaries.append(summary)
        print(pd.DataFrame(result).loc[:, ["statistic", "estimate", "ci_low", "ci_high"]].to_string(index=False), flush=True)
    for name, records in [("estimates", rows), ("direction_counts", counts), ("null_and_coverage", summaries)]:
        pd.DataFrame(records).to_csv(args.output / (name + ".tsv"), sep="\t", index=False, float_format="%.8g", na_rep="NA")


if __name__ == "__main__":
    main()
