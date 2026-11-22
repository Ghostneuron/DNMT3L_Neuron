"""Exploratory pooled donor-panel power under a log-OR location shift."""
import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t


BOOT_SEED = 20261120
PERM_SEED = 2026112002
PANEL_SEED = 2026112101
N_PERM = 9999
N_BOOT = 2500
N_PANELS = 50
SIZES = (28, 56, 112, 224)
ALTERNATIVES = (1.6, 2., 2.4)


def load_module(path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def signature(path):
    return pd.read_csv(path, sep="\t", index_col="probe_id").delta_beta


def create_cells(meta):
    grouping = meta.reset_index(drop=True).groupby("donor", sort=False)
    records, cells = [], {}
    for _, block in grouping:
        assert block.dx.nunique() == block.sex.nunique() == 1 and block.region.is_unique
        dx, sex = int(block.dx.iloc[0]), block.sex.iloc[0]
        rows = block.sort_values("region").index.to_numpy()
        pattern = tuple(sorted(block.region))
        key = (dx, sex, pattern)
        cells.setdefault(key, []).append(rows)
        records.append(dx)
    return {key: np.asarray(cells[key]) for key in sorted(cells)}, np.mean(records), len(records)


def sizes_by_cell(cells, prevalence, n):
    counts = {}
    ds = int(np.floor(n * prevalence + .5))
    for dx, total in ((0, n - ds), (1, ds)):
        keys = [key for key in cells if key[0] == dx]
        raw = np.asarray([len(cells[key]) for key in keys], dtype=float)
        expected = total * raw / raw.sum()
        rounded = np.maximum(1, np.floor(expected).astype(int))
        while rounded.sum() < total:
            rounded[np.argmax(expected - rounded)] += 1
        while rounded.sum() > total:
            rounded[np.argmax(np.where(rounded > 1, rounded - expected, -np.inf))] -= 1
        counts.update(zip(keys, rounded))
    assert sum(counts.values()) == n
    return counts


def expand(cells, counts, rng):
    source, permutations, boots, member = [], {}, {}, []
    n_rows, n_people = 0, 0
    for (dx, sex, pattern), choices in cells.items():
        count = counts[(dx, sex, pattern)]
        indices = rng.integers(len(choices), size=count)
        selected = choices[indices]
        source.extend(selected.ravel())
        indices = np.arange(n_rows, n_rows + selected.size).reshape(selected.shape)
        permutations.setdefault((sex, pattern), []).append(indices)
        boots.setdefault((dx, pattern), []).extend(range(n_people, n_people + count))
        member.extend(np.repeat(np.arange(n_people, n_people + count), selected.shape[1]))
        n_rows += selected.size
        n_people += count
    blocks = [np.concatenate(permutations[key]) for key in sorted(permutations)]
    boot_blocks = [np.asarray(boots[key]) for key in sorted(boots)]
    return np.asarray(source), blocks, boot_blocks, np.asarray(member)


def collapse(weights, source, n_rows):
    result = np.zeros((n_rows, weights.shape[1]))
    np.add.at(result, source, weights)
    return np.ascontiguousarray(result)


def expanded_weights(x, source, blocks, bootstrap_blocks, member, n_perm, n_boot, p_rng, b_rng):
    x_pseudo = x[source]
    z = np.delete(x_pseudo, 1, axis=1)
    observed = np.linalg.pinv(x_pseudo)[1]
    inverse_z = np.linalg.pinv(z)
    shuffled = np.empty((len(source), n_perm))
    for positions in blocks:
        for p in range(n_perm):
            order = p_rng.permutation(len(positions))
            shuffled[positions[order].ravel(), p] = observed[positions.ravel()]
    fl = shuffled + z @ (inverse_z @ (observed[:, None] - shuffled))
    assert np.max(np.abs(z.T @ fl)) < 1e-8
    n_donor = int(member.max() + 1)
    multiplicity = np.empty((n_donor, n_boot))
    for inds in bootstrap_blocks:
        multiplicity[inds] = b_rng.multinomial(len(inds), np.full(len(inds), 1/len(inds)), size=n_boot).T
    boot = np.full((len(source), n_boot), np.nan)
    for b, counts in enumerate(multiplicity.T):
        root = np.sqrt(counts[member])
        matrix = x_pseudo * root[:, None]
        if np.linalg.matrix_rank(matrix) == x.shape[1]:
            boot[:, b] = np.linalg.pinv(matrix)[1] * root
    return collapse(fl, source, len(x)), collapse(boot, source, len(x))


def prepare_panel(original, beta, ann, effects):
    chosen = beta.index[beta.index.isin(effects.index)]
    mapping = ann.loc[chosen, ["CpG_chrm", "island", "designType"]].copy()
    mapping["mean_bin"] = np.digitize(beta.loc[chosen].mean(axis=1), [.2, .4, .6, .8])
    strata = original.stratum_codes(mapping)
    group = 2 * strata + effects.loc[chosen].lt(0).to_numpy(dtype=int)
    values = np.ascontiguousarray(beta.loc[chosen], dtype=np.float64)
    return values, group, int(strata.max() + 1)


def assay(panel, weights):
    y, groups, n = panel
    result = np.empty((weights.shape[1], 2))
    for start in range(0, weights.shape[1], 128):
        block = y @ weights[:, start:start + 128]
        for j, delta in enumerate(block.T):
            if not np.isfinite(delta).all():
                result[start+j] = np.nan
                continue
            up = np.bincount(groups[delta > 1e-7], minlength=2*n).reshape(n, 2)
            down = np.bincount(groups[delta < -1e-7], minlength=2*n).reshape(n, 2)
            h_plus, l_plus = up.T
            h_minus, l_minus = down.T
            denom = np.maximum((up + down).sum(axis=1), 1)
            terms = ((h_plus * l_minus / denom).sum(), (h_minus * l_plus / denom).sum())
            crude = (h_plus.sum()*l_minus.sum(), h_minus.sum()*l_plus.sum())
            result[start+j] = [a/b if b else np.nan for a,b in (crude, terms)]
    return result


def slice_draws(n, k, panels):
    return n // panels + int(k < n % panels)


def interval(values):
    return np.quantile(values, (.025, .975))


def compute(original, sensitivity, cache, cohort, signatures, sizes, panels, out):
    region = "cerebrum" if cohort == "GSE74486" else None
    meta,beta = original.load_cohort(cache, cohort, region, False)
    x, _ = original.designs(meta)
    ann = pd.read_pickle(cache / "eligible_annotations.pkl")
    arrays = {name: prepare_panel(original, beta, ann, effects) for name,effects in signatures.items()}
    print(cohort, {key: arr[0].shape for key,arr in arrays.items()}, flush=True)
    cells, prevalence, n_observed = create_cells(meta)
    # These are the same donor draws used by the frozen-size analysis.
    obs_boot, _, _ = original.fit_weights(meta, np.random.default_rng(BOOT_SEED), N_BOOT)
    blocks = sensitivity.donor_blocks(meta)
    obs_perm = sensitivity.residual_weights(x, np.delete(x,1,axis=1), blocks)
    observed = []
    for name, a in arrays.items():
        b = assay(a, obs_boot[:, :1+N_BOOT])
        p = assay(a, obs_perm)
        for j, stat in enumerate(("Crude", "Mantel–Haenszel")):
            null = p[:, j]
            assert np.isfinite(null).all() and (null > 0).all()
            draws = b[1:,j]
            draws = draws[np.isfinite(draws) & (draws > 0)]
            observed.append([cohort,name,stat,n_observed,len(a[0]),b[0,j],*interval(draws),
                             np.std(np.log(draws),ddof=1),
                             (1+(np.abs(np.log(null))>=abs(np.log(b[0,j]))-1e-10).sum())/(1+N_PERM)])
    pd.DataFrame(observed,columns=["cohort","signature","statistic","donors","probes","or",
                                  "ci_low","ci_high","bootstrap_se","permutation_p"]).to_csv(
        out/f"{cohort}_observed.tsv",sep="\t",index=False,float_format="%.9g")
    summary = []
    for size in sizes:
        alloc = sizes_by_cell(cells, prevalence, size)
        p_rng = np.random.default_rng(PERM_SEED)
        b_rng = np.random.default_rng(BOOT_SEED)
        outer = np.random.default_rng(np.random.SeedSequence([PANEL_SEED,size]))
        metrics = {(name,stat): [] for name in arrays for stat in ("Crude", "Mantel–Haenszel")}
        for k in range(panels):
            n_perm = slice_draws(N_PERM,k,panels)
            n_boot = slice_draws(N_BOOT,k,panels)
            for _ in range(100):
                source,p_blocks,b_blocks,member = expand(cells,alloc,outer)
                if np.linalg.matrix_rank(x[source])==x.shape[1]:
                    break
            else:
                raise RuntimeError("Unable to generate a full-rank donor panel")
            p, b = expanded_weights(x,source,p_blocks,b_blocks,member,n_perm,n_boot,p_rng,b_rng)
            for name,a in arrays.items():
                pn,bn = assay(a,p),assay(a,b)
                for j,stat in enumerate(("Crude","Mantel–Haenszel")):
                    null = pn[:,j]
                    b_values = bn[:,j]
                    valid_b = b_values[np.isfinite(b_values) & (b_values>0)]
                    if not np.isfinite(null).all() or np.any(null<=0) or len(valid_b)<2:
                        raise RuntimeError("Undefined ratio in simulated inference")
                    null = np.log(null)
                    cutoff = np.quantile(np.abs(null), .95, method="higher")
                    centered = null - np.median(null)
                    power = [np.mean(np.abs(centered + np.log(effect)) > cutoff) for effect in (1., *ALTERNATIVES)]
                    metrics[(name,stat)].append((null,len(valid_b),float(np.std(np.log(valid_b),ddof=1)),cutoff,power))
            if (k+1) % 10 == 0:
                print(cohort,size,k+1,flush=True)
        for (name,stat), replicates in metrics.items():
            null = np.concatenate([row[0] for row in replicates])
            boot_se = np.mean([row[2] for row in replicates])
            crit = np.median([row[3] for row in replicates])
            powers = np.asarray([row[4] for row in replicates])
            average = powers.mean(axis=0)
            margin = t.ppf(.975,len(replicates)-1)*powers.std(axis=0,ddof=1)/np.sqrt(len(replicates))
            row = [cohort,name,stat,size,sum(num for key,num in alloc.items() if key[0]),len(replicates),
                   len(null),sum(rep[1] for rep in replicates),np.std(null,ddof=1),boot_se,crit,average[0]]
            for j in range(1,4):
                row.extend((average[j],max(0.,average[j]-margin[j]),min(1.,average[j]+margin[j])))
            summary.append(row)
        cols = ["cohort","signature","statistic","donors","ds_donors","outer_panels",
                "permutations","bootstrap_valid","null_log_or_sd","bootstrap_log_or_se",
                "median_abs_log_or_critical","power_or1"]
        cols += [f"or{e:g}_{field}" for e in ALTERNATIVES for field in ("power","low","high")]
        key = out/f"{cohort}_power.tsv"
        work = key.with_suffix(".tmp")
        pd.DataFrame(summary,columns=cols).to_csv(work,sep="\t",index=False,float_format="%.9g")
        work.replace(key)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("compendium","cache","corrected","out"):
        p.add_argument(name,type=Path)
    p.add_argument("--cohort",choices=("GSE74486","GSE73747"),required=True)
    args = p.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    original = load_module(args.compendium/"analyze.py")
    sensitivity = load_module(args.compendium/"sensitivity"/"run_sensitivity.py")
    signatures = {"Published":signature(args.compendium/"GSE157078_candidate_dmcs.tsv"),
                  "ssNoob":signature(args.corrected)}
    compute(original,sensitivity,args.cache,args.cohort,signatures,SIZES,N_PANELS,args.out)


if __name__ == "__main__":
    main()
