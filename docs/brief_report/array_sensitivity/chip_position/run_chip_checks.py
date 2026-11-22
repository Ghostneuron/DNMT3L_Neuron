"""Leave-one-array-out and within-group row diagnostics for GSE157078."""
import argparse
from itertools import permutations
from pathlib import Path

import numpy as np
import pandas as pd

import importlib.util


def read_module(path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def omissions(beta, frozen, out):
    signatures = {"effect_size": frozen}
    for i in range(6):
        control = [j for j in range(3) if j != i]
        oe = [j for j in range(3, 6) if j != i]
        delta = beta.iloc[:, oe].mean(axis=1) - beta.iloc[:, control].mean(axis=1)
        signatures[f"omit_{i+1}"] = delta.loc[delta.abs().ge(.05)]
    rows = []
    for label, sig in signatures.items():
        common = frozen.index.intersection(sig.index)
        n_gain = sig.gt(0).sum()
        rows.append({"signature":label, "discovery_n":len(sig), "gain_n":n_gain,
                     "gain_fraction":n_gain / len(sig), "overlap_original":len(common),
                     "same_sign_overlap":int((np.sign(frozen.loc[common])==np.sign(sig.loc[common])).sum())})
    pd.DataFrame(rows).to_csv(out / "omission_signature_counts.tsv", sep="\t", index=False, float_format="%.9g")
    return signatures


def row_diagnostics(beta, frozen, out):
    delta = beta.iloc[:, 3:].mean(axis=1) - beta.iloc[:, :3].mean(axis=1)
    np.testing.assert_allclose(delta.loc[frozen.index], frozen, atol=5e-7, rtol=5e-6)
    assert len(delta.loc[delta.abs().ge(.05)]) == len(frozen)
    orders = list(permutations(range(3)))
    masks = {"all_discovery_probes":delta.ne(0), "frozen_all":delta.index.isin(frozen.index),
             "frozen_gains":delta.index.isin(frozen.index[frozen>0]),
             "frozen_losses":delta.index.isin(frozen.index[frozen<0])}
    rows = []
    contrasts = []
    for label, keep in masks.items():
        values = beta.loc[keep].to_numpy()
        sign = np.sign(delta.loc[keep].to_numpy())
        signed = values * sign[:,None]
        steps = np.diff(signed, axis=1)
        control = signed[:,2] - signed[:,0]
        oe = signed[:,5] - signed[:,3]
        left = [np.all(np.diff(signed[:,:3][:,p], axis=1)>0,axis=1) for p in orders]
        right = [np.all(np.diff(signed[:,3:][:,p], axis=1)>0,axis=1) for p in orders]
        distribution = np.array([np.count_nonzero(a & b) for a in left for b in right])
        strict = distribution[0]
        rows.append({"probe_set":label, "n":len(signed), "control_endpoint_aligned":np.mean(control>0),
                     "oe_endpoint_aligned":np.mean(oe>0), "both_endpoints_aligned":np.mean((control>0)&(oe>0)),
                     "control_strict_aligned":np.mean(left[0]), "oe_strict_aligned":np.mean(right[0]),
                     "both_strict_n":strict,"both_strict_aligned":strict/len(signed),
                     "orderings_at_least_observed":np.count_nonzero(distribution>=strict),
                     "row_orders":len(distribution), "row_order_median_fraction":np.median(distribution)/len(signed)})
        data = {
            **dict(zip(("R03_minus_R02","R04_minus_R03","R05_minus_R04","R06_minus_R05","R07_minus_R06"),
                       steps.T, strict=True)),
            "R04_minus_R02":control, "R07_minus_R05":oe,
            "full_group_difference":delta.loc[keep].abs().to_numpy(),
            "R02_beta_minus_other_control_mean":values[:,0]-values[:,1:3].mean(axis=1),
        }
        for key,y in data.items():
            lo,med,hi=np.quantile(y,[.25,.5,.75])
            contrasts.append({"probe_set":label,"contrast":key,"median":med,"q25":lo,"q75":hi,"positive_fraction":np.mean(y>0)})
    pd.DataFrame(rows).to_csv(out / "within_group_row_trends.tsv",sep="\t",index=False,float_format="%.9g")
    pd.DataFrame(contrasts).to_csv(out / "row_beta_contrasts.tsv",sep="\t",index=False,float_format="%.9g")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("signals","manifest","signature_module","compendium","fetal_cache","out"):
        p.add_argument(name,type=Path)
    a=p.parse_args()
    a.out.mkdir(exist_ok=True,parents=True)
    base=read_module(a.signature_module)
    frozen=pd.read_csv(a.compendium/"GSE157078_candidate_dmcs.tsv",sep="\t",index_col="probe_id").delta_beta
    beta=base.load_discovery(a.signals,a.manifest)
    signs=omissions(beta,frozen,a.out)
    row_diagnostics(beta,frozen,a.out)
    base.donor_checks(a.compendium,a.fetal_cache,signs).to_csv(
        a.out/"omission_fetal_results.tsv",sep="\t",index=False,float_format="%.9g",na_rep="NA")


if __name__=="__main__":
    main()
