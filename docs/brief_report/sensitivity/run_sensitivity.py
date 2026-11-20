"""Run the 2026-11-20 sensitivity plan against the frozen v1.2.0 analysis."""
import argparse
import gzip
import importlib.util
from itertools import combinations, product
from math import comb, prod
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

WINDOWS = {"overlap_70_112": (70, 112), "offset_84_126": (84, 126), "full_23_184": (23, 184)}
N_PERM = 9999
TIE = 1e-7


def original_code(path):
    spec = importlib.util.spec_from_file_location("brief_v12", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def spiers_slopes(matrix):
    fields = {}
    with gzip.open(matrix, "rt") as handle:
        for index, line in enumerate(handle):
            if line.startswith("!series_matrix_table_begin"):
                skiprows = index + 1
                break
            items = line.strip().replace('"', '').split("\t")
            if items[0] == "!Sample_geo_accession":
                fields["sample"] = items[1:]
            if items[0] == "!Sample_characteristics_ch1":
                key = items[1].split(": ", 1)[0]
                fields[key] = [item.split(": ", 1)[1] for item in items[1:]]
        else:
            raise ValueError("The GEO series matrix lacks its data-table marker")
    meta = pd.DataFrame(fields).set_index("sample")
    days = meta["days post-conception"].to_numpy(dtype=float)
    sex = meta.gender.map({"M": 1, "F": 0}).to_numpy(dtype=float)
    assert len(meta) == 179 and meta.index.is_unique and np.isfinite(days).all() and np.isfinite(sex).all()
    windows = {}
    counts = []
    for name, (lo, hi) in WINDOWS.items():
        keep = np.flatnonzero((days >= lo) & (days <= hi))
        assert len(keep) >= 20 and np.ptp(days[keep]) >= 28
        x = np.column_stack([np.ones(len(keep)), sex[keep], (days[keep]-days[keep].mean())/28.])
        assert np.linalg.matrix_rank(x) == x.shape[1]
        windows[name] = keep, np.linalg.pinv(x)[-1]
        counts.append({"window": name, "reference_donors": len(keep), "days_min": int(days[keep].min()),
                       "days_max": int(days[keep].max()), "female": int(sum(sex[keep]==0)),
                       "male": int(sum(sex[keep]==1))})
    source = pd.read_csv(matrix, sep="\t", skiprows=skiprows, index_col="ID_REF", comment="!",
                         chunksize=12000, dtype={key: "float64" for key in meta.index})
    chunks = []
    for table in source:
        assert table.columns.equals(meta.index)
        vals = table.to_numpy()
        slopes = {}
        for name, (keep, weights) in windows.items():
            y = vals[:, keep]
            valid = np.isfinite(y).all(axis=1) & (y >= 0).all(axis=1) & (y <= 1).all(axis=1)
            slopes[name] = np.where(valid, y @ weights, np.nan)
        chunks.append(pd.DataFrame(slopes, index=table.index))
    slope_table = pd.concat(chunks)
    assert slope_table.index.is_unique
    return slope_table, pd.DataFrame(counts)


def probe_frame(original, beta, annotation, signature):
    probes = annotation.loc[beta.index].copy()
    probes["mean_bin"] = np.digitize(beta.mean(axis=1), [.2, .4, .6, .8])
    probes["stratum"] = original.stratum_codes(probes)
    probes["signature"] = signature.reindex(probes.index)
    return probes


def signature_panel(probes):
    sig = probes.loc[probes.signature.notna()]
    strata, uniques = pd.factorize(sig.stratum, sort=True)
    n = len(sig)
    rows = strata*2 + sig.signature.lt(0).to_numpy()
    panel = csr_matrix((np.ones(n), (rows, np.arange(n))), shape=(2*len(uniques), n))
    return sig, panel, len(uniques)


def values_to_or(original, delta, panel, n_strata):
    pos = (panel @ (delta > TIE)).reshape(n_strata, 2, -1)
    neg = (panel @ (delta < -TIE)).reshape(n_strata, 2, -1)
    return original.mh_or(pos[:,0], neg[:,0], pos[:,1], neg[:,1])


def or_from_weights(original, y, weights, panel, n_strata):
    output=[]
    for start in range(0,weights.shape[1],128):
        result=y @ weights[:,start:start+128]
        output.append(values_to_or(original,result,panel,n_strata))
    return np.concatenate(output)


def donor_blocks(meta):
    groups=[]
    for _, sub in meta.reset_index(drop=True).groupby("donor", sort=False):
        assert sub.sex.nunique()==sub.dx.nunique()==1 and sub.region.is_unique
        ordered=sub.sort_values("region")
        key=(ordered.sex.iloc[0], tuple(ordered.region))
        groups.append((key, ordered.index.to_numpy(), int(ordered.dx.iloc[0])))
    blocks=[]
    for key in sorted({group[0] for group in groups}):
        members=[(i,dx) for block,i,dx in groups if block==key]
        positions=np.vstack([i for i,_ in members])
        diagnoses=np.array([dx for _,dx in members])
        blocks.append((key,positions,diagnoses))
    return blocks


def residual_weights(x, z, blocks):
    rng=np.random.default_rng(2026112002)
    h=z@np.linalg.pinv(z)
    rr=np.eye(len(x))-h
    observed=np.linalg.pinv(x)[1]
    shuffled=np.empty((len(x),N_PERM))
    for _,positions,_ in blocks:
        for i in range(N_PERM):
            perm=positions[rng.permutation(len(positions))]
            shuffled[perm.ravel(),i]=observed[positions.ravel()]
    return h @ observed[:,None] + rr @ shuffled


def label_weights(z, blocks):
    rng=np.random.default_rng(2026112003)
    total=prod(comb(len(dx), int(dx.sum())) for _,_,dx in blocks)
    exact=total <= 10000
    def enumerated():
        pools=[list(combinations(range(len(dx)), int(dx.sum()))) for _,_,dx in blocks]
        yield from product(*pools)
    def sampled():
        for _ in range(N_PERM):
            yield tuple(rng.choice(len(dx), int(dx.sum()), replace=False) for _,_,dx in blocks)
    draw_count=total if exact else N_PERM
    labels=np.zeros((len(z),draw_count))
    for col,selections in enumerate(enumerated() if exact else sampled()):
        for (_,positions,_),chosen in zip(blocks,selections,strict=True):
            labels[positions[np.array(chosen,dtype=int)].ravel(),col]=1
    adjusted=labels-z@(np.linalg.pinv(z)@labels)
    ss=(adjusted**2).sum(axis=0)
    weights=np.divide(adjusted,ss,out=np.full_like(adjusted,np.nan),where=ss>1e-10)
    return weights,total,exact


def interval(values):
    return np.nanquantile(values,[.025,.5,.975])


def run_permutation(original, cohort, meta, beta, probes, frozen):
    x,_=original.designs(meta)
    z=np.delete(x,1,axis=1)
    blocks=donor_blocks(meta)
    sig,panel,n_strata=signature_panel(probes)
    y=np.ascontiguousarray(beta.loc[sig.index],dtype="float64")
    obs=or_from_weights(original,y,np.linalg.pinv(x)[1,:,None],panel,n_strata)[0]
    assert np.isclose(obs,frozen.estimate,rtol=2e-5), (cohort,obs,frozen.estimate)
    label,total,exact=label_weights(z,blocks)
    rows=[]
    for method,weights,enumerated in (
        ("Freedman–Lane",residual_weights(x,z,blocks),False),
        ("Blocked diagnosis labels",label,exact),
    ):
        scores=or_from_weights(original,y,weights,panel,n_strata)
        valid_weight = np.isfinite(weights).all(axis=0)
        scores=scores[valid_weight & np.isfinite(scores) & (scores > 0)]
        assert len(scores)>0
        tail=int(np.sum(np.abs(np.log(scores))>=abs(np.log(obs))-1e-10))
        add=0 if enumerated else 1
        low,median,high=interval(scores)
        rows.append({"cohort":cohort,"method":method,"enumerated":enumerated,
                     "label_assignments":total if method=="Blocked diagnosis labels" else None,
                     "draws":weights.shape[1],"valid":len(scores),"tail":tail,"p_two_sided":(tail+add)/(len(scores)+add),
                     "or":obs,"null_low":low,"null_median":median,"null_high":high,
                     "leave_one_out_min":frozen.jackknife_min,"leave_one_out_max":frozen.jackknife_max,
                     "blocks":" | ".join(f"{key[0]}:{'+'.join(key[1])} {int(dx.sum())}/{len(dx)} DS" for key,_,dx in blocks)})
    return rows


def run_maturation(original, cohort, meta, beta, probes, slopes, counts):
    weights,_,_=original.fit_weights(meta,np.random.default_rng(original.SEED),original.N_BOOT)
    weights=weights[:,:original.N_BOOT+1]
    results=[]
    for win in WINDOWS:
        sl=slopes[win].reindex(probes.index).dropna()
        panel_probes=probes.loc[sl.index]
        sig,panel,n_strata=signature_panel(panel_probes)
        center=sl-sl.groupby(panel_probes.stratum,sort=False).transform("mean")
        center_values=center.to_numpy()
        projector=center_values @ beta.loc[sl.index].to_numpy() / (center_values@center_values)
        base=np.ascontiguousarray(beta.loc[sig.index],dtype="float64")
        adjusted=base-center.loc[sig.index].to_numpy()[:,None]*projector
        baseline=or_from_weights(original,base,weights,panel,n_strata)
        corrected=or_from_weights(original,adjusted,weights,panel,n_strata)
        external=values_to_or(original,sl.loc[sig.index].to_numpy()[:,None],panel,n_strata)[0]
        ratio=corrected/baseline
        vals={}
        for label,v in [("baseline",baseline),("age_adjusted",corrected),("adjusted_to_baseline_ratio",ratio)]:
            lo,median,hi=interval(v[1:])
            vals.update({label:v[0],label+"_ci_low":lo,label+"_bootstrap_median":median,label+"_ci_high":hi})
        info=counts.loc[counts.window.eq(win)].iloc[0].to_dict()
        results.append({"cohort":cohort,**info,"shared_probes":len(sl),"signature_gains":int(sig.signature.gt(0).sum()),
                        "signature_losses":int(sig.signature.lt(0).sum()),"valid_paired_bootstrap":int(np.isfinite(ratio[1:]).sum()),
                        "dnmt3l_vs_spiers_sign_or":external,
                        "ds_slope_projection_scale":float(projector @ weights[:,0]),**vals})
    return results


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("frozen",type=Path,help="Extracted v1.2.0 compendium with analyze.py and signature")
    parser.add_argument("cache",type=Path,help="Output cache from original prepare_inputs.py")
    parser.add_argument("spiers",type=Path,help="GSE58885_series_matrix.txt.gz")
    parser.add_argument("out",type=Path)
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    original=original_code(args.frozen/"analyze.py")
    result=pd.read_csv(args.frozen/"estimates.tsv",sep="\t")
    sig=pd.read_csv(args.frozen/"GSE157078_candidate_dmcs.tsv",sep="\t",index_col="probe_id").delta_beta
    annotation=pd.read_pickle(args.cache/"eligible_annotations.pkl")
    perm=[]
    for cohort in ("GSE74486","GSE73747"):
        meta,beta=original.load_cohort(args.cache,cohort,None,False)
        probes=probe_frame(original,beta,annotation,sig)
        frozen=result.query('context==@cohort and statistic=="DNMT3L_OR"').iloc[0]
        perm.extend(run_permutation(original,cohort,meta,beta,probes,frozen))
        print("Completed donor permutations",cohort,flush=True)
    pd.DataFrame(perm).to_csv(args.out/"donor_permutation.tsv",sep="\t",index=False,na_rep="NA",float_format="%.9g")
    slopes,counts=spiers_slopes(args.spiers)
    print("Estimated Spiers slopes in reference windows",counts.to_dict("records"),flush=True)
    mat=[]
    for cohort in ("GSE74486","GSE73747"):
        meta,beta=original.load_cohort(args.cache,cohort,None,False)
        probes=probe_frame(original,beta,annotation,sig)
        mat.extend(run_maturation(original,cohort,meta,beta,probes,slopes,counts))
        print("Completed Spiers maturation analysis",cohort,flush=True)
    pd.DataFrame(mat).to_csv(args.out/"spiers_maturation.tsv",sep="\t",index=False,na_rep="NA",float_format="%.9g")


if __name__=="__main__":
    main()
