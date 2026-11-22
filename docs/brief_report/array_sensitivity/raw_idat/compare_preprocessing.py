"""Compare GSE157078 IDAT-derived signatures with the frozen fetal donor models."""
import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


WORKFLOWS = ("IDAT_raw","Noob_background","ssNoob_dye")


def load_module(path):
    spec=importlib.util.spec_from_file_location(path.stem,path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def signature(matrix, omit):
    ref=[i for i in range(3) if i != omit]
    oe=[i for i in range(3,6) if i != omit]
    delta=matrix.iloc[:,oe].mean(axis=1)-matrix.iloc[:,ref].mean(axis=1)
    return delta[delta.abs().ge(.05)]


def pca_scores(matrix, mode, support):
    # Eigenvectors of the 6x6 sample Gram matrix. Common per-probe level is
    # removed in the centered analysis. Fix arbitrary PC signs toward row.
    values=matrix.to_numpy(dtype=float)
    centered=values-values.mean(axis=1,keepdims=True)
    rows=[]
    for name,y in (("centered",centered),("uncentered",values)):
        eigen,vectors=np.linalg.eigh(y.T@y)
        idx=np.argsort(eigen)[::-1]
        for pc in range(2):
            score=vectors[:,idx[pc]]
            rho=spearmanr(np.arange(2,8),score).statistic
            if rho<0:
                score=-score
            for i,value in enumerate(score,1):
                rows.append({"support":support,"workflow":mode,"centering":name,
                             "pc":pc+1,"sample":i,"score":value,"rc_fraction":eigen[idx[pc]]/eigen.sum(),
                             "abs_row_spearman":abs(rho)})
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ("idat_beta","geo_signals","epic_annotation","signature_module","compendium","fetal_cache","out"):
        p.add_argument(n,type=Path)
    a=p.parse_args()
    base=load_module(a.signature_module)
    geo=base.load_discovery(a.geo_signals,a.epic_annotation)
    frozen=pd.read_csv(a.compendium/"GSE157078_candidate_dmcs.tsv",sep="\t",index_col="probe_id").delta_beta
    b=pd.read_csv(a.idat_beta,sep="\t",index_col="probe_id",low_memory=False)
    cols=[f"{mode}_{i}" for mode in WORKFLOWS for i in range(1,7)]
    valid=np.isfinite(b[cols]).all(axis=1) & (b[cols]>=0).all(axis=1) & (b[cols]<=1).all(axis=1)
    shared=geo.index.intersection(b.index[valid])
    b=b.loc[shared]
    archived=geo.loc[shared].to_numpy()
    raw=b[cols[:6]].to_numpy()
    differences=np.abs(archived-raw)
    result=[]
    for i in range(6):
        d=differences[:,i]
        result.append({"sample":i+1,"matched_probes":len(d),"fraction_close_1e_8":np.mean(d<1e-8),
                       "median_absolute_difference":np.median(d),"p99_absolute_difference":np.quantile(d,.99),
                       "max_absolute_difference":d.max()})
    pd.DataFrame(result).to_csv(a.out/"raw_beta_comparison.tsv",sep="\t",index=False,float_format="%.9g")
    sigs={"effect_size":frozen}
    summary=[]
    scores=[]
    for support,keep in (("GEO_QC", b.index),("IDAT_QC", b.index[b.pass_ewastools.eq(True)])):
        block=b.loc[keep]
        for mode in WORKFLOWS:
            values=block[[f"{mode}_{i}" for i in range(1,7)]]
            scores.extend(pca_scores(values,mode,support))
            for omit in (None,*range(6)):
                label="all" if omit is None else f"omit_{omit+1}"
                selected=signature(values,omit)
                common=frozen.index.intersection(selected.index)
                key=f"{support}:{mode}:{label}"
                sigs[key]=selected
                up=int(selected.gt(0).sum())
                summary.append({"signature":key,"support":support,"workflow":mode,"omission":label,
                                "callable":len(block),"discovery_n":len(selected),"gain_n":up,
                                "gain_fraction":up/len(selected),"overlap_frozen":len(common),
                                "same_sign_overlap":np.count_nonzero(np.sign(frozen.loc[common])==np.sign(selected.loc[common]))})
    pd.DataFrame(summary).to_csv(a.out/"preprocessed_signature_counts.tsv",sep="\t",index=False,float_format="%.9g")
    pd.DataFrame(scores).to_csv(a.out/"preprocessing_pca.tsv",sep="\t",index=False,float_format="%.9g")
    base.donor_checks(a.compendium,a.fetal_cache,sigs).to_csv(
        a.out/"preprocessed_fetal_results.tsv",sep="\t",index=False,float_format="%.9g",na_rep="NA")


if __name__=="__main__":
    main()
