"""Assess enrichment in a public EPIC row reference and a published slide-effect set."""
import argparse
import gzip
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from run_chip_checks import read_module


AUTOSOMES = [str(i) for i in range(1, 23)]
BOOTSTRAPS = 10000


def epi_features(path):
    with gzip.open(path, "rt") as stream:
        for offset, line in enumerate(stream):
            if line.startswith("IlmnID,"):
                break
        else:
            raise ValueError("Missing EPIC header")
    out = pd.read_csv(path, skiprows=offset, usecols=["Name","CHR","Infinium_Design_Type"],
                      dtype=str, index_col="Name")
    out = out[out.CHR.isin(AUTOSOMES) & out.Infinium_Design_Type.isin(["I","II"])]
    assert out.index.is_unique
    return out.rename(columns={"CHR":"chromosome","Infinium_Design_Type":"design"})


def technical_samples(path):
    samples = []
    for rec in path.read_text().split("^SAMPLE = ")[1:]:
        pooled = re.search(r"^!Sample_title = subject([A-D])_pooled_replicate(\d+)$",rec,re.M)
        if pooled:
            chip, pos = re.search(r"_(\d+)_(R\d+C\d+)_Grn\.idat\.gz",rec).groups()
            samples.append({"donor":pooled[1],"chip":chip,"position":pos,"row":int(pos[1:3]),
                            "geo_accession":rec.splitlines()[0],"array":f"{chip}_{pos}"})
    notes = pd.DataFrame(samples)
    assert len(notes)==32 and notes.groupby("donor").row.nunique().eq(8).all()
    # The primary source dropped this pooled outlier under its prespecified data QC.
    keep = notes.row.between(2,7) & notes.array.ne("205832330169_R03C01")
    return notes.loc[keep].set_index("array")


def fitted_row_effects(path, meta, mode):
    other = pd.get_dummies(meta[["donor","chip"]],drop_first=True,dtype=float)
    design = np.column_stack((np.ones(len(meta)), meta.row.ge(5), other))
    assert np.linalg.matrix_rank(design)==design.shape[1]
    weights = np.linalg.pinv(design)[1]
    if mode == "raw":
        meth = [k+".M" for k in meta.index]
        unmeth = [k+".U" for k in meta.index]
        columns = meth + unmeth
    else:
        columns = [f"X{k}_betas" for k in meta.index]
    effects = []
    for block in pd.read_csv(path,usecols=["ID_REF",*columns],index_col="ID_REF",chunksize=40000):
        vals = block[columns].to_numpy(dtype=float)
        if mode == "raw":
            halves = len(columns)//2
            m,u = vals[:,:halves], vals[:,halves:]
            with np.errstate(divide="ignore",invalid="ignore"):
                vals=m/(m+u+100)
        valid = np.isfinite(vals).all(axis=1) & (vals>=0).all(axis=1) & (vals<=1).all(axis=1)
        effects.append(pd.Series(vals[valid] @ weights,index=block.index[valid],name="row_effect"))
    result = pd.concat(effects)
    assert result.index.is_unique
    return result


def slide_sets(directory):
    with tempfile.TemporaryDirectory(prefix="dnmt3l-seff-") as tmp:
        export = Path(tmp)/"slide.tsv"
        command = '''
args <- commandArgs(TRUE)
labels <- c("S_0_20","S_20_40","S_40_60","S_60_80","S_80_100")
for (x in c(labels, "S_high")) load(file.path(args[1], paste0(x, ".rda")))
all_probes <- unique(unlist(mget(labels), use.names=FALSE))
stopifnot(all(S_high %in% all_probes))
write.table(data.frame(probe=all_probes, high=all_probes %in% S_high),
            args[2],sep="\t",quote=FALSE,row.names=FALSE)
'''
        subprocess.run(["Rscript","-e",command,str(directory),str(export)],check=True)
        return pd.read_csv(export,sep="\t",index_col="probe").high


def compare(frame,target,exposure,rng):
    # A = target and exposed; D = neither. Resample whole chromosomes.
    flags=frame.assign(a=target&exposure,b=target&~exposure,c=~target&exposure,d=~target&~exposure)
    counts=flags.groupby(["chromosome","design","beta_bin"],observed=True)[["a","b","c","d"]].sum()
    a,b,c,d=counts.T.to_numpy()
    n=a+b+c+d
    numer=a*d/n
    denom=b*c/n
    numden=pd.DataFrame({"numer":numer,"denom":denom},index=counts.index).groupby(level=0).sum()
    numden=numden.reindex(AUTOSOMES,fill_value=0)
    boot = rng.multinomial(len(AUTOSOMES), np.full(len(AUTOSOMES),1/len(AUTOSOMES)),size=BOOTSTRAPS)
    q = (boot@numden.numer)/(boot@numden.denom)
    chrom=numden.loc[(numden.numer+numden.denom)>0]
    overlap_expect=np.divide((a+b)*c,c+d,out=np.zeros(len(c),dtype=float),where=c+d>0)
    matched=(c+d)>0
    return {"n":int(n.sum()),"target_n":int((a+b).sum()),"exposed_n":int((a+c).sum()),"overlap":int(a.sum()),
            "matched_target_n":int((a[matched]+b[matched]).sum()),
            "matched_overlap":int(a[matched].sum()),"matched_expected_overlap":float(overlap_expect.sum()),
            "or":float(numer.sum()/denom.sum()),"chromosome_ci_low":np.quantile(q,.025),"chromosome_ci_high":np.quantile(q,.975),
            "informative_chromosomes":len(chrom),"chromosomes_or_gt_one":int((chrom.numer>chrom.denom).sum())}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("signals","manifest","signature_module","compendium","technical_soft","raw","sesame","hecker_data","out"):
        p.add_argument(name,type=Path)
    a=p.parse_args()
    orig=read_module(a.signature_module)
    beta=orig.load_discovery(a.signals,a.manifest)
    frozen=pd.read_csv(a.compendium/"GSE157078_candidate_dmcs.tsv",sep="\t",index_col="probe_id").delta_beta
    probes=epi_features(a.manifest).join(beta.mean(axis=1).rename("mean_beta"),how="inner")
    probes["beta_bin"]=np.digitize(probes.mean_beta,[.2,.4,.6,.8])
    probes["signature"]=frozen.reindex(probes.index)
    meta=technical_samples(a.technical_soft)
    meta.to_csv(a.out/"technical_reference_arrays.tsv",sep="\t")
    comparisons=[]
    directions=[]
    fitted={}
    rng=np.random.default_rng(2026112006)
    for mode,path in (("raw",a.raw),("sesame",a.sesame)):
        row=fitted_row_effects(path,meta,mode)
        fitted[mode]=row
        frame=probes.join(row,how="inner")
        threshold=frame.row_effect.abs().quantile(.9)
        s=frame.signature.notna()
        comparisons.append({"reference":f"Butler_{mode}_top_decile","reference_threshold":threshold,
                            **compare(frame,s,frame.row_effect.abs().ge(threshold),rng)})
        selected=frame.loc[s & frame.row_effect.ne(0)]
        up=selected.signature.gt(0)
        row_up=selected.row_effect.gt(0)
        directions.append({"reference":f"Butler_{mode}","agreement":(up==row_up).mean(),
                           "expected_raw_agreement":up.mean()*row_up.mean()+(1-up.mean())*(1-row_up.mean()),
                           "spearman":spearmanr(selected.signature,selected.row_effect).statistic,
                           **compare(selected,up,row_up,rng)})
        print(f"Completed Butler {mode} reference.",flush=True)
    common = probes.join(pd.DataFrame(fitted).dropna(),how="inner")
    for mode in ("raw","sesame"):
        threshold=common[mode].abs().quantile(.9)
        comparisons.append({"reference":f"Butler_{mode}_common_top_decile","reference_threshold":threshold,
                            **compare(common,common.signature.notna(),common[mode].abs().ge(threshold),rng)})
    s_high=slide_sets(a.hecker_data)
    frame=probes.join(s_high,how="inner")
    comparisons.append({"reference":"Hecker_S_high","reference_threshold":np.nan,
                        **compare(frame,frame.signature.notna(),frame.high,rng)})
    for rows,name in ((comparisons,"position_susceptibility_enrichment"),(directions,"external_chamber_direction")):
        pd.DataFrame(rows).to_csv(a.out/(name+".tsv"),sep="\t",index=False,float_format="%.9g",na_rep="NA")


if __name__=="__main__":
    main()
