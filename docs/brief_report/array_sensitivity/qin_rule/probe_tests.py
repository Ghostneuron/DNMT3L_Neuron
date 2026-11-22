import sys,gzip
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
d=Path(sys.argv[1])
geo=d/'GSE157078_series_matrix.txt.gz'
with gzip.open(geo,'rt') as f:
 skip=next(i for i,line in enumerate(f) if line.startswith('"ID_REF"'))
geo=pd.read_csv(geo,sep='\t',skiprows=skip,index_col='ID_REF',comment='!').dropna()
geo.columns=[f'GEO_{i}' for i in range(1,7)]
cols=['chr','pos','Name','CpG_maf','SBE_maf','UCSC_RefGene_Name','UCSC_RefGene_Group','Relation_to_Island','Islands_Name','Type','Color']
ann=pd.read_csv(d/'epic_annotation.tsv.gz',sep='\t',usecols=cols,low_memory=False).set_index('Name')
zcols=['probeID','CpG_chrm','CpG_beg','MASK_general','MASK_mapping','MASK_snp5_common','MASK_sub30_copy','MASK_sub35_copy','MASK_sub40_copy','MASK_typeINextBaseSwitch']
zhou=pd.read_csv(d/'zhou_epic.tsv.gz',sep='\t',usecols=zcols,index_col='probeID')
geop=Path(sys.argv[2])
detcols=[f'SAMPLE {i} Detection Pval' for i in range(1,7)]
det=pd.read_csv(geop,sep='\t',usecols=['ID_REF',*detcols],index_col='ID_REF').le(.01).all(axis=1).rename('detect')
rerun=pd.read_csv(d/'swan_betas.tsv.gz',sep='\t',index_col='probe_id')
pre='ssNoob_dye'
nob=pd.read_csv(Path(sys.argv[3]),sep='\t',index_col='probe_id',usecols=['probe_id',*[f'{pre}_{i}' for i in range(1,7)]])
common=geo.join([ann,det,rerun,nob],how="inner").join(zhou,how="left")
common=common.loc[common.chr.isin([f'chr{i}' for i in range(1,23)]) & common.index.str.startswith('cg')]
for field in ['MASK_mapping','MASK_snp5_common','MASK_general','MASK_sub30_copy','MASK_sub35_copy','MASK_sub40_copy','MASK_typeINextBaseSwitch']:
 common[field]=common[field].eq(True) | common.CpG_chrm.isna()
methods=['GEO','SWAN','ssNoob_SWAN','ssNoob_dye','ssNoob_0']
for method in methods:
 vals=common[[f'{method}_{i}' for i in range(1,7)]].to_numpy()
 common[f'{method}_delta']=vals[:,3:].mean(axis=1)-vals[:,:3].mean(axis=1)
 common[f'{method}_p']=ttest_ind(vals[:,3:],vals[:,:3],axis=1,equal_var=True).pvalue
w4=pd.ExcelFile(d/'mmc4.xlsx')
pub=w4.parse(0,index_col='Name')
print("Missing author probes",len(pub.index.difference(common.index)),list(pub.index.difference(common.index)),flush=True)
assert len(pub)==347
comp=pub.loc[pub.index.intersection(common.index)]
assert np.max(abs(common.loc[comp.index,'GEO_p']-comp.pValue))<2e-9
assert np.max(abs(common.loc[comp.index,'GEO_delta']-comp['D-Z.Beta']))<2e-9
masks={
'Detection only':True,
'Mapping + common SNP':~(common.MASK_mapping|common.MASK_snp5_common),
'Mapping + common SNP + sub30':~(common.MASK_mapping|common.MASK_snp5_common|common.MASK_sub30_copy),
'Mapping + common SNP + sub35':~(common.MASK_mapping|common.MASK_snp5_common|common.MASK_sub35_copy),
'Mapping + common SNP + sub40':~(common.MASK_mapping|common.MASK_snp5_common|common.MASK_sub40_copy),
'Zhou general mask':~common.MASK_general,
}
rows=[]
for name,mask in masks.items():
 mask=common.detect & mask
 for method in methods:
  m=common.loc[common[f'{method}_p'].le(.05) & mask]
  case=m.loc[m[f'{method}_delta'].abs().ge(.1)]
  ov=case.index.intersection(pub.index)
  old=pub.loc[ov,'D-Z.Beta'].gt(0)
  new=case.loc[ov,f'{method}_delta'].gt(0)
  ch21=m[m.chr.eq('chr21')][f'{method}_delta']
  rows.append([name, method,sum(np.broadcast_to(mask,len(common))),len(case),case[f'{method}_delta'].gt(0).sum(),case[f'{method}_delta'].lt(0).sum(),len(ov),sum(old!=new),sum(ch21>0),sum(ch21<0)])
pd.DataFrame(rows,columns=['mask','method','eligible','total','gains','losses','overlap347','reversed_among_overlap','chr21_nominal_gains','chr21_nominal_losses']).to_csv(d/'probe_counts.tsv',sep='\t',index=False)
common.to_pickle(d/'probe_data.pkl')
