from pathlib import Path
import sys,re
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind,pearsonr
from scipy.sparse import csr_matrix
d=Path(sys.argv[1])
p=pd.read_pickle(d/'probe_data.pkl')
mask=p.detect & ~(p.MASK_mapping|p.MASK_snp5_common|p.MASK_sub30_copy)
data=p.loc[mask]
w4=pd.ExcelFile(d/'mmc4.xlsx')
w5=pd.ExcelFile(d/'mmc5.xlsx')
pub=w4.parse(0,index_col='Name')
summary=w4.parse(11)
pub_sig=set(summary.Merged.dropna())
rna_short=set(w5.parse(9).RNAseq.dropna()); rna_long=set(w5.parse(10).RNAseq.dropna())
expected=w5.parse(12).set_index('GeneSymbol')
pub_bump_names=set(summary.Bump.dropna())
def gene_tuple(names):
 return [x.strip() for x in names.split(';') if x.strip()]
def genes_called(frame):
 return {item for names in frame.UCSC_RefGene_Name.dropna() for item in gene_tuple(names)}
body,utr3,prom={},{},{}
for i,(name,group) in enumerate(data[['UCSC_RefGene_Name','UCSC_RefGene_Group']].fillna('').itertuples(index=False,name=None)):
 names=gene_tuple(name); regions=gene_tuple(group)
 for n,r in set(zip(names,regions)):
  if r=='Body': body.setdefault(n,[]).append(i)
  if r=="3'UTR": utr3.setdefault(n,[]).append(i)
  if r in ("TSS200","5'UTR"): prom.setdefault(n,[]).append(i)
def make_matrix(g):
 names=list(g)
 obs=np.concatenate(list(g.values()))
 n=np.fromiter((len(g[x]) for x in names),int)
 mat=csr_matrix((np.repeat(1/n,n),(np.repeat(np.arange(len(n)),n),obs)),shape=(len(n),len(data)))
 return names,mat
groups={r:make_matrix(g) for r,g in [('Body',body),("3'UTR",utr3),('Promoter',prom)]}
# Re-test the listed significant regional features at their original hg19 coordinates.
# Isoforms collapse at gene/region, and shore/shelf labels are believed to refer
# to annotated CpG islands whose coordinates are given by the authors.
a=pd.DataFrame({'Island':data.Islands_Name,'Group':data.Relation_to_Island,'Name':data.UCSC_RefGene_Name}).reset_index(drop=True)
island_groups=a.dropna(subset=['Island','Group']).groupby(['Island','Group']).indices
# groupby indices are relative to filtered frame so map through its original index.
valid_index=a.dropna(subset=['Island','Group']).index.to_numpy()
regional_rows={}
for k,label in zip(range(6,11),['Island','N_Shelf','N_Shore','S_Shelf','S_Shore']):
 rows=w4.parse(k)
 rows=rows.loc[(rows['P-Value']<=.05)&(rows['beta-Difference'].abs()>=.05)]
 for row in rows.itertuples(index=False,name=None):
  identifier=row[0]
  if not isinstance(identifier,str): continue
  hits=island_groups.get((identifier,label),[])
  if len(hits):
   regional_rows[(identifier,label)]=valid_index[hits].tolist()
# Refit exact intervals with a mappable published Bump gene rather than
# inventing the original smoothing and bootstrap arguments.
bump=w4.parse(3)
bump_groups={}
old_bump=[]
for row in bump.itertuples(index=False,name=None):
 identity,chrom,start,end,value,area,clust,idx1,idx2,n,n_total,pval,_,_,_,genes=row
 if pd.isna(genes) or not(set(gene_tuple(str(genes))) & pub_bump_names) or not np.isfinite(start) or not np.isfinite(end): continue
 chrom=str(chrom).replace('ch','chr',1) if str(chrom).startswith('ch') and not str(chrom).startswith('chr') else str(chrom)
 if pval>.05: continue
 keep=np.flatnonzero(data.chr.eq(chrom).to_numpy() & data.pos.between(start,end).to_numpy())
 if len(keep)<2: continue
 bump_groups[(str(identity),str(genes))]=keep.tolist()
groups['Published CGI bins']=make_matrix(regional_rows)
groups['Published bump intervals']=make_matrix(bump_groups)
reggenes={x:set(gene_tuple(';'.join(a.loc[regional_rows[x],'Name'].dropna()))) for x in regional_rows}
agg={}
coefs_rows=[]
overlap_rows=[]
regional_detail=[]
pub_gene_call=genes_called(data.loc[data.index.intersection(pub.index)])
groupsets={}
for meth in ['GEO','SWAN','ssNoob_SWAN','ssNoob_dye']:
 cases=data.loc[data[f'{meth}_p'].le(.05)&data[f'{meth}_delta'].abs().ge(.1)]
 vals=data[[f'{meth}_{i}' for i in range(1,7)]].to_numpy()
 features={}
 groupcalls={}
 for region,(ids,A) in groups.items():
  B=A @ vals
  effects=B[:,3:].mean(axis=1)-B[:,:3].mean(axis=1)
  prob=ttest_ind(B[:,3:],B[:,:3],axis=1,equal_var=True).pvalue
  features[region]=pd.Series(effects,index=pd.Index(ids,tupleize_cols=False))
  cutoff=.1 if region=='Published bump intervals' else .05
  ind=(prob<=.05)&(abs(effects)>=cutoff)
  selected=features[region].index[ind]
  if region=='Published bump intervals':
   called={g for (_,gene) in selected for g in gene_tuple(gene) if g in pub_bump_names}
  elif region=='Published CGI bins':
   called={g for key in selected for g in reggenes[key]}
  else:
   called=set(selected)
  groupcalls[region]=called
  regional_detail.append([meth,region,A.shape[0],int(ind.sum()),len(called)])
 cg_genes=genes_called(cases)
 region_genes=groupcalls['Body']|groupcalls["3'UTR"]|groupcalls['Published CGI bins']
 merged=cg_genes|region_genes|groupcalls['Published bump intervals']
 groupsets[meth]=merged
 for desc,g in [('Single-site only',cg_genes),('Reconstructed union',merged)]:
  overlap_rows.append([meth,desc,len(g),len(g&rna_short),len(g&rna_long),len(g&pub_sig),';'.join(sorted(g&rna_short))])
 for desc,genes in [('Published 225 gene panel',set(expected.index)),('Reconstructed RNA–methylation intersection',merged&rna_long)]:
  if desc.startswith('Published'):
   expr=expected.log2FoldChange
  else:
   expr=w5.parse(0,usecols=['GeneSymbol','log2FoldChange','qvalue']).dropna().sort_values('qvalue').drop_duplicates('GeneSymbol').set_index('GeneSymbol').log2FoldChange
  for reg in ['Body','Promoter']:
   vals=pd.concat([expr.rename('expression'),features[reg].rename('beta_delta')],axis=1).reindex(sorted(genes)).dropna()
   vals.to_csv(d/f'corr_{meth}_{reg}_{desc[:3]}.tsv',sep='\t')
   if len(vals)>2:
    result=pearsonr(vals.expression,vals.beta_delta)
    coefs_rows.append([meth,desc,reg,len(vals),result.statistic,result.pvalue])
 agg[meth]=features
for sample_cols in ['Methylation difference in gene body','Methylation difference in promoter']:
 complete=expected[['log2FoldChange',sample_cols]].dropna()
 for mode,adjust in [('as stored', complete),('excluding abs beta >1',complete[complete[sample_cols].abs()<=1])]:
  result=pearsonr(adjust.iloc[:,0],adjust.iloc[:,1])
  coefs_rows.append(['Published cells',mode,sample_cols,len(adjust),result.statistic,result.pvalue])
for display,g in [('Published merged',pub_sig),('Published single-site workbook',set(summary['Single site'].dropna()))]:
 overlap_rows.append(['Published',display,len(g),len(g&rna_short),len(g&rna_long),len(g&pub_sig),';'.join(sorted(g&rna_short))])
pd.DataFrame(regional_detail,columns=['method','region','tested','calls','genes']).to_csv(d/'regions.tsv',sep='\t',index=False)
pd.DataFrame(coefs_rows,columns=['method','contrast','region','n','r','p']).to_csv(d/'correlations.tsv',sep='\t',index=False)
pd.DataFrame(overlap_rows,columns=['method','set','methylation_genes','rna_strict_overlap','rna_broad_overlap','overlap_published_union','rna_strict_genes']).to_csv(d/'overlap.tsv',sep='\t',index=False)
