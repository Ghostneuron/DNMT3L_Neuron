import sys,numpy as np,pandas as pd
from scipy.stats import pearsonr,ttest_ind
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, fcluster, leaves_list
from pathlib import Path
d=Path(sys.argv[1])
a=pd.read_pickle(d/'probe_data.pkl')
w4=pd.ExcelFile(d/'mmc4.xlsx')
w5=pd.ExcelFile(d/'mmc5.xlsx')
w7=pd.ExcelFile(d/'mmc7.xlsx')
pub=w4.parse(0,index_col='Name')
a['Qin_selected']=a.index.isin(pub.index)
a['Zhou_selected']=a.detect & ~(a.MASK_mapping|a.MASK_snp5_common|a.MASK_sub30_copy)
a['All_detected']=a.detect
methods=['GEO','SWAN','ssNoob_SWAN','ssNoob_dye','ssNoob_0']
if (d/'background_swan_betas.tsv.gz').exists():
 bg=pd.read_csv(d/'background_swan_betas.tsv.gz',sep='\t',index_col='probe_id')
 m='Noob_background_SWAN'
 bg=bg.reindex(a.index)
 dat=bg.to_numpy()
 a[bg.columns]=bg
 a[f'{m}_delta']=dat[:,3:].mean(axis=1)-dat[:,:3].mean(axis=1)
 a[f'{m}_p']=ttest_ind(dat[:,3:],dat[:,:3],axis=1,equal_var=True).pvalue
 methods.insert(2,m)
chunks=[]; stats=[]; locations=[]; f1f=[]; fig1b=[]; clusters=[]; fig1e=[]
target_motif=w7.parse(12,index_col='Name')
for label,subset in [('Published 347 (fixed)',a.Qin_selected),('Zhou filter',a.Zhou_selected),('Detection only',a.All_detected),('Full autosomal joint support',np.ones(len(a),bool))]:
 for method in methods:
  beta=a[f'{method}_delta']
  prob=a[f'{method}_p']
  c=a.loc[subset & prob.le(.05) & beta.abs().ge(.1)].copy()
  c['delta']=beta.loc[c.index]
  c['p']=prob.loc[c.index]
  c['gain']=c.delta.gt(0)
  q21=beta[subset & a.chr.eq('chr21') & prob.le(.05)]
  published_found=c.index.intersection(pub.index)
  stats.append([label,method,int(np.sum(subset)),len(c),c.gain.sum(),sum(~c.gain),len(published_found),sum(pub.loc[published_found,'D-Z.Beta']>0),len(c.index.intersection(target_motif.index)),sum(q21>0),sum(q21<0)])
  for chr_,t in c.groupby('chr'):
   fig1b.append([label,method,chr_,sum(t.gain),sum(~t.gain)])
  for col,categories in [('Relation_to_Island',['N_Shelf','N_Shore','Island','S_Shore','S_Shelf','OpenSea']),('UCSC_RefGene_Group',['TSS1500','TSS200',"5'UTR",'1stExon','Body',"3'UTR",'Intergene'])]:
   for region in categories:
    annotations=c[col].fillna('OpenSea' if col=='Relation_to_Island' else 'Intergene')
    hits=annotations.str.contains(r'(?:^|;)'+region+r'(?:;|$)',regex=True)
    locations.append([label,method,region,int((hits&c.gain).sum()),int((hits&~c.gain).sum()),int(hits.sum()),len(c)])
  if label not in ('Published 347 (fixed)','Zhou filter') or method not in ('GEO','SWAN','ssNoob_SWAN','ssNoob_dye'): continue
  chunks.append(c[['chr','pos','UCSC_RefGene_Name','Type','Color','delta','p','gain',*[f'{method}_{i}' for i in range(1,7)]]].assign(method=method,set=label).rename(columns={f'{method}_{i}':sample for i,sample in enumerate(['Z1','Z2','Z3','D1','D2','D3'],1)}))
  for s_label,keys in [('Fixed published surviving',a.index[a.Qin_selected]),('Newly selected',c.index)]:
   v=a.loc[keys,[f'{method}_{i}' for i in range(1,7)]].to_numpy()
   l=linkage(pdist(v.T),method='complete')
   order=leaves_list(l)
   split=fcluster(l,2,criterion='maxclust')
   is_split=np.unique(split[:3]).size==1 and np.unique(split[3:]).size==1 and split[0]!=split[3]
   dists=pdist(v.T)
   fig1e.append([method,label,s_label,len(v),bool(is_split),','.join(np.array(['Z1','Z2','Z3','D1','D2','D3'])[order])])
  for category in ['cell junction probes','synapse probes']:
   f=w4.parse(category).rename(columns={'TargetID.y':'probe'})
   combined=f[['UCSC_REFGENE_NAME','probe']].drop_duplicates().join(a[[f'{method}_delta',f'{method}_p']],on='probe')
   combined['sig']=combined[f'{method}_p'].le(.05)
   combined['strict']=combined.sig & combined[f'{method}_delta'].abs().ge(.1)
   mask1=combined.sig & combined[f'{method}_delta'].gt(0)
   mask2=combined.sig & combined[f'{method}_delta'].lt(0)
   f1f.append([label,method,category,len(combined),sum(mask1),sum(mask2),sum(combined.strict)])
   combined['method']=method; combined['category']=category
   combined=combined.rename(columns={f'{method}_delta':'delta',f'{method}_p':'p'})
   if label=='Published 347 (fixed)': clusters.append(combined)
res=pd.DataFrame(stats,columns=['universe','method','eligible','called','gains','losses','published_sites_retained','published_gains_retained','T5_STAT3_probes_passing','chr21_nominal_gains','chr21_nominal_losses'])
res.to_csv(d/'figure1_counts.tsv',sep='\t',index=False)
pd.DataFrame(locations,columns=['universe','method','annotation','gains','losses','annotated','total_calls']).to_csv(d/'figure1D.tsv',sep='\t',index=False)
pd.DataFrame(fig1b,columns=['universe','method','chr','gains','losses']).to_csv(d/'figure1B.tsv',sep='\t',index=False)
pd.DataFrame(fig1e,columns=['method','universe','cpg_set','n','clean_group_split','sample_order']).to_csv(d/'figure1E.tsv',sep='\t',index=False)
pd.DataFrame(f1f,columns=['universe','method','category','probe_gene_pairs','nominal_gains','nominal_losses','at_strict_cutoff']).to_csv(d/'figure1F.tsv',sep='\t',index=False)
pd.concat(clusters).to_csv(d/'figure1F_probes.tsv',sep='\t',index=False)
pd.concat(chunks).to_csv(d/'figure1A_probes.tsv',sep='\t',index_label='probe')
corr=w5.parse(12,index_col='GeneSymbol')
out=[]
for reg,col in [('Body','Methylation difference in gene body'),('Promoter','Methylation difference in promoter')]:
 ids=corr.index[corr[col].abs().le(.01)&corr.log2FoldChange.abs().le(1)]
 for meth in ['Published cells','GEO','SWAN','ssNoob_SWAN','ssNoob_dye']:
  v=corr[['log2FoldChange',col]].rename(columns={'log2FoldChange':'expression',col:'beta_delta'}) if meth=='Published cells' else pd.read_csv(d/f'corr_{meth}_{reg}_Pub.tsv',sep='\t',index_col=0)
  modes={'Fixed published display-range genes':v.reindex(ids).dropna(),'Full 225-gene panel':v.dropna(),'Reapply plotted range':v[v.beta_delta.abs().le(.01)&v.expression.abs().le(1)]}
  for mode,f in modes.items():
   t=pearsonr(f.expression,f.beta_delta)
   out.append([reg,meth,mode,len(f),t.statistic,t.pvalue])
pd.DataFrame(out,columns=['region','method','scope','n','r','p']).to_csv(d/'figure2E.tsv',sep='\t',index=False)
g=w5.parse(14); h=w5.parse(16)
TF=set(g.TF.dropna()); ds=set(h['DS microarray'].dropna())
short=set(w5.parse(9).RNAseq.dropna()); broad=set(w5.parse(10).RNAseq.dropna())
pd.DataFrame([['TF',len(TF)],['DS microarray',len(ds)],['Strict RNA-seq',len(short)],['Broad RNA-seq',len(broad)],['Strict RNA-seq and TF',len(short&TF)],['Broad RNA-seq and TF',len(broad&TF)],['Broad RNA-seq and DS',len(broad&ds)],['TF and DS',len(TF&ds)],['All three',len(TF&ds&broad)]],columns=['set','genes']).to_csv(d/'figure2F.tsv',sep='\t',index=False)
