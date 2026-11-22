import sys,numpy as np,pandas as pd
from scipy.stats import pearsonr
from pathlib import Path
d=Path(sys.argv[1]); a=pd.read_pickle(d/'probe_data.pkl')
expected=pd.read_excel(d/'mmc5.xlsx',sheet_name=12,index_col='GeneSymbol')
by_gene={}
for i, names in enumerate(a.UCSC_RefGene_Name.fillna('')):
 for gene in set(names.split(';')):
  if gene in expected.index: by_gene.setdefault(gene,[]).append(i)
rows=[]; outputs={}; equal=[]
for (region, col, label) in [({'Body'},'Methylation difference in gene body','Body'),({"TSS200","5'UTR"},'Methylation difference in promoter','Promoter')]:
 for universe,u in [('All mapped',np.ones(len(a),bool)),('Detection',a.detect.to_numpy()),('Zhou', (a.detect & ~(a.MASK_mapping|a.MASK_snp5_common|a.MASK_sub30_copy)).to_numpy())]:
  for match in ['Gene and region paired','Any region at matching gene']:
   results=[]
   for gene,indices in by_gene.items():
    subsets=a.iloc[indices]
    chrom='chr'+str(expected.loc[gene,'Chr'])
    masks=[(bool(set(g.split(';'))&region) if match.startswith('Any') else any(n==gene and tag in region for n,tag in zip(names.split(';'),g.split(';')))) for names,g in subsets[['UCSC_RefGene_Name','UCSC_RefGene_Group']].fillna('').itertuples(index=False,name=None)]
    take=np.array(masks)&u[indices]&subsets.chr.eq(chrom).to_numpy()
    if take.any():
     results.append([gene,take.sum(), *[subsets[f'{method}_delta'].to_numpy()[take].mean() for method in ['GEO','SWAN','ssNoob_SWAN','ssNoob_dye']]])
   frame=pd.DataFrame(results,columns=['gene','probes','GEO','SWAN','ssNoob_SWAN','ssNoob_dye']).set_index('gene').join(expected[[col,'log2FoldChange']])
   can=frame[[col,'GEO']].dropna()
   err=np.abs(can[col]-can.GEO)
   rows.append([label,universe,match,len(err),np.median(err),np.sum(err<1e-8),pearsonr(can[col],can.GEO).statistic])
   outputs[(label,universe,match)]=frame
pd.DataFrame(rows,columns=['region','universe','match','n','med_abs_err','within_1e_8','r_of_effects']).to_csv(d/'correlation_mapping.tsv',sep='\t',index=False)
for (label,universe,match),frame in outputs.items():
 col='Methylation difference in gene body' if label=='Body' else 'Methylation difference in promoter'
 for method in ['GEO','SWAN','ssNoob_SWAN','ssNoob_dye']:
  complete=frame[[col,'log2FoldChange',method]].dropna()
  for scope,t in [('Published range fixed',complete.loc[complete[col].abs().le(.01)&complete.log2FoldChange.abs().le(1)]),('Full',complete),('Reapply ranges',complete.loc[complete[method].abs().le(.01)&complete.log2FoldChange.abs().le(1)])]:
   result=pearsonr(t[method],t.log2FoldChange)
   equal.append([label,universe,match,method,scope,len(t),result.statistic,result.pvalue])
 if universe=='All mapped' and match=='Any region at matching gene': frame.to_csv(d/f'figure2E_{label}_full.tsv',sep='\t')
pd.DataFrame(equal,columns=['region','universe','match','method','scope','N','r','P']).to_csv(d/'correlation_scope.tsv',sep='\t',index=False)
