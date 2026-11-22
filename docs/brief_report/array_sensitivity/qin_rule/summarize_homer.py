import math,re,sys
import pandas as pd
from pathlib import Path
d=Path(sys.argv[1])
rows=[]
for f in sorted(d.glob('motif/*/knownResults.txt')):
 t=pd.read_csv(f,sep='\t')
 name=t.columns[0]
 q=next(c for c in t if c.startswith('q-value'))
 target=next(c for c in t if c.startswith('# of Target'))
 n=int(re.search(r'\(of (\d+)\)',target).group(1))
 h=(t[q]<.05).sum()
 for i,r in t.iterrows():
  consensus=str(r['Consensus'])
  motif=r[name].split('(')[0]
  if motif.lower()=='stat3' and consensus=='CTTCCGGGAA': motif='STAT3'
  if motif not in ['STAT3','Stat3+il21','GABPA','ERG','ELF1','ETS1','Lhx2'] and i>=8: continue
  rows.append([f.parent.name,i+1,motif,consensus,math.exp(float(r['Log P-value'])),r[q],r[target],n,int(h)])
frame=pd.DataFrame(rows,columns=['set','rank','motif','consensus','p','q','targets_with_motif','targets','FDR_hits'])
frame.to_csv(d/'homer_summary.tsv',sep='\t',index=False)
# De novo pages include the ranked motif IDs, best matches and percentages.
for p in sorted(d.glob('motif/*/homerResults.html')):
 t=pd.read_html(p)[0]
 print(p.parent.name,t.head(4).to_string(index=False),sep='\n')
