from pathlib import Path
import sys,re,math,pandas as pd
d=Path(sys.argv[1]); rows=[]
for f in d.glob("motif/*/targeted_results.motifs"):
 for line in f.read_text().splitlines():
  if not line.startswith('>'): continue
  fields=line.lstrip('>').split('\t')
  mat=re.search(r'T:([0-9.]+)\(([0-9.]+)%\),B:([0-9.]+)\(([0-9.]+)%\)',fields[-1])
  if not mat: continue
  x,t,y,b=map(float,mat.groups())
  motif=fields[1].split('(')[0]
  if motif=='Stat3':motif='STAT3'
  rows.append([f.parent.name,motif,fields[0],int(x),round(100*x/t) if t else 0,t,b,math.exp(float(fields[3]))])
stats=pd.DataFrame(rows,columns=['set','motif','consensus','targets_with_motif','targets','target_pct','background_pct','p'])
stats.to_csv(d/'targeted_summary.tsv',sep='\t',index=False)
print(stats[(stats.p<.05)|stats.motif.eq('STAT3')].sort_values(['set','p']).to_string(index=False))
