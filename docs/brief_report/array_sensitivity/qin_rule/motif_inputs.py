import sys
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import SeqIO
d=Path(sys.argv[1])
fetal=Path(sys.argv[2])
data=pd.read_pickle(d/'probe_data.pkl')
mask=data.detect & ~(data.MASK_mapping|data.MASK_snp5_common|data.MASK_sub30_copy)
a=data.loc[mask,['chr','pos']].copy()
w4=pd.ExcelFile(d/'mmc4.xlsx')
published=w4.parse(0,index_col='Name')
ds=w4.parse(1,index_col='TARGETID')
pub_xy=published[["CHR","MAPINFO"]].rename(columns={"CHR":"chr","MAPINFO":"pos"})
pub_xy.chr="chr"+pub_xy.chr.astype(str)
ds_xy=ds[['CHR','MAPINFO']].rename(columns={'CHR':'chr','MAPINFO':'pos'})
ds_xy.chr='chr'+ds_xy.chr.astype(str)
manifestcols=['probeID','CpG_chrm','CpG_beg','MASK_snp5_common','MASK_mapping']
hm=pd.read_csv(fetal/'HM450_Zhou_manifest.tsv.gz',sep='\t',usecols=manifestcols,index_col='probeID')
ep=pd.read_csv(d/'zhou_epic.tsv.gz',sep='\t',usecols=manifestcols,index_col='probeID')
ep=ep.loc[ep.index.str.startswith('cg') & ep.CpG_chrm.isin([f'chr{i}' for i in range(1,23)]) & ~(ep.MASK_snp5_common.eq(True)|ep.MASK_mapping.eq(True))]
hm=hm.loc[hm.index.str.startswith('cg') & hm.CpG_chrm.isin([f'chr{i}' for i in range(1,23)]) & ~(hm.MASK_snp5_common.eq(True)|hm.MASK_mapping.eq(True))]
bg=pd.concat([ep,hm.loc[~hm.index.isin(ep.index)]])
bg=bg.rename(columns={'CpG_chrm':'chr','CpG_beg':'pos'})[['chr','pos']]
bg['pos']+=1
sets={}
for name,df,delta in [('Qin_Table_S2',pub_xy,published['D-Z.Beta']),('DS_Table_S2',ds_xy,ds['T2118B.CON18B'])]:
 for direction,sel in [('gain',delta.gt(0)),('loss',delta.lt(0))]: sets[f'{name}_{direction}']=df.loc[sel,['chr','pos']]
for method in ['GEO','SWAN','ssNoob_SWAN','ssNoob_dye']:
 f=data.loc[mask & data[f'{method}_p'].le(.05)]
 for direction,sel in [('gain',f[f'{method}_delta'].ge(.1)),('loss',f[f'{method}_delta'].le(-.1))]:
  sets[f'{method}_{direction}']=f.loc[sel,['chr','pos']]
# Reserve all foreground probes across runs to hold background fixed.
excluded=set().union(*(set(s.index) for s in sets.values()))
bg=bg.loc[~bg.index.isin(excluded)].sort_index()
genome=SeqIO.parse(d/'hg19.2bit','twobit')
def write_fasta(path,coords):
 valid=0
 with path.open('w') as out:
  for chrom, grp in coords.groupby('chr'):
   reference=genome[chrom]
   for probe,pos in zip(grp.index,grp.pos):
    center=int(pos)-1
    seq=str(reference.seq[center-100:center+100]).upper()
    if len(seq)!=200 or set(seq)-set('ACGT'): continue
    out.write(f'>{probe}\n{seq}\n')
    valid+=1
 return valid
motif=d/'motif'
motif.mkdir(exist_ok=True)
rows=[]
for name,coords in sets.items():
 rows.append([name,len(coords),write_fasta(motif/f'{name}.fa',coords)])
# Independent reproducible subset of the array-wide union, comparable in size
# to the archived HOMER background. Full bg is retained for a known-motif check.
counts=[len(bg),write_fasta(motif/'background.fa',bg.sample(n=60000,random_state=20261120))]
rows.append(['background',*counts])
pd.DataFrame(rows,columns=['set','input','valid_sequences']).to_csv(motif/'input_counts.tsv',sep='\t',index=False)
