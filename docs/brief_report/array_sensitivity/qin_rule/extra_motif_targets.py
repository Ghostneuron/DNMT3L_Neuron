from pathlib import Path
import pandas as pd
import sys
from Bio import SeqIO
d=Path(sys.argv[1])
a=pd.read_pickle(d/'probe_data.pkl')
s5=pd.ExcelFile(d/'mmc7.xlsx')
gain=s5.parse(3).PositionID.unique()
loss=a.index[a.GEO_p.le(.05)&a.GEO_delta.le(-.1)]
assert len(gain)==385 and len(loss)==193 and set(gain)<=set(a.index)
print('385 gain archive exact, 193 loss candidates reconstructed, missing original hit IDs',len(set(s5.parse(1).PositionID.unique())-set(loss)))
out={'Archive_gain':a.loc[gain],'Archive_loss_reconstructed':a.loc[loss]}
for method in ['ssNoob_SWAN','ssNoob_dye']:
 for sign in ['gain','loss']:
  take=a[f'{method}_delta']>=.1 if sign=='gain' else a[f'{method}_delta']<=-.1
  out[f'{method}_broad_{sign}']=a.loc[take&a[f'{method}_p'].le(.05)]
genome=SeqIO.parse(d/'hg19.2bit','twobit')
for name,f in out.items():
 with (d/'motif'/f'{name}.fa').open('w') as o:
  for chrom,g in f.groupby('chr'):
   ch=genome[chrom]
   for probe,pos in g.pos.items():
    center=int(pos)-1
    seq=str(ch.seq[center-100:center+100]).upper()
    if len(seq)==200 and set(seq)<=set('ACGT'): o.write(f'>{probe}\n{seq}\n')
