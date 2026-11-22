import sys,shutil
from pathlib import Path
import pandas as pd
d=Path(sys.argv[1]); motif_dir=Path(shutil.which('findMotifs.pl')).resolve().parents[1]/'data'/'knownTFs'
f=motif_dir/'vertebrates'/'known.motifs'
if not f.exists(): f=motif_dir/'known.motifs'
names=['Stat3','Stat3+il21','GABPA','ERG','ELF1','ETS1','Lhx2']
selected={}
record=False
for l in f.read_text().splitlines(keepends=True):
 if l.startswith('>'):
  motif=l.split('\t')[1].split('(')[0]
  record=motif in names
  if record: selected[motif]=[]
 if record: selected[motif].append(l)
assert len(selected)>=6,selected.keys()
(d/'motif'/'targeted.motifs').write_text(''.join(line for block in selected.values() for line in block))
