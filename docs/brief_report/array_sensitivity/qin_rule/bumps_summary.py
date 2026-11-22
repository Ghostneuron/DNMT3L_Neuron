from pathlib import Path
import pandas as pd,numpy as np,sys
d=Path(sys.argv[1])
pub=pd.read_excel(d/'mmc4.xlsx',sheet_name=3)
pub=pub.assign(chr=pub.chr.str.replace('ch','',regex=False).astype(int))
rows=[]
for m,label in [('GEO','GEO'),('ssNoob_SWAN','ssNoob + SWAN'),('ssNoob_dye','ssNoob alone')]:
  a=pd.read_csv(d/f'bumps_{m}.tsv',sep='\t')
  shared=a.merge(pub,how='inner',on=['chr','start','end'],suffixes=('','_pub'))
  e=np.abs(shared.value-shared.value_pub)
  sig=a.loc[(a.value.abs()>=.1)&(a['p.value']<=.05)&(a.L>=2)]
  rows.append([label,len(a),len(sig),sum(sig.value>0),sum(sig.value<0),sum(sig.fwer<=.05),len(shared),sum(e<1e-6)])
out=pd.DataFrame(rows,columns=['workflow','total','calls','gains','losses','fwer_05','coordinates_shared','difference_within_1e_6'])
out.to_csv(d/'bumps_summary.tsv',sep='\t',index=False)
print(out.to_string(index=False))
