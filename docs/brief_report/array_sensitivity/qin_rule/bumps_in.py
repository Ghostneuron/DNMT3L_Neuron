import sys
from pathlib import Path
import pandas as pd
d=Path(sys.argv[1])
a=pd.read_pickle(d/'probe_data.pkl')
keep=a.detect & ~(a.MASK_mapping|a.MASK_snp5_common|a.MASK_sub30_copy)
a=a.loc[keep].copy()
a['chrnum']=a.chr.str[3:].astype(int)
c=['chrnum','pos','UCSC_RefGene_Name',*[f'{m}_{i}' for m in ['GEO','ssNoob_SWAN','ssNoob_dye'] for i in range(1,7)]]
a.sort_values(['chrnum','pos'],kind='stable')[c].to_csv(d/'bumps_input.tsv.gz',sep='\t',index_label='probe')
