from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
from pathlib import Path
import gzip, re, subprocess, sys, time

soft, out=map(Path,sys.argv[1:])
urls=re.findall(r'^!Sample_supplementary_file = (ftp://[^\s]+\.idat\.gz)',soft.read_text(),re.M)
assert len(set(urls))==12, len(urls)
def fetch(url):
    name=url.split('/')[-1]
    f=out/name
    u=url.replace('ftp://','https://',1)
    for address in (u,u+'?download=1'):
        proc=subprocess.run(['curl','-fLsS','--connect-timeout','15','--max-time','70',address,'-o',str(f)])
        if not proc.returncode and f.stat().st_size>1000000:
            with gzip.open(f,'rb') as stream:
                assert stream.read(4)==b'IDAT',name
            return f"{sha256(f.read_bytes()).hexdigest()}  {name}"
    raise RuntimeError(name)
with ThreadPoolExecutor(max_workers=6) as runner:
    records=list(runner.map(fetch,urls))
(out/'SHA256SUMS.txt').write_text('\n'.join(sorted(records))+'\n')
print(len(records),'verified IDAT channel files')
