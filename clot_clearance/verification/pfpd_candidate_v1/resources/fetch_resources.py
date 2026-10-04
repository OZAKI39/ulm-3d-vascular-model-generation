from pathlib import Path
import urllib.request,json,hashlib,subprocess,tarfile,datetime,re
O=Path(__file__).resolve().parents[1]; E=O/'external';R=O/'resources'
def get(url):
 return urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'clot-pfpd-reproducibility-audit'}),timeout=90).read()
versions=json.loads(get('https://julialang-s3.julialang.org/bin/versions.json'))
v=sorted([v for v,d in versions.items() if v.startswith('1.10.') and d['stable']],key=lambda x:tuple(map(int,x.split('.'))))[-1]
f=next(f for f in versions[v]['files'] if f['triplet']=='x86_64-linux-gnu' and f['extension']=='tar.gz')
(R/'JULIA_DOWNLOAD.json').write_text(json.dumps(dict(f, selected_version=v),indent=2))
p=E/f'julia-{v}.tar.gz';p.write_bytes(get(f['url']));assert hashlib.sha256(p.read_bytes()).hexdigest()==f['sha256']
with tarfile.open(p) as t:t.extractall(E,filter='data')
print('JULIA READY',v,flush=True)
