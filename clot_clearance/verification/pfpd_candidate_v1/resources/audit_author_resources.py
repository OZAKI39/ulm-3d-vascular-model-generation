from pathlib import Path
import urllib.request,json,subprocess,datetime,re,hashlib
O=Path(__file__).resolve().parents[1];R=O/'resources';repo=O/'external/Peridynamics.jl'
def git(*args):return subprocess.check_output(['git','-C',str(repo),*args],text=True)
sha=git('rev-parse','HEAD').strip();subprocess.run(['git','-C',str(repo),'checkout','--detach',sha],check=True)
meta=dict(checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),repository='https://github.com/kaipartmann/Peridynamics.jl',sha=sha,version='0.6.0-DEV',license_sha256=hashlib.sha256((repo/'LICENSE').read_bytes()).hexdigest(),tags_at_sha=git('tag','--points-at','HEAD').splitlines())
for name,query in [('branches',''),('tags',''),('releases',''),('pulls','&state=all')]:
 data=[]
 for page in range(1,20):
  url=f'https://api.github.com/repos/kaipartmann/Peridynamics.jl/{name}?per_page=100&page={page}{query}'
  req=urllib.request.Request(url,headers={'User-Agent':'clot-pfpd-audit'})
  batch=json.load(urllib.request.urlopen(req,timeout=60));data+=batch
  if len(batch)<100:break
 (R/f'github_{name}.json').write_text(json.dumps(data,indent=2));meta[name+'_count']=len(data)
refs=git('for-each-ref','--format=%(refname)','refs/remotes/origin','refs/tags').splitlines()
pattern=r'PFPD|PhaseField|phase_field|critical_phase|\bsc\b|\bYc\b|crack_driving|kinematic_degradation'
hits=[]
for ref in refs:
 p=subprocess.run(['git','-c','grep.threads=1','-C',str(repo),'grep','-n','-i','-E',pattern,ref,'--','src','test','docs','examples'],text=True,capture_output=True)
 if p.returncode not in [0,1]:raise RuntimeError(p.stderr)
 if p.stdout:hits.append(p.stdout)
(R/'OFFICIAL_ALL_REFS_KEYWORD_SEARCH.txt').write_text('\n'.join(hits))
meta['refs_scanned']=refs;meta['keyword_pattern']=pattern
(R/'OFFICIAL_REPOSITORY.json').write_text(json.dumps(meta,indent=2))
try:
 html=urllib.request.urlopen('https://arxiv.org/abs/2604.09215',timeout=60).read().decode()
 (R/'arxiv_current.html').write_text(html)
except Exception as e:(R/'arxiv_fetch_error.txt').write_text(str(e))
print(json.dumps({k:v for k,v in meta.items() if k!='refs_scanned'},indent=2))
