from pathlib import Path
import urllib.request,urllib.parse,json,subprocess,re,datetime,hashlib
R=Path(__file__).resolve().parent;O=R.parent;repo=R/'external/Peridynamics.jl'
def fetch(url):return urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'clot-at1-reproducibility-audit'}),timeout=90).read()
def save_json(name,value):(R/name).write_text(json.dumps(value,indent=2))
def paged(path):
 result=[]
 for i in range(1,30):
  sep='&' if '?' in path else '?';data=json.loads(fetch('https://api.github.com/'+path+f'{sep}per_page=100&page={i}'));result+=data
  if len(data)<100:return result
 raise RuntimeError('Pagination incomplete')
metadata={}
for name,path in [('branches','repos/kaipartmann/Peridynamics.jl/branches'),('tags','repos/kaipartmann/Peridynamics.jl/tags'),('releases','repos/kaipartmann/Peridynamics.jl/releases'),('pulls','repos/kaipartmann/Peridynamics.jl/pulls?state=all'),('author_repositories','users/kaipartmann/repos')]:
 data=paged(path);save_json('current_'+name+'.json',data);metadata[name]=len(data)
sha=subprocess.check_output(['git','-C',str(repo),'rev-parse','origin/main'],text=True).strip();subprocess.run(['git','-C',str(repo),'checkout','--detach',sha],check=True)
refs=subprocess.check_output(['git','-C',str(repo),'for-each-ref','--format=%(refname)','refs/remotes','refs/tags'],text=True).splitlines()
pattern=r'PFPD|phase[ _-]?field|kinematic[ _-]?degradation|crack[ _-]?driving|\bYc\b|critical[ _-]?phase|\bsc\b'
hits=[]
for ref in refs:
 p=subprocess.run(['git','-c','grep.threads=1','-C',str(repo),'grep','-I','-n','-i','-E',pattern,ref,'--','src','docs','test','examples'],capture_output=True,text=True)
 if p.returncode not in [0,1]:raise RuntimeError(p.stderr)
 if p.stdout:hits.append(p.stdout)
(R/'AUTHOR_ALL_REFS_SEARCH.txt').write_text('\n'.join(hits))
html=fetch('https://arxiv.org/abs/2604.09215').decode();(R/'arxiv_current.html').write_text(html)
for name,url in [('zenodo_author','https://zenodo.org/api/records?q='+urllib.parse.quote('creators.name:"Partmann, Kai"')+'&size=100'),('zenodo_peridynamics','https://zenodo.org/api/records/8187477')]:
 try:(R/(name+'.json')).write_bytes(fetch(url))
 except Exception as e:save_json(name+'_fetch_error.json',dict(url=url,error=str(e)))
# Full peer-reviewed paper, hosted by Scuola Normale Superiore course materials.
url='https://crm.sns.it/media/course/3060/miehe%2Bhofacker%2Bwelschinger10.pdf'
try:
 raw=fetch(url);assert raw.startswith(b'%PDF');(R/'Miehe_Hofacker_Welschinger_2010.pdf').write_bytes(raw)
 save_json('MIEHE_PAPER_IDENTITY.json',dict(url=url,doi='10.1016/j.cma.2010.04.011',sha256=hashlib.sha256(raw).hexdigest()))
except Exception as e:save_json('MIEHE_FETCH_ERROR.json',dict(url=url,error=str(e)))
repos={}
for name in ['Peridynamics.jl','phasefieldx','fenicsx-fracture']:
 p=R/'external'/name;license=next(iter(p.glob('LICENSE*')),None)
 repos[name]=dict(commit=subprocess.check_output(['git','-C',str(p),'rev-parse','HEAD'],text=True).strip(),url=subprocess.check_output(['git','-C',str(p),'remote','get-url','origin'],text=True).strip(),license_file=str(license.relative_to(R)) if license else None,license_sha256=hashlib.sha256(license.read_bytes()).hexdigest() if license else None)
 if license:(R/(name+'_LICENSE.txt')).write_bytes(license.read_bytes())
save_json('AUTHOR_RESOURCE_SNAPSHOT.json',dict(checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),metadata_counts=metadata,author_main_sha=sha,refs=refs,repositories=repos,arxiv_revision_markers=re.findall(r'\[v\d+\]',html)))
print(json.dumps(dict(metadata_counts=metadata,repositories=repos,refs=len(refs)),indent=2))
