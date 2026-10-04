from pathlib import Path
import urllib.request,urllib.parse,json,hashlib,subprocess,concurrent.futures
R=Path(__file__).resolve().parent
repos=['BBVV.jl','BBVV1D.jl','FiniteElements.jl','mfo2024_xwave','paper_gamm_mitteilungen_2025','PeriFast.jl','SimZoo']
def one(name):
 p=R/'external'/name;url='https://github.com/kaipartmann/'+name+'.git'
 q=subprocess.run(['git','clone','--depth','1',url,str(p)],capture_output=True,text=True)
 if q.returncode:return dict(repo=name,error=q.stderr)
 sha=subprocess.check_output(['git','-C',str(p),'rev-parse','HEAD'],text=True).strip()
 q=subprocess.run(['git','-C',str(p),'grep','-I','-n','-i','-E',r'PFPD|phase[ _-]?field|kinematic[ _-]?degradation|crack[ _-]?driving|critical[ _-]?phase|\bYc\b'],capture_output=True,text=True)
 (R/(name+'_keyword_search.txt')).write_text(q.stdout)
 return dict(repo=name,url=url,sha=sha,keyword_hits=len(q.stdout.splitlines()))
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:rows=list(ex.map(one,repos))
(R/'AUTHOR_OTHER_REPO_CHECK.json').write_text(json.dumps(rows,indent=2))
fetches=[]
for url in ['https://www.crm.sns.it/media/course/3060/miehe+hofacker+welschinger10.pdf','http://crm.sns.it/media/course/3060/miehe+hofacker+welschinger10.pdf','https://zenodo.org/api/records?q=Partmann&size=100','https://zenodo.org/api/records/8187477/versions?size=100']:
 try:
  raw=urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'clot-at1-audit'}),timeout=30).read()
  if raw.startswith(b'%PDF'):
   name='Miehe_Hofacker_Welschinger_2010.pdf'
   (R/'MIEHE_PAPER_IDENTITY.json').write_text(json.dumps(dict(url=url,doi='10.1016/j.cma.2010.04.011',sha256=hashlib.sha256(raw).hexdigest()),indent=2))
  elif 'zenodo' in url:name='zenodo_author_current.json' if 'q=' in url else 'zenodo_all_versions.json'
  else:raise ValueError('Not a PDF')
  (R/name).write_bytes(raw);fetches.append(dict(url=url,status='DOWNLOADED',file=name))
 except Exception as e:fetches.append(dict(url=url,error=str(e)))
(R/'EXTRA_FETCH_LOG.json').write_text(json.dumps(fetches,indent=2));print(json.dumps(rows+fetches,indent=2))
