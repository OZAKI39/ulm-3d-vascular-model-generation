import json,urllib.request,concurrent.futures,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'evidence/official'; OUT.mkdir(exist_ok=True)
repos=['tryolabs/norfair','bonsai-rx/bonsai','facebookresearch/co-tracker','DeepLabCut/DeepLabCut-live','quantumjot/btrack','opencv/opencv','krontech/chronos-examples','bonsai-rx/sleap','bonsai-rx/zeromq']
def get(url):
 r=urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'microbubble-vision-research'}),timeout=25); return r.read()
def work(repo):
 tag=repo.replace('/','__'); result={'repo':repo,'url':'https://github.com/'+repo,'checked_client_date':'2026-10-16'}
 for name,endpoint in [('repo',''),('commits','/commits?per_page=1'),('release','/releases/latest')]:
  try:
   data=get('https://api.github.com/repos/'+repo+endpoint); (OUT/f'{tag}_{name}.json').write_bytes(data); obj=json.loads(data)
   if name=='repo':result.update(default_branch=obj['default_branch'],archived=obj['archived'],pushed_at=obj['pushed_at'],license=(obj.get('license') or {}).get('spdx_id'))
   if name=='commits':result.update(head=obj[0]['sha'],commit_date=obj[0]['commit']['committer']['date'])
   if name=='release':result.update(release=obj['tag_name'],release_date=obj['published_at'])
  except Exception as e:result[name+'_error']=str(e)
 for file in ['README.md','pyproject.toml','setup.py','LICENSE']:
  try:
   data=get(f'https://raw.githubusercontent.com/{repo}/{result.get("head",result.get("default_branch","main"))}/{file}'); (OUT/f'{tag}_{file.replace("/","_")}').write_bytes(data)
  except Exception:pass
 return result
with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex: rows=list(ex.map(work,repos))
(OUT/'repository_metadata.json').write_text(json.dumps(rows,indent=2)); print(json.dumps(rows,indent=2))
