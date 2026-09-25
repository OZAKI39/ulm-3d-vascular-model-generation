from pathlib import Path
import json,csv,hashlib,shutil,gzip,os,subprocess
A=Path(__file__).resolve().parents[1];C=json.loads((A/'context.json').read_text());D=Path(C['destination']);M=D/'sync_metadata'/C['snapshot'];M.mkdir(parents=True,exist_ok=True)
SOURCE=Path('/home/lzy/projects/ulm_particle_population_inlet_p9a4');rows=[];packed=[];inherited=[]
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def copy(src,target,group):
 src=Path(src);dst=D/target;st=src.stat();h=sha(src)
 dst.parent.mkdir(parents=True,exist_ok=True)
 if not dst.exists() or sha(dst)!=h:shutil.copy2(src,dst)
 assert sha(dst)==h and src.stat().st_mtime_ns==st.st_mtime_ns
 rows.append(dict(source=str(src),target=str(target),group=group,bytes=st.st_size,sha256=h,representation='IDENTICAL_BYTES'))
old=D/'sync_metadata/network_h0_particle_20260925';packaging=set(json.loads((old/'packaging_edits.json').read_text()))
# Recheck every previously selected original working file; copy actual updates.
with (old/'local_file_inventory.csv').open() as f:
 for r in csv.DictReader(f):
  if r['status']!='INCLUDED' or not r.get('sha256'):continue
  p=Path(r['source']);entry=dict(source=str(p),target=r['target'],previous_sha256=r['sha256'])
  if not p.is_file():entry['status']='SOURCE_NO_LONGER_PRESENT_RETAINED_HISTORICAL_SNAPSHOT'
  else:
   h=sha(p);entry['current_sha256']=h
   if h==r['sha256']:entry['status']='UNCHANGED'
   elif r['target'] in packaging:entry['status']='NAVIGATION_CHANGED_SAVED_AS_PROVENANCE';copy(p,M.relative_to(D)/'original_docs'/r['target'],'UPDATED_SOURCE_NAVIGATION')
   else:entry['status']='CURRENT_SOURCE_UPDATE';copy(p,r['target'],'UPDATED_PRIOR_SELECTION')
  inherited.append(entry)
# P9-A.4 is an additive scientific stage. Preserve all outputs, including failures.
paths=[SOURCE/'particle_3d/contracts/P9A4_CONTINUOUS_INFUSION_V1.json',SOURCE/'particle_3d/src/particle_3d/continuous_infusion.py',SOURCE/'particle_3d/src/particle_3d/population_inlet_p9a4.py']
paths+=list((SOURCE/'particle_3d/tests/particle9a4_population_inlet').rglob('*'))+list((SOURCE/'particle_3d/reports/particle9a4_population_inlet').rglob('*'))
compress={'particle_3d/reports/particle9a4_population_inlet/data/inlet100k/proposal_ledger.jsonl','particle_3d/reports/particle9a4_population_inlet/data/inlet100k/accepted_births.json','particle_3d/reports/particle9a4_population_inlet/data/legacy_b_inlet_500.json'}
for p in sorted(paths):
 if not p.is_file() or '__pycache__' in p.parts:continue
 rel=str(p.relative_to(SOURCE))
 if rel not in compress:copy(p,rel,'P9A4');continue
 h=sha(p);q=D/(rel+'.gz');q.parent.mkdir(parents=True,exist_ok=True)
 with p.open('rb') as src,q.open('xb') as raw:
  with gzip.GzipFile(fileobj=raw,mode='wb',filename='',mtime=0,compresslevel=6) as z:shutil.copyfileobj(src,z,1024**2)
 with gzip.open(q,'rb') as z:assert hashlib.file_digest(z,'sha256').hexdigest()==h
 spec=dict(source=str(p),original_path=rel,packed_path=rel+'.gz',original_bytes=p.stat().st_size,packed_bytes=q.stat().st_size,original_sha256=h,packed_sha256=sha(q),codec='gzip')
 packed.append(spec);rows.append(dict(source=str(p),target=rel+'.gz',group='P9A4',bytes=p.stat().st_size,sha256=h,packed_sha256=spec['packed_sha256'],representation='GZIP_OF_IDENTICAL_BYTES'))
for p in [Path('/home/lzy/projects/WORKFLOW_PATHS_20260925_ZH.md'),*Path('/home/lzy/projects/workflow_path_inventory_20260925').glob('*')]:
 if p.is_file():copy(p,M.relative_to(D)/'path_inventory'/p.relative_to('/home/lzy/projects'),'PATH_INVENTORY')
copy('/home/lzy/projects/ACTIVE_VASCULAR_WORKFLOW.md',M.relative_to(D)/'source_workspace_index.md','SOURCE_NAVIGATION')
(M/'packed_artifacts.json').write_text(json.dumps(packed,indent=2)+'\n')
(M/'local_added_inventory.json').write_text(json.dumps(rows,indent=2)+'\n')
(M/'inherited_source_recheck.json').write_text(json.dumps(inherited,indent=2)+'\n')
# All inherited file hashes are valid lookup targets; new large data use lossless representations.
known={}
for line in (old/'SNAPSHOT_SHA256.txt').read_text().splitlines():
 h,p=line.split('  ',1);known.setdefault(h,dict(path=p,representation='IDENTICAL_BYTES'))
for r in rows:known[r['sha256']]=dict(path=r['target'],representation=r['representation'])
(A/'audit/known_content.json').write_text(json.dumps(known)+'\n')
from collections import Counter
summary=dict(inherited_source_checks=len(inherited),inherited_status=dict(Counter(r['status'] for r in inherited)),new_or_updated_files=len(rows),original_bytes=sum(r['bytes'] for r in rows),compressed_artifacts=packed)
(A/'audit/local_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)
