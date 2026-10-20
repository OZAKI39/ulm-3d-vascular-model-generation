"""Prepare, hash and list a compact local Git delivery; does not stage or commit."""
from pathlib import Path
import json,hashlib,gzip,os,datetime
R=Path(__file__).resolve().parents[3];TASK='clot_irregular_lps_discretization_v4';RUN='20261019T221000Z'
D=R/'analysis'/TASK;O=R/'verification'/TASK/RUN;F=R/'results'/TASK/RUN
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
def write(p,obj):p.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n')
def main():
 assert json.loads((D/'DELIVERY_VALIDATION.json').read_text())['status']=='PASS'
 assert json.loads((D/'NOTEBOOK_UPDATE.json').read_text())['status']=='UPDATED_WITH_HASH_GUARD_AND_BACKUP'
 spec=json.loads((O/'GIT_COMPRESSED_EVIDENCE.json').read_text());records={q['source']:q for q in spec['records']};transforms=[]
 for p in sorted(O.iterdir()):
  if p.is_file() and p.suffix in ['.json','.csv','.log'] and p.stat().st_size>2**20:
   if p.name in records:
    previous=records[p.name]
    if sha(p)==previous['source_sha256'] and sha(O/previous['compressed'])==previous['compressed_sha256']:continue
   raw=p.read_bytes();dest=p.with_name(p.name+'.gz');encoded=gzip.compress(raw,compresslevel=6,mtime=0)
   assert gzip.decompress(encoded)==raw;dest.write_bytes(encoded)
   records[p.name]=dict(source=p.name,source_bytes=len(raw),source_sha256=sha(p),compressed=dest.name,compressed_bytes=len(encoded),compressed_sha256=sha(dest))
 write(O/'GIT_COMPRESSED_EVIDENCE.json',dict(lossless=True,records=list(records.values())))
 # Identical tables in the friendly analysis entrypoint and native case need no duplicate blobs.
 pairs=[(D/name,O/name) for name in ['ENERGY_LEDGER.csv','FRAGMENT_HISTORY.csv','FRAGMENT_LINEAGE.json','FRAGMENT_MORPHOLOGY.csv']]
 pairs += [(O/'cases/H_IRREGULAR_LPS_V4/history.csv',O/'FRAGMENT_HISTORY.csv'),(O/'cases/H_IRREGULAR_LPS_V4/LINEAGE.json',O/'FRAGMENT_LINEAGE.json')]
 for alias,target in pairs:
  assert sha(alias)==sha(target),f'Duplicate contents differ: {alias}'
  if not alias.is_symlink():alias.unlink();alias.symlink_to(os.path.relpath(target,alias.parent))
  transforms.append(dict(path=str(alias.relative_to(R)),symlink=os.readlink(alias),content_sha256=sha(target)))
 selected={};omitted={};compressed_sources={str((O/q['source']).relative_to(R)) for q in records.values()}
 def add(p):
  key=str(p.relative_to(R))
  if '__pycache__' in p.parts or p.suffix in ['.pyc','.o','.so']:return
  if p.is_symlink():selected[key]=dict(symlink=os.readlink(p));return
  if not p.is_file():return
  if key in compressed_sources:omitted[key]=dict(bytes=p.stat().st_size,sha256=sha(p),reason='Losslessly compressed in GIT_COMPRESSED_EVIDENCE.json');return
  assert p.stat().st_size<32*2**20,f'Unexpected large Git payload: {p}'
  selected[key]=dict(bytes=p.stat().st_size,sha256=sha(p))
 for p in (D/'code').rglob('*'):add(p)
 for p in D.iterdir():
  if p.name in ['GIT_EXPORT_MANIFEST.json','GIT_STAGE_PATHS.txt']:continue
  if p.is_file() or p.is_symlink():add(p)
 for p in F.rglob('*'):add(p)
 for p in O.iterdir():
  if p.is_file() or p.is_symlink():add(p)
 for p in (O/'provenance').rglob('*'):
  if p.name not in ['git_before.diff','git_before_status.txt']:add(p)
 case=O/'cases/H_IRREGULAR_LPS_V4'
 for p in case.iterdir():
  if p.suffix in ['.json','.csv','.txt']:add(p)
 keep={'initial.npz','final.npz','first_target_observed.npz','max_target_count.npz'}
 for p in (case/'keys').iterdir():
  if p.name in keep:add(p)
  elif p.is_file():omitted[str(p.relative_to(R))]=dict(bytes=p.stat().st_size,sha256=sha(p),reason='Additional event/pilot snapshot retained locally and on server; selected initial, first-target, peak-target and final states are in Git')
 for p in (case/'last_accepted.npz',):
  if p.exists():omitted[str(p.relative_to(R))]=dict(bytes=p.stat().st_size,sha256=sha(p),reason='Duplicate final accepted state retained locally/server')
 for p in (O/'p0').rglob('*'):
  if p.is_file():omitted[str(p.relative_to(R))]=dict(bytes=p.stat().st_size,reason='Raw profiling evidence retained outside compact Git; full timing/hash provenance in PERFORMANCE_BASELINE_V3_COMPACT.json')
 manifest=dict(status='PREPARED_AND_HASHED',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),selected=selected,omitted=omitted,local_export_alias_transforms=transforms,restore='Run code/restore_export.py to recover all compressed original bytes before analysis.',server_archive=str(O).replace(str(R),'/workspace/y_clot/project'),no_github_push=True)
 write(D/'GIT_EXPORT_MANIFEST.json',manifest)
 paths=sorted(selected)+[str((D/'GIT_EXPORT_MANIFEST.json').relative_to(R)),str((D/'GIT_STAGE_PATHS.txt').relative_to(R))]
 (D/'GIT_STAGE_PATHS.txt').write_text('\n'.join(paths)+'\n')
 print('COMPACT_GIT_PREPARED',len(paths),'paths',sum(q.get('bytes',0) for q in selected.values()),'bytes')
if __name__=='__main__':main()
