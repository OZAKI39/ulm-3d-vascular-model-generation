"""WSL-only selective terminal export, with source hashes and explicit omissions."""
from pathlib import Path
import hashlib,json,subprocess,tempfile,datetime
from monitor_remote import SSH
R=Path(__file__).resolve().parents[3]
D=Path(__file__).resolve().parents[1]
REMOTE=r'''
from pathlib import Path
import json,hashlib
r=Path('/workspace/y_clot/project');task='clot_irregular_lps_discretization_v4';run='20261019T221000Z'
o=r/'verification'/task/run;d=r/'analysis'/task;f=r/'results'/task/run;c=o/'cases/H_IRREGULAR_LPS_V4'
assert json.loads((o/'POSTPROCESS_STATUS.json').read_text())['status']=='COMPLETE'
assert json.loads((o/'DELIVERY_VALIDATION.json').read_text())['status']=='PASS'
selected={};omitted={}
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as stream:
  for b in iter(lambda:stream.read(2**20),b''):h.update(b)
 return h.hexdigest()
def add(p,limit=True):
 key=str(p.relative_to(r))
 if p.is_symlink():
  import os
  target=os.readlink(p);assert not target.startswith('/');selected[key]=dict(symlink=target);return
 if not p.is_file():return
 size=p.stat().st_size
 if limit and size>=32*2**20:omitted[key]=dict(bytes=size,reason='Large raw evidence retained on server; compact evidence exported');return
 selected[key]=dict(bytes=size,sha256=digest(p))
for base in [o,d]:
 for p in base.iterdir():
  if p.suffix in ['.json','.csv','.md','.log','.txt','.html'] or p.is_symlink():add(p)
for p in f.rglob('*'):
 if p.is_file():add(p,limit=False)
for p in c.iterdir():
 if p.suffix in ['.json','.csv','.txt']:add(p)
for name in ['initial','pilot','final','first_softening','first_failure','first_free_component','max_target_count']:
 add(c/'keys'/(name+'.npz'))
for name in ['last_accepted.npz','failed_native_trial.npz']:add(c/name)
for base in [c/'checkpoints',c/'keys',o/'build']:
 for p in base.rglob('*'):
  if p.is_file() and str(p.relative_to(r)) not in selected:
   omitted.setdefault(str(p.relative_to(r)),dict(bytes=p.stat().st_size,reason='Full native/build or dense accepted-state archive retained on server'))
for name in ['damage_events.jsonl.gz','CONNECTED_COMPONENT_HISTORY.jsonl.gz']:
 p=c/name
 if p.exists():omitted[str(p.relative_to(r))]=dict(bytes=p.stat().st_size,reason='Dense FP64 damage/member history retained on server; scalar curves and exact lineage exported')
print(json.dumps(dict(server_root=str(r),selected=selected,omitted=omitted)))
'''
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(2**20),b''):h.update(block)
 return h.hexdigest()
def main():
 q=subprocess.run(SSH+['python3','-'],input=REMOTE,capture_output=True,text=True,timeout=60)
 if q.returncode:raise RuntimeError(q.stderr)
 manifest=json.loads(q.stdout)
 with tempfile.NamedTemporaryFile() as listing:
  listing.write(b'\0'.join(s.encode() for s in manifest['selected'])+b'\0');listing.flush()
  subprocess.run(['rsync','-a','--from0','--files-from='+listing.name,'-e','ssh -o BatchMode=yes -o ConnectTimeout=15 -p 15283 -i /home/lzy/.ssh/vast_step3b_ed25519','root@95.3.33.46:/workspace/y_clot/project/',str(R)+'/'],check=True)
 for name,record in manifest['selected'].items():
  p=R/name
  if 'symlink' in record:assert p.is_symlink() and str(p.readlink())==record['symlink'],name
  else:assert p.stat().st_size==record['bytes'] and digest(p)==record['sha256'],name
 manifest.update(status='SELECTED_FILES_VERIFIED',verified_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),git_policy='Local export includes requested videos; Git selection is a separate reviewed operation, without native libraries or dense archives.')
 (D/'LOCAL_EXPORT_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
 print('LOCAL_EXPORT_VERIFIED',len(manifest['selected']),'files;',len(manifest['omitted']),'explicit server archive entries')
if __name__=='__main__':main()
