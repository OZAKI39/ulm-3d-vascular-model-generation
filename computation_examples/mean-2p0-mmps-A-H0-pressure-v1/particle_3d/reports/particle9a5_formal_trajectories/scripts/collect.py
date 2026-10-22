"""Incremental collection without deletion; final collection verifies remote bytes."""
from pathlib import Path
import argparse,json,subprocess,sys,shlex
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT/'particle_3d/src'))
from particle_3d.formal_cohort_p9a5 import *
R=ROOT/REL;c=json.loads((R/'data/run_context.json').read_text());remote=c['remote']
p=argparse.ArgumentParser();p.add_argument('--progress',action='store_true');a=p.parse_args()
(R/'logs/server').mkdir(exist_ok=True)
for source,target in [(remote+'/logs/',str(R/'logs/server')+'/'),(remote+'/'+REL+'/data/',str(R/'data')+'/')]:
 subprocess.run(['rsync','-a','--partial','--delay-updates','-e','ssh -o BatchMode=yes','vast4090:'+source,target],check=True)
if a.progress:
 subprocess.run(['rsync','-a','--partial','--delay-updates','--include=*/','--include=batch_result.json','--include=BATCH_COMPLETE.json','--exclude=*','-e','ssh -o BatchMode=yes','vast4090:'+remote+'/'+REL+'/outputs/',str(R/'outputs')+'/'],check=True)
else:
 code='''from pathlib import Path
import json,hashlib
p=Path(REMOTE)/REL
assert (p/'data/production_complete.json').is_file() and (p/'data/point_complete.json').is_file()
files={}
for f in sorted((p/'outputs').rglob('*')):
 if f.is_file():
  with f.open('rb') as stream:files[str(f.relative_to(p))]=hashlib.file_digest(stream,'sha256').hexdigest()
print(json.dumps(files))
'''.replace('REMOTE',repr(remote)).replace('REL',repr(REL))
 manifest=json.loads(subprocess.check_output(['ssh','-o','BatchMode=yes','vast4090','python3','-'],input=code,text=True))
 subprocess.run(['rsync','-a','--partial','--delay-updates','-e','ssh -o BatchMode=yes','vast4090:'+remote+'/'+REL+'/outputs/',str(R/'outputs')+'/'],check=True)
 for name,expected in manifest.items():
  if digest(R/name)!=expected:raise ValueError('Remote/local mismatch: '+name)
 write_new(R/'data/server_output_manifest.json',manifest)
 write_new(R/'data/server_collection_verification.json',dict(files=len(manifest),all_sha256_match=True,remote=remote,deleted_server_files=0))
 print('ALL_SERVER_OUTPUTS_SHA_VERIFIED',len(manifest),flush=True)
