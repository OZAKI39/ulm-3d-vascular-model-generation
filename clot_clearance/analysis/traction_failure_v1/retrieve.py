"""Fetch only explicitly inventoried current-run files; no broad project sync."""
from pathlib import Path
import json,sys,subprocess,shlex,hashlib,tarfile,io
C=json.loads(Path('/home/lzy/projects/.traction_failure_current.json').read_text());ROOT=Path(C['root']);COMMON=ROOT/C['common']
SSH=['ssh','-p','15283','-i','/home/lzy/.ssh/vast_step3b_ed25519','-o','IdentitiesOnly=yes','-o','BatchMode=yes','root@95.3.33.46']
mode=sys.argv[1]
if mode=='preview':files=[C['FA'],C['FB']]
else:
 manifest=json.loads(subprocess.check_output(SSH+[shlex.join(['cat','/workspace/y_clot/project/'+C['common']+'/COMPACT_MANIFEST.json'])]))
 files=[d['path'] for d in manifest['files']]+[C['common']+'/'+s for s in ['COMPACT_MANIFEST.json','SERVER_ARTIFACT_MANIFEST.json']]
assert all(not Path(p).is_absolute() and '..' not in Path(p).parts for p in files)
proc=subprocess.Popen(SSH+[shlex.join(['tar','-cf','-','-C','/workspace/y_clot/project','--',*files])],stdout=subprocess.PIPE)
with tarfile.open(fileobj=proc.stdout,mode='r|') as archive:
 for member in archive:
  target=ROOT/member.name
  assert not member.name.startswith('/') and '..' not in Path(member.name).parts and not member.issym() and not member.islnk()
  if member.isdir():target.mkdir(parents=True,exist_ok=True)
  elif member.isfile():
   target.parent.mkdir(parents=True,exist_ok=True)
   with archive.extractfile(member) as src,target.open('wb') as dst:
    import shutil;shutil.copyfileobj(src,dst)
assert proc.wait()==0
if mode!='preview':
 checked=[]
 for d in manifest['files']:
  p=ROOT/d['path'];actual=hashlib.sha256(p.read_bytes()).hexdigest();assert p.stat().st_size==d['bytes'] and actual==d['sha256'],d['path'];checked.append(d['path'])
 result=dict(all_sha256_match=True,files=len(checked),bytes=manifest['total_bytes'],manifest_sha256=hashlib.sha256((COMMON/'COMPACT_MANIFEST.json').read_bytes()).hexdigest(),large_previous_CFD_transferred=False,scope='Current run only; existing literature duplicates verified by SHA. Receipts excluded from self inventory.')
 (COMMON/'LOCAL_RETRIEVAL_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
else:print('downloaded seven figure groups for visual inspection')
