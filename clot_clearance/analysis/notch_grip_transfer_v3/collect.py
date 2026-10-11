"""Only this task's review evidence; preserve heavy originals on server."""
import json,subprocess,hashlib,tarfile,shlex,time
from pathlib import Path
C=json.loads(Path('/home/lzy/projects/.notch_grip_v3_current.json').read_text());root=Path(C['root'])
ssh=['ssh','-p','15283','-i','/home/lzy/.ssh/vast_step3b_ed25519','-o','IdentitiesOnly=yes','-o','BatchMode=yes','root@95.3.33.46']
prefixes=[C['relative'],'results/lammps_notch_fracture_grip_transfer_v3/'+C['run_id']]
script="""import json,hashlib
from pathlib import Path
root=Path('/workspace/y_clot/project');files={}
for prefix in PREFIXES:
 for p in (root/prefix).rglob('*'):
  if not p.is_file() or p.name in ['body.data','history.csv','native_energy_density.dat','boundary_dynamics_audit.csv','COMPACT_DOWNLOAD.json','LOCAL_FINAL_VALIDATION.json'] or '__pycache__' in p.parts or (p.parent.name=='M5C' and p.name.startswith('frame_')):continue
  files[str(p.relative_to(root))]={'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
print(json.dumps(files))
""".replace('PREFIXES',repr(prefixes))
files=json.loads(subprocess.check_output(ssh+['python3 -c '+shlex.quote(script)],text=True));need=[p for p,v in files.items() if not (root/p).exists() or hashlib.sha256((root/p).read_bytes()).hexdigest()!=v['sha256']]
if need:
 with subprocess.Popen(ssh+[shlex.join(['tar','-cf','-','-C','/workspace/y_clot/project','--',*need])],stdout=subprocess.PIPE) as proc:
  with tarfile.open(fileobj=proc.stdout,mode='r|') as tf:
   for m in tf:
    assert m.name in need and not m.name.startswith('/') and '..' not in Path(m.name).parts
    tf.extract(m,root,filter='data')
  assert proc.wait()==0
assert all(hashlib.sha256((root/p).read_bytes()).hexdigest()==v['sha256'] for p,v in files.items())
receipt=dict(server_root='/workspace/y_clot/project',local_root=str(root),files=files,total_bytes=sum(v['bytes'] for v in files.values()),newly_transferred=len(need),excluded_names=['body.data','history.csv','native_energy_density.dat','boundary_dynamics_audit.csv','M5C/frame_*'],scope='Current V3 only; review_curves retain exact scalar values and selected accepted states retain full precision. Heavy original files and interval native states stay on server.',retrieval_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
(root/C['relative']/'provenance/COMPACT_DOWNLOAD.json').write_text(json.dumps(receipt,indent=2));print(json.dumps({k:v for k,v in receipt.items() if k!='files'},indent=2))
