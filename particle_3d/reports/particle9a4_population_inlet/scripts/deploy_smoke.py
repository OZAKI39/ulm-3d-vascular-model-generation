"""Package only current smoke dependencies into a newly created Vast directory."""
from pathlib import Path
import datetime,hashlib,io,json,shlex,subprocess,tarfile
ROOT=Path(__file__).resolve().parents[4];R=ROOT/'particle_3d/reports/particle9a4_population_inlet'
assert (R/'data/smoke_gate.json').is_file()
stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ');remote='/workspace/particle9a4_population_inlet_'+stamp
archive=Path('/home/lzy/archives')/('p9a4_smoke_'+stamp);archive.mkdir(exist_ok=False)
ref='particle_3d/reports/network_derived_flow_mb_validation_v1';new='particle_3d/reports/particle9a4_population_inlet'
paths=['particle_3d/src','particle_3d/contracts','sonovue_size_distribution_v0',
 'formal_3D_flow_solver/FEM_SimVascular/frozen_reference',
 ref+'/data/old_new_flow_contract.json',ref+'/data/code_snapshot_manifest.json',ref+'/scripts/runner.py',
 ref+'/server_bundle/formal_3D_flow_solver/FEM_SimVascular',ref+'/server_bundle/data/old_new_flow_contract.json',ref+'/server_bundle/inputs/NEW.vtu',
 new+'/data/smoke_gate.json',new+'/data/inlet100k/smoke30_births.json',new+'/scripts/run_smoke30.py']
files={}
for rel in paths:
 p=ROOT/rel
 for q in p.rglob('*') if p.is_dir() else [p]:
  if q.is_file() and not q.is_symlink() and '__pycache__' not in q.parts and q.suffix not in ['.pyc','.pyo']:
   files[str(q.relative_to(ROOT))]=q
manifest={name:hashlib.sha256(p.read_bytes()).hexdigest() for name,p in sorted(files.items())}
package=archive/'deployment.tar.gz'
with tarfile.open(package,'x:gz') as tar:
 for name,p in sorted(files.items()):tar.add(p,arcname=name,recursive=False)
 data=(json.dumps(manifest,indent=2)+'\n').encode();info=tarfile.TarInfo('deployment_manifest.json');info.size=len(data);tar.addfile(info,io.BytesIO(data))
command='test ! -e '+shlex.quote(remote)+' && mkdir '+shlex.quote(remote)
subprocess.run(['ssh','-o','BatchMode=yes','vast4090',command],check=True)
subprocess.run(['scp','-q',str(package),'vast4090:'+remote+'/deployment.tar.gz'],check=True)
subprocess.run(['ssh','-o','BatchMode=yes','vast4090','tar -xzf '+shlex.quote(remote+'/deployment.tar.gz')+' -C '+shlex.quote(remote)],check=True)
record=dict(server_root=remote,local_archive=str(archive),files=len(files),archive_sha256=hashlib.sha256(package.read_bytes()).hexdigest(),archive_bytes=package.stat().st_size,manifest=manifest)
with (R/'data/server_deployment.json').open('x') as f:json.dump(record,f,indent=2)
print(json.dumps({k:v for k,v in record.items() if k!='manifest'}),flush=True)
