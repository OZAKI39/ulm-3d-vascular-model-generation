"""One new server deployment; never reuse or overwrite any historical directory."""
from pathlib import Path
import json,tarfile,io,subprocess,shlex,sys
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT/'particle_3d/src'))
from particle_3d.formal_cohort_p9a5 import digest,write_new,REL
R=ROOT/REL;c=json.loads((R/'data/run_context.json').read_text());remote=c['remote'];A=Path(c['archive'])
old='particle_3d/reports/network_derived_flow_mb_validation_v1';p9='particle_3d/reports/particle9a4_population_inlet'
paths=['particle_3d/src','particle_3d/contracts','sonovue_size_distribution_v0',
 'formal_3D_flow_solver/FEM_SimVascular/frozen_reference',
 old+'/data/old_new_flow_contract.json',old+'/data/code_snapshot_manifest.json',old+'/scripts/runner.py',
 old+'/server_bundle/formal_3D_flow_solver/FEM_SimVascular',old+'/server_bundle/data/old_new_flow_contract.json',old+'/server_bundle/inputs/NEW.vtu',
 p9+'/data/inlet100k/accepted_births.json.gz',REL+'/data/CORE500_COHORT.json',REL+'/data/BENCHMARK24_COHORT.json',
 REL+'/data/historical_horizon_audit.json',REL+'/data/input_audit.json',REL+'/data/run_context.json',REL+'/scripts/run_server.py']
paths += [p9+f'/outputs/smoke30/trajectories/mb_{i:06d}.npz' for i in range(1,25)]
files={}
for rel in paths:
 p=ROOT/rel
 for q in p.rglob('*') if p.is_dir() else [p]:
  if q.is_file() and not q.is_symlink() and '__pycache__' not in q.parts and q.suffix not in ['.pyc','.pyo']:
   files[str(q.relative_to(ROOT))]=q
manifest={name:digest(p) for name,p in sorted(files.items())};package=A/'deployment.tar.gz'
with tarfile.open(package,'x:gz') as tar:
 for name,p in sorted(files.items()):tar.add(p,arcname=name,recursive=False)
 data=(json.dumps(manifest,indent=2)+'\n').encode();info=tarfile.TarInfo('deployment_manifest.json');info.size=len(data);tar.addfile(info,io.BytesIO(data))
subprocess.run(['scp','-q',str(package),'vast4090:'+remote+'/deployment.tar.gz'],check=True)
subprocess.run(['ssh','-o','BatchMode=yes','vast4090','test ! -e '+shlex.quote(remote+'/deployment_manifest.json')+' && tar -xzf '+shlex.quote(remote+'/deployment.tar.gz')+' -C '+shlex.quote(remote)],check=True)
write_new(R/'data/deployment.json',dict(remote=remote,files=len(files),archive_sha256=digest(package),archive_bytes=package.stat().st_size,manifest=manifest))
launch='''import subprocess,os,json
from pathlib import Path
p=Path(REMOTE);env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
f=(p/'logs/benchmark.txt').open('xb')
proc=subprocess.Popen(['/root/particle8_2_runs/env/bin/python',str(p/REL/'scripts/run_server.py'),'--stage','benchmark'],cwd=p,env=env,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
(p/'logs/benchmark_pid.json').write_text(json.dumps(dict(pid=proc.pid,command='P9A5 fixed24 worker benchmark'))+'\\n');print(proc.pid)
'''.replace('REMOTE',repr(remote)).replace('REL',repr(REL))
subprocess.run(['ssh','-o','BatchMode=yes','vast4090','python3','-'],input=launch,text=True,check=True)
print('DEPLOYED_ONCE_BENCHMARK_STARTED',remote,flush=True)
