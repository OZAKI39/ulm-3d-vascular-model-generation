from pathlib import Path
import subprocess,json
S=Path(__file__).resolve().parents[1]
d=json.loads((S/'provenance/REMOTE_INPUTS_TO_PRESERVE.json').read_text())
new='/workspace/microbubble_lammps/results/microbubble_vascular_workflow_v0_'+S.name
work=new.replace('/results/','/work/')
donor='/workspace/microbubble_lammps/results/microbubble_vascular_workflow_v0_adaptive_flux_20260916_222104'
script='import json,hashlib,shutil,socket\nfrom pathlib import Path\n'+f'D={d!r}\nnew=Path({new!r})\nwork=Path({work!r})\ndonor=Path({donor!r})\n'+'''
assert socket.gethostname()==D['hostname']
for n,item in D['inputs'].items():assert hashlib.sha256(Path(n).read_bytes()).hexdigest()==item['sha256'],n
new.mkdir();work.mkdir()
for name in ['geometry','fields','inputs','contracts']:shutil.copytree(donor/name,new/name)
cache=Path('/workspace/microbubble_lammps/work/stable_rigid_sphere_near_field_v0_20260916_101617/build/CMakeCache.txt').read_text()
engine=[l.split('=',1)[1] for l in cache.splitlines() if l.startswith('ENGINE_ROOT:')][0]
p=Path('/workspace/microbubble_lammps/work/stable_rigid_sphere_near_field_v0_20260916_101617/build/rigid_lmp_cpu')
assert hashlib.sha256(p.read_bytes()).hexdigest()=='830d3a1bf41f4944614d6899d8ef74199ff4cc4b496f6037dbe72a2004007789'
print(json.dumps(dict(status='PASS',hostname=socket.gethostname(),protected_files=len(D['inputs']),engine_root=engine,new_results=str(new),new_work=str(work))))
'''
r=subprocess.run(['ssh','vast4090','python3 -'],input=script,text=True,capture_output=True)
(S/'provenance/REMOTE_STAGE_INIT.stderr').write_text(r.stderr);print(r.stdout);assert r.returncode==0,r.stderr
(S/'provenance/REMOTE_STAGE_INIT.json').write_text(r.stdout)
