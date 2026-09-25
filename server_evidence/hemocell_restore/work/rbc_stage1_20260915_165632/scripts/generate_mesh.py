from pathlib import Path
import subprocess,os,time,json,hashlib
R=Path(__file__).resolve().parents[1]
with (R/'geometry/GENERATION_STARTED.json').open('x') as f:json.dump(dict(unix=time.time(),solver_steps=0,MPI1_SANITY_performed=False),f)
env=dict(os.environ,HWLOC_COMPONENTS='-opencl',OMP_NUM_THREADS='1')
cmd=[str(R/'build/rbc_geometry_probe'),str(R/'contracts/NEW_MEDIUM_NUMERICS_CONTRACT.json'),str(R/'geometry')]
t=time.monotonic()
with (R/'logs/geometry_generation.log').open('w') as f:p=subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=600)
(R/'geometry/GENERATION_TERMINAL.json').write_text(json.dumps(dict(status='PASS' if p.returncode==0 else 'FAIL',returncode=p.returncode,seconds=time.monotonic()-t,command=cmd,solver_steps=0),indent=2))
raise SystemExit(p.returncode)
