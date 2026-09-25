from pathlib import Path
import os,subprocess,time,json
W=Path('/workspace/microbubble_lammps/work/fixed_multiblob_wall_audit_20260916_143602');R=Path('/workspace/microbubble_lammps/results/fixed_multiblob_wall_audit_20260916_143602');name='patch_0p125'
start=time.monotonic()
env=dict(os.environ,AUDIT_ROOT=str(W/'bundle'),AUDIT_RESULT=str(R),PECNUT_CODE=str(W/'bundle/stokesian_dynamics'),OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4',NUMBA_NUM_THREADS='4')
try:
 p=subprocess.run([str(W/'env/bin/python'),'-u',str(W/'bundle/scripts/run_remote_audit.py'),'0.125','patch'],env=env,timeout=1984.8710042780149)
 receipt={'exit_code':p.returncode,'wall_seconds':time.monotonic()-start,'state':'COMPLETED' if p.returncode==0 else 'FAILED'}
except subprocess.TimeoutExpired:
 receipt={'exit_code':124,'wall_seconds':time.monotonic()-start,'state':'RESOURCE_TIME_BOUND'}
(R/(name+'_RECEIPT.json')).write_text(json.dumps(receipt,indent=2))
