import copy,json
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1];R=ROOT/'reports/sv1_3m'
def read(name):return json.loads((R/(name+'.json')).read_text())
def accepted(name):
 d=read(name)
 if d.get('status')=='NOT_RUN':pytest.skip(d.get('reason','Upstream execution gate did not pass'))
 assert d['status']=='PASS',d.get('reason',d.get('acceptance_errors',d['status']))
 return d
def ghost(variant='CUDA_MPI'):
 return copy.deepcopy(next(r for r in read('ghost_probe_after')['runs'] if r['variant']==variant))
def good_flow():
 return dict(exit_code=0,linear_solves=3,linear_failures=0,nonlinear_failures=0,VTU_count=1,velocity_finite=True,pressure_finite=True,reload_pass=True,mat_type='seqaijcusparse',vec_type='seqcuda',KSP='gmres',PC='asm',steps_completed=20,initial_state='t=0',wall_noslip_pass=True,flows_finite=True,PETSc_error_detected=False)
