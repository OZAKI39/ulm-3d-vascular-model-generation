import copy,json,hashlib
from pathlib import Path
import pytest
from sv_validation.sv13n import *
ROOT=Path(__file__).resolve().parents[1];R=ROOT/'reports/sv1_3n'
def read(name):return json.loads((R/(name+'.json')).read_text())
def accepted(name):
 d=read(name)
 if d.get('status')=='NOT_RUN':pytest.skip(d.get('reason','Blocked by documented preceding gate'))
 assert d['status']=='PASS',d.get('reason',d.get('acceptance_errors',d['status']))
 return d
def ghost(variant='CUDA_MPI'):
 d=read('ghost_probe_new')
 if d.get('status')=='NOT_RUN':pytest.skip(d['reason'])
 return copy.deepcopy(next(r for r in d['runs'] if r['variant']==variant))
def good_flow():
 return dict(exit_code=0,linear_solves=40,linear_failures=0,nonlinear_failures=0,VTU_count=2,velocity_finite=True,pressure_finite=True,reload_pass=True,mat_type='seqaijcusparse',vec_type='seqcuda',KSP='gmres',PC='asm',steps_completed=20,requested_steps=20,initial_state='t=0',wall_noslip_pass=True,flows_finite=True,PETSc_error_detected=False,history={'linear_solves':[{'step':s} for s in range(1,21)]},runtime_semantics=[dict(SOLVER_SEMANTICS,ordering='natural',zero_pivot=2.22045e-14) for _ in range(40)])
def good_timing(t=100):
 return dict(good_flow(),accepted=True,wall_time_s=t,xml_sha256='input',solver_sha256='solver',PETSc_library_sha256='petsc',PETSC_OPTIONS='fixed',mpi_ranks=1,OMP_NUM_THREADS=1,profiling=False)
def synthetic_ghost():
 # Explicit small arithmetic fixture, never used as an experimental result.
 values=[]
 for phase,arrays in [('forward',[[1,2,3,4,5,6],[5,6,7,8,1,2]]),('reverse',[[111,113,3,4,100,101],[105,107,7,8,110,111]])]:
  for rank,array in enumerate(arrays):
   for i,v in enumerate(array):values.append([phase,str(rank),str(i),str(v),str(v)])
 return dict(exit_code=0,timeout=False,hard_error=False,ranks=2,variant='CPU_MPI',types=[['0','mpi','mpi','4','1','2'],['1','mpi','mpi','4','1','0']],locals=[['0','seq','6','6','1'],['1','seq','6','6','1']],checks=[['0','0','0','6','6','0'],['1','0','0','6','6','0']],values=values)
