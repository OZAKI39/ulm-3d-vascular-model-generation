"""Close successful standalone provenance without executing another solve."""
from runner_remote import *
h=load('petsc_hypre_build');s=load('hypre_ilu_standalone');assert h['status']==s['status']=='PASS'
write('petsc_hypre_build_before_selftest_close',h)
h.update(selftest='PASS',selftest_artifact='hypre_ilu_standalone.json',selftest_scope='Only one MPI rank and CUDA sparse solve; no new CPU or 2-rank test')
write('petsc_hypre_build',h)
for name in ['hypre_ilu_sparse_profile.txt','hypre_boomeramg_sparse_profile.txt']:
 p=B/name
 if p.exists():
  import shutil
  shutil.copyfile(p,R/name)
print('Standalone provenance closed; no extra solver run',flush=True)
