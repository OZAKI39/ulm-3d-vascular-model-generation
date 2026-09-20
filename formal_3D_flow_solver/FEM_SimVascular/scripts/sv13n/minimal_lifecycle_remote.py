"""One tiny CUDA object lifetime test, no vascular solve."""
from runner_remote import *
w=load('petsc_gpu13_build');prefix=Path(w['prefix']);mpi=Path(w['MPI_prefix'])
extra={'LD_LIBRARY_PATH':str(prefix/'lib'),'PETSC_OPTIONS':'-skip_petscrc -use_gpu_aware_mpi 0'}
binary=B/'petsc_lifecycle_probe';source=B/'petsc_lifecycle_probe.c'
c=run([mpi/'bin/mpicc',source,'-I'+str(prefix/'include'),'-L'+str(prefix/'lib'),'-Wl,-rpath,'+str(prefix/'lib'),'-lpetsc','-o',binary],'minimal_lifecycle_compile',cuda=w['candidate_wrapper'],extra_env=extra);assert okay(c)
r=run([w['launcher'],'-n','1',w['candidate_wrapper'],binary,'-vec_type','cuda'],'minimal_lifecycle_run',cuda=w['candidate_wrapper'],extra_env=extra)
passed=okay(r) and 'LIFECYCLE object_destroyed=1 petsc_finalize_calls=1 mpi_finalized_last=1 PASS' in text(r) and not any(s in text(r) for s in ('PETSC ERROR','MPI_ABORT','SEGV'))
write('minimal_lifecycle',dict(status='PASS' if passed else 'FAIL',compile=c,execution=r,source_sha256=digest(source),binary_sha256=digest(binary),scope='Public PETSc CUDA lifecycle contract; actual svMP lifecycle checked separately in official GPU smoke'))
assert passed
