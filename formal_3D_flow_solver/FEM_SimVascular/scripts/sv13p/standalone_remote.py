"""One GPU sparse solve with live PETSc help before choosing the CFD profile."""
import sys,shlex,re
from runner_remote import *
kind=sys.argv[1];pc=sys.argv[2];assert kind in ('hypre','amgx');assert pc in ('ilu','boomeramg','amgx')
w=load('petsc_'+kind+'_build');prefix=Path(w['prefix']);extra={'PKG_CONFIG_PATH':str(prefix/'lib/pkgconfig'),'LD_LIBRARY_PATH':str(prefix/'lib')}
binary=BASE/('external/gpu_sparse_'+kind)
if not binary.exists():
 d=run(['pkg-config','--cflags','--libs','PETSc'],kind+'_pkgconfig',cuda=w['candidate_wrapper'],extra_env=extra);assert okay(d)
 d=run([Path(w['MPI_prefix'])/'bin/mpicc',BASE/'gpu_sparse.c','-o',binary,*shlex.split(text(d))],kind+'_sparse_compile',cuda=w['candidate_wrapper'],extra_env=extra);assert okay(d)
options=['-skip_petscrc','-use_gpu_aware_mpi','0','-mat_type','aijcusparse','-vec_type','cuda','-ksp_type','gmres','-ksp_pc_side','right','-ksp_rtol','1e-10','-ksp_atol','1e-24','-ksp_max_it','2000','-ksp_gmres_restart','100','-ksp_converged_reason','-ksp_view','-help','-log_view',':'+str(BASE/'benchmarks'/(kind+'_'+pc+'_sparse_profile.txt')),'-log_view_gpu_time','-pc_type','amgx' if pc=='amgx' else 'hypre']
if pc!='amgx':options+=['-pc_hypre_type',pc]
d=run([w['launcher'],'-n','1',w['candidate_wrapper'],binary,*options],kind+'_'+pc+'_standalone_help',timeout=120,cuda=w['candidate_wrapper'],extra_env=extra)
raw=text(d);m=re.search(r'SV13P_GPU_SPARSE reason=(\d+) error_inf=(\S+)',raw)
valid=okay(d) and m and 'type: seqaijcusparse' in raw and 'type: seqcuda' in raw and int(m[1])>0 and float(m[2])<1e-7
write(kind+'_'+pc+'_standalone',dict(status='PASS' if valid else 'FAIL',execution=d,options=options,binary_sha256=digest(binary),reason=int(m[1]) if m else None,error_inf=float(m[2]) if m else None,GPU_device_build=w.get('hypre_features'),CPU_only_build=False if valid else None))
raise SystemExit(0 if valid else 1)
