import sys
from runner_remote import *
w=load('petsc_'+sys.argv[1]+'_build');prefix=Path(w['prefix']);mpi=Path(w['MPI_prefix']);wrapper=w['candidate_wrapper'];launcher=json.loads((BASE/'configs/baseline_L_mpi_application_gate.json').read_text())['working_launcher']
binary=B/'petsc_cuda_ghost_coherence';extra={'LD_LIBRARY_PATH':str(prefix/'lib'),'PETSC_OPTIONS':'-skip_petscrc -use_gpu_aware_mpi 0'}
r=run([mpi/'bin/mpicc',B/'petsc_cuda_ghost_coherence.c','-I'+str(prefix/'include'),'-L'+str(prefix/'lib'),'-Wl,-rpath,'+str(prefix/'lib'),'-lpetsc','-lm','-o',binary],'coherence_compile',cuda=wrapper,extra_env=extra);assert okay(r)
runs=[]
for backend in ('standard','cuda'):
 for i in range(1,2):
  r=run([launcher,'-n','2',wrapper,binary,'-vec_type',backend],f'coherence_{backend}_{i}',timeout=60,cuda=wrapper,extra_env=extra)
  r['accepted']=okay(r) and 'COHERENCE device_arithmetic=PASS owned_writeback=PASS reverse=PASS duplicate=PASS' in text(r) and 'PETSC ERROR' not in text(r);r['backend']=backend;runs.append(r)
write('ghost_coherence',{'status':'PASS' if all(r['accepted'] for r in runs) else 'FAIL','runs':runs,'source_sha256':digest(B/'petsc_cuda_ghost_coherence.c')})
