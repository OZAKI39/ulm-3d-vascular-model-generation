from runner_remote import *
w=load('compatibility_winner');p=Path(w['prefix']);mpi=Path(w['MPI_prefix']);wrapper=w['candidate_wrapper'];launcher=json.loads((BASE/'configs/mpi_resolution.json').read_text())['working_launcher'];extra={'LD_LIBRARY_PATH':str(p/'lib'),'PETSC_OPTIONS':'-skip_petscrc -use_gpu_aware_mpi 0'}
binary=B/'petsc_finalize_probe';c=run([mpi/'bin/mpicc',B/'petsc_finalize_probe.c','-I'+str(p/'include'),'-L'+str(p/'lib'),'-Wl,-rpath,'+str(p/'lib'),'-lpetsc','-lm','-o',binary],'finalize_probe_compile',cuda=wrapper,extra_env=extra);assert okay(c)
runs=[]
for backend,proper in [('standard',False),('cuda',False),('cuda',True)]:
 r=run([launcher,'-n','1',wrapper,binary,'-vec_type',backend,'-proper_finalize',str(proper).lower()],f'finalize_probe_{backend}_{proper}',timeout=60,cuda=wrapper,extra_env=extra);r.update(backend=backend,proper_finalize=proper,diagnostic_only=True);runs.append(r)
write('finalize_probe',{'status':'EXECUTED','runs':runs,'source_sha256':digest(B/'petsc_finalize_probe.c'),'solver_source_changed':False,'compatibility_patch_applied':False})
