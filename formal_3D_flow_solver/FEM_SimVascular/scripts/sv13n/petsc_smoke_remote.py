"""One necessary GPU sparse solve under the convergence-driven strategy."""
import json, math, re, shutil, sys
from pathlib import Path
from runner_remote import *
candidate=sys.argv[1];w=load('petsc_'+candidate+'_build');assert w['status']=='PASS'
prefix=Path(w['prefix']);source=Path(w['source']);wrapper=w['candidate_wrapper']
mpi=json.loads((BASE/'configs/baseline_L_mpi_application_gate.json').read_text());mpiprefix=Path(w['MPI_prefix'])
extra={'OMPI_CC':w.get('CC','/usr/bin/gcc-12'),'OMPI_CXX':w.get('CXX','/usr/bin/g++-12'),'PETSC_OPTIONS':'-skip_petscrc -use_gpu_aware_mpi 0','LD_LIBRARY_PATH':str(prefix/'lib')}
def execute(args,name,cwd=source,timeout=120):return run(args,name,cwd=cwd,timeout=timeout,cuda=wrapper,extra_env=extra)
binary=B/'petsc_cuda_smoke'
compile=execute([mpiprefix/'bin/mpicc',B/'petsc_cuda_smoke.c','-I'+str(prefix/'include'),'-L'+str(prefix/'lib'),'-Wl,-rpath,'+str(prefix/'lib'),'-lpetsc','-lm','-o',binary],'petsc_gpu_smoke_compile')
assert okay(compile),'PETSC_GPU_SMOKE_COMPILE_FAIL'
runs=[]
for i in range(1,2):
    r=execute([mpi['working_launcher'],'-n','1',wrapper,binary,'-mat_type','aijcusparse','-vec_type','cuda','-ksp_type','gmres','-ksp_converged_reason','-ksp_monitor_true_residual','-mat_view','::ascii_info','-vec_view','::ascii_info','-log_view','-log_view_gpu_time'],'petsc_gpu_smoke_'+str(i))
    match=re.search(r'SMOKE Mat=(\S+) Vec=(\S+) KSP=(\S+) PC=(\S+) reason=(-?\d+) iterations=(\d+) relative_residual=(\S+) error_inf=(\S+)',text(r))
    if match:
        mt,vt,kt,pt,reason,its,resid,error=match.groups()
        r.update(mat_type=mt,vec_type=vt,KSP=kt,PC=pt,converged_reason=int(reason),iterations=int(its),true_relative_residual=float(resid),solution_error_inf=float(error))
        r['accepted']=okay(r) and 'cusparse' in mt and 'cuda' in vt and kt=='gmres' and int(reason)>0 and math.isfinite(float(resid)) and float(resid)<=1e-10 and math.isfinite(float(error)) and float(error)<=1e-10
    else:r.update(accepted=False,parse_error='Missing runtime SMOKE record')
    runs.append(r)
    write('petsc_gpu_smoke',{'status':'RUNNING','runs':runs,'compile':compile})
passed=len(runs)==1 and all(r['accepted'] for r in runs)
write('petsc_gpu_smoke',{'status':'PASS' if passed else 'FAIL','runs':runs,'compile':compile,'binary':str(binary),'binary_sha256':digest(binary),'source_sha256':digest(B/'petsc_cuda_smoke.c')})
write('petsc_cuda_types',{'status':'PASS' if passed else 'FAIL','mat_type':runs[0].get('mat_type'),'vec_type':runs[0].get('vec_type'),'KSP':runs[0].get('KSP'),'PC':runs[0].get('PC'),'evidence':'One runtime MatGetType/VecGetType/KSPView/PCView/OptionsView solve'})
assert passed,'PETSC_GPU_RUNTIME_FAIL'
print('PETSc GPU sparse smoke 1/1 PASS',flush=True)
