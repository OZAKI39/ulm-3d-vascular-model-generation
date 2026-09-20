"""Three real sparse GMRES CUDA solves using the pinned standalone C test."""
import json,re
from pathlib import Path
from runner_remote import *
build=load('petsc_cuda12_build');assert build['status']=='PASS'
prefix=Path(build['prefix']);mpi=json.loads((BASE/'configs/mpi_resolution.json').read_text());mpiprefix=Path(mpi['prefix'])
binary=B/'petsc_cuda_smoke'
extra={'LD_LIBRARY_PATH':str(prefix/'lib'),'PETSC_OPTIONS':'-skip_petscrc -use_gpu_aware_mpi 0'}
compile=run([mpiprefix/'bin/mpicc',B/'petsc_cuda_smoke.c','-I'+str(prefix/'include'),'-L'+str(prefix/'lib'),'-Wl,-rpath,'+str(prefix/'lib'),'-lpetsc','-lm','-o',binary],'petsc_gpu_smoke_compile',timeout=120,extra_env=extra)
assert okay(compile),'PETSC_GPU_SMOKE_FAIL: compilation'
runs=[]
for i in range(1,4):
    r=run([mpi['working_launcher'],'-n','1',BASE/'scripts/use_cuda12_gpu_env.sh',binary,'-mat_type','aijcusparse','-vec_type','cuda','-ksp_type','gmres','-ksp_converged_reason','-ksp_monitor_true_residual','-log_view','-log_view_gpu_time'],'petsc_gpu_smoke_'+str(i),timeout=120,extra_env=extra)
    match=re.search(r'SMOKE Mat=(\S+) Vec=(\S+) KSP=(\S+) PC=(\S+) reason=(\d+) iterations=(\d+) relative_residual=(\S+) error_inf=(\S+)',text(r))
    if match:
        mt,vt,kt,pt,reason,its,resid,error=match.groups()
        r.update(mat_type=mt,vec_type=vt,KSP=kt,PC=pt,converged_reason=int(reason),iterations=int(its),true_relative_residual=float(resid),solution_error_inf=float(error),tolerance=1e-10,linear_failures=0 if int(reason)>0 else 1,log_view_present='Event' in text(r) and 'Time (sec)' in text(r))
    else:r['parse_error']='Missing runtime SMOKE record'
    runs.append(r)
    passed=okay(r) and bool(match) and r['mat_type'] in ('seqaijcusparse','mpiaijcusparse','aijcusparse') and r['vec_type'] in ('seqcuda','mpicuda','cuda') and r['converged_reason']>0 and r['true_relative_residual']<=1e-10 and r['solution_error_inf']<=1e-10
    write('petsc_gpu_smoke',{'status':'RUNNING' if passed else 'FAIL','runs':runs,'compile':compile})
    if not passed:raise SystemExit('PETSC_GPU_SMOKE_FAIL')
snapshot=run(['nvidia-smi'],'petsc_gpu_smoke_snapshot')
linked=run(['ldd',binary],'petsc_gpu_smoke_ldd',extra_env=extra)
write('petsc_gpu_smoke',{'status':'PASS','runs':runs,'compile':compile,'binary_sha256':digest(binary),'ldd':text(linked),'GPU_snapshot':snapshot,'source_sha256':digest(B/'petsc_cuda_smoke.c'),'source':'Unchanged SV1.3G standalone smoke C source'})
write('petsc_cuda_types',{'status':'PASS','source':'Actual MatGetType/VecGetType/KSPView/PCView from three native executions','mat_types':[r['mat_type'] for r in runs],'vec_types':[r['vec_type'] for r in runs],'KSP':[r['KSP'] for r in runs],'PC':[r['PC'] for r in runs]})
print('PETSc GPU standalone 3/3 PASS',flush=True)
