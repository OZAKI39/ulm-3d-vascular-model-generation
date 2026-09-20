#!/usr/bin/env python3
"""Classify the observed compiler errors; never convert an unrun GPU stage to pass."""
import json,shutil,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import sha256,write_json,now
R=ROOT/'reports/sv1_3g';remote=R/'remote'
for name in ('cuda_configure','preconfigure_smoke','gpu_hardware','mpi_runtime_libraries','mpi_native_artifacts','remote_completion'):
    shutil.copyfile(remote/(name+'.json'),R/(name+'.json'))
build=json.loads((remote/'cuda_build.json').read_text());log=ROOT/'logs/sv1_3g/remote/petsc_make.log';txt=log.read_text()
assert json.loads((R/'cuda_configure.json').read_text())['exit_code']==0
assert build['failed_step']=='petsc_make' and build['steps'][-1]['exit_code']!=0
assert all(x in txt for x in ("has no member named 'clockRate'","has no member named 'memoryClockRate'",'namespace "thrust" has no member class "unary_function"'))
inventory=json.loads((R/'cuda_toolkit_inventory.json').read_text());assert not inventory['usable_cuda12_found'] and not inventory['errors']
diagnostic=[{'line':i,'text':line} for i,line in enumerate(txt.splitlines(),1) if 'error:' in line]
build.update(status='BLOCKED',reason='CUDA_TOOLKIT_COMPATIBILITY',failure_layer='CUDA toolkit API / PETSc source compatibility',failure_class='CUDA_API_VERSION_FAIL',diagnostic_evidence={'path':str(log.relative_to(ROOT)),'sha256':sha256(log),'errors':diagnostic},configure='PASS',make='FAIL',install='NOT_RUN',self_test='NOT_RUN',library_sha256=None,usable_cuda12_found=False,toolkit_fallback_attempts=0,source_patch_applied=False,upgraded_PETSc=False)
write_json(R/'cuda_build.json',build)
write_json(R/'failure_classification.json',{'status':'BLOCKED','reason':'CUDA_TOOLKIT_COMPATIBILITY','layer':build['failure_layer'],'configure':'PASS','make':'FAIL','actual_compiler_errors':diagnostic,'MPI_repaired':True,'CUDA_arch_rejected':False,'CUDA_library_rejected':False,'GPU_performance_inference':None,'next_authorized_action':'Stop this stage. Recommend a future side-by-side CUDA 12.x toolkit trial, with driver/PETSc/svMultiPhysics unchanged; installation is not performed in this stage.','official_sources':[{'url':'https://docs.nvidia.com/cuda/archive/13.0.3/pdf/CUDA_Toolkit_Release_Notes.pdf','supports':'CUDA 13 removed clockRate and memoryClockRate fields'},{'url':'https://nvidia.github.io/cccl/unstable/cccl/3.0_migration_guide.html','supports':'CCCL 3 removed thrust::unary_function'},{'url':'https://www.open-mpi.org/software/ompi/v4.1/','supports':'Open MPI 4.1.6 official archive SHA256'},{'url':'https://www.open-mpi.org/projects/hwloc/doc/v2.14.0/plugins.html','supports':'GL discovery may hang connecting to X server; agrees with observed trace, not proof of exact stack frame'}]})
names=['petsc_gpu_smoke','cuda_types','svmp_gpu_build','svmp_gpu_smoke','gpu_proof','cpu_proof','gpu_cpu_equivalence','gpu_residency','gpu_transfer','gpu_memory','benchmark']
for name in names:
    d={'status':'NOT_RUN','executed':False,'reason':'CUDA_TOOLKIT_COMPATIBILITY: PETSc library build failed; prerequisite hard gate not met','measurements':None}
    if name=='cuda_types':d.update(mat_type=None,vec_type=None,KSP=None,PC=None)
    if name=='gpu_proof':d.update(steps=None,linear_failures=None,nonlinear_failures=None,velocity_finite=None,pressure_finite=None,mass_error=None,reload_pass=None)
    if name=='gpu_residency':d.update(measured=False,matrix_resident=None,vectors_resident=None,pc_apply_device=None)
    if name=='gpu_transfer':d.update(measured=False,H2D_count=None,D2H_count=None,large_transfer_per_iteration=None)
    if name=='gpu_memory':d.update(measured=False,peak_bytes=None)
    if name=='benchmark':d.update(CPU_1R=None,CPU_4R=None,GPU_1R=None,speedup=None)
    write_json(R/(name+'.json'),d)
write_json(ROOT/'benchmarks/sv1_3g/status.json',json.loads((R/'benchmark.json').read_text()))
write_json(R/'stage_result.json',{'timestamp':now(),'status':'BLOCKED','reason':'CUDA_TOOLKIT_COMPATIBILITY','MPI':'PASS','PETSc_configure':'PASS','PETSc_build':'FAIL','GPU_runtime':'NOT_RUN','GPU_recommendation':'Retain CPU_EARLY_STOP_PRODUCTION. Future authorization: side-by-side CUDA 12.x toolkit compatibility validation.','production_changed':False,'mesh_convergence_started':False,'GPU_full_production_started':False})
readme=ROOT/'external/gpu_mpi/README.md';readme.parent.mkdir(parents=True,exist_ok=True)
readme.write_text('# Stage SV1.3G isolated MPI\n\nOpen MPI 4.1.6 source archive is mirrored here. The verified native installation is on the RTX 4090 server at `/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3g/external/gpu_mpi`. Its archive is `outputs/sv1_3g/remote_mpi/gpu_mpi_install.tar.gz`; wrapper paths refer to that server. No system MPI was overwritten. See `reports/sv1_3g/mpi_fallback_build.json` and `mpi_resolution.json`.\n')
write_json(R/'native_artifact_mirror.json',{'files':[{'path':p,'sha256':sha256(ROOT/p),'size':(ROOT/p).stat().st_size} for p in ('external/gpu_mpi/openmpi-4.1.6.tar.gz','outputs/sv1_3g/remote_mpi/gpu_mpi_install.tar.gz','outputs/sv1_3g/remote_mpi/mpi_hello_local')],'purpose':'Preserve native MPI installation/source off the ephemeral remote filesystem; not a WSL production solver'})
print('BLOCKED: CUDA_TOOLKIT_COMPATIBILITY (MPI and CUDA configure PASS)')
