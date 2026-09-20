"""Retest the preserved MPI stack under CUDA12, then run native sm_89 kernels."""
import json,os,subprocess
from pathlib import Path
from runner_remote import *
assert load('cuda12_install')['status']=='PASS' and load('cuda13_preservation')['status']=='PASS'
prefix=Path(load('cuda12_install')['prefix'])
prior=json.loads((BASE/'configs/mpi_resolution.json').read_text())
mpiprefix=Path(prior['prefix']);mpi=prior['working_launcher']
assert digest(mpi)==prior['wrapper_sha256']
artifacts=json.loads((BASE/'configs/mpi_native_artifacts.json').read_text())
for a in artifacts:assert digest(a['realpath'])==a['sha256']
hello=next(a['realpath'] for a in artifacts if a['path'].endswith('mpi_hello_local'))
wrapper=BASE/'scripts/use_cuda12_gpu_env.sh';wrapper.parent.mkdir(exist_ok=True)
wrapper.write_text('#!/usr/bin/env bash\nset -euo pipefail\n# Stage SV1.3H: scoped CUDA12 environment; no permanent shell or system edits.\nexport CUDA_HOME="'+str(prefix)+'"\nexport PATH="'+str(prefix/'bin')+':'+str(mpiprefix/'bin')+':/usr/bin:/bin${PATH:+:$PATH}"\nexport LD_LIBRARY_PATH="'+str(prefix/'lib64')+':'+str(mpiprefix/'lib')+'${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"\nif [[ "$#" -gt 0 ]]; then exec "$@"; fi\n')
wrapper.chmod(0o755)
repeats=[]
for i in range(1,6):
    r=run([mpi,'-n','1',wrapper,hello],'cuda12_mpi_rank1_'+str(i),timeout=10);r['stdout']=text(r);repeats.append(r)
rank2=run([mpi,'-n','2',wrapper,hello],'cuda12_mpi_rank2',timeout=10);rank2['stdout']=text(rank2)
passed=all(okay(r) and r['stdout'].strip()=='rank=0 size=1' for r in repeats) and okay(rank2) and set(rank2['stdout'].splitlines())=={'rank=0 size=2','rank=1 size=2'}
write('mpi_with_cuda12',{'status':'PASS' if passed else 'FAIL','rank1':repeats,'rank2':rank2,'MPI_rebuilt':False,'MPI_wrapper_modified':False,'wrapper_sha256':digest(wrapper),'MPI_prefix':str(mpiprefix)})
assert passed,'MPI_CUDA12_ENVIRONMENT_FAIL'
nvcc=run([prefix/'bin/nvcc','--version'],'cuda12_nvcc_version');assert okay(nvcc) and 'release 12.6' in text(nvcc)
old=Path(load('pre_install_environment')['cuda13_prefix'])
oldver=run([old/'bin/nvcc','--version'],'cuda13_nvcc_after_install');assert okay(oldver) and 'release 13.2' in text(oldver)
binary=B/'cuda_kernel_smoke'
compile=run([prefix/'bin/nvcc','-std=c++17','-arch=sm_89','--cudart','shared',B/'cuda_kernel_smoke.cu','-o',binary],'cuda12_kernel_compile',timeout=120)
assert okay(compile),'CUDA12_RUNTIME_FAIL: compilation'
runs=[]
for i in range(1,4):
    r=run([binary],'cuda12_kernel_run_'+str(i));r['stdout']=text(r);r['correct']='correct=1 last_error=0' in r['stdout'];runs.append(r)
gpu=run(['nvidia-smi'],'cuda12_gpu_snapshot')
ldd=run(['ldd',binary],'cuda12_kernel_ldd');assert okay(ldd)
assert 'not found' not in text(ldd) and str(prefix) in text(ldd)
runtime={'status':'PASS' if all(okay(r) and r['correct'] for r in runs) else 'FAIL','selected_release':json.loads((BASE/'configs/cuda_selection.json').read_text())['version'],'prefix':str(prefix),'nvcc_path':str(prefix/'bin/nvcc'),'nvcc_version':text(nvcc),'libcudart_path':str((prefix/'lib64/libcudart.so').resolve()),'cusparse_path':str((prefix/'lib64/libcusparse.so').resolve()),'curand_path':str((prefix/'lib64/libcurand.so').resolve()),'include_path':str((prefix/'include').resolve()),'wrapper_enabled':True,'wrapper_sha256':digest(wrapper),'compile':compile,'runs':runs,'binary':str(binary),'binary_sha256':digest(binary),'ldd':text(ldd),'GPU_snapshot':gpu,'cuda_version_manifest':json.loads((prefix/'version.json').read_text())}
write('cuda12_runtime',runtime);assert runtime['status']=='PASS','CUDA12_RUNTIME_FAIL'
