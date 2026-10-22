import hashlib,json,os,subprocess,traceback
from datetime import datetime,timezone
from pathlib import Path

ROOT=Path('/workspace/brava_cfd');SOURCE=ROOT/'sources/petsc-3.25.5';ARCH='arch-brava-cuda126'
WRAP=ROOT/'bin/brava-env';steps=[]
def run(name,args,cwd=ROOT,env=None):
    print(name,flush=True)
    with (ROOT/'logs'/f'{name}.log').open('w') as f:
        p=subprocess.run([str(WRAP),*args],cwd=cwd,stdout=f,stderr=subprocess.STDOUT,env=env)
    steps.append(dict(name=name,command=[str(WRAP),*args],cwd=str(cwd),returncode=p.returncode))
    (ROOT/'reports/NATIVE_BUILD_STEPS.json').write_text(json.dumps(steps,indent=2)+'\n')
    if p.returncode:raise RuntimeError(f'{name} failed: {p.returncode}')

result=dict(start_utc=datetime.now(timezone.utc).isoformat(),status='RUNNING',purpose='Environment compatibility rebuild, not a new discretization or CFD campaign')
try:
    assert json.loads((ROOT/'reports/ENVIRONMENT_INSTALL.json').read_text())['status']=='PASS'
    run('mpi_version',[str(ROOT/'deps/openmpi-4.1.6/bin/mpiexec'),'--version'])
    run('mpi_wrapper_show',[str(ROOT/'deps/openmpi-4.1.6/bin/mpicc'),'--showme'])
    run('mpi_launch',[str(ROOT/'bin/brava-mpi'),'-n','2','/bin/hostname'])
    configure=['/usr/bin/python3','-B','configure',f'--prefix={ROOT}/deps/petsc-cuda12',f'PETSC_ARCH={ARCH}',f'--with-cc={ROOT}/deps/openmpi-4.1.6/bin/mpicc',f'--with-cxx={ROOT}/deps/openmpi-4.1.6/bin/mpicxx',f'--with-mpiexec={ROOT}/bin/brava-mpi','--with-fc=0','--with-fortran-bindings=0','--with-debugging=0','--with-shared-libraries=1','--with-scalar-type=real','--with-precision=double','--with-64-bit-indices=0','--with-blaslapack-lib=[/usr/lib/x86_64-linux-gnu/liblapack.so.3,/usr/lib/x86_64-linux-gnu/libblas.so.3]','--with-x=0','--with-cxx-dialect=C++17','--with-hip=0','--with-opencl=0','--with-cuda=1','--with-cuda-dir=/usr/local/cuda','--with-cudac=/usr/local/cuda/bin/nvcc','--with-cuda-arch=89','--with-cuda-dialect=C++17','CUDAFLAGS=-ccbin /usr/bin/g++']
    # The wrapper exports the install prefix; configure/make explicitly select the source.
    run('petsc_configure',['env',f'PETSC_DIR={SOURCE}',f'PETSC_ARCH={ARCH}',*configure],SOURCE)
    run('petsc_make',['make','-j8',f'PETSC_DIR={SOURCE}',f'PETSC_ARCH={ARCH}','all'],SOURCE)
    run('petsc_install',['make',f'PETSC_DIR={SOURCE}',f'PETSC_ARCH={ARCH}','install'],SOURCE)
    run('petsc_check',['make',f'PETSC_DIR={ROOT}/deps/petsc-cuda12','PETSC_ARCH=','check'],SOURCE)
    args=f'-DSV_USE_TETGEN:BOOL=OFF;-DTETGEN_LIBRARY_NAME:FILEPATH={ROOT}/deps/tetgen/lib/lib_simvascular_thirdparty_tetgen.a;-DCMAKE_CXX_FLAGS:STRING=-I{ROOT}/deps/tetgen/include;-DBLAS_LIBRARIES:FILEPATH=/usr/lib/x86_64-linux-gnu/libblas.so.3;-DLAPACK_LIBRARIES:FILEPATH=/usr/lib/x86_64-linux-gnu/liblapack.so.3;-DMPI_C_COMPILER:FILEPATH={ROOT}/deps/openmpi-4.1.6/bin/mpicc;-DMPI_CXX_COMPILER:FILEPATH={ROOT}/deps/openmpi-4.1.6/bin/mpicxx'
    run('svmp_configure',['cmake','-S',str(ROOT/'project/gpu_solver_fix/source'),'-B',str(ROOT/'build/svmp-cuda12'),'-DCMAKE_C_COMPILER=/usr/bin/gcc','-DCMAKE_CXX_COMPILER=/usr/bin/g++','-DCMAKE_BUILD_TYPE=Release','-DCMAKE_CXX_STANDARD=17',f'-DSV_LOCAL_VTK_PATH={ROOT}/deps/vtk-9.3',f'-DSV_PETSC_DIR={ROOT}/deps/petsc-cuda12',f'-DSV_ADDITIONAL_CMAKE_ARGS={args}'])
    run('svmp_make',['cmake','--build',str(ROOT/'build/svmp-cuda12'),'--parallel','8'])
    binary=ROOT/'build/svmp-cuda12/svMultiPhysics-build/bin/svmultiphysics'
    run('svmp_ldd',['ldd',str(binary)])
    assert 'not found' not in (ROOT/'logs/svmp_ldd.log').read_text()
    result.update(status='PASS',binary=str(binary),binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),petsc_library_sha256=hashlib.sha256((ROOT/'deps/petsc-cuda12/lib/libpetsc.so').read_bytes()).hexdigest(),PETSc_version='3.25.5',CUDA_toolkit='12.6',GPU_arch='89',compiler_parallelism=8)
except Exception:
    result.update(status='FAIL',exception=traceback.format_exc());print(result['exception'],flush=True)
finally:
    result['end_utc']=datetime.now(timezone.utc).isoformat()
    (ROOT/'reports/NATIVE_BUILD.json').write_text(json.dumps(result,indent=2)+'\n')
raise SystemExit(result['status']!='PASS')
