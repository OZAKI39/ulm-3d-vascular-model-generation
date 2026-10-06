"""Read-only final CUDA runtime identity check, without launching GPU computation."""
from common import *
import ctypes

def main():
    old=read(BASE/'production_workflow/GPU_MANUFACTURED_RUN_IDENTITY.json')['loaded_GPU_libraries']
    libraries={p:sha(p) for p in old['loaded_libraries']}
    assert libraries==old['loaded_libraries']
    runtime=next(p for p in libraries if 'libcudart.so' in p)
    lib=ctypes.CDLL(runtime);version=ctypes.c_int()
    lib.cudaRuntimeGetVersion.argtypes=[ctypes.POINTER(ctypes.c_int)]
    lib.cudaRuntimeGetVersion.restype=ctypes.c_int
    code=lib.cudaRuntimeGetVersion(ctypes.byref(version));assert code==0
    assert version.value==old['CUDA_runtime_version']
    put('CUDA_LIBRARY_IDENTITY.json',dict(status='CONFIRMED',current_libraries_sha256=libraries,
        current_runtime_version=version.value,version_query_exit_code=code,
        baseline='production_workflow/GPU_MANUFACTURED_RUN_IDENTITY.json',GPU_computation_launched=False,
        method='Rehash the actual retained production CUDA driver/runtime/cuSPARSE/PETSc paths and query runtime version; no benchmark or CFD.'))
    print('Current CUDA/PETSc identities unchanged; CUDA runtime',version.value)

if __name__=='__main__':main()
