from sv13h_support import *
def test_actual_official_fluid_gpu_smoke():
    d=actual('svmp_gpu_smoke');gpu_smoke_gate(d)
def test_downstream_work_respects_failed_build_gate():
    assert load('stage_result')['reason']=='PETSC_CUDA12_BUILD_FAIL'
    for name in ('petsc_self_test','petsc_cuda_types','petsc_gpu_smoke','svmp_gpu_build',
                 'svmp_cuda_link','svmp_gpu_smoke','gpu_proof','cpu_proof','benchmark'):
        d=load(name)
        assert d['status']=='NOT_RUN' and not d['executed'] and d['measurements'] is None

