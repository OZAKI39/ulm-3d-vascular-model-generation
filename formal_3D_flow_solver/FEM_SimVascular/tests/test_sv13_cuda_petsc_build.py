from sv13_support import gpu_evidence
def test_actual_cuda_build_not_library_presence():
    b=gpu_evidence('cuda_petsc_build')
    assert b['status']=='PASS' and b['version']=='3.19.6'
    assert b['cuda_support'] and b['binary_sha256'] and b['configure_options']
