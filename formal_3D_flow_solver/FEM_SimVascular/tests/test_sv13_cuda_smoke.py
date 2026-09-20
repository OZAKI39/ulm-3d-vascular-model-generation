from sv13_support import gpu_evidence
def test_petsc_and_official_fluid_smoke():
    s=gpu_evidence('cuda_smoke')
    assert s['petsc']=='PASS' and s['official_fluid']=='PASS'
    assert s['real_gpu_kernel'] and s['solution_correct']
