from sv13l_support import *
def test_actual_official_gpu_smoke_acceptance():official_flow_gate(actual('svmp_gpu_smoke'))
def test_ksp_convergence_does_not_hide_prior_ghost_error():
    d=load('svmp_gpu_smoke');d['exit_code']=0
    assert d['KSP_reasons'][0]['reason']=='CONVERGED_RTOL'
    with pytest.raises(GateError,match='PETSC_ERROR'):official_flow_gate(d)
