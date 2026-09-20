from sv13h_support import *
def test_actual_three_sparse_gpu_solves():
    d=actual('petsc_gpu_smoke')
    assert len(d['runs'])==3
    for r in d['runs']:gpu_smoke_gate(r)
@pytest.mark.parametrize('change',[{'converged_reason':-3},{'linear_failures':1},{'true_relative_residual':float('nan')}])
def test_exit_zero_with_bad_ksp_or_residual_is_rejected(change):
    d=smoke_fixture();d.update(change)
    with pytest.raises(GateError):gpu_smoke_gate(d)

