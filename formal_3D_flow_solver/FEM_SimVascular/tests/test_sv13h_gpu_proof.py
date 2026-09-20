from sv13h_support import *
def test_actual_vascular_gpu_twenty_steps():
    proof_gate(actual('gpu_proof'))
@pytest.mark.parametrize('change',[{'velocity_finite':False},{'pressure_finite':False},
    {'mass_error':float('nan')},{'Qin':float('nan')},{'linear_failures':1},{'nonlinear_failures':1}])
def test_nonfinite_or_unconverged_field_rejected(change):
    d=proof_fixture();d.update(change)
    with pytest.raises(GateError):proof_gate(d)

