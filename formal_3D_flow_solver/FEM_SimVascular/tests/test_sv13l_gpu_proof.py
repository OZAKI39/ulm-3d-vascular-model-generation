from sv13l_support import *
def test_actual_gpu_twenty_step_proof():
    proof_gate(actual('gpu_proof'))
@pytest.mark.parametrize('key',['velocity_finite','pressure_finite'])
def test_nonfinite_field_rejected(key):
    d=proof_fixture();d[key]=False
    with pytest.raises(GateError,match='NONFINITE'):proof_gate(d)
def test_nan_flux_rejected():
    d=proof_fixture();d['Qin']=float('nan')
    with pytest.raises(GateError,match='NONFINITE'):proof_gate(d)

