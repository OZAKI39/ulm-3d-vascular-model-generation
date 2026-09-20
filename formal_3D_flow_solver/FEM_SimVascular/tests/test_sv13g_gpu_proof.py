from sv13g_support import *
from sv_validation.sv13g import *
def test_actual_gpu_vascular_proof():proof_gate(actual('gpu_proof'))
def test_four_rank_checkpoint_cannot_restart_one_rank():
    with pytest.raises(GateError,match='INCOMPATIBLE'):restart_gate('native',4,1)
    restart_gate('t=0',None,1)
def test_nonfinite_fields_rejected():
    d=dict(executed=True,steps=20,linear_failures=0,nonlinear_failures=0,velocity_finite=False,pressure_finite=True)
    with pytest.raises(GateError,match='NONFINITE'):proof_gate(d)
