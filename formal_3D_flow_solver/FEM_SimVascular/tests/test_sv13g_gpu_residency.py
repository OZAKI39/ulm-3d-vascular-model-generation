from sv13g_support import *
from sv_validation.sv13g import *
def test_actual_residency():residency_gate(actual('gpu_residency'))
def test_unmeasured_residency_is_not_a_pass():
    with pytest.raises(GateError,match='NOT_MEASURED'):residency_gate({'measured':False})
