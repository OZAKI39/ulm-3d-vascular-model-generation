from sv13m_support import *
from sv_validation.sv13m import *
def test_actual_science():science_gate(accepted('science_equivalence'))
@pytest.mark.parametrize('field',list(SCIENCE_LIMITS))
def test_outside_frozen_tolerance_rejected(field):
 d=dict(pressure_shift_applied=False,errors={k:0 for k in SCIENCE_LIMITS});d['errors'][field]=1.01*SCIENCE_LIMITS[field]
 with pytest.raises(GateError):science_gate(d)
def test_pressure_alignment_forbidden():
 with pytest.raises(GateError):science_gate(dict(pressure_shift_applied=True,errors={k:0 for k in SCIENCE_LIMITS}))
