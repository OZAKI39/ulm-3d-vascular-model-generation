from sv13m_support import *
from sv_validation.sv13m import *
def test_real_twenty_steps():flow_gate(accepted('GPU_PROOF_20_acceptance'),steps=20)
def test_nan_fields_rejected():
 d=good_flow();d['velocity_finite']=False
 with pytest.raises(GateError):flow_gate(d,steps=20)
