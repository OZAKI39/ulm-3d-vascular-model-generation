from sv13m_support import *
from sv_validation.sv13m import *
def test_official_complete_flow():flow_gate(accepted('svmp_gpu_smoke'))
def test_missing_vtu_rejected():
 d=good_flow();d['VTU_count']=0
 with pytest.raises(GateError):flow_gate(d)
