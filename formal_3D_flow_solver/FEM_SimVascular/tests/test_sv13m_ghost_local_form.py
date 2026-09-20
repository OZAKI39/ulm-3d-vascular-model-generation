from sv13m_support import *
from sv_validation.sv13m import *
def test_actual_device_arithmetic_and_local_writeback():
 d=accepted('ghost_coherence');assert len(d['runs'])==6 and all(r['accepted'] for r in d['runs'])
def test_wrong_layout_rejected():
 d=ghost();d['locals'][0][2]='4'
 with pytest.raises(GateError):ghost_run_gate(d)
