from sv13m_support import *
from sv_validation.sv13m import *
def test_all_forward_data():
 ghost_gate(read('ghost_probe_after'))
def test_silent_corruption_rejected():
 d=ghost();d['values'][0][-2]='999'
 with pytest.raises(GateError):ghost_run_gate(d)
def test_falsified_expected_values_rejected():
 d=ghost();d['values'][0][-2:]=['999','999']
 with pytest.raises(GateError):ghost_run_gate(d)
