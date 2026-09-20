from sv_validation.sv13 import gpu_gate
from test_sv13_gpu_transfer_audit import good
from sv13_support import gpu_evidence
def test_changed_gpu_field_is_rejected():
    d=good();d['scientific_equivalence']=False;assert gpu_gate(d)=='GPU_REJECTED_NUMERICAL_DIFFERENCE'
def test_fixed_window_uses_same_native_checkpoint():
    d=gpu_evidence('fixed_window_equivalence')
    assert d['status']=='PASS' and d['window_steps']==20 and d['same_native_checkpoint']
    assert all(v=='PASS' for v in d['comparisons'].values())
