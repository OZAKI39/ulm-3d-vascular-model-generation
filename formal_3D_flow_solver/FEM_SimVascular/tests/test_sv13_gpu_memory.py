import pytest
from sv_validation.sv13 import gpu_gate
from test_sv13_gpu_transfer_audit import good
from sv13_support import gpu_evidence
@pytest.mark.parametrize('update',[{'oom':True},{'peak_memory_fraction':.9},{'peak_memory_fraction':1.}])
def test_oom_and_small_margin_rejected(update):
    d=good();d.update(update);assert gpu_gate(d)=='DEVICE_MEMORY_MARGIN_TOO_SMALL'
def test_actual_gpu_memory():
    d=gpu_evidence('gpu_memory');assert d['peak_bytes']<.9*d['physical_bytes']
