import pytest
from sv_validation.sv13 import gpu_gate
from test_sv13_gpu_transfer_audit import good
from sv13_support import artifact
@pytest.mark.parametrize('speed,expected',[(.8,'GPU_SLOWER_THAN_CPU'),(1.,'GPU_SLOWER_THAN_CPU'),(1.24,'GPU_WORKS_BUT_NOT_WORTH_COMPLEXITY'),(1.25,'PASS')])
def test_gpu_adoption_threshold(speed,expected):
    d=good();d['speedup']=speed;assert gpu_gate(d)==expected
def test_actual_cpu_saves_time():
    d=artifact('cpu_performance')
    assert d['same_initial_checkpoint'] and d['measured_continuation_speedup']>1
    assert d['measured_continuation_speedup']==d['reference_measured_wall_s']/d['candidate_measured_wall_s']
