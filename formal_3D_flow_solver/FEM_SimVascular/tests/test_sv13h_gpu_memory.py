from sv13h_support import *
def test_actual_solver_vram_peak():
    d=actual('gpu_memory');assert d['measured'] and d['peak_bytes']>0
def test_device_capacity_does_not_substitute_for_solver_peak():
    d=load('gpu_memory')
    assert d['device_capacity_MiB']==24564 and not d['device_capacity_is_solver_peak']
    assert d['peak_bytes'] is None and not d['measured']

