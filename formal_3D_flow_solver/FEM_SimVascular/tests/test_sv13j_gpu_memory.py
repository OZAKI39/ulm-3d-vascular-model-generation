from sv13j_support import *
def test_actual_vascular_memory_peak():
    d=actual('gpu_memory');assert d['measured'] and d['peak_bytes']>0
def test_capacity_is_not_a_measured_peak():
    d=load('gpu_memory');assert d['device_capacity_MiB']==24564 and d['peak_bytes'] is None

