from sv13g_support import *
def test_actual_peak_gpu_memory():
    d=actual('gpu_memory');assert d['measured'] and 0<d['peak_bytes']<=d['capacity_bytes']
