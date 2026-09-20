from sv13h_support import *
def test_actual_profiled_transfer_cost():
    d=actual('gpu_transfer')
    assert d['measured'] and d['H2D_bytes']>=0 and d['D2H_bytes']>=0
def test_transfer_bound_retains_cpu_production():
    c=[run_record(4.),run_record(4.)];g=[run_record(1.),run_record(1.)]
    d=performance_classification(c,c,g,True,'GPU_TRANSFER_BOUND')
    assert d['classification']=='GPU_TRANSFER_BOUND' and d['production']=='CPU_EARLY_STOP_PRODUCTION'

