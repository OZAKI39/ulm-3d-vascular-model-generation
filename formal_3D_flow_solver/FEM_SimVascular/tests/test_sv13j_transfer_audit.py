from sv13j_support import *
def test_actual_measured_transfer():
    d=actual('gpu_transfer');assert d['measured']
def test_unrun_transfer_is_not_zero():
    d=load('gpu_transfer')
    for k in ('H2D_count','D2H_count','H2D_bytes','D2H_bytes','transfer_bound'):assert d[k] is None

