from sv13h_support import *
def test_actual_measured_residency():
    d=actual('gpu_residency');assert d['measured'] and d['pc_apply_device'] in ('CPU','GPU','MIXED')
def test_unmeasured_residency_prevents_promising_classification():
    runs=[run_record(1.),run_record(1.)]
    d=performance_classification([run_record(2.),run_record(2.)],runs,runs,True,'NOT_MEASURED')
    assert d['classification']=='BLOCKED'

