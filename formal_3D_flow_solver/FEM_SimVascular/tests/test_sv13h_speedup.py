from sv13h_support import *
def test_actual_gpu_speedup():
    d=actual('speedup');assert d['speedup']>0
@pytest.mark.parametrize('gpu_time,label',[(2.,'GPU_A_PROMISING'),(3.5,'GPU_WORKS_BUT_NOT_WORTH_COMPLEXITY'),(5.,'GPU_SLOWER_THAN_CPU')])
def test_classification_uses_faster_cpu_baseline(gpu_time,label):
    c1=[run_record(6.),run_record(6.)];c4=[run_record(4.),run_record(4.)]
    d=performance_classification(c1,c4,[run_record(gpu_time),run_record(gpu_time)],True,'GPU_RESIDENCY_PASS')
    assert d['speedup']==pytest.approx(4/gpu_time) and d['classification']==label
def test_fast_but_incorrect_gpu_rejected():
    r=[run_record(1.),run_record(1.)]
    with pytest.raises(GateError,match='GPU_SCIENCE_FAIL'):performance_classification(r,r,r,False,'GPU_RESIDENCY_PASS')

