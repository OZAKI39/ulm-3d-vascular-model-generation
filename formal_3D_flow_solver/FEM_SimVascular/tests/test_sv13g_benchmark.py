from sv13g_support import *
from sv_validation.sv13g import *
def test_actual_repeated_benchmark():
    d=actual('benchmark');assert d['same_host'] and d['same_science']
    gpu_performance(d['CPU_1R'],d['CPU_4R'],d['GPU_1R'],True,True)
def test_one_measurement_cannot_form_performance_conclusion():
    with pytest.raises(GateError,match='REPEATABILITY'):benchmark_time([1.0])
def test_third_run_required_for_timing_variation():
    with pytest.raises(GateError,match='THIRD_REPEAT'):benchmark_time([10,8])
    assert benchmark_time([10,8,9])==9
    assert benchmark_time([10,9.5])==9.5
