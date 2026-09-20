from sv13j_support import *
def test_actual_repeated_benchmarks():
    d=actual('benchmark')
    for k in ('CPU_1R','CPU_4R','GPU_1R'):median_benchmark(d[k])
def test_single_sample_has_no_performance_conclusion():
    with pytest.raises(GateError,match='REPEATABILITY'):median_benchmark([benchmark_sample(1.)])
def test_profiling_not_formal_timing():
    runs=[benchmark_sample(1.),benchmark_sample(1.)];runs[0]['profiling']=True
    with pytest.raises(GateError,match='PROFILING'):median_benchmark(runs)
def test_median_and_required_third_repeat():
    assert median_benchmark([benchmark_sample(10.),benchmark_sample(10.8)])==pytest.approx(10.4)
    with pytest.raises(GateError,match='THIRD_REPEAT'):median_benchmark([benchmark_sample(10.),benchmark_sample(12.)])
    assert median_benchmark([benchmark_sample(10.),benchmark_sample(12.),benchmark_sample(11.)])==11.

