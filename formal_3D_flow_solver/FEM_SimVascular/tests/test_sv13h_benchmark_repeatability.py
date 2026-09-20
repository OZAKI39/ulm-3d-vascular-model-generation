from sv13h_support import *
def test_actual_repeated_benchmarks():
    d=actual('benchmark')
    for key in ('CPU_1R','CPU_4R','GPU_1R'):median_benchmark(d[key])
def test_single_benchmark_rejected():
    with pytest.raises(GateError,match='REPEATABILITY'):median_benchmark([run_record(1.)])
def test_profiling_cannot_count_as_benchmark():
    runs=[run_record(1.),run_record(1.)];runs[1]['profiling']=True
    with pytest.raises(GateError,match='PROFILING'):median_benchmark(runs)
def test_median_uses_all_repeats_not_best_or_second():
    assert median_benchmark([run_record(10.),run_record(10.8)])==pytest.approx(10.4)
    assert median_benchmark([run_record(10.),run_record(20.),run_record(15.)])==15.
def test_more_than_ten_percent_requires_third_run():
    with pytest.raises(GateError,match='THIRD_REPEAT'):median_benchmark([run_record(10.),run_record(12.)])

