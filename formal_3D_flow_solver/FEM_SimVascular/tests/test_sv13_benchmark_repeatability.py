import pytest
from sv_validation.sv13 import repeated_timing
from sv_validation.validation import ValidationError
from sv13_support import gpu_evidence
def test_second_run_used_not_best():
    assert repeated_timing([100.,105.])==105.
def test_third_required_and_median_used():
    with pytest.raises(ValidationError):repeated_timing([100.,150.])
    assert repeated_timing([100.,150.,120.])==120.
def test_fixed_window_real_repeated_measurements():
    d=gpu_evidence('benchmark_repeatability')
    for candidate in ('CPU_1R','CPU_4R','GPU_1R'):
        run=d[candidate];assert repeated_timing(run['wall_seconds'])==run['reported_wall_seconds']
        assert not run['gpu_timing_profiler_enabled']
