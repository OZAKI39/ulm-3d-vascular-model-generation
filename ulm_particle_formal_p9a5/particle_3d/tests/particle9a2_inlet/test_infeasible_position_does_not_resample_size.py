import pytest
from method_c_test_helpers import source,FixedDistribution
from particle_3d.injection_method_c import PositionSamplingProgressFailure

def test_guard_is_progress_failure_not_size_rejection():
    d=FixedDistribution(1e-6);s=source(d)
    class Reject:
        def check(self,*args):return None,'FORCED_TEST_REJECTION',{}
    s.checker=Reject()
    with pytest.raises(PositionSamplingProgressFailure) as err:s.event(1)
    assert d.calls==1
    assert err.value.record['size_resampled_after_position_failure'] is False
    assert err.value.record['diameter_m']==1e-6
