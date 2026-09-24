from copy import deepcopy
import pytest
from particle_3d.routing_stationary_audit import verify_events

def test_original_admission_identity_and_tamper_rejected(candidates,events):
    result=verify_events(candidates,events)
    assert result['accepted_count']==len(events)
    altered=deepcopy(events);altered[0]['birth_center_m'][0]+=1e-9
    with pytest.raises(AssertionError):verify_events(candidates,altered)
