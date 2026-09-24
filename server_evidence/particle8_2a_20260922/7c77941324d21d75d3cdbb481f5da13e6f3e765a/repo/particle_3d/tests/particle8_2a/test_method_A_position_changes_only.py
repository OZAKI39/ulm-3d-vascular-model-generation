from copy import deepcopy
from particle_3d.particle82a_admission import method_a


def test_position_retries_keep_id_and_ledger(simple_context,common_inputs):
    original=deepcopy(common_inputs);changed=False
    for event in common_inputs:
        result=method_a(event,simple_context,guard=64)
        assert result['event_id']==event['event_id']
        assert result['attempts'][0]['position_m']==event['anchor_m']
        changed |= len(result['attempts'])>1
    assert changed and common_inputs==original
