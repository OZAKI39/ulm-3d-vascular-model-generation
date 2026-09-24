from copy import deepcopy
from particle_3d.particle82a_admission import common_event,method_a,method_b


def test_first_proposal_and_size_are_paired(simple_context,common_inputs):
    for event in common_inputs:
        a=method_a(event,simple_context,guard=16);b=method_b(event,simple_context,guard=16)
        assert a['attempts'][0]['position_m']==b['attempts'][0]['position_m']==event['anchor_m']
        assert a['attempts'][0]['diameter_um']==b['attempts'][0]['diameter_um']==event['first_diameter_um']
        assert event==common_event(event['event_id'],ctx=simple_context)
