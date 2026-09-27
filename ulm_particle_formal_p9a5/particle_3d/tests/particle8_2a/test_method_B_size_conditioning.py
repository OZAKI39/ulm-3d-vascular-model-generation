import numpy as np
from particle_3d.particle82a_admission import method_b


def test_accepted_sizes_are_conditioned_and_can_differ(simple_context,common_inputs):
    changed=0
    for event in common_inputs:
        event=dict(event,anchor_m=[1e-6,0,0],first_diameter_um=2.,first_radius_m=1e-6)
        result=method_b(event,simple_context,guard=512)
        if result['accepted']:
            assert result['diameter_um']<=.99600001
            assert result['diameter_um']!=event['first_diameter_um'];changed+=1
        assert 'NOT_ORIGINAL_SONOVUE' in result['size_distribution_role']
    assert changed>0
