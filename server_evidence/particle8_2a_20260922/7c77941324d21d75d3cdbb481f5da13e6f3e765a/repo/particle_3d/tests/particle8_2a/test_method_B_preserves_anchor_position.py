from particle_3d.particle82a_admission import method_b


def test_every_size_retry_stays_at_anchor(simple_context,common_inputs):
    for event in common_inputs:
        result=method_b(event,simple_context,guard=64)
        assert all(x['position_m']==event['anchor_m'] for x in result['attempts'])
        assert result['anchor_m']==event['anchor_m']
