from particle_3d.particle82a_admission import method_a


def test_all_retries_retain_original_size(simple_context,common_inputs):
    for event in common_inputs:
        result=method_a(event,simple_context,guard=64)
        assert result['radius_m']==event['first_radius_m']
        assert all(x['diameter_um']==event['first_diameter_um'] for x in result['attempts'])
