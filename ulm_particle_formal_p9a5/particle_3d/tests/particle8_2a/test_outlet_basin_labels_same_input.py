from copy import deepcopy
from particle_3d.particle82a_admission import method_a,method_b


def test_strategy_does_not_use_basin_to_choose_samples(simple_context,common_inputs):
    event=common_inputs[0]
    for method in [method_a,method_b]:
        results=[method(dict(event,point_tracer_basin=basin),simple_context,guard=16)
                 for basin in ['OUTLET_01','OUTLET_02','OUTLET_03']]
        assert results[0]==results[1]==results[2]
