from copy import deepcopy
import numpy as np
from particle_3d.particle82a_admission import method_c


def test_inward_search_preserves_inputs_and_does_not_certify_swept_entry(simple_context,common_inputs):
    event=dict(common_inputs[0],anchor_m=[0.,0.,0.],first_diameter_um=1.,first_radius_m=.5e-6)
    before=deepcopy(event);path=np.array([[0,0,0,0],[.01,0,0,6e-6]])
    result=method_c(event,path,simple_context)
    assert result['accepted'] and event==before
    assert result['anchor_m']==event['anchor_m'] and result['radius_m']==event['first_radius_m']
    assert not result['normal_fallback'] and not result['physical_crossing_certified']
