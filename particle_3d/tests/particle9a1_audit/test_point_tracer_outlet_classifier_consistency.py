import numpy as np
from particle_3d.routing_stationary_audit import read
from particle_3d.particle81_simulation import environment

def test_every_saved_point_exit_reclassifies_on_original_triangles(audit):
    env=environment()
    for r in read(audit/'data/point_completed.json')['results']:
        if r['point_outlet']=='NO_EXIT':continue
        a=np.load(audit/'diagnostic_outputs/point'/f"point_{r['particle_id']:06d}.npz")['path']
        hit=env.classifier.first_event(a[-2,1:],a[-1,1:])
        assert hit is not None and hit.role==r['point_outlet']
