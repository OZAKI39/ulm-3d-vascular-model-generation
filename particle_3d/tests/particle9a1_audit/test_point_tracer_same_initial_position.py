import numpy as np
from particle_3d.routing_stationary_audit import read,verify_events

def test_all_saved_point_initial_positions_exact(audit,candidates,events):
    assert verify_events(candidates,events)['same_initial_position']
    rows=read(audit/'data/point_completed.json')['results']
    assert len(rows)==len(candidates)
    for r,c in zip(rows,candidates):
        a=np.load(audit/'diagnostic_outputs/point'/f"point_{r['particle_id']:06d}.npz")['path']
        assert np.array_equal(a[0,1:],c['common_event']['anchor_m'])
        assert a[0,0]==0 and np.all(np.diff(a[:,0])>0)
