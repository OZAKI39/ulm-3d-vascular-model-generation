import numpy as np,hashlib
from particle_3d.routing_stationary_audit import read
from particle_3d.stationary_radius_audit import DIAGNOSTIC_ROLE

def test_all_virtual_trials_are_bounded_and_keep_checkpoint(audit,root):
    for r in read(audit/'data/stationary_audit.json'):
        assert r['role']==DIAGNOSTIC_ROLE
        pid=r['particle_id'];s=np.load(root/'particle_3d/outputs/particle9a1_2mmps/P9A1/trajectories'/f'mb_{pid:06d}.npz')['samples'];k=r['checkpoint']['index']
        prefix=np.load(audit/'diagnostic_outputs/virtual_radius'/f'prefix_{pid:06d}.npz')['samples']
        assert np.array_equal(prefix,s[:k+1]);assert hashlib.sha256(prefix.tobytes()).hexdigest()==r['checkpoint']['prefix_sha256']
        passed=False
        for t,expected in zip(r['trials'],[1.,.99,.95,.90]):
            assert not passed and t['radius_ratio']==expected
            assert t['role']==DIAGNOSTIC_ROLE and t['production_result'] is False and t['criteria']['local_only']
            a=np.load(audit/'diagnostic_outputs/virtual_radius'/f'mb_{pid:06d}_ratio_{round(expected*100):03d}.npz')['samples']
            assert np.array_equal(a[0,1:4],s[k,1:4]);assert np.array_equal(a[0,10:14],s[k,10:14])
            assert a[0,0]==s[k,0] and a[-1,0]-a[0,0]<=t['criteria']['local_time_budget_s']+1e-15
            passed=t['status']=='PASSED_LOCAL_HOTSPOT'
