import numpy as np
from particle_3d.routing_stationary_audit import read,sha
from particle_3d.particle8_replay import canonical_hash

def test_paired_every_birth_and_nominal_settings(audit,root,events):
    assert len(events)==500
    for e in events:
        name=f"mb_{e['particle_id']:06d}"
        f=audit/'diagnostic_outputs/P65_NEW/trajectories'/(name+'.json')
        g=root/'particle_3d/outputs/particle9a1_2mmps/P9A1/trajectories'/(name+'.json')
        p,q=read(f),read(g)
        assert p['birth_metadata']==q['birth_metadata']==e
        assert p['birth_metadata_sha256']==q['birth_metadata_sha256']==canonical_hash(e)
        assert p['integration_config']==q['integration_config']
        a,b=np.load(f.with_suffix('.npz'))['samples'],np.load(g.with_suffix('.npz'))['samples']
        assert np.array_equal(a[0,1:4],b[0,1:4]) and np.array_equal(a[0,10:14],b[0,10:14])
        assert sha(f.with_suffix('.npz'))==p['samples_sha256']
        assert p['audit_identity']['production_source_sha256']==q['p9a1_identity']['source_sha256']
