import gzip,json
import numpy as np
def validate(c):
    if c['proof']=='CONTINUOUS_ORIGINAL_SUPPORT_PLANE_MINUS_H_LOWER':
        assert c['minimum_g_nf_bound_m']>=-c['roundoff_m']
    else:
        assert c['proof']=='UNION_OF_ORIGINAL_CONTINUOUS_HANDOFF_CERTIFICATES'
        assert c['proof_partition_only'] and c['held_velocity_path_unchanged']
        endpoint=0.
        for part in c['subcertificates']:
            assert part['t0_fraction']==endpoint and part['t1_fraction']>endpoint
            endpoint=part['t1_fraction'];validate(part['certificate'])
        assert endpoint==1.
def test_no_handoff_violation_new_smoke(report,newpaths,newenv):
    for e,a in newpaths:
        assert np.min(a[:,15])>=-newenv.wall.roundoff_m
        p=report/'outputs/NEW/trajectories'/f"mb_{e['particle_id']:06d}.audit.jsonl.gz"
        with gzip.open(p,'rt') as f:
            for line in f:
                r=json.loads(line)
                if r['accepted']:
                    for c in r['continuous_certificates']:
                        validate(c)
