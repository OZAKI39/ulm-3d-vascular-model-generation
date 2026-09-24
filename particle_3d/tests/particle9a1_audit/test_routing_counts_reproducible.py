from collections import Counter
import numpy as np
from particle_3d.routing_stationary_audit import read,ROLES

def test_saved_counts_and_transitions_are_reproducible(audit):
    s=read(audit/'data/audit_summary.json');raw=read(audit/'data/point_completed.json')['results'];p65=read(audit/'data/p65_completed.json')['results'];pair=read(audit/'data/paired_routing.json')
    for rows,key,col in [(raw,'raw_candidate_point_split','point_outlet'),([r for r in raw if r['accepted']],'accepted500_point_split','point_outlet'),(p65,'p65_500_split','outlet')]:
        c=Counter(r[col] for r in rows);assert s[key]=={role:c[role] for role in ROLES}
    for a,b,key in [('point_outlet','p65_outlet','point_to_p65_transition'),('p65_outlet','p9a1_outlet','p65_to_p9a1_transition')]:
        m=np.zeros((4,4),int)
        for r in pair:m[ROLES.index(r[a]),ROLES.index(r[b])]+=1
        assert m.tolist()==s[key]['counts'];assert m.sum()==500
