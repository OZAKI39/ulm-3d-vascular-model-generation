import numpy as np
from particle_3d.routing_stationary_audit import read,downstream_feasibility
from particle_3d.contact_redundancy import independent_contact_rows

def test_contact_audit_real_and_known_cones(audit):
    for r in read(audit/'data/stationary_audit.json'):
        assert r['contact_rank']==r['contact_redundancy']['rank_after']
        assert r['retained_rows_independent']
        d=r['downstream'];assert d['status']=='RESOLVED'
        if not d['feasible_downstream_direction']:
            N=np.array(r['contact_normals']);u=np.array(r['free_fem_velocity_m_s']);u/=np.linalg.norm(u)
            assert np.linalg.norm(N.T@d['certificate_multipliers']+u)<=d['verification_tolerance']
    N=np.array([[1.,0,0],[1.,0,0],[0,1,0],[0,0,1]])
    kept,d=independent_contact_rows(N,list(range(4)));assert len(kept)==3 and d['rank_after']==3
    assert downstream_feasibility(N,[1,1,1])['feasible_downstream_direction']
    assert not downstream_feasibility(N,[-1,-1,-1])['feasible_downstream_direction']
