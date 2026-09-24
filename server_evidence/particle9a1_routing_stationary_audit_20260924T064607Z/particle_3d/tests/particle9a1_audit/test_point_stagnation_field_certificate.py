import numpy as np
from particle_3d.routing_stationary_audit import read

def test_saved_field_certificate_is_an_attracting_wall_face(audit):
    evidence=read(audit/'data/point_stagnation_field.json')
    assert evidence['all_confirmed'] and evidence['count']>0
    for r in evidence['records']:
        v=np.array(r['nodal_velocities_m_s']);G=np.array(r['velocity_gradient_s_inv']);b=np.array(r['limit_barycentric']);j=r['interior_node_local_index']
        assert np.sum(np.all(v==0,axis=1))==3
        assert np.trace(G)<0 and np.isclose(np.trace(G),r['divergence_s_inv'])
        assert abs(b[j])<=4*r['barycentric_tolerance'] and b.min()>=-4*r['barycentric_tolerance']
        assert r['wall_triangle_id'] is not None and r['gradient_rank_one_relative_error']<=r['rank_one_verification_budget']
        assert np.linalg.norm(b@v)<=4*r['barycentric_tolerance']*np.linalg.norm(v)
