"""Reproduce pinned native higher-order cap normals and Flat/rim-zero shape.

The source updates edge normals inside the face loop. Preserving that order is
necessary even though the original cap is nearly planar. This changes only a
single amplitude when normalizing the physical signed flux to the target.
"""
import numpy as np
import math
from p2 import TRI_EDGES,triangle_rule7,tri_basis

def native_profile(measure):
    tri,av=measure.boundary['INLET'];s=np.zeros_like(measure.points)
    for face,area_vector in zip(tri,av):
        s[face[:3]]+=area_vector/3
        for a,(i,j) in enumerate(TRI_EDGES):s[face[3+a]]=(s[face[i]]+s[face[j]])/2
    ids=np.unique(tri);norm=s[ids]/np.linalg.norm(s[ids],axis=1)[:,None]
    shape=np.ones(len(ids));shape[np.isin(ids,measure.wall_nodes)]=0
    velocity=np.zeros_like(s);velocity[ids]=-norm*shape[:,None]
    fq,fw=triangle_rule7();N=tri_basis(fq);uq=np.einsum('qa,tai->tqi',N,velocity[tri])
    directed=-math.fsum(np.einsum('tqi,ti->tq',uq,av)@fw)
    scalar=np.zeros(len(s));scalar[ids]=shape
    integrated_scalar=math.fsum(np.linalg.norm(av,axis=1)*((scalar[tri]@N.T)@fw))
    return dict(node_ids=ids,flat_shape=shape,normals_outward=norm,velocity_m_s=velocity[ids]*(measure.Q/directed),
        unit_directed_flux_m2=directed,unit_scalar_integral_m2=integrated_scalar,
        alpha_m_s=measure.Q/directed,native_direction_factor=directed/integrated_scalar,
        XML_value_for_exact_requested_Q=-measure.Q*integrated_scalar/directed)

if __name__=='__main__':
    from common import CASE,REPORT,dump
    from p2_measure import P2Measurements
    m=P2Measurements(CASE);a=native_profile(m)
    np.savez_compressed(REPORT/'data/native_tet10_inlet_profile.npz',**{k:v for k,v in a.items() if isinstance(v,np.ndarray)})
    j={k:v for k,v in a.items() if not isinstance(v,np.ndarray)}
    j.update(source='baf_ini.cpp:face_ini higher-order normals and bc_ini Flat/perimeter-zero normalization',
        target_Q=m.Q,expected_corrected_Q=a['alpha_m_s']*a['unit_directed_flux_m2'],
        expected_unadjusted_actual_Q=m.Q*a['native_direction_factor'],XML_single_scalar_normalization=1/a['native_direction_factor'],
        face_quadrature='Native TRI6 seven-point rule, pinned nn_elem_gip.h:740-766',runtime_profile_check_pending=True)
    dump(REPORT/'data/native_inlet_direction_audit.json',j)
    print(j)
