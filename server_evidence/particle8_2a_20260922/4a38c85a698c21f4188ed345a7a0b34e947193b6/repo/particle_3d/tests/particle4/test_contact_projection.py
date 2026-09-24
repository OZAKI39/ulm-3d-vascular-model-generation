from itertools import permutations
import numpy as np
import pytest
from particle_3d.particle_shapes import Sphere,Ellipsoid,Capsule,EPS,unit
from particle_3d.pair_geometry import pair_gap
from particle_3d.kinematic_contact import ContactConstraint,project_contacts,MultiContactInfeasible
from particle_3d.particle4_cases import shapes_for_geometry,place_pair
from particle_3d.particle3_cases import plane_triangle
from particle_3d.wall_geometry import WallGeometry
from particle_3d.wall_gap import touching_contacts


def solve(shapes,v,w=None,wall=None):
    ids=sorted(shapes);contacts=[]
    for k,i in enumerate(ids):
        for j in ids[k+1:]:
            g=pair_gap(shapes[i],shapes[j],i,j)
            if g.state=='TOUCHING':contacts.append(ContactConstraint.pair(g))
    if wall is not None:
        for i,s in shapes.items():contacts.extend(ContactConstraint.wall(i,g) for g in touching_contacts(s,wall))
    return project_contacts(shapes,v,{i:np.zeros(3) for i in shapes} if w is None else w,contacts)


@pytest.mark.parametrize('tangent',[0.,.3])
def test_mb_mb_head_on_glancing_symmetry_tangent_no_restitution_or_spin(tangent):
    shapes={7:Sphere([-1e-6,0,0],1e-6),11:Sphere([1e-6,0,0],1e-6)}
    v={7:np.array([1.,tangent,0])*1e-6,11:np.array([-1.,-tangent,0])*1e-6};w={7:np.array([1.,2,3]),11:np.array([-1.,2,5])}
    p=solve(shapes,v,w)
    assert np.linalg.norm(p.translation_corrections[7]+p.translation_corrections[11])<=p.record['velocity_budget_m_s']
    for i in shapes:
        assert abs(p.velocities[i][0])<=p.record['velocity_budget_m_s']
        assert p.velocities[i][1]==v[i][1]
        np.testing.assert_array_equal(p.omegas[i],w[i])


@pytest.mark.parametrize('separation',[0.,1e-6])
def test_separating_or_noncontact_does_not_change_velocity(separation):
    a,b=place_pair(*shapes_for_geometry()[:2],[1,0,0],separation)
    v={1:np.array([-1.,0,0])*1e-6,2:np.array([1.,0,0])*1e-6};p=solve({1:a,2:b},v)
    for i in v:np.testing.assert_array_equal(p.velocities[i],v[i])


def test_offcenter_rbc_angular_correction_and_contact_jacobian_identity():
    a,b=place_pair(shapes_for_geometry()[1],shapes_for_geometry()[0],[1,.4,.3])
    gap=pair_gap(a,b,1,2);n=gap.normal_j_to_i;v={1:-n*1e-6,2:n*1e-6};w={1:np.array([.1,.2,.3]),2:np.zeros(3)}
    p=solve({1:a,2:b},v,w);r=gap.point_i_m-a.center_m
    assert np.linalg.norm(p.angular_corrections[1])>0
    expected=p.record['multipliers'][0]*np.cross(r,n)/a.bounding_radius_m**2
    np.testing.assert_allclose(p.angular_corrections[1],expected,rtol=64*EPS,atol=0)
    free=n@(v[1]+np.cross(w[1],r)-v[2])
    corrected=n@(p.velocities[1]+np.cross(p.omegas[1],r)-p.velocities[2])
    assert abs(free-p.record['free_normal_speeds_m_s'][0])<=p.record['velocity_budget_m_s']
    assert corrected>=-p.record['velocity_budget_m_s']
    assert np.linalg.norm(p.translation_corrections[1]+p.translation_corrections[2])<=p.record['velocity_budget_m_s']


@pytest.mark.parametrize('other',[0,1,2])
def test_capsule_contact_never_changes_axis_spin_or_dimensions(other):
    base=shapes_for_geometry();a,b=place_pair(base[2],base[other],[.3,1,.2]);axis=a.axis_world.copy()
    g=pair_gap(a,b,1,2);n=g.normal_j_to_i;w={1:np.array([1.,2,3]),2:np.array([3.,1,2])}
    p=solve({1:a,2:b},{1:-n*1e-6,2:n*1e-6},w)
    np.testing.assert_array_equal(p.omegas[1],w[1]);np.testing.assert_array_equal(a.axis_world,axis)
    assert np.linalg.norm(p.angular_corrections[1])==0 and a.radius_m==base[2].radius_m
    assert a.cylindrical_length_m==base[2].cylindrical_length_m


@pytest.mark.parametrize('scene',['chain','triangle','four'])
def test_simultaneous_contacts_kkt_and_all_input_permutations(scene):
    centers=[[-2,0,0],[0,0,0],[2,0,0]] if scene=='chain' else [[-1,0,0],[1,0,0],[0,np.sqrt(3),0]]
    if scene=='four':centers=[[-3,0,0],[-1,0,0],[1,0,0],[3,0,0]]
    ids=[17,3,91,42][:len(centers)];shapes={i:Sphere(np.array(c)*1e-6,1e-6) for i,c in zip(ids,centers)}
    centroid=np.mean(list(s.center_m for s in shapes.values()),axis=0);v={i:centroid-s.center_m for i,s in shapes.items()}
    baseline=solve(shapes,v)
    for order in permutations(ids):
        p=solve({i:shapes[i] for i in order},v)
        for i in ids:
            np.testing.assert_array_equal(p.velocities[i],baseline.velocities[i]);np.testing.assert_array_equal(p.omegas[i],baseline.omegas[i])
        assert p.record['canonical_ids']==baseline.record['canonical_ids']
        assert p.record['max_primal_violation_m_s']<=p.record['velocity_budget_m_s']
        assert p.record['max_complementarity_residual']<=p.record['complementarity_budget']
        assert np.linalg.norm(sum(p.translation_corrections.values()))<=p.record['velocity_budget_m_s']


@pytest.mark.parametrize('rbc',[False,True])
def test_wall_and_pair_enter_same_solve(rbc):
    wall=WallGeometry(plane_triangle()[None,:,:]);a=Sphere([0,0,1e-6],1e-6)
    b=shapes_for_geometry()[1] if rbc else Sphere([0,0,0],1e-6)
    a,b=place_pair(a,b,[0,0,1]);shapes={1:a,2:b};v={1:np.array([0,0,-1e-6]),2:np.array([0,0,-2e-6])}
    p=solve(shapes,v,wall=wall)
    assert {c[0] for c in p.record['canonical_ids']}=={'WALL','PAIR'}
    assert p.record['max_primal_violation_m_s']<=p.record['velocity_budget_m_s']
    assert p.record['max_complementarity_residual']<=p.record['complementarity_budget']


def test_rank_incompatible_constraints_report_explicit_infeasibility():
    shape=Sphere([0,0,0],1e-6);constraints=[ContactConstraint(('WALL',1,k),1,None,np.array([sign,0,0]),np.zeros(3),offset_m_s=-1e-6) for k,sign in enumerate([-1,1])]
    with pytest.raises(MultiContactInfeasible,match='MULTI_CONTACT_INFEASIBLE'):
        project_contacts({1:shape},{1:np.zeros(3)},{1:np.zeros(3)},constraints)
