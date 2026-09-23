"""Analytical and geometric guards for the new derived WSS diagnostics."""
from pathlib import Path
import sys
import numpy as np
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'scripts/flow_2mmps'))
from compute_field_diagnostics import p1_gradients, tangential_traction, boundary_owners, wall_geometry, nodal_average


def sample_tetra():
    return np.array([[2,3,1],[4,3,1],[2.2,6,1],[2.5,3.8,5]], dtype=float)*1e-6, np.array([[0,1,2,3]])


def test_exact_affine_gradient_in_si_coordinates():
    points, tetra = sample_tetra()
    expected = np.array([[2,3,-4],[1,7,2],[-8,1,4]], dtype=float)*1e3
    velocity = points@expected.T + [0.01,-0.02,0.03]
    np.testing.assert_allclose(p1_gradients(points,tetra,velocity)[0],expected,rtol=1e-13,atol=1e-10)


def test_couette_wall_stress_has_known_magnitude():
    g = np.zeros((1,3,3));g[0,0,1]=1234.
    result = tangential_traction(g,np.array([[0.,1.,0.]]),.00345312)
    np.testing.assert_allclose(result,[[.00345312*1234,0,0]],rtol=1e-14)


def test_rigid_body_rotation_produces_no_viscous_stress():
    g = np.array([[[0,-7,3],[7,0,-4],[-3,4,0.]]])
    n = np.array([[1.,2.,3.]]);n/=np.linalg.norm(n)
    np.testing.assert_allclose(tangential_traction(g,n,.003),0,atol=1e-15)


def test_normal_pressure_does_not_change_tangential_traction():
    rng=np.random.default_rng(23);g=rng.normal(size=(10,3,3))*1000
    n=rng.normal(size=(10,3));n/=np.linalg.norm(n,axis=1)[:,None]
    expected=tangential_traction(g,n,.003)
    total=.003*np.einsum('nij,nj->ni',g+g.swapaxes(1,2),n)-1234*n
    tangential=total-n*np.einsum('ij,ij->i',total,n)[:,None]
    np.testing.assert_allclose(tangential,expected,atol=1e-12)
    np.testing.assert_allclose(np.einsum('ij,ij->i',expected,n),0,atol=1e-14)


def test_wss_rotates_covariantly_and_magnitude_ignores_normal_sign():
    g=np.array([[[2.,3.,1.],[4.,-2.,5.],[3.,6.,0.]]]);n=np.array([[1.,2.,3.]])/np.sqrt(14)
    angle=.63;r=np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1.]])
    t=tangential_traction(g,n,.003)
    np.testing.assert_allclose(tangential_traction(r@g@r.T,n@r.T,.003),t@r.T,atol=1e-16)
    np.testing.assert_allclose(tangential_traction(g,-n,.003),-t,atol=1e-16)


def test_boundary_owner_and_outward_normal():
    p,t=sample_tetra();faces=np.array([[2,1,0],[0,1,3],[0,3,2],[3,1,2]])
    owners=boundary_owners(t,faces)
    np.testing.assert_array_equal(owners,[0,0,0,0])
    c,a,n=wall_geometry(p,t,faces,owners)
    assert np.all(a>0) and np.all(np.sum((c-p.mean(0))*n,axis=1)>0)
    np.testing.assert_allclose(np.sum(n*a[:,None],axis=0),0,atol=1e-26)


def test_interior_face_cannot_be_treated_as_wall():
    tetra=np.array([[0,1,2,3],[0,1,2,4]])
    with pytest.raises(AssertionError,match='exactly one'):
        boundary_owners(tetra,np.array([[0,1,2]]))


def test_display_reconstruction_preserves_area_integral():
    faces=np.array([[0,1,2],[1,2,3]]);area=np.array([1.,3.]);values=np.array([2.,6.])
    avg,weights=nodal_average(faces,values,area,4)
    np.testing.assert_allclose(avg,[2,5,5,6])
    np.testing.assert_allclose(np.sum(avg*weights)/3,np.sum(area*values))


def test_fixed_axis_rotation_preserves_axis_and_distances():
    from render_field_diagnostics import axis_rotation
    axis=np.array([1.,-2.,3.]);axis/=np.linalg.norm(axis)
    q=np.array([[1.,4.,5.],[2.,7.,-3.]])
    for angle in [0,37,90,180,271,360]:
        rotation=axis_rotation(axis,angle)
        np.testing.assert_allclose(rotation@axis,axis,atol=1e-14)
        np.testing.assert_allclose(rotation.T@rotation,np.eye(3),atol=1e-14)
        assert np.linalg.det(rotation)==pytest.approx(1.)
        np.testing.assert_allclose(np.linalg.norm(q@rotation.T,axis=1),np.linalg.norm(q,axis=1))
    np.testing.assert_allclose(axis_rotation(axis,360),np.eye(3),atol=1e-14)


def test_turntable_never_changes_camera_center_or_scale():
    from render_field_diagnostics import FixedAxisRotation,FRAMES
    rng=np.random.default_rng(42)
    points=rng.normal(size=(120,3))*[10,5,30]+[50,80,90]
    turntable=FixedAxisRotation(points)
    for key in ['camera_position_um','focal_point_um','rotation_center_um','rotation_axis_unit','parallel_scale_um']:
        values=np.array([r[key] for r in turntable.trace])
        np.testing.assert_array_equal(values,np.repeat(values[:1],FRAMES,axis=0))
    np.testing.assert_allclose(np.diff([r['rotation_angle_deg'] for r in turntable.trace]),360/FRAMES)
    for matrix in turntable.transforms:
        for p in [turntable.center,turntable.center+17*turntable.axis]:
            np.testing.assert_allclose(matrix[:3,:3]@p+matrix[:3,3],p,atol=1e-12)


def test_continuous_turntable_bound_keeps_full_geometry_visible():
    from render_field_diagnostics import FixedAxisRotation,axis_rotation,ASPECT
    rng=np.random.default_rng(7);points=rng.normal(size=(150,3))*[10,5,30]
    turntable=FixedAxisRotation(points)
    for angle in [2.37,37.123,107.4,248.913,359.7]:
        q=(points-turntable.center)@axis_rotation(turntable.axis,angle).T
        assert np.max(abs(q@turntable.right))/(turntable.scale*ASPECT)<.946
        assert np.max(abs(q@turntable.up))/turntable.scale<.946
