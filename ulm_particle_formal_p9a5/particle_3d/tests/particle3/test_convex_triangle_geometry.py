import numpy as np
import pytest
from particle_3d.particle_shapes import Sphere,Ellipsoid,Capsule
from particle_3d.convex_triangle import triangle_gap
from particle_3d.rbc_orientation import rotation_matrix,quaternion_from_short_axis


@pytest.mark.parametrize('height',[2.,1.,.6])
@pytest.mark.parametrize('xy,kind',[(np.array([2.,2.]),'FACE'),(np.array([4.,-.3]),'EDGE'),(np.array([-.3,-.4]),'VERTEX')])
def test_sphere_finite_triangle_face_edge_vertex_and_penetration(height,xy,kind):
    tri=np.array([[0,0,0],[0,10,0],[10,0,0]],dtype=float)*1e-6
    center=np.r_[xy,height]*1e-6
    projected=np.r_[np.maximum(xy,0),0]*1e-6
    expected=np.linalg.norm(center-projected)-1e-6
    g=triangle_gap(Sphere(center,1e-6),tri)
    assert abs(g.gap_m-expected)<64*np.finfo(float).eps*1e-5
    np.testing.assert_allclose(g.wall_point_m,projected,rtol=0,atol=1e-19)
    assert g.wall_feature==kind
    assert abs(np.linalg.norm(g.normal_inward)-1)<8*np.finfo(float).eps


@pytest.mark.parametrize('angle',[0.,.2,.8,np.pi/2])
@pytest.mark.parametrize('height',[3.,1.,.1])
def test_rotated_ellipsoid_plane_signed_gap(plane_triangle,angle,height):
    r=np.array([[np.cos(angle),0,np.sin(angle)],[0,1,0],[-np.sin(angle),0,np.cos(angle)]])
    shape=Ellipsoid([0,0,height*1e-6],np.array([2,2,.5])*1e-6,r)
    support=np.sqrt((2*np.sin(angle))**2+(.5*np.cos(angle))**2)*1e-6
    g=triangle_gap(shape,plane_triangle)
    assert abs(g.gap_m-(height*1e-6-support))<2e-18
    assert abs(np.linalg.norm((r.T@(g.particle_point_m-shape.center_m))/shape.axes_m)-1)<1e-12


@pytest.mark.parametrize('center',[[2,2,1.5],[2,2,.2],[4,-.4,.3],[-.3,-.4,.2],[-2,-2,2]])
def test_ellipsoid_sphere_limit_including_edge_vertex_overlap(center):
    tri=np.array([[0,0,0],[0,10,0],[10,0,0]],dtype=float)*1e-6
    c=np.asarray(center)*1e-6
    e=triangle_gap(Ellipsoid(c,np.ones(3)*1e-6,np.eye(3)),tri)
    s=triangle_gap(Sphere(c,1e-6),tri)
    assert abs(e.gap_m-s.gap_m)<2e-18


@pytest.mark.parametrize('length',[0.,1.,4.])
@pytest.mark.parametrize('height',[3.,.7,.1])
def test_capsule_plane_gap_and_overlap(plane_triangle,length,height):
    axis=np.array([1,0,2])/np.sqrt(5)
    cap=Capsule([0,0,height*1e-6],axis,.5e-6,length*1e-6)
    g=triangle_gap(cap,plane_triangle)
    expected=height*1e-6-(.5e-6+length*1e-6/2*axis[2])
    assert abs(g.gap_m-expected)<3e-18


def test_support_mapping_independent_constraints():
    r=rotation_matrix(quaternion_from_short_axis([1,2,3]));axes=np.array([3,2,.5])*1e-6
    e=Ellipsoid([1e-5,2e-5,0],axes,r)
    rng=np.random.default_rng(3103)
    for d in rng.normal(size=(20,3)):
        p=e.support(d);body=r.T@(p-e.center_m)
        assert abs(np.sum((body/axes)**2)-1)<1e-13
        expected=np.linalg.norm(axes*(r.T@d))
        assert abs(d@(p-e.center_m)-expected)<1e-19
        cap=Capsule(e.center_m,[1,2,3],1e-6,4e-6)
        assert abs(d@(cap.support(d)-cap.center_m)-(2e-6*abs(d@cap.axis_world)+1e-6*np.linalg.norm(d)))<1e-19
