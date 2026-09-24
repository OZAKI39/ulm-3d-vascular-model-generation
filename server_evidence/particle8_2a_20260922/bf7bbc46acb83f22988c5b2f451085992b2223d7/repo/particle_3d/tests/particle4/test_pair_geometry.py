import numpy as np
import pytest
from scipy.optimize import minimize, differential_evolution
from particle_3d.particle_shapes import Sphere,Ellipsoid,Capsule,EPS,unit
from particle_3d.pair_geometry import pair_gap
from particle_3d.particle4_cases import shapes_for_geometry,place_pair


@pytest.mark.parametrize('i,j',[(0,0),(0,1),(0,2),(1,1),(1,2),(2,2)])
@pytest.mark.parametrize('target',[1e-6,0.,-.02e-6])
def test_six_shape_pairs_signed_gap_reciprocity_and_support_witness(i,j,target):
    base=shapes_for_geometry();a,b=place_pair(base[i],base[j],[1,.3,.7],target)
    g=pair_gap(a,b,17,42);r=pair_gap(b,a,42,17)
    assert abs(g.gap_m-target)<=g.roundoff_budget_m
    assert g.gap_m==r.gap_m and g.canonical_pair_id==(17,42)
    np.testing.assert_array_equal(g.normal_j_to_i,-r.normal_j_to_i)
    np.testing.assert_array_equal(g.point_i_m,r.point_j_m)
    assert np.linalg.norm(g.point_i_m-g.point_j_m-g.gap_m*g.normal_j_to_i)<=g.roundoff_budget_m
    assert abs(np.linalg.norm(g.normal_j_to_i)-1)<16*EPS
    assert g.state==('SEPARATED' if target>0 else 'PENETRATING' if target<0 else 'TOUCHING')


def _segment(shape):
    if isinstance(shape,Sphere):return shape.center_m,shape.center_m
    v=shape.axis_world*shape.cylindrical_length_m/2
    return shape.center_m-v,shape.center_m+v


def segment_distance(a,b):
    p,q=_segment(a);r,s=_segment(b);u=q-p;v=s-r;w=p-r
    candidates=[]
    for t in [0.,1.]:
        z=np.clip(((p+t*u-r)@v)/(v@v),0,1) if v@v else 0.
        candidates.append(np.linalg.norm(p+t*u-r-z*v))
    for z in [0.,1.]:
        t=np.clip(((r+z*v-p)@u)/(u@u),0,1) if u@u else 0.
        candidates.append(np.linalg.norm(p+t*u-r-z*v))
    mat=np.column_stack([u,-v]);x=np.linalg.lstsq(mat,-w,rcond=None)[0]
    if np.all(x>=0) and np.all(x<=1):candidates.append(np.linalg.norm(w+mat@x))
    return min(candidates)


@pytest.mark.parametrize('kind',['sphere-sphere','sphere-capsule','capsule-capsule'])
def test_random_analytic_segment_references(kind):
    rng=np.random.default_rng(2026092004)
    for _ in range(24):
        c=rng.normal(size=(2,3))*1.5e-6;r=rng.uniform(.2,1.,2)*1e-6
        a=Sphere(c[0],r[0]) if kind!='capsule-capsule' else Capsule(c[0],rng.normal(size=3),r[0],1.2e-6)
        b=Sphere(c[1],r[1]) if kind=='sphere-sphere' else Capsule(c[1],rng.normal(size=3),r[1],1.4e-6)
        result=pair_gap(a,b);exact=segment_distance(a,b)-sum(r)
        assert abs(result.gap_m-exact)<=result.roundoff_budget_m


def primal_reference(a,b):
    """Independent convex point-in-shape optimization, with support dual bound."""
    scale=max(a.bounding_radius_m,b.bounding_radius_m);delta=(a.center_m-b.center_m)/scale
    def point(s,z):
        if isinstance(s,Ellipsoid):return s.rotation@(s.axes_m/scale*z[:3])
        return s.radius_m/scale*z[:3]+(s.axis_world*s.cylindrical_length_m/(2*scale)*z[3] if isinstance(s,Capsule) else 0)
    def difference(z):return delta+point(a,z[:4])-point(b,z[4:])
    cons=[{'type':'ineq','fun':lambda z:1-z[:3]@z[:3]}, {'type':'ineq','fun':lambda z:1-z[4:7]@z[4:7]}]
    opt=minimize(lambda z:.5*np.dot(difference(z),difference(z)),np.zeros(8),method='SLSQP',bounds=[(-1,1)]*8,constraints=cons,options={'ftol':16*EPS,'maxiter':1000})
    vector=difference(opt.x);upper=np.linalg.norm(vector)*scale;n=unit(vector)
    lower=float(n@(a.support(-n)-b.support(n)))
    # Feasible primal endpoints and a support plane bound are independent of GJK.
    assert max(np.linalg.norm(opt.x[:3]),np.linalg.norm(opt.x[4:7]))<=1+256*EPS
    assert upper-lower<=16384*EPS*scale
    return lower,upper


@pytest.mark.parametrize('j',[0,1,2])
def test_ellipsoid_pair_independent_convex_primal_dual_reference(j):
    base=shapes_for_geometry(1)
    for n in [[1,2,-1],[-2,1,3],[1,.1,.2]]:
        a,b=place_pair(base[1],base[j],n,.5e-6)
        result=pair_gap(a,b);lower,upper=primal_reference(a,b)
        assert lower-result.roundoff_budget_m<=result.gap_m<=upper+result.roundoff_budget_m


@pytest.mark.parametrize('j',[0,1,2])
def test_ellipsoid_overlap_independent_global_angular_reference(j):
    base=shapes_for_geometry(3);a,b=place_pair(base[1],base[j],[1,.3,.7],-.15e-6)
    result=pair_gap(a,b);scale=max(a.bounding_radius_m,b.bounding_radius_m)
    def objective(angles):
        t,p=angles;n=np.array([np.sin(t)*np.cos(p),np.sin(t)*np.sin(p),np.cos(t)])
        return -float(n@(a.support(-n)-b.support(n)))/scale
    ref=differential_evolution(objective,[(0,np.pi),(-np.pi,np.pi)],seed=20260920,tol=64*EPS,atol=64*EPS,popsize=20,maxiter=1000,polish=True)
    assert abs(result.gap_m+ref.fun*scale)<=result.roundoff_budget_m


def test_capsule_kink_witness_and_equal_center_penetration():
    for a,b in [(Capsule([0,0,0],[0,0,1],.4e-6,2e-6),Capsule([0,0,0],[0,1,0],.6e-6,3e-6)),
                (Sphere([0,0,0],.4e-6),Sphere([0,0,0],.6e-6))]:
        g=pair_gap(a,b);assert abs(g.gap_m+1e-6)<=g.roundoff_budget_m
        assert np.linalg.norm(g.point_i_m-g.point_j_m-g.gap_m*g.normal_j_to_i)<=g.roundoff_budget_m


def deep_overlap_reference_records():
    from particle_3d.rbc_orientation import rotation_matrix,quaternion_from_short_axis
    rng=np.random.default_rng(2026092404)
    records=[]
    for index in range(8):
        g=shapes_for_geometry(index%5)[1]
        a=Ellipsoid([0,0,0],g.axes_m,rotation_matrix(quaternion_from_short_axis(rng.normal(size=3))))
        b=Ellipsoid(rng.normal(size=3)*.2e-6,shapes_for_geometry((index+2)%5)[1].axes_m,rotation_matrix(quaternion_from_short_axis(rng.normal(size=3))))
        result=pair_gap(a,b);scale=max(a.bounding_radius_m,b.bounding_radius_m)
        def objective(angles):
            t,p=angles;n=np.array([np.sin(t)*np.cos(p),np.sin(t)*np.sin(p),np.cos(t)])
            return -float(n@(a.support(-n)-b.support(n)))/scale
        # A single DE population can collapse onto the wrong antipodal basin
        # despite reporting success. Independently search all eight angular
        # sectors; never seed the reference with the production normal.
        solutions=[]
        for hemisphere in range(2):
            for quadrant in range(4):
                bounds=[(hemisphere*np.pi/2,(hemisphere+1)*np.pi/2),(-np.pi+quadrant*np.pi/2,-np.pi+(quadrant+1)*np.pi/2)]
                reference=differential_evolution(objective,bounds,seed=100+index*8+hemisphere*4+quadrant,popsize=20,tol=64*EPS,atol=64*EPS,maxiter=1000,polish=True)
                solutions.append(float(reference.fun))
        exact=-min(solutions)*scale
        records.append(dict(pair_type='Ellipsoid-Ellipsoid',reference='INDEPENDENT_EIGHT_SECTOR_GLOBAL_ANGULAR_SEARCH',case=index,
            gap_m=result.gap_m,state=result.state,reference_gap_m=exact,sector_depths_m=np.array(solutions)*scale,
            error_m=abs(result.gap_m-exact),roundoff_m=result.roundoff_budget_m))
    return records


def test_noncoaxial_deep_ellipsoids_against_independent_global_reference():
    for row in deep_overlap_reference_records():
        assert row['state']=='PENETRATING'
        assert row['error_m']<=row['roundoff_m']


def test_pair_geometry_rigid_rotation_and_translation_covariance():
    from particle_3d.rbc_orientation import rotation_matrix,quaternion_from_short_axis
    q=rotation_matrix(quaternion_from_short_axis([1,-2,3]));offset=np.array([1.,2.,-1.])*1e-4
    a,b=place_pair(shapes_for_geometry()[1],shapes_for_geometry()[2],[1,.2,.5],.3e-6)
    c=Ellipsoid(q@a.center_m+offset,a.axes_m,q@a.rotation)
    d=Capsule(q@b.center_m+offset,q@b.axis_world,b.radius_m,b.cylindrical_length_m)
    before=pair_gap(a,b);after=pair_gap(c,d)
    assert abs(before.gap_m-after.gap_m)<=after.roundoff_budget_m
    assert np.linalg.norm(q@before.point_i_m+offset-after.point_i_m)<=after.roundoff_budget_m


@pytest.mark.parametrize('axis_a,axis_b,normal',[
    ([0,0,1],[0,0,1],[0,0,1]),  # both broad sides
    ([1,2,3],[-1,1,2],[1,.2,.3]),  # different oblique orientations
    ([0,0,1],[0,0,1],[1,0,0]),  # equatorial ends
])
def test_different_rbc_samples_broadside_oblique_and_end_on_contact(axis_a,axis_b,normal):
    from particle_3d.rbc_orientation import rotation_matrix,quaternion_from_short_axis
    a=shapes_for_geometry(1)[1];b=shapes_for_geometry(3)[1]
    a=Ellipsoid(a.center_m,a.axes_m,rotation_matrix(quaternion_from_short_axis(axis_a)))
    b=Ellipsoid(b.center_m,b.axes_m,rotation_matrix(quaternion_from_short_axis(axis_b)))
    a,b=place_pair(a,b,normal);g=pair_gap(a,b)
    assert g.state=='TOUCHING'
    assert np.linalg.norm(g.normal_j_to_i+unit(normal))<=4096*EPS
