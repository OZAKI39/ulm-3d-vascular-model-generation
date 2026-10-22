import numpy as np
from scipy.optimize import minimize_scalar
from particle_3d.particle_shapes import Ellipsoid,Capsule
from particle_3d.convex_triangle import triangle_gap,capsule_triangle_many
from particle_3d.rbc_orientation import rotation_matrix,quaternion_from_short_axis


def test_nonspherical_edge_vertex_distance_against_independent_surface_minimization():
    # An independent angle-parameterized smooth ellipse boundary provides the
    # edge/vertex reference. No production secular roots or support solver used.
    tri=np.array([[0.,0.,0.],[0.,10.,0.],[10.,0.,0.]])
    for center in [[4.,-.7,.2],[-.4,-.5,.3],[4.,-2.5,.2],[-2.,-2.,.5]]:
        shape=Ellipsoid(center,[2.,2.,.6],np.eye(3))
        actual=triangle_gap(shape,tri)
        radial=-center[1] if center[0]>0 else np.linalg.norm(center[:2])
        target=np.array([radial,-center[2]])
        def objective(angle):return float(np.sum((np.array([2*np.cos(angle),.6*np.sin(angle)])-target)**2))
        angles=np.linspace(-np.pi,np.pi,33)
        squared=min(minimize_scalar(objective,bounds=(a,b),method='bounded',options={'xatol':1e-14}).fun for a,b in zip(angles[:-1],angles[1:]))
        expected=np.sqrt(squared)*(-1 if np.sum((target/[2,.6])**2)<1 else 1)
        assert abs(actual.gap_m-expected)<5e-12
        # Independent rigid-coordinate covariance catches world/body and finite-edge errors.
        rot=rotation_matrix(quaternion_from_short_axis([1,2,3]));shift=np.array([5.,-4.,2.])
        moved=Ellipsoid(rot@shape.center_m+shift,shape.axes_m,rot)
        other=triangle_gap(moved,tri@rot.T+shift)
        assert abs(actual.gap_m-other.gap_m)<5e-12


def test_vectorized_capsule_is_scalar_finite_feature_geometry():
    rng=np.random.default_rng(303)
    triangles=rng.normal(size=(40,3,3))*2e-6
    for axis in [[0,0,1],[1,2,3],[1,0,0]]:
        cap=Capsule([.2e-6,-.3e-6,.4e-6],axis,.8e-6,3e-6)
        values=capsule_triangle_many(cap,triangles)[0]
        for i,tri in enumerate(triangles):
            assert abs(values[i]-triangle_gap(cap,tri).gap_m)<3e-18


def test_random_ellipsoid_separation_against_independent_convex_primal():
    from scipy.optimize import minimize,least_squares
    rng=np.random.default_rng(30303)
    for _ in range(24):
        rotation=rotation_matrix(quaternion_from_short_axis(rng.normal(size=3)))
        shape=Ellipsoid(rng.normal(size=3)*.2,[2.,1.3,.45],rotation)
        tri=rng.normal(size=(3,3))*.3+np.array([4.,1.,1.])
        transform=rotation@np.diag(shape.axes_m)
        def residual(x):return shape.center_m+transform@x[:3]-(x[3:]@tri)
        matrix=np.column_stack((transform,-tri.T))
        result=minimize(lambda x:residual(x)@residual(x),np.r_[np.zeros(3),[1/3]*3],method='SLSQP',
            jac=lambda x:2*matrix.T@residual(x),
            constraints=[{'type':'ineq','fun':lambda x:1-x[:3]@x[:3],'jac':lambda x:np.r_[-2*x[:3],np.zeros(3)]},
                         {'type':'eq','fun':lambda x:x[3:].sum()-1,'jac':lambda x:np.r_[np.zeros(3),np.ones(3)]}],
            bounds=[(None,None)]*3+[(0,1)]*3,options={'ftol':1e-13,'maxiter':400})
        # Polish the active-feature KKT system independently of production
        # secular roots; SLSQP's objective stopping criterion alone is too weak.
        ids=np.flatnonzero(result.x[3:]>1e-7);selected=tri[ids]
        edges=(selected[1:]-selected[0]).T
        def kkt(x):
            z=x[:3];weights=x[3:-1];delta=shape.center_m+transform@z-selected[0]-edges@weights
            return np.r_[transform.T@delta+x[-1]*z,edges.T@delta,z@z-1]
        z=result.x[:3];lam=-z@(transform.T@residual(result.x))
        refined=least_squares(kkt,np.r_[z,result.x[3+ids[1:]],lam],ftol=2e-15,xtol=2e-15,gtol=2e-15,max_nfev=200)
        assert np.linalg.norm(kkt(refined.x))<1e-11
        z=refined.x[:3]/max(1.,np.linalg.norm(refined.x[:3]))
        weights=np.r_[1-refined.x[3:-1].sum(),refined.x[3:-1]]
        weights=np.maximum(weights,0);weights/=weights.sum()
        # Certify the reference by matching independent convex primal/dual
        # bounds; an optimizer status string is not a numerical certificate.
        delta=shape.center_m+transform@z-weights@selected;upper=np.linalg.norm(delta);n=delta/upper
        lower=n@shape.center_m-np.linalg.norm(transform.T@n)-np.max(tri@n)
        assert upper-lower<2e-10,(result.message,upper-lower)
        assert abs(triangle_gap(shape,tri).gap_m-(upper+lower)/2)<2e-10


def test_bvh_and_enclosing_ball_pruning_match_every_actual_wall_triangle(real_wall,p3_repo):
    import json
    from particle_3d.wall_gap import wall_gap
    records=json.loads((p3_repo/'particle_3d/reports/particle3/data/11_real_rbc_g0_dt0_states.json').read_text())
    # Exhaustive 45,221-triangle checks independently validate both levels of
    # candidate pruning, including a deliberately overlapping larger capsule.
    for index,scale in [(0,1.),(len(records)//2,1.),(len(records)-1,1.1)]:
        row=records[index];shape=Capsule(row['center_m'],row['capsule_axis_world'],row['R_cap_m']*scale,row['L_cap_m'])
        all_gaps=capsule_triangle_many(shape,real_wall.triangles)[0]
        actual=wall_gap(shape,real_wall)
        assert abs(actual.gap_m-float(all_gaps.min()))<=actual.roundoff_m
