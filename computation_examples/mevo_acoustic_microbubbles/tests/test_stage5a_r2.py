"""R2 real-matrix regressions, native mapping, geometry and full integrator checks."""
import json
from functools import partial
from pathlib import Path
import numpy as np
import pytest
from scipy import sparse
from scipy.integrate import solve_ivp
import test_stage5a_motion as old
from mevo_mb.bubbles.motion_qp_r2 import NativeQP,certificate
from mevo_mb.bubbles.motion_r2 import MotionR2,SemiImplicit
from mevo_mb.bubbles.motion_sweep_r2 import pair_sweep,segment_triangle2,triangle_minimum
from mevo_mb.bubbles.constrained_motion import AcceptedCache, SubdivideInterval


@pytest.mark.parametrize('fixture',[
 '20261006T200738743816Z_failure_replay','20261006T210101208569Z_B2_osqp','20261006T210102848863Z_B3_osqp'])
def test_real_failed_qp(fixture):
    p=Path('results/stage5/motion_solver_r2/fixtures')/fixture
    d=dict(np.load(next(z for z in p.glob('*.npz') if z.name not in ['A.npz','G.npz'])))
    a,g=sparse.load_npz(p/'A.npz'),sparse.load_npz(p/'G.npz')
    qp=NativeQP();u,l,diag=qp.solve(a,d['b'],g,d['lower'],feasible=d.get('feasible'))
    assert diag['pass_gate'] and diag['kkt']<=1e-9
    v,m,w=qp.solve(a,d['b'],g,d['lower'],warm_primal=u,warm_dual=l)
    np.testing.assert_allclose(u,v,atol=1e-8,rtol=1e-8)
    np.testing.assert_allclose(g.T@l,g.T@m,atol=1e-8,rtol=1e-7)
    assert abs(w['objective']-diag['objective'])<1e-9*max(1,abs(w['objective']))


@pytest.mark.parametrize('solver',['prox','prox_sparse','osqp','clarabel'])
def test_native_lower_sign_and_bound(solver):
    a=sparse.eye(2,format='csc');g=sparse.csc_matrix([[1.,0.]])
    u,l,d=NativeQP(solver,1e-10).candidate(a,np.array([-1.,2.]),g,np.array([.5]))
    np.testing.assert_allclose(u,[.5,2.],atol=1e-8)
    np.testing.assert_allclose(l,[1.5],atol=1e-8)
    assert certificate(a,np.array([-1.,2.]),g,np.array([.5]),u,-l)['kkt']>1


@pytest.mark.parametrize('epsilon',[0.,1e-12,1e-8])
def test_parallel_different_bounds_and_nonunique_dual(epsilon):
    a=sparse.eye(2,format='csc');g=sparse.csc_matrix([[1.,0.],[1.,epsilon],[1.,0.]])
    lo=np.array([.1,.2,.2]);b=np.array([-1.,1.])
    u,l,d=NativeQP().solve(a,b,g,lo)
    assert d['pass_gate'] and np.min(g@u-lo)>-1e-9
    np.testing.assert_allclose(u,[.2,1.],atol=1e-7)


@pytest.mark.parametrize('objective_scale',[.01,1.,100.])
@pytest.mark.parametrize('row_scale',[.1,1.,10.])
def test_objective_units_row_scaling(objective_scale,row_scale):
    a=sparse.csc_matrix([[2.,.2],[.2,3.]])
    b=np.array([-1.,2.]);g=sparse.eye(2,format='csc');lo=np.array([0.,0.])
    u,l,d=NativeQP().solve(a*objective_scale,b*objective_scale,g*row_scale,lo*row_scale)
    audit=certificate(a,b,g,lo,u,l*row_scale/objective_scale)
    assert audit['pass_gate']
    np.testing.assert_allclose(u,[0.,2/3],atol=1e-9)


@pytest.mark.parametrize('name',['test_T1_stokes_force_and_flow','test_T2_common_translation',
 'test_T3_squeeze_large_dt','test_T4_separating_inactive','test_T5_tangential_not_locked',
 'test_T6_wall_large_dt','test_T10_components_global','test_T11_warm_cold_and_cache_ids'])
def test_physical_limits(monkeypatch,name):
    # Exercise the new adapter with the original model-limit reference tests.
    class Adapter(MotionR2):
        def constrained(self,*args,**kw):
            return super().constrained(*args,**kw)
    monkeypatch.setattr(old,'ConstrainedLubrication',Adapter)
    getattr(old,name)()


@pytest.mark.parametrize('count',[3,30])
def test_chains(monkeypatch,count):
    monkeypatch.setattr(old,'ConstrainedLubrication',MotionR2)
    old.test_contact_chain(count)


def test_endpoint_safe_midpath_wall():
    tri=np.array([[-2.,-2.,0.],[2.,-2.,0.],[0.,2.,0.]])
    p=np.array([0.,0.,-1.]);q=np.array([0.,0.,1.])
    assert segment_triangle2(p,p,*tri)==pytest.approx(1.)
    assert segment_triangle2(q,q,*tri)==pytest.approx(1.)
    assert segment_triangle2(p,q,*tri)==0.
    # Finite radius detects an edge graze even without centerline intersection.
    p=np.array([2.1,-2.,-1.]);q=np.array([2.1,-2.,1.])
    assert triangle_minimum(p,q,tri[None])==pytest.approx(.1)


def analytic_sweep(geometry):
    def sweep(x,y,r,tol):
        pg=pair_sweep(x,y,r);w=min(geometry.wall_gap_and_normal(x,r)[0].min(),geometry.wall_gap_and_normal(y,r)[0].min())
        return dict(pass_gate=min(pg,w)>=-tol,minimum_pair_swept_gap_m=pg,
                    minimum_wall_swept_gap_bound_m=w)
    return sweep


@pytest.mark.parametrize('wall',[False,True])
def test_whole_integrator_squeeze_wall_refinement(wall):
    R=old.R;v=.001;model=MotionR2(old.MU)
    geometry=old.Plane() if wall else old.NoWall()
    x0=np.array([[0.,0.,-1.05*R]]) if wall else old.pair(.05)
    flow=np.array([[0.,0.,v]]) if wall else np.array([[v,0.,0.],[-v,0.,0.]])
    def forcing(x,r,t):
        return dict(primary_force_N=np.zeros_like(x),secondary_force_N=np.zeros_like(x),flow_velocity_m_s=flow)
    def rhs(t,h):
        c=model.coefficient(h,R,'wall' if wall else 'pair')[0]
        return [-v/(1+c/old.DRAG)] if wall else [-2*v/(1+2*c/old.DRAG)]
    duration=20e-6
    ref=solve_ivp(rhs,[0,duration],[.05*R],rtol=1e-12,atol=1e-20).y[0,-1]
    errors=[]
    for dt in [5e-6,2.5e-6,1.25e-6]:
        step=SemiImplicit(MotionR2(old.MU),geometry,forcing,analytic_sweep(geometry))
        x=x0.copy();r=np.full(len(x),R)
        for i in range(round(duration/dt)):
            x,_,_,_=step.advance(x,r,np.arange(len(x)),i*dt,dt,'test')
        gap=-x[0,2]-R if wall else np.linalg.norm(x[1]-x[0])-2*R
        errors.append(abs(gap-ref))
    assert errors[1]<errors[0]*.65 and errors[2]<errors[1]*.65


def test_negative_initial_band_release_and_cache_rollback():
    geometry=old.NoWall();x=old.pair();x[1,0]-=5e-12;r=np.full(2,old.R)
    def forcing(x,r,t):
        return dict(primary_force_N=np.zeros_like(x),secondary_force_N=np.zeros_like(x),
                    flow_velocity_m_s=np.array([[-.001,0.,0.],[.001,0.,0.]]))
    model=MotionR2(old.MU);prior=AcceptedCache((5,9),{}, {},(), 'accepted');model.cache=prior
    step=SemiImplicit(model,geometry,forcing,analytic_sweep(geometry))
    y,u,cache,d=step.advance(x,r,np.array([5,9]),0,1e-8,'test')
    assert np.linalg.norm(y[1]-y[0])>np.linalg.norm(x[1]-x[0])
    assert model.cache is prior


def test_atomic_checkpoint_restore(tmp_path):
    from scripts.stage5a_r2.run import checkpoint,restore
    c=AcceptedCache((3,),{3:np.ones(3)},{('wall',3):2.},(('wall',3),),'test')
    p=tmp_path/'checkpoint.npz';x=np.zeros((1,3));r=np.array([old.R]);ids=np.array([3]);mask=np.ones(1,bool)
    checkpoint(p,x,x,r,ids,mask,3.,dict(next_dt=1e-6),c,'identity')
    d,m,recovered=restore(p)
    np.testing.assert_array_equal(d['positions_m'],x)
    assert recovered.dual==c.dual and m['controller']['next_dt']==1e-6
    assert not p.with_suffix('.tmp').exists()
