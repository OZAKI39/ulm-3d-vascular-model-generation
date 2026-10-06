"""Independent manufactured motion, QP and frozen-physics regressions T1–T15."""
import ast
import hashlib
import json
from pathlib import Path
import numpy as np
import pytest
from scipy import sparse
from scipy.optimize import minimize, LinearConstraint
from mevo_mb.bubbles.constrained_motion import (
    ConstrainedLubrication, horizon_constraints, PredictorCorrector,
    SparseActiveSetQP, AcceptedCache, SubdivideInterval,
)
from mevo_mb.bubbles.translation_backends import LubricationCoupledTranslation
R=1.78e-6;MU=.001;DRAG=6*np.pi*MU*R

class NoWall:
    def wall_gap_and_normal(self,x,r):return np.full(len(x),np.inf),np.zeros_like(x),np.ones(len(x),bool)

class Plane:
    def wall_gap_and_normal(self,x,r):return -x[:,2]-r,np.tile([0.,0.,1.],(len(x),1)),np.ones(len(x),bool)

def solve(x,flow,dt=.001,geometry=None,model=None,ids=None):
    m=model or ConstrainedLubrication(MU);r=np.full(len(x),R);wall=(geometry or NoWall()).wall_gap_and_normal(x,r)
    result=m.solve(x,r,np.zeros_like(x),flow,wall);ids=np.arange(len(x)) if ids is None else ids
    c=horizon_constraints(x,r,wall,result.velocity,dt,1e-11,ids)
    v,cache,d=m.constrained(result,c,ids)
    return v,cache,d,result,c

def pair(gap=0):return np.array([[0.,0.,0.],[(2+gap)*R,0.,0.]])

def test_T1_stokes_force_and_flow():
    m=ConstrainedLubrication(MU);x=np.zeros((1,3));f=np.array([[1e-11,-2e-11,3e-11]]);flow=np.array([[.001,.002,-.001]])
    r=np.array([R]);result=m.solve(x,r,f,flow,NoWall().wall_gap_and_normal(x,r))
    c=horizon_constraints(x,r,result.wall,result.velocity,.01,1e-11,[37]);u,_,d=m.constrained(result,c,[37])
    np.testing.assert_allclose(u,flow+f/DRAG,rtol=1e-14,atol=1e-17);assert d['constraint_count']==0

def test_T2_common_translation():
    flow=np.tile([.01,-.02,.003],(2,1));u,_,_,result,_=solve(pair(.001),flow)
    np.testing.assert_allclose(u,flow,rtol=1e-12,atol=1e-15)
    np.testing.assert_allclose(result.lubrication_forces(u)[0],0,atol=1e-20)

def test_T3_squeeze_large_dt():
    x=pair(.05);flow=np.array([[.01,0,0],[-.01,0,0]]);dt=.001
    u,_,d,result,_=solve(x,flow,dt)
    legacy=LubricationCoupledTranslation(MU).timestep(x,np.full(2,R),result.velocity,result.wall,dt)
    assert dt>10*legacy
    assert np.linalg.norm((x+dt*u)[1]-(x+dt*u)[0])-2*R>=-1e-17
    assert d['kkt']<=1e-9

def test_T4_separating_inactive():
    flow=np.array([[-.01,0,0],[.01,0,0]]);u,cache,_,result,_=solve(pair(0),flow)
    np.testing.assert_allclose(u,result.velocity,atol=1e-16);assert not cache.active

def test_T5_tangential_not_locked():
    flow=np.array([[0,.01,0],[0,-.01,0]]);u,*_=solve(pair(),flow)
    np.testing.assert_allclose(u,flow,rtol=1e-14,atol=1e-16)

def test_T6_wall_large_dt():
    x=np.array([[0.,0.,-1.05*R]]);u,_,d,_,_=solve(x,np.array([[.001,.002,.01]]),.01,Plane())
    assert (x+.01*u)[0,2]<=-R+1e-17
    np.testing.assert_allclose(u[0,:2],[.001,.002],rtol=1e-13);assert d['kkt']<=1e-9

@pytest.mark.parametrize('count',[3,30],ids=['T7_three','T8_thirty'])
def test_contact_chain(count):
    x=np.zeros((count,3));x[:,0]=2*R*np.arange(count)
    flow=np.zeros_like(x);flow[:,0]=np.linspace(.02,-.01,count);u,_,d,_,_=solve(x,flow)
    np.testing.assert_allclose(u[:,0],.005,atol=2e-11,rtol=2e-8);assert d['kkt']<=1e-9

@pytest.mark.parametrize('seed',range(6))
def test_T9_random_sparse_dense(seed):
    rng=np.random.default_rng(seed);n=9;m=11
    w=rng.normal(size=(n,n));a=w.T@w+np.eye(n);b=rng.normal(size=n);g=rng.normal(size=(m,n));lower=-rng.uniform(.1,1,m)
    u,lam,d=SparseActiveSetQP().solve(sparse.csc_matrix(a),b,sparse.csc_matrix(g),lower)
    # Independent primal SLSQP reference, only small-N test path.
    fit=minimize(lambda z:(.5*z@a@z-b@z,a@z-b),np.zeros(n),jac=True,method='SLSQP',
                 constraints=LinearConstraint(g,lower,np.inf),options=dict(ftol=1e-13,maxiter=1000))
    assert fit.success
    np.testing.assert_allclose(u,fit.x,rtol=2e-6,atol=2e-7)
    assert d['kkt']<=1e-9;assert np.min(g@u-lower)>=-1e-9
    assert abs((.5*u@a@u-b@u)-fit.fun)<1e-10
    assert np.max(abs(lam*(g@u-lower)))<1e-8


def test_T10_components_global():
    x=np.vstack([pair(),pair(.005)+[10*R,0,0],[[100*R,0,0]]]);flow=np.array([[.01,.003,0],[-.01,0,0],[.02,0,0],[-.01,0,0],[.001,0,0]])
    a=solve(x,flow,1e-5,model=ConstrainedLubrication(MU,components=True))
    b=solve(x,flow,1e-5,model=ConstrainedLubrication(MU,components=False))
    np.testing.assert_allclose(a[0],b[0],rtol=1e-9,atol=1e-12);assert a[2]['component_count']>=3


def test_T11_warm_cold_and_cache_ids():
    x=pair(.005);flow=np.array([[.01,.003,0],[-.01,0,0]]);m=ConstrainedLubrication(MU)
    first=solve(x,flow,model=m,ids=np.array([7,29]));m.cache=first[1]
    warm=solve(x,flow*1.01,model=m,ids=np.array([7,29]));cold=solve(x,flow*1.01,model=ConstrainedLubrication(MU,warm_start=False),ids=np.array([7,29]))
    np.testing.assert_allclose(warm[0],cold[0],rtol=1e-10,atol=1e-12)
    assert m.cache.particle_ids==(7,29);assert m.cache.sparsity_pattern
    assert m.cache is first[1]  # trial solve does not commit


def test_T12_predictor_corrector_refinement():
    def force(x,r,t):return dict(primary_force_N=np.zeros_like(x),secondary_force_N=np.zeros_like(x),flow_velocity_m_s=x)
    errors=[]
    for dt in [.1,.05,.025]:
        step=PredictorCorrector(ConstrainedLubrication(MU),NoWall(),force,position_atol_R=1e6,gap_atol_R=1e6)
        x=np.array([[R,0.,0.]]);r=np.array([R]);t=0.
        for _ in range(round(.4/dt)):
            x,_,cache,_=step.trial(x,r,[0],t,dt);step.backend.cache=cache;t+=dt
        errors.append(abs(x[0,0]-R*np.exp(.4)))
    assert errors[1]<errors[0]/3.5 and errors[2]<errors[1]/3.5


def test_T13_frozen_classes_and_physics_bytes():
    root=Path('src/mevo_mb/bubbles');expected=json.load(open('tests/data/stage5a_frozen_hashes.json'))
    text=(root/'translation_backends.py').read_text()
    for node in ast.parse(text).body:
        if isinstance(node,ast.ClassDef):assert hashlib.sha256(ast.get_source_segment(text,node).encode()).hexdigest()==expected[node.name]
    for name in ['forces.py','contact.py','secondary_ablation.py']:assert hashlib.sha256((root/name).read_bytes()).hexdigest()==expected[name]


def test_T14_secondary_disabled_force_path():
    import yaml
    from mevo_mb.bubbles.tracker import ParticleTracker
    class Field:
        frequency_Hz=1e6
        def interpolate_complex_pressure(self,x):return np.full(len(x),1000+0j)
        def interpolate_complex_gradient(self,x):return np.tile([1e5+0j,2e5+0j,0j],(len(x),1))
    class Flow:
        def velocity(self,x,t):return np.tile([.001,0.,0.],(len(x),1))
    values=[]
    for name in ['lubrication_coupled','constrained_lubrication_v1']:
        c=yaml.safe_load(Path('configs/demo.yaml').read_text());c['particles'].update(translation_backend=name,secondary_enabled=False)
        tracker=ParticleTracker(c,Field(),Flow())
        tracker.pair.forces=lambda *args:(_ for _ in ()).throw(AssertionError('Disabled secondary queried'))
        values.append(tracker.evaluate(pair(.05),np.full(2,R),0.))
    for key in values[1]:np.testing.assert_array_equal(values[0][key],values[1][key])


@pytest.mark.gpu
def test_T15_cpu_gpu_frozen_pair_force():
    cp=pytest.importorskip('cupy')
    from mevo_mb.bubbles.forces import ClassicalFarFieldBjerknes
    rng=np.random.default_rng(512);x=rng.normal(size=(186,3))*50*R;r=np.full(186,R);vol=(rng.normal(size=(186,4))+1j*rng.normal(size=(186,4)))*1e-18
    args=(x,r,vol,1e6,1000.)
    cpu=ClassicalFarFieldBjerknes('numpy',enforce_far_field=False).forces(*args)
    gpu=ClassicalFarFieldBjerknes('cupy',enforce_far_field=False).forces(*args).get()
    np.testing.assert_allclose(gpu,cpu,rtol=2e-12,atol=1e-26)


def test_rejection_does_not_commit_cache():
    def force(x,r,t):return dict(primary_force_N=np.zeros_like(x),secondary_force_N=np.zeros_like(x),flow_velocity_m_s=x*1000)
    m=ConstrainedLubrication(MU);prior=AcceptedCache((3,),{3:np.zeros(3)},{},(),'prior');m.cache=prior
    step=PredictorCorrector(m,NoWall(),force)
    with pytest.raises(SubdivideInterval):step.trial(np.array([[R,0.,0.]]),np.array([R]),[3],0,.1)
    assert m.cache is prior
