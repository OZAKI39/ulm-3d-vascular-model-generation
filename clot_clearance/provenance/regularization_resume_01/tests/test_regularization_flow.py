import json
from pathlib import Path
import numpy as np
from pd_clot.regularization.flow import ManufacturedStreamingProvider
from pd_clot.streaming import fluid_traction


def fixture():
    c=json.loads((Path(__file__).parents[1]/'configs/streaming_regularized_demo_coarse.json').read_text())
    return c,ManufacturedStreamingProvider(c)


def test_manufactured_analytic_gradient_divergence_and_time_derivative():
    c,p=fixture();rng=np.random.default_rng(250);x=p.center+rng.normal(size=(40,3))*p.L;t=.000173
    actual=p.sample(x,t);numerical=np.zeros_like(actual.gradient);eps=1e-8
    for j in range(3):
        dx=np.eye(3)[j]*eps
        numerical[:,:,j]=(p.sample(x+dx,t).velocity-p.sample(x-dx,t).velocity)/(2*eps)
    np.testing.assert_allclose(actual.gradient,numerical,atol=3e-7,rtol=2e-7)
    assert np.abs(np.trace(actual.gradient,axis1=1,axis2=2)).max()<1e-10
    dt=1e-8
    np.testing.assert_allclose(actual.du_dt,(p.sample(x,t+dt).velocity-p.sample(x,t-dt).velocity)/(2*dt),atol=1e-7,rtol=1e-7)


def test_manufactured_spatial_decay_and_recirculation():
    _,p=fixture();near=p.center+np.array([[p.L,0,0],[-p.L,0,0],[0,0,p.L],[0,0,-p.L]])
    u=p.localized(near,0).velocity
    assert u[0,2]>0 and u[1,2]<0 and u[2,0]<0 and u[3,0]>0
    far=p.center+np.array([[8*p.L,0,0]])
    assert np.linalg.norm(p.localized(far,0).velocity)<1e-10*np.linalg.norm(u)


def test_traction_sign_and_consistency():
    c,p=fixture();x=p.center[None,:];normal=np.array([[1.,0,0]])
    sample=p.sample(x,0)
    expected=fluid_traction(sample.gradient,sample.pressure,normal,p.mu)
    np.testing.assert_allclose(p.traction(x,normal,0),expected,rtol=1e-14)
    pure=fluid_traction(np.zeros((1,3,3)),np.array([2.]),normal,p.mu)
    np.testing.assert_array_equal(pure,[[-2.,0,0]])
    assert p.center[0]<c['clot']['origin_m'][0]
    assert np.isclose(p.center[2],c['clot']['origin_m'][2]+.5*c['clot']['cells'][2]*c['clot']['particle_spacing_m'])
