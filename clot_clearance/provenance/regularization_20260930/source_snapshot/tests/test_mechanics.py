import json
from pathlib import Path
import numpy as np
import pytest
from pd_clot.geometry import make_cloud
from pd_clot.mechanics import evaluate_cloud, prepare_shape
from pd_clot.damage import instantaneous, accumulate, fragments

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def model():
    c = json.loads((ROOT/'configs/straight_pipe.json').read_text())
    c['clot']['cells'] = [5, 4, 4]
    g = make_cloud(c['clot'])
    return c, g, np.ones(len(g.pairs))


def test_affine_and_rigid_objectivity(model):
    c,g,b = model
    A=np.array([[1.03,.06,0],[0,.99,.02],[0,0,1.01]])
    f,F,J,e,hg=evaluate_cloud(g,g.X@A.T,b,c['material'],c['safety'])
    np.testing.assert_allclose(F,np.broadcast_to(A,F.shape),atol=3e-14)
    assert hg.max()<1e-20
    theta=.4;R=np.array([[np.cos(theta),-np.sin(theta),0],[np.sin(theta),np.cos(theta),0],[0,0,1]])
    fr,*_=evaluate_cloud(g,(g.X@A.T)@R.T+0.001,b,c['material'],c['safety'])
    np.testing.assert_allclose(fr,f@R.T,atol=1e-15,rtol=1e-9)
    assert np.linalg.norm(f.sum(axis=0))<1e-15
    assert np.linalg.norm(np.cross(g.X@A.T,f).sum(axis=0))<1e-17


def test_energy_gradient_and_stabilization(model):
    c,g,b=model
    rng=np.random.default_rng(100)
    x=g.X+rng.normal(size=g.X.shape)*1e-7
    b=rng.uniform(.65,1,len(b))
    f,F,J,e,hg=evaluate_cloud(g,x,b,c['material'],c['safety'])
    assert hg.sum()>0
    d=rng.normal(size=x.shape);d/=np.linalg.norm(d)
    eps=1e-9
    def energy(y):
        z=evaluate_cloud(g,y,b,c['material'],c['safety'])
        return np.dot(g.volume,z[3]+z[4])
    fd=(energy(x+eps*d)-energy(x-eps*d))/(2*eps)
    np.testing.assert_allclose(fd,-np.sum(f*d),rtol=2e-5,atol=1e-11)


def test_no_load_and_single_fragment(model):
    c,g,b=model
    f,F,J,e,hg=evaluate_cloud(g,g.X,b,c['material'],c['safety'])
    assert np.max(np.abs(f))<1e-14
    assert fragments(g,b,0)[2]['number_of_fragments']==1


def test_shape_and_jacobian_guards(model):
    c,g,b=model
    x=g.X.copy();x[:,2]*=-1
    with pytest.raises(ValueError,match='Jacobian'):
        evaluate_cloud(g,x,b,c['material'],c['safety'])
    b[:]=0;b[0]=1
    with pytest.raises(FloatingPointError,match='shape support'):
        prepare_shape(g,b,c['safety'])


def test_instantaneous_irreversibility():
    b=instantaneous(np.ones(3),np.array([1.,1.2,1.4]),1.1,1.3)
    np.testing.assert_allclose(b,[1,.5,0])
    np.testing.assert_array_equal(instantaneous(b,np.ones(3),1.1,1.3),b)


def test_cyclic_rate_and_cycle_monotonicity(model):
    c,_,b=model;p=c['damage'].copy();q=np.full(len(b),.05)
    a=accumulate(b,q,p);p['DeltaN']*=2;z=accumulate(b,q,p)
    assert np.all(z<=a)
    p['C_damage']=0
    np.testing.assert_array_equal(accumulate(b,q,p),b)


def test_artificial_fragmentation_and_supported_shapes(model):
    c,g,b=model
    x=g.X.copy();split=g.X[:,2].mean();x[g.X[:,2]>split,2]+=5*g.spacing
    stretch=np.linalg.norm(x[g.pairs[:,1]]-x[g.pairs[:,0]],axis=1)/g.length
    broken=instantaneous(b,stretch,1.01,1.02)
    lab,_,info=fragments(g,broken,0)
    assert info['number_of_fragments']==2
    assert sorted(info['fragment_sizes'])==[40,40]
    prepare_shape(g,broken,c['safety'])
