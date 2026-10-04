"""Independent candidate gates: boundaries, mechanics, variational force, author oracle."""
import json
from pathlib import Path
import numpy as np
import pytest
from scipy.spatial import cKDTree
from pd_clot.pfpd_candidate import Family,evaluate,Material

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'verification/pfpd_candidate_v1'


def family(n=4,ratio=2.01,unequal=False):
    X=(np.indices((n,n,n)).reshape(3,-1).T+.5)*.000125
    V=np.full(len(X),.000125**3)
    if unequal:V*=1+.2*np.sin(np.arange(len(X)))
    return Family.from_arrays(X,V,cKDTree(X).query_pairs(.000125*ratio,output_type='ndarray'),.000125*ratio)


def deformed(f):
    X=f.X;q=X/.000125
    return X@np.array([[1.02,.01,0],[0,.99,.004],[0,0,1.01]]).T+2e-7*np.sin(q[:,[1,2,0]])


def test_omega0_endpoint_sum_positive_and_boundary():
    f=family(9);i,j=f.pairs.T
    direct=np.zeros(len(f.X))
    for a,b,w in zip(i,j,f.omega):direct[a]+=w*f.volume[b];direct[b]+=w*f.volume[a]
    np.testing.assert_allclose(f.omega0,direct,rtol=2e-15)
    assert np.min(f.omega0)>0 and np.ptp(f.omega0)>0


def test_orientation_single_physical_identity_and_wbar_symmetry():
    f=family(unequal=True);b=f.oriented_bond
    assert np.array_equal(np.bincount(b),np.full(len(f.pairs),2))
    assert np.array_equal(f.source[:len(f.pairs)],f.target[len(f.pairs):])
    np.testing.assert_array_equal(f.wbar[b[:len(f.pairs)]],f.wbar[b[len(f.pairs):]])
    assert np.all(f.sign[:len(f.pairs)]==-f.sign[len(f.pairs):])


def test_translation_zero_energy_force():
    f=family();s=evaluate(f,f.X+[.0007,-.0005,.001])
    np.testing.assert_allclose(s.Fbar,np.broadcast_to(np.eye(3),s.Fbar.shape),atol=3e-14)
    assert s.energy<1e-25 and np.max(abs(s.force))<1e-15


@pytest.mark.parametrize('model',['svk','neo_hookean'])
def test_rotation_objectivity_and_force_covariance(model):
    f=family();x=deformed(f);a=.51;R=np.array([[1,0,0],[0,np.cos(a),-np.sin(a)],[0,np.sin(a),np.cos(a)]])
    s=evaluate(f,x,Material(model));r=evaluate(f,x@R.T,Material(model))
    assert abs(s.energy-r.energy)/s.energy<1e-10
    np.testing.assert_allclose(r.force,s.force@R.T,rtol=1e-9,atol=1e-15)


@pytest.mark.parametrize('F',[np.diag([1.1,1,1]),np.diag([.9,1,1]),np.eye(3)*1.02,np.array([[1,.03,0],[0,1,0],[0,0,1]])])
def test_affine_bulk_and_truncated_boundary(F):
    f=family(9);s=evaluate(f,f.X@F.T)
    interior=np.all((f.X-f.X.min(0)>f.horizon)&(f.X.max(0)-f.X>f.horizon),axis=1)
    assert interior.any() and (~interior).any()
    for mask in [interior,~interior]:assert np.max(abs(s.Fbar[mask]-F))<5e-12
    assert np.max(abs(s.Ftilde-F))<5e-12


def test_bond_direction_exact_for_nonaffine_deformation():
    f=family();x=deformed(f);s=evaluate(f,x);i,j=f.pairs.T
    np.testing.assert_allclose(np.einsum('bij,bj->bi',s.Ftilde,f.xi),x[j]-x[i],rtol=1e-13,atol=1e-18)


def test_force_and_moment_with_unequal_volumes():
    f=family(unequal=True);x=deformed(f);s=evaluate(f,x)
    scale=np.linalg.norm(s.force,axis=1).sum()
    assert np.linalg.norm(s.force.sum(0))/scale<5e-12
    assert np.linalg.norm(np.cross(x-x.mean(0),s.force).sum(0))/(scale*.0005)<5e-12


@pytest.mark.parametrize('model',['svk','neo_hookean'])
def test_multi_epsilon_energy_force_with_unequal_volumes(model):
    f=family(3,unequal=True);x=deformed(f);m=Material(model);s=evaluate(f,x,m)
    direction=np.random.default_rng(581).normal(size=x.shape);direction/=np.linalg.norm(direction)
    exact=-np.sum(s.force*direction);errors=[]
    for e in [1e-8,3e-9,1e-9,3e-10]:
        fd=(evaluate(f,x+e*direction,m).energy-evaluate(f,x-e*direction,m).energy)/(2*e)
        errors.append(abs(fd-exact)/abs(exact))
    assert max(errors)<2e-6


def test_endpoint_energy_volume_counting():
    f=family();x=deformed(f);s=evaluate(f,x)
    scaled=Family.from_arrays(f.X,2*f.volume,f.pairs,f.horizon,f.omega);r=evaluate(scaled,x)
    np.testing.assert_allclose(r.force,2*s.force,rtol=3e-14,atol=1e-18)
    assert abs(r.energy-2*s.energy)<1e-24
    assert abs(s.energy-s.point_energy.sum())<1e-24


def test_inputs_read_only_and_unchanged():
    f=family();x=deformed(f);g=np.linspace(.2,1,len(f.pairs));h=np.ones_like(g)
    originals=[a.copy() for a in [x,g,h]]
    for a in [x,g,h]:a.setflags(write=False)
    evaluate(f,x,g=g,h=h)
    for a,b in zip([x,g,h],originals):np.testing.assert_array_equal(a,b)
    assert all(not getattr(f,k).flags.writeable for k in ['X','volume','pairs','omega','wbar','source'])


def test_rank_deficiency_is_explicitly_out_of_scope():
    f=family()
    with pytest.raises(ValueError,match='Rank deficient'):evaluate(f,f.X,h=np.zeros(len(f.pairs)))


@pytest.mark.parametrize('name',['affine','nonaffine'])
def test_against_executed_official_julia_oracle(name):
    co=OUT/'comparisons';X=np.loadtxt(co/'input_X.csv',delimiter=',');V=np.loadtxt(co/'input_volume.csv',delimiter=',');delta=np.loadtxt(co/'input_parameters.csv',delimiter=',')[0]
    f=Family.from_arrays(X,V,cKDTree(X).query_pairs(delta,output_type='ndarray'),delta)
    s=evaluate(f,np.loadtxt(co/f'input_x_{name}.csv',delimiter=','));ids=np.loadtxt(co/f'julia_{name}_ids.csv',delimiter=',',dtype=int)-1
    force=np.loadtxt(co/f'julia_{name}_force.csv',delimiter=',');E=float(np.loadtxt(co/f'julia_{name}_energy.csv',delimiter=','))
    assert np.linalg.norm(s.force[ids]-force)/np.linalg.norm(force)<1e-9
    assert abs(s.energy-E)/E<1e-9
