"""Targeted spectral constitutive checks; no repetition of old BAQP campaigns."""
import numpy as np
import pytest
from scripts.spectral_split_audit import (tensor_cases,F_from_E,svk_split_response,svk_split_energy,
    svk_split_second_pk,degraded_svk_split_response,Material,constitutive,OUT,read,SpectralProblem,baqp,v2)

@pytest.mark.parametrize('name',list(tensor_cases()))
@pytest.mark.parametrize('bulk',[666.6666666666666,2300.])
def test_energy_and_stress_decomposition(name,bulk):
    mat=Material('svk',1000.,bulk);E=tensor_cases()[name];F=F_from_E(E);r=svk_split_response(F,mat);W,P=constitutive(F,mat)
    np.testing.assert_allclose(r.Psi_total,W,rtol=2e-12,atol=1e-25)
    np.testing.assert_allclose(r.S_total,mat.lame*np.trace(r.E)*np.eye(3)+2*mat.mu*r.E,rtol=2e-12,atol=1e-12)
    np.testing.assert_allclose(r.P_total,P,rtol=2e-12,atol=1e-12)

def test_three_epsilon_energy_derivative():
    r=read(OUT/'constitutive/VERIFICATION.json')
    assert {d['epsilon'] for d in r['derivatives']}=={1e-5,3e-6,1e-6}
    assert r['maxima']['energy_derivative']<5e-6

def test_128_random_intact_F_states():
    rows=read(OUT/'constitutive/VERIFICATION.json')['intact_random_states']
    assert len(rows)>=100
    assert max(max(r[k] for k in ['energy','stress','PK']) for r in rows)<2e-12

def test_pure_compression_has_zero_active_energy_and_no_phase_increase():
    r=read(OUT/'constitutive/BAQP_SMALL_FIXTURES.json')['fixtures'][0]
    assert r['driver_max']<1e-20 and r['maximum_phase_increase']==0

def test_mixed_state_contains_both_energies():
    r=svk_split_response(F_from_E(tensor_cases()['mixed']),Material('svk',1000.,666.6666666666666))
    assert r.Psi_plus>0 and r.Psi_minus>0

def test_shear_fixture_has_both_energies():
    r=read(OUT/'constitutive/BAQP_SMALL_FIXTURES.json')['fixtures'][1]
    assert r['Psi_plus_min']>0 and r['Psi_minus_min']>0 and r['mixed_state_admitted']

def test_compression_energy_and_stress_are_not_degraded():
    mat=Material('svk',1000.,2300.);F=np.eye(3)*.98
    a=degraded_svk_split_response(F,0,mat);b=degraded_svk_split_response(F,.6,mat)
    np.testing.assert_allclose(a.Psi,b.Psi,rtol=1e-15)
    np.testing.assert_allclose(a.P,b.P,rtol=1e-15)

def test_repeated_eigenspaces_are_basis_independent():
    mat=Material('svk',1000.,2300.);rng=np.random.default_rng(481);Q,_=np.linalg.qr(rng.normal(size=(3,3)))
    E=tensor_cases()['repeated_eigenvalues'];a,b=svk_split_energy(E,mat.lame,mat.mu);S,T=svk_split_second_pk(E,mat.lame,mat.mu)
    c,d=svk_split_energy(Q@E@Q.T,mat.lame,mat.mu);U,V=svk_split_second_pk(Q@E@Q.T,mat.lame,mat.mu)
    np.testing.assert_allclose([c,d],[a,b],rtol=2e-12)
    np.testing.assert_allclose(U,Q@S@Q.T,rtol=2e-12,atol=1e-12)

def test_degraded_PK_is_energy_gradient_wrt_F():
    mat=Material('svk',1000.,2300.);F=np.array([[1.03,.015,0],[0,.98,.01],[0,0,1.02]])
    H=np.array([[.3,-.2,.1],[.1,.4,-.3],[.2,0,-.2]]);H/=np.linalg.norm(H)
    a=degraded_svk_split_response(F,.3,mat);exact=np.sum(a.P*H)
    for e in [1e-5,3e-6,1e-6]:
        fd=(degraded_svk_split_response(F+e*H,.3,mat).Psi-degraded_svk_split_response(F-e*H,.3,mat).Psi)/(2*e)
        assert abs(fd-exact)/max(abs(exact),np.linalg.norm(a.P))<5e-6

def test_new_BAQP_callback_restores_existing_constitutive():
    original=baqp.constitutive;p=SpectralProblem(3,0);p.state(p.family.X*.99,p.lower)
    assert baqp.constitutive is original

def test_split_phase_gradient_uses_Psi_plus():
    p=SpectralProblem(3,0);x=p.family.X.copy();x[p.free,0]*=.995;s=np.full(p.ns,.03)
    state=p.state(x,s);z=p.pack(x,s);_,g=p.objective_gradient(z)
    np.testing.assert_allclose(g[p.nx:]*p.energy_scale/p.weights,-2*(1-s)*state.split_response.split.Psi_plus+2*p.yc,rtol=2e-12)

def test_C3_pure_tension_limit_reproduced():
    rows=read(OUT/'constitutive/BAQP_SMALL_FIXTURES.json')['C3_pure_tension_limit']
    for r in rows:
        assert r['energy_relative_error']<2e-12 and r['reaction_relative_error']<2e-12
        assert r['driver_relative_error']<2e-12 and r['Psi_minus_max']<1e-20

def test_nonpositive_J_is_still_rejected():
    with pytest.raises(ValueError,match='Positive determinant'):svk_split_response(np.diag([-1.,1.,1.]),Material())
