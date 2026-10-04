import numpy as np
from scipy.linalg import cho_factor,cho_solve,eigh
from scripts import c5_stability_audit as a


def test_schur_equals_relaxed_coupled_curvature():
    for r in a.read(a.OUT/'ALL_STATES.json')['states']:
        for t in r['schur_tests']:
            assert t['algebraic_relative_error']<1e-10
            assert t['direct']['passed']


def test_condensation_solve_residual():
    for r in a.read(a.OUT/'ALL_STATES.json')['states']:
        assert all(t.get('condensation_residual',0)<1e-10 for t in r['eigenvalues'])


def test_positive_schur_certifies_all_cone_subsets():
    S=np.array([[2.,.4],[.4,1.]])
    groups={key:np.zeros(2,bool) for key in ['WEAK_LOWER','WEAK_UPPER']};groups['WEAK_LOWER'][:]=True
    witness,record=a.select_cone(S,groups,np.arange(2),np.ones(2))
    assert record['method']=='POSITIVE_UNCONSTRAINED_CERTIFICATE' and not record['optimizer_run']
    assert eigh(S,eigvals_only=True)[0]>0 and witness[0]>0


def test_no_negative_witness_is_not_a_stable_certificate():
    S=np.array([[1.,2.],[2.,1.]])
    groups={'WEAK_LOWER':np.ones(2,bool),'WEAK_UPPER':np.zeros(2,bool)}
    witness,record=a.select_cone(S,groups,np.arange(2),np.ones(2))
    assert record['lambda_min']<0 and witness[0]>0
    assert record['method']=='CONSTRAINED_WITNESS_SEARCH'
