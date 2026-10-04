import numpy as np
from scripts import c7_early_branch as a


def test_saved_geometry_hash_and_quadrature():
    r=a.read(a.OUT/'C7_GEOMETRY_CHECK.json');assert r['passed']
    assert r['geometry_SHA256']==a.digest(a.old.v2.OUT/'geometry/C7.npz')
    assert (r['particles'],r['bonds'],r['mechanical_DOF'])==(343,46059,735)


def test_q_classification_at_declared_tolerance():
    tol=a.cfg()['q_tol'];q=np.array([2*tol,tol,0,-tol,-2*tol]);c=a.classify(q)
    assert c['STRONG_LOWER'].tolist()==[True,False,False,False,False]
    assert c['WEAK_LOWER'].tolist()==[False,True,True,True,False]
    assert c['VIOLATED_LOWER'].tolist()==[False,False,False,False,True]


def test_activation_uses_equilibrated_full_bond_states():
    b=a.read(a.OUT/'activation/BRACKET.json')
    for row in [b['lower'],b['upper']]:
        assert row['intact_accepted'] and row['gates']['force'] and row['phase_bonds']==46059
        assert row['metrics']['force_residual']<=1e-9
        state=a.arrays(a.OUT/f'intact/{row["label"]}/STATE.npz')
        assert len(state['q'])==46059 and np.count_nonzero(state['s'])==0
    assert b['lower']['min_q']>=0 and b['upper']['min_q']<0


def test_bisection_contracts_and_meets_resolution():
    b=a.read(a.OUT/'activation/BRACKET.json');assert b['width']<=1e-7
    for r in b['history']:
        assert r['before'][0]<=r['after'][0]<r['after'][1]<=r['before'][1]
        assert np.isclose(np.ptp(r['after']),np.ptp(r['before'])/2,atol=1e-16,rtol=0)


def test_c7_intact_resume_checks_history_and_reuses_exact_state(monkeypatch):
    from scripts import c7_reduced_vi as reduced
    b=a.read(a.OUT/'activation/BRACKET.json');upper=b['upper'];lower=b['lower']
    previous=a.arrays(a.OUT/f'intact/{lower["label"]}/STATE.npz')
    # The last midpoint used this lower endpoint as its mechanical history.
    assert b['history'][-1]['label']==upper['label']
    def forbidden(*args,**kwargs):raise AssertionError('Cached C7 point was re-solved')
    monkeypatch.setattr(reduced,'solve_once',forbidden)
    row,data=a.intact(upper['strain'],upper['label'],previous,None)
    expected=a.arrays(a.OUT/f'intact/{upper["label"]}/STATE.npz')
    assert row['checkpoint_key']==upper['checkpoint_key']
    assert np.array_equal(data['x'],expected['x']) and np.array_equal(data['q'],expected['q'])
