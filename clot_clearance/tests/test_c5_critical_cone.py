import numpy as np
from scripts import c5_stability_audit as a


def test_strong_lower_excluded_from_reduced_phase():
    r=a.read(a.OUT/'states/S0_LAST/STABILITY.json');assert r['m']==0 and r['counts']['STRONG_LOWER']==5392
    for r in a.read(a.OUT/'ALL_STATES.json')['states']:
        assert r['m']==sum(r['counts'][k] for k in ['WEAK_LOWER','INTERIOR','WEAK_UPPER'])


def test_weak_lower_cone_nonnegative():
    groups={'WEAK_LOWER':np.array([True,False]),'WEAK_UPPER':np.zeros(2,bool)}
    assert a.cone_feasible(np.array([.1,-3.]),groups,np.arange(2))
    assert not a.cone_feasible(np.array([-.1,3.]),groups,np.arange(2))


def test_interior_phase_unrestricted_and_upper_nonpositive():
    groups={'WEAK_LOWER':np.zeros(2,bool),'WEAK_UPPER':np.array([False,True])}
    assert a.cone_feasible(np.array([-4.,-.1]),groups,np.arange(2))
    assert a.cone_feasible(np.array([4.,-.1]),groups,np.arange(2))
    assert not a.cone_feasible(np.array([4.,.1]),groups,np.arange(2))


def test_negative_cone_witness_is_feasible():
    groups={'WEAK_LOWER':np.ones(2,bool),'WEAK_UPPER':np.zeros(2,bool)}
    witness,record=a.select_cone(np.diag([-2.,1.]),groups,np.arange(2),np.ones(2))
    assert witness[0]<0 and a.cone_feasible(witness[2],groups,np.arange(2))


def test_every_final_classification_has_direct_energy_sign_confirmation():
    for r in a.read(a.OUT/'ALL_STATES.json')['states']:
        if r['status'] in ['STATIONARY_STABLE','STATIONARY_UNSTABLE']:
            d=a.read(a.OUT/f'states/{r["label"]}/MODE_ENERGY.json')
            assert d['cone_feasible'] and d['passed'] and all(x['sign_consistent'] for x in d['rows'])
            expected=-1 if r['status']=='STATIONARY_UNSTABLE' else 1
            assert all(expected*c>0 for row in d['rows'] for c in row['selected_curvatures'])
