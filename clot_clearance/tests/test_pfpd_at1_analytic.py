import numpy as np
import pytest
from pd_clot.pfpd_candidate.at1_evolution import (
    at1_trial_state, at1_constrained_update, critical_driver, audit_update, HighDamageGate)
from pd_clot.pfpd_candidate.normalization import c0


@pytest.mark.parametrize('ratio, expected', [(0.,0.),(.01,0.),(.9,0.),(1.,0.),
                                              (1.0001,1-1/1.0001),(1.5,1/3),(2.,.5),(100.,.99)])
def test_virgin_threshold_and_stationary_minimum(ratio, expected):
    update = at1_constrained_update(0., ratio*3., 3.)
    assert float(update.s_new) == pytest.approx(expected, abs=2e-16)
    assert audit_update(update, ratio*3., 3.)['passed']


def test_vector_broadcast_bounds_and_no_input_mutation():
    old = np.array([0., .5, .9]); y = np.array([[0.], [2.], [100.]])
    old.setflags(write=False); y.setflags(write=False)
    old_copy, y_copy = old.copy(), y.copy()
    update = at1_constrained_update(old, y, 1.)
    assert update.s_new.shape == (3,3)
    assert np.all((update.s_new >= old) & (update.s_new <= 1))
    assert np.array_equal(old, old_copy) and np.array_equal(y, y_copy)
    assert not np.shares_memory(update.s_new, old)


@pytest.mark.parametrize('old,y,yc', [(np.nan,1,1),(0,np.nan,1),(0,1,np.inf),
                                    (-.1,1,1),(1.1,1,1),(0,-1,1),(0,1,0),(0,1,-1)])
def test_reject_invalid_state(old,y,yc):
    with pytest.raises(ValueError): at1_constrained_update(old,y,yc)


def test_lower_bound_kkt_and_total_energy():
    update = at1_constrained_update([.9,.5,0.], [1.,1.,2.], 1.)
    assert list(update.active_set_status) == ['LOWER_BOUND','LOWER_BOUND','INTERIOR']
    assert np.all(update.derivative_at_solution[:2] >= 0)
    assert abs(update.derivative_at_solution[2]) < 1e-15
    assert np.all(update.total_after <= update.total_before)


def test_singleton_upper_bound_is_valid_only_in_unrestricted_analytic_case():
    update = at1_constrained_update(1., 0., 1.)
    assert update.s_new == 1 and str(update.active_set_status) == 'FIXED_SINGLETON'
    assert audit_update(update,0,1)['passed']
    with pytest.raises(HighDamageGate): at1_constrained_update(1.,0.,1.,.95)


def test_high_damage_gate_preserves_unclipped_proposal_and_old_state():
    old = np.array([0.,.8]); before = old.copy()
    with pytest.raises(HighDamageGate) as exc:
        at1_constrained_update(old, [100.,3.], 1., .95)
    assert np.array_equal(old,before)
    assert exc.value.update.s_new[0] == .99
    assert exc.value.update.s_old[0] == 0
    assert 'HIGH_DAMAGE_KINEMATIC_COUPLING_NOT_IMPLEMENTED' in str(exc.value)


def test_exact_sc_reaches_gate():
    with pytest.raises(HighDamageGate): at1_constrained_update(.95,1.,1.,.95)


def test_critical_driver_uses_existing_kernel_normalization():
    gc, delta = .013,.007
    for kernel in ['constant','linear','quadratic']:
        assert critical_driver(gc,delta,kernel) == gc/(2*c0(kernel)['c0']*delta)


def test_eq325_is_distinguished_by_threshold_and_supercritical_solution():
    assert at1_trial_state(.5,1.) == 0.
    assert at1_trial_state(2.,1.) == .5


def test_unloading_and_reload_beyond_historical_peak():
    phase = []; s = 0.
    for y in [0.,1.,2.,3.,2.,0.,1.,2.5,4.]:
        s = float(at1_constrained_update(s,y,1.).s_new); phase.append(s)
    assert np.all(np.diff(phase) >= 0)
    assert np.all(np.array(phase[3:8]) == phase[3])
    assert phase[-1] > phase[3]


def test_energy_derivative_against_finite_difference():
    y,yc,s = 7.,2.,.3
    h = 1e-6
    energy = lambda a:(1-a)**2*y+2*yc*a
    assert (energy(s+h)-energy(s-h))/(2*h) == pytest.approx(-2*(1-s)*y+2*yc,rel=1e-10)
