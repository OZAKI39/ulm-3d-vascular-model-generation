"""Preserve the difference between a converged optimizer and an accepted state."""
from scripts.spectral_split_audit import OUT,read

def test_C3_spectral_stability_regression():
    r=read(OUT/'comparisons/C3_STABILITY_REGRESSION.json')
    assert r['passed'] and r['interval']==[.0099921875,.01]
    endpoint=read(OUT/'critical_brackets/C3_A.json')['first_unstable']
    assert endpoint['stationary'] and endpoint['negative_witness']['verified_negative']

def test_xtol_never_overrides_force_or_phase_KKT_gate():
    r=read(OUT/'branch_tracking/C5_W_005/STATIONARY.json')
    assert r['solver']['status']==3 and r['gates']['least_squares_converged']
    assert not r['accepted'] and not r['gates']['force'] and not r['gates']['phase_KKT']
    assert r['metrics']['force_residual']>1e-9 and r['metrics']['projected_phase_KKT']>1e-8
    # A bound-feasible AT1 predictor improves the final residual but cannot
    # turn optimizer termination into physical acceptance.
    p=read(OUT/'branch_tracking/C5_P_005/STATIONARY.json')
    assert p['solver']['status']==4 and p['gates']['least_squares_converged']
    assert not p['accepted'] and not p['gates']['force'] and not p['gates']['phase_KKT']
    assert p['gates']['bounds'] and p['gates']['irreversibility']

def test_C7_heavy_work_stays_stopped_after_C5_activation_failure():
    gate=read(OUT/'GATES.json');assert gate['C7_heavy']=='NOT_RUN'
    assert not list((OUT/'branch_tracking').glob('C7_*/STATE.npz'))
    r=read(OUT/'critical_brackets/C5_A.json');assert r['interval'] is None and not r['complete']
