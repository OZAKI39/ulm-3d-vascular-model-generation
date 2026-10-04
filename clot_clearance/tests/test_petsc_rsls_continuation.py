import numpy as np
from scripts import petsc_vi_audit as audit


def test_continuation_uses_only_previous_petsc_state():
    x,s,eps=audit.initial()
    for r in audit.targets():
        a=audit.prev.ca.read_npz(audit.OUT/f'continuation/PETSC_RSLS/{r["label"]}/STATE.npz')
        p=audit.Problem(5,r['strain'],s)
        expected=p.pack(audit.prev.tracker.displacement_predictor(p,x,eps),s)
        np.testing.assert_array_equal(a['z_initial'],expected)
        np.testing.assert_array_equal(a['lower'],s)
        assert np.all(a['s']>=s)
        x,s,eps=a['x'],a['s'],r['strain']


def test_original_load_is_in_chronological_order():
    values=[t['strain'] for t in audit.targets()]
    assert values==sorted(values) and len(values)==7


def test_ssls_initial_conditions_match_rsls():
    for label in ['first_A1','same_load_as_original']:
        a=audit.prev.ca.read_npz(audit.OUT/f'continuation/PETSC_RSLS/{label}/STATE.npz')
        b=audit.prev.ca.read_npz(audit.OUT/f'ssls/{label}/STATE.npz')
        np.testing.assert_array_equal(a['z_initial'],b['z_initial'])
        np.testing.assert_array_equal(a['lower'],b['lower'])
