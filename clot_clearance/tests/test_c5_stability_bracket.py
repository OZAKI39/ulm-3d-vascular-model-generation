import hashlib
import numpy as np
from scripts import c5_stability_audit as a


def test_continuation_and_bracket_preserve_irreversibility():
    for p in (a.OUT/'continuation').glob('*/SUMMARY.json'):
        r=a.read(p);d=a.old.ca.read_npz(p.with_name('STATE.npz'));old=a.old.ca.read_npz(a.OUT/f'states/{r["previous_label"]}/STATE.npz')
        assert r['scientifically_accepted'] and r['strain']>r['previous_strain']
        np.testing.assert_array_equal(d['lower'],old['s']);assert np.all(d['s']>=old['s'])
        np.testing.assert_array_equal(d['z_initial'][225:],old['s'])


def test_no_bracket_without_verified_stable_and_unstable_endpoints():
    b=a.read(a.OUT/'STABILITY_BRACKET.json')
    if b['interval'] is None:
        assert not b['resolved']
        assert not any(r['status']=='STATIONARY_UNSTABLE' for r in a.read(a.OUT/'ALL_STATES.json')['states'])
    else:
        lo=a.read(a.OUT/f'states/{b["lower_label"]}/STABILITY.json');hi=a.read(a.OUT/f'states/{b["upper_label"]}/STABILITY.json')
        assert lo['status']=='STATIONARY_STABLE' and hi['status']=='STATIONARY_UNSTABLE'
        assert b['width']<=1e-6 or b['kind']=='NONSMOOTH_STABILITY_BRACKET'


def test_existing_states_all_checked_before_extension_and_no_overrun():
    r=a.read(a.OUT/'ALL_STATES.json');assert r['existing_roots_resolved_before_extension']
    assert len(a.read(a.OUT/'EXISTING_STATES.json')['states'])==9
    extended=[x for x in r['states'] if x['label'].startswith('E')]
    assert all(x['strain']<=.006 for x in extended)
    assert all(x['metrics']['maximum_s']<.1 for x in extended[:-1])


def test_saved_roots_production_and_existing_tests_unchanged():
    m=a.read(a.OUT/'PRE_RUN_MANIFEST.json')
    for name,digest in m['files'].items():
        if name.startswith(('pd_clot/','tests/','verification/pfpd_at1_petsc_vi_solver_v1/continuation/')):
            assert hashlib.sha256((a.ROOT/name).read_bytes()).hexdigest()==digest
