import hashlib
from pathlib import Path
import numpy as np
from pd_clot.pfpd_candidate.coupling_stability_audit import read_json,read_npz
from pd_clot.pfpd_candidate.branch_stability_tracking import OUT,settings,align_modes


def test_both_full_predeclared_grids_are_saved_without_dropping_points():
    for k,folder in [('A','coarse_tracking'),('B','fine_tracking')]:
        r=read_json(OUT/folder/'TRACKING.json')
        assert r['status']=='TRACKED'
        np.testing.assert_array_equal([s['strain'] for s in r['states']],settings()['schedules'][k])


def test_no_continuation_state_heals_or_uses_a_future_load_parent():
    for f in (OUT/'branch_states').glob('*/STATIONARITY.json'):
        r=read_json(f);d=read_npz(f.with_name('STATE.npz'));assert np.all(d['s']>=d['lower'])
        if r['initial_kind']=='continuation_predictor':
            parent=read_json(OUT/'branch_states'/r['previous_accepted_state']/'STATIONARITY.json')
            assert parent['strain']<r['strain'] and parent['accepted']


def test_active_set_mass_switch_records_all_bond_IDs():
    row=read_json(OUT/'active_sets/B_grid_04_EVENT.json')
    assert row['status']=='BRANCH_CONTINUITY_EVENT'
    assert len(row['switched_bond_IDs'])==185 and row['switched_bonds']==185
    assert row['active_set_changes']['LOWER_WEAK_NEAR_SWITCH']['entered']==list(range(185))


def test_eigenvector_sign_alignment_and_missing_phase_component():
    a=np.array([1.,2.,3.,4.]);v,c=align_modes(a,-a,2)
    np.testing.assert_array_equal(v,a);assert all(abs(x-1)<1e-14 for x in c.values())
    _,c=align_modes(np.array([1.,2.,0.,0.]),a,2);assert c['phase'] is None


def test_production_and_all_existing_tests_unchanged():
    m=read_json(OUT/'provenance/PRE_RUN_MANIFEST.json');root=Path(m['repo_root']);checked=0
    for name,r in m['files'].items():
        if name.startswith(('clot_clearance/pd_clot/','clot_clearance/tests/')) and name.endswith('.py'):
            assert hashlib.sha256((root/name).read_bytes()).hexdigest()==r['sha256'],name;checked+=1
    assert checked>20


def test_historical_results_and_previous_spectrum_status_unchanged():
    m=read_json(OUT/'provenance/PRE_RUN_MANIFEST.json');root=Path(m['repo_root']);checked=0
    for name,r in m['files'].items():
        if name.startswith(('clot_clearance/results/pfpd_at1_evolution_v1/','clot_clearance/verification/pfpd_at1_coupling_stability_audit_v1/')) and not ('runtime_cache/' in name or '_tmp/' in name):
            assert hashlib.sha256((root/name).read_bytes()).hexdigest()==r['sha256'],name;checked+=1
    assert checked>100
    prior=read_json(OUT.parent/'pfpd_at1_coupling_stability_audit_v1/STATUS.json')
    assert prior['FIXED_POINT_SPECTRAL_RADIUS']=='UNRESOLVED'
