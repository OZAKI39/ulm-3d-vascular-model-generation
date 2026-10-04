import numpy as np
from scripts import c7_early_branch as a


def test_every_accepted_state_passes_full_bond_kkt():
    rows=a.read(a.OUT/'branch/ALL.json')['states'];assert len(rows)>=5
    for r in rows:
        assert r['scientifically_accepted'] and r['full_KKT_bond_count']==46059
        assert r['metrics']['force_residual']<=1e-9 and r['metrics']['projected_phase_KKT']<=1e-8
        assert r['q_counts']['WEAK_LOWER']==r['weak_count']
        assert r['q_counts']['INTERIOR']==r['interior_count']
        data=a.arrays(a.OUT/r['label']/'STATE.npz');assert len(data['q_scaled'])==46059
        p=a.problem(7,r['strain'],data['lower']);k=data['s']-np.clip(data['s']-data['q_scaled'],p.lower,p.sc)
        assert abs(k).max()<=1e-8


def test_early_branch_irreversibility_and_predecessor_history():
    rows=a.read(a.OUT/'branch/ALL.json')['states'];previous=None
    for r in rows:
        data=a.arrays(a.OUT/r['label']/'STATE.npz');assert np.all(data['s']>=data['lower'])
        if previous is not None:assert np.array_equal(data['lower'],previous)
        previous=data['s']


def test_stop_rules_and_no_stability_claim():
    b=a.read(a.OUT/'branch/ALL.json');assert 5<=len(b['states'])<=20 and b['stability']=='NOT_ASSESSED'
    for r in b['states']:assert r['pool_size']<=.25*46059
    assert b['stop_reason']!='RUNNING'


def test_spatial_diagnostics_have_physical_units_and_all_bonds():
    s=a.read(a.OUT/'activation/SPATIAL.json');assert s['count']==len(s['bonds'])==len(s['IDs'])
    assert s['RMS_radius_m']>=0 and len(s['nearest_clamp_percentiles_m'])==5
    assert sum(s['orientation_counts'])==s['count']
    assert 0<s['equivalent_volume_m3']<=a.old.v2.L**3
    assert len(s['particle_aggregate'])==343


def test_production_and_historical_geometry_unchanged():
    pre=a.read(a.OUT/'PRE_RUN_MANIFEST.json')
    for name,h in pre['files'].items():
        if name.startswith('pd_clot/') or name.endswith('/geometry/C7.npz'):
            assert a.digest(a.ROOT/name)==h
