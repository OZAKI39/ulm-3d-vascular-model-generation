"""Source geometry and conservation gates, before a continuum domain is chosen."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scripts.audit_spatial_geometry import geometry_facts, conservative_hull_volumes, OUT, ROOT


def test_legacy_represented_volume_matches_cell_box_not_point_hull():
    facts, cloud, family, _, _ = geometry_facts()
    np.testing.assert_allclose(family.X.min(axis=0), -0.000125, rtol=0, atol=1e-19)
    np.testing.assert_allclose(family.X.max(axis=0), 0.000125, rtol=0, atol=1e-19)
    assert np.isclose(cloud.volume.sum(), (0.000375)**3, rtol=1e-15)
    assert np.isclose(facts['represented_over_hull_volume'], 27/8, rtol=1e-15)


def test_horizon_threshold_and_clamp_coordinates_are_original():
    facts, cloud, family, _, top = geometry_facts()
    assert facts['physical_horizon_m'] == cloud.spacing*2.01
    assert facts['Yc_Pa'] == 0.1010025
    assert len(family.pairs) == 185
    assert np.all(cloud.X[top, 2] == .000125)
    assert np.all(cloud.X[cloud.fixed & ~top, 2] == -.000125)


def test_proposed_hull_rule_integrates_volume_and_first_moments():
    _, _, family, _, _ = geometry_facts()
    before = family.volume.copy()
    V = conservative_hull_volumes(family.X)
    assert np.isclose(V.sum(), .00025**3, rtol=1e-15)
    np.testing.assert_allclose(V@family.X, 0, rtol=0, atol=1e-30)
    assert np.all(V > 0)
    assert np.array_equal(before, family.volume)


def test_same_coarse_coordinates_compare_affine_energy_density_independently():
    rows = json.loads((OUT/'comparisons/ISOLATED_COARSE_QUADRATURE_COMPARISON.json').read_text())
    # SVK with lambda=0: W=mu*(Green_zz)^2 at uniaxial stretch 1.01.
    expected = 1000 * (((1.01**2)-1)/2)**2
    for row in rows:
        assert row['coordinates_unchanged']
        assert np.isclose(row['energy_per_represented_volume_Pa'], expected, rtol=1e-12)
        assert row['maximum_Fbar_error'] < 1e-12
        assert row['maximum_Ftilde_error'] < 1e-12
    assert np.isclose(rows[0]['elastic_energy_J']/rows[1]['elastic_energy_J'], 27/8, rtol=1e-12)


def test_legacy_interval_and_cone_are_reproduced():
    data = json.loads((OUT/'level_0/COARSE_REFERENCE_REPRODUCTION.json').read_text())
    assert data['status'] == 'VERIFIED'
    assert data['reproduced_interval'] == [.0099921875, .01]
    assert data['width'] <= 1e-5
    assert data['method'] == 'NONSMOOTH_STABILITY_BRACKET'
    for row in data['comparisons'].values():
        assert row['active_counts_equal'] and row['classification_equal']
        assert row['witness_absolute_cosine'] > 1-1e-6


def test_geometry_audit_did_not_change_legacy_mechanics_sources():
    records = json.loads((OUT/'provenance/BEFORE_MANIFEST.json').read_text())['files']
    names = ['pd_clot/geometry.py', 'pd_clot/pfpd_candidate/adr_adapter_at1.py',
             'pd_clot/pfpd_candidate/kernel.py', 'pd_clot/pfpd_candidate/baqp.py',
             'pd_clot/pfpd_candidate/at1_evolution.py',
             'pd_clot/pfpd_candidate/branch_stability_tracking.py']
    for name in names:
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == records[name]['sha256']
