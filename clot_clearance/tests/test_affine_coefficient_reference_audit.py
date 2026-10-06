"""Exactly six focused acceptance tests for the saved local coefficient audit."""
import numpy as np
import pytest

from pd_clot.pfpd_candidate.affine_reference_audit import (
    OLD, OUT, PATHS, LD, closed_form, frob, historical_integrity, load,
)
from pd_clot.pfpd_candidate.energy import constitutive
from pd_clot.pfpd_candidate.force_stencil_audit import material


@pytest.fixture(scope='module')
def data():
    with np.load(OUT/'LOCAL_COEFFICIENT_ARRAYS.npz', allow_pickle=False) as source:
        arrays = {k: source[k] for k in source.files}
    return load(OUT/'SCIENTIFIC_RESULTS.json'), arrays


def test_closed_form_affine_svk_matches_current_constitutive(data):
    _, a = data
    W, P = constitutive(a['F_exact'], material())
    decimal = closed_form('0.003')
    stored = closed_form(LD(a['F_exact'][2, 2])-LD(1))
    assert float(frob(P.astype(LD)-stored['P'])/frob(P)) < 1e-12
    assert float(frob(P.astype(LD)-decimal['P'])/frob(P)) < 1e-12
    assert abs(W-float(decimal['W']))/W < 1e-12
    assert np.array_equal(P, a['P_exact'])
    assert np.count_nonzero(decimal['P']) == 1


def test_raw_n27_reproduces_previous_audit(data):
    s, a = data
    old = load(OLD/'STATUS.json')
    error = np.max(frob(a[PATHS[0]]-a['Zbulk']))/frob(a['Zbulk'])
    assert abs(error-old['Z_BULK_MATCH']) <= 1e-16
    assert error > 1e-12  # The historical exceedance must remain visible.
    n27 = load(OLD/'N27.json')
    center = next(row for row in n27['observations'] if row['label'] == 'center')
    for row in center['coefficients']['representative_coefficients']:
        i = np.where(a['Z_node_ids'] == row['ID'])[0][0]
        np.testing.assert_array_equal(a[PATHS[0]][i], row['Z'])
    assert s['scope']['Z_nodes'] == 257 and s['scope']['Fbar_nodes'] == 1797


def test_exact_p_z_matches_saved_diagnostic(data):
    _, a = data
    old = load(OLD/'NOISE_BASELINE.json')['N27_stress_reconstruction']
    error = np.max(frob(a[PATHS[1]]-a['Zbulk']))/frob(a['Zbulk'])
    assert abs(error-old['constant_P_Z_max_bulk_relative_error']) <= 1e-16
    assert error <= 1e-12


def test_delta_p_to_delta_z_linear_identity(data):
    _, a = data
    np.testing.assert_array_equal(a['delta_P'], a['Praw']-a['P_exact'])
    observed = a[PATHS[0]]-a[PATHS[1]]
    np.testing.assert_array_equal(observed, a['delta_Z_observed'])
    error = np.max(frob(a['delta_Z_predicted']-observed))/frob(a['Zbulk'])
    assert error <= 1e-13


def test_pairwise_and_longdouble_preserve_coefficients(data):
    s, a = data
    scale = frob(a['Zbulk'])
    for key in ['Z_pairwise', 'Z_longdouble_same_terms']:
        error = np.max(frob(a[key].astype(LD)-a[PATHS[0]].astype(LD)))/scale
        assert error <= 1e-13
    matched = np.max(frob(a[PATHS[3]]-a['Zbulk_high_precision']))/scale
    assert matched <= 1e-13
    assert np.finfo(LD).eps < np.finfo(float).eps
    assert s['summation']['convention'].startswith('Same float64 P')


def test_historical_review_results_and_production_unchanged():
    integrity = historical_integrity()
    assert integrity['protected_files'] >= 51
    assert integrity['PASS'], integrity['differences']
    assert load(OLD/'STATUS.json')['N27_COEFFICIENT_GATE'] == 'FAIL'
    assert load(OLD/'PREDECLARED_CONFIG.json')['coefficient_match_tolerance'] == 1e-12
