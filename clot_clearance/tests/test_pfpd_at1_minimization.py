import json
import numpy as np
from scripts.verify_at1_evolution import OUT, config, scalar_oracle


def test_all_64_independent_bounded_minimizations_and_kkt():
    rows = json.loads((OUT/'minimization/CASES.json').read_text())
    assert len(rows) == 64
    for row in rows:
        assert row['scipy_success'] and row['audit']['passed']
        assert row['absolute_s_error'] <= config()['tolerances']['analytic_numeric_absolute_s']


def test_external_phasefieldx_local_energy_after_normalization_conversion():
    rows = json.loads((OUT/'minimization/CASES.json').read_text())
    assert all(r['external_success'] and r['external_absolute_s_error'] <= 5e-8 for r in rows)
    assert all(abs(r['external_local_coefficient']-2*r['Yc_Pa']) < 1e-15 for r in rows)


def test_live_scipy_oracle_does_not_receive_analytic_answer():
    c = config(); delta = c['geometry']['spacing_m']*c['geometry']['horizon_ratio']
    row = scalar_oracle((1.5,.1,c['phase']['Yc_design_Pa'],delta,c['phase']['Gc_demo_J_per_m2']))
    assert row['absolute_s_error'] < 5e-8
    assert len(row['endpoint_candidates']) == 3


def test_saved_bond_path_has_no_healing_and_has_second_damage_increment():
    rows = json.loads((OUT/'analytic/SINGLE_BOND_PATH.json').read_text())
    assert np.all(np.diff([r['s_new'] for r in rows]) >= 0)
    assert rows[-1]['s_new'] > rows[-2]['s_new']
    assert all(r['total_after'] <= r['total_before']+1e-15 for r in rows)
