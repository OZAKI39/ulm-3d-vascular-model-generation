"""Guard the actual nonlinear stopping rule and residual normalization."""
from pathlib import Path
import csv
import json
import sys
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts/flow_2mmps'))
from residual_figures import nonlinear_stop,orders_reduced,OUT


def test_order_reduction_is_invariant_to_common_reference():
    assert orders_reduced(1e-2,1e-12)==pytest.approx(10.)
    assert orders_reduced(1e-2/.0057,1e-12/.0057)==pytest.approx(10.)
    assert orders_reduced(2.,2.)==pytest.approx(0.)
    with pytest.raises(ValueError):orders_reduced(1.,0.)


def test_nonlinear_stop_uses_either_global_or_step_relative_reference():
    assert nonlinear_stop(2,8.586e-11,1.835e-6,1e-10,2,12)
    assert nonlinear_stop(3,1.2e-10,8e-11,1e-10,2,12)
    assert not nonlinear_stop(2,1.594e-10,1.824e-6,1e-10,2,12)


def test_minimum_iterations_and_iteration_cap_have_distinct_roles():
    assert not nonlinear_stop(1,1e-14,1e-14,1e-10,2,12)
    assert nonlinear_stop(2,1e-14,1e-14,1e-10,2,12)
    # The native implementation stops at the cap; that alone is NOT convergence.
    assert nonlinear_stop(12,1.,1.,1e-10,2,12)


def test_logged_transitions_and_fixed_reference_are_audited():
    a=json.loads((OUT/'NONLINEAR_RESIDUAL_AUDIT.json').read_text())
    assert a['all_pass'] and a['maximum_iteration_cap_never_reached']
    assert a['reference_norm']==pytest.approx(.005714650454390)
    assert not a['numerical_floor_established'] and not a['near_zero_substitution_triggered']
    assert [(r['step'],r['previous_iterations'],r['current_iterations']) for r in a['transitions']]==[(4,4,3),(23,3,2)]
    assert all(r['no_preconditioner_rebuild_in_current_step'] for r in a['transitions'])
    assert a['transitions'][1]['global_tolerance_met'] and not a['transitions'][1]['step_relative_tolerance_met']
    with (OUT/'data/nonlinear_step_audit.csv').open() as f:rows=list(csv.DictReader(f))
    assert len(rows)==71
    for row in rows:
        start,end=float(row['start_norm']),float(row['end_norm'])
        assert float(row['start_global'])==pytest.approx(start/a['reference_norm'])
        assert float(row['end_global'])==pytest.approx(end/a['reference_norm'])
        assert float(row['reduction_orders'])==pytest.approx(np.log10(start/end))
