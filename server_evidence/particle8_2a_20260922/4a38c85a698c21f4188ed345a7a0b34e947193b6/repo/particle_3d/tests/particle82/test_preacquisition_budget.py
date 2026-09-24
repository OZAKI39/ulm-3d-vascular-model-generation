import pytest
from particle_3d.particle82_budget import choose_natural_guard

def test_choose_largest_guard_that_fits_without_new_outcome_information():
    scale={'parallel_tracks_hour':3000}
    diagnostics={'factors':[dict(guard_factor=4,performance={'wall_seconds':700},all_original_prefixes_exact=True),dict(guard_factor=16,performance={'wall_seconds':2500},all_original_prefixes_exact=True)]}
    result=choose_natural_guard(scale,diagnostics,dict(natural_initial_scheduled=5000,natural_batch_wall_budget_hours=3))
    assert result['chosen_guard_factor']==4
    assert result['selection_before_new_birth_ledger'] and result['no_new_sample_outcome_or_outlet_conditioning']
    assert result['estimates'][-1]['predicted_seconds_with_margin']>10800

def test_budget_selection_rejects_invalid_prefix_and_infeasible_budget():
    with pytest.raises(ValueError,match='baseline forecast'):
        choose_natural_guard({'parallel_tracks_hour':1},{'factors':[]},dict(natural_initial_scheduled=5000,natural_batch_wall_budget_hours=3))
    with pytest.raises(ValueError,match='invalid continuation'):
        choose_natural_guard({'parallel_tracks_hour':3000},{'factors':[{'all_original_prefixes_exact':False}]},dict(natural_initial_scheduled=5000,natural_batch_wall_budget_hours=3))
