import pytest
@pytest.mark.parametrize('index',[0,1])
def test_actual_40_destroy_restore_60_matches_continuous(restart_cases,index):
    d=restart_cases[index]
    assert d['checkpoint_step']==40 and d['steps']==100
    assert d['physical_time_equal'] and d['step_index_equal']
    assert all(e['exact_equal'] and e['neighbor_mismatch_count']==0 for e in d['errors'])
