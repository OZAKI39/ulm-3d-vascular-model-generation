import pytest
@pytest.mark.parametrize('index',[0,1])
def test_one_step_exact(parity_cases,index):
    d=parity_cases[index];e=d['errors'][0]
    assert e['exact_equal'] and e['neighbor_mismatch_count']==0
    assert all(v==0 for v in e['field_errors'].values())
