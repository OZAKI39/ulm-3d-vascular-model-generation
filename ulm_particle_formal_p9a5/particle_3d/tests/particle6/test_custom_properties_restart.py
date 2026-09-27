import pytest
@pytest.mark.parametrize('index',[0,1])
def test_custom_property_binary_restart(restart_cases,index):
    d=restart_cases[index];m=d['metadata']
    assert m['errors']['exact_equal'] and m['destroyed_before_new_instance']
    assert m['global_before']==m['global_after']
    assert all(e==0 for e in m['errors']['field_errors'].values())
    if index==1:assert {p['mode_code'] for p in m['after']}=={1,2,3}
