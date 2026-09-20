from sv13j_support import *
def test_no_later_candidate_executed():
    d=load('compatibility_matrix')
    assert [r['status'] for r in d['candidates']]==['PASS','NOT_REQUIRED','NOT_REQUIRED']
    assert load('compatibility_winner')['frozen_on_first_make_pass']
def test_continuing_b_after_a_success_is_rejected():
    d=load('compatibility_matrix')['candidates'];d[1].update(status='FAIL',executed=True,make='FAIL')
    with pytest.raises(GateError,match='CONTINUED_AFTER_SUCCESS'):matrix_order_gate(d)

