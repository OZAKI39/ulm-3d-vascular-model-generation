from sv13p_support import *
def test_fastest_healthy_complete_window_only():
 w=read('winner');expected=choose_winner(read('reference_freeze')['window_baseline']['wall_time_s'],cases())
 assert all(w[k]==v for k,v in expected.items())
def test_fast_failed_or_incomplete_is_not_winner():
 ds=[dict(candidate='bad',status='FAIL',mode='window',steps=10,wall_time_s=1),dict(candidate='partial',status='PASS',mode='window',steps=3,wall_time_s=2),dict(candidate='good',status='PASS',mode='window',steps=10,wall_time_s=89)]
 assert choose_winner(100,ds)==dict(fastest='good',improvement=pytest.approx(.11),full_required=True)
 assert gain_classification(100,ds[-1])=='USEFUL_GAIN'
