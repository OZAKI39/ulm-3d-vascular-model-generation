from sv13_support import artifact
from sv_validation.sv13 import load
def test_replay_matches_original_actual_first_crossing():
    r=artifact('steady_stop_replay');ref=artifact('reference_freeze')
    assert r['status']=='PASS' and r['first_stop_step']==ref['baseline_steady']['first_five_joint_intervals_end_step']
    assert r['actual_field_count']==len(r['rows'])==len(load('saved_state_qc','sv1_2')['states'])
    assert not any(x['stop'] for x in r['rows'] if x['step']<r['first_stop_step'])
