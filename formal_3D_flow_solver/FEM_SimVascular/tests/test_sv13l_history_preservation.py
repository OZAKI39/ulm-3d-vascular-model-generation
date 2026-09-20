from sv13l_support import *
def test_actual_history_preservation():
    d=actual('preservation_audit');assert d['historical']['changes']==[] and d['reference_changed']==[]
    assert d['old_FEM']['status']=='PASS'
    assert all(actual('environment_preservation')['checks'].values())
