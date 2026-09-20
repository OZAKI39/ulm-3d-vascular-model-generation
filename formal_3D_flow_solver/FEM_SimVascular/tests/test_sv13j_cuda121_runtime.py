from sv13j_support import *
def test_candidate_c_explicitly_not_required_after_a_success():
    d=load('cuda121_runtime')
    assert d['status']=='NOT_REQUIRED' and not d['executed']
    assert d['kernel']=='NOT_RUN' and d['Thrust_version'] is None
    assert candidate()['make']=='PASS'

