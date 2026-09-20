from sv13m_support import *
from sv_validation.sv13m import *
def test_history_including_stage_l_immutable():
 d=accepted('preservation_audit');assert d['historical']['entries']>=2000 and not d['historical']['changes']
 assert not d['reference_changed'] and d['old_FEM']['status']=='PASS'
def test_remote_historical_stacks_immutable():accepted('remote_preservation')
