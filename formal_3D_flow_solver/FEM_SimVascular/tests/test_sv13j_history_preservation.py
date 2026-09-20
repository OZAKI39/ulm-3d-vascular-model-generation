from sv13j_support import *
def test_all_historical_evidence_and_old_fem_preserved():
    d=actual('preservation_audit')
    assert d['historical']['entries']==1245 and not d['historical']['changes']
    assert d['old_FEM']['entries']==7440 and d['old_FEM']['status']=='PASS' and d['old_FEM']['git_unchanged']
    assert d['official_source']['status']=='PASS' and not d['reference_changed']

