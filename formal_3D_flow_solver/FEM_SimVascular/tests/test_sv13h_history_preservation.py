from sv13h_support import *
def test_all_historical_entries_old_fem_and_official_source_preserved():
    d=actual('preservation_audit')
    assert d['historical']['entries']==1017 and not d['historical']['changes']
    assert d['old_FEM']['entries']==7440 and d['old_FEM']['status']=='PASS'
    assert d['old_FEM']['git_unchanged'] and d['official_source']['status']=='PASS'
    assert not d['reference_changed']

