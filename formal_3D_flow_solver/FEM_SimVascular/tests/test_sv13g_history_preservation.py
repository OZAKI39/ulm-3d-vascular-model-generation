from sv13g_support import *
def test_historical_and_old_fem_preserved():
    d=load('preservation_audit');assert d['status']=='PASS'
    assert d['historical']['status']=='PASS' and d['old_FEM']['status']=='PASS'
    assert d['official_source']['status']=='PASS'
