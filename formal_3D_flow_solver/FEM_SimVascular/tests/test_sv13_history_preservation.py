from sv13_support import artifact
def test_all_historical_stages_and_old_fem_preserved():
    d=artifact('preservation_audit')
    assert d['status']=='PASS' and not d['history_changes'] and d['old_fem_git_unchanged']
    assert d['frozen_inputs_unchanged']
