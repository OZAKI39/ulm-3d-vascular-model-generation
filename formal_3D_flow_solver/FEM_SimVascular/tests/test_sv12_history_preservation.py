from sv12_support import artifact
def test_history_and_old_fem_unchanged():
    data=artifact('preservation_audit')
    assert data['status']=='PASS' and data['old_fem_git_unchanged']
    assert not data['history_changes']
