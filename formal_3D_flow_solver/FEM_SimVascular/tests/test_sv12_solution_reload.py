from sv12_support import artifact
def test_new_process_reload():
    r=artifact('solution_reload');assert r['status']=='PASS' and r['fresh_process']
    assert all(r['checks'].values())
