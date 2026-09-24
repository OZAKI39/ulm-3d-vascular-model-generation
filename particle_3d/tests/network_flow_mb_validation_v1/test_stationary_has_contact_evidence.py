import json
def test_stationary_has_contact_evidence(report):
    rows=json.loads((report/'data/new_paired_metrics.json').read_text())
    for row in rows:
        assert row['status']!='SOLVER_FAIL'
        if row['status']=='STATIONARY':assert row['stationary_three_contact_supported']
