import json,csv
from collections import Counter

def test_all_original_proposals_and_distinct_partitions(evidence_root):
    root=evidence_root/'admission';s=json.loads((root/'ADMISSION_BASIN_AUDIT.json').read_text())
    recheck=json.loads((root/'ORIGINAL_ADMISSION_RECHECK.json').read_text());assert recheck['all_exact']
    rows=list(csv.DictReader((root/'proposal_by_basin.csv').open()))
    assert len(rows)==s['totals']['proposals']==recheck['checked_proposals']==149716
    assert s['totals']==dict(proposals=149716,scheduled=2200,admitted=1969,guard_exhausted=231)
    for basin,c in s['by_basin'].items():
        r=[x for x in rows if x['basin']==basin]
        assert len(r)==c['proposal_count']
        assert sum(x['status']=='ACCEPTED' for x in r)==c.get('admitted_count',0)
        assert sum(x['draw']=='0' for x in r)==c.get('scheduled_first_proposal_count',0)
    assert len(s['per_id'])==2200 and s['no_radius_resampling']
