import json,csv

def test_sizes_are_original_draws_at_the_correct_partition(evidence_root):
    root=evidence_root/'admission';s=json.loads((root/'ADMISSION_BASIN_AUDIT.json').read_text())
    rows=list(csv.DictReader((root/'proposal_by_basin.csv').open()))
    for basin,sizes in s['diameters'].items():
        part=[r for r in rows if r['basin']==basin]
        for key,predicate in [('admitted_diameter_um',lambda r:r['status']=='ACCEPTED'),('scheduled_first_proposal_diameter_um',lambda r:r['draw']=='0'),('guard_exhausted_diameter_um',lambda r:r['draw']=='0' and r['guard_exhausted']=='True')]:
            assert sizes.get(key,[])==[float(r['diameter_um']) for r in part if predicate(r)]
    by_id={}
    for r in rows:by_id.setdefault(r['stable_id'],set()).add(float(r['diameter_um']))
    assert all(len(v)==1 for v in by_id.values())
