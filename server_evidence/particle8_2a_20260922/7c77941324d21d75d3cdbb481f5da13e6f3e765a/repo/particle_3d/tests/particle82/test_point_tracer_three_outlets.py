import json
from collections import Counter
from particle_3d.particle82_diagnostics import basin_rows

def test_direct_actual_point_basin_counts(point_root):
    s=json.loads((point_root/'POINT_TRACER_SUMMARY.json').read_text());rows=basin_rows(point_root)
    assert len(rows)==s['count']>=100000
    c=Counter(r['outlet'] for r in rows)
    for outlet in ['OUTLET_01','OUTLET_02','OUTLET_03']:assert c[outlet]==s['outlet_counts'][outlet]>0
    assert c[None]==s['unresolved']
    assert all(r['end_reason']==r['outlet'] for r in rows if r['outlet'])
