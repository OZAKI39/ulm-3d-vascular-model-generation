import json,csv
from collections import Counter
def test_transition_table_consistency(report,events):
    rows=json.loads((report/'data/paired_outcome_transition.json').read_text())
    summary=json.loads((report/'data/analysis_summary.json').read_text())
    assert len(rows)==len(events)==30 and len({r['bubble_id'] for r in rows})==30
    assert len(list(csv.DictReader((report/'data/paired_outcome_transition.csv').open())))==30
    matrix=Counter()
    for r in rows:
        for label in ['OLD','NEW']:
            m=json.loads((report/'outputs'/label/'trajectories'/f"mb_{r['bubble_id']:06d}.json").read_text())
            if m['completed']:assert r[label.lower()+'_outlet']=='O'+str(int(m['exit_outlet'].split('_')[-1]))
            else:assert not r[label.lower()+'_outlet']
        matrix[r['old_outcome']+' -> '+r['new_outcome']]+=1
    assert dict(matrix)==summary['transition_matrix']
    assert sum(r['outlet_changed'] for r in rows)==summary['outlet_changed_count']
