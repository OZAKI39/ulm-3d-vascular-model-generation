import json
from particle_3d.nearfield_regularization import contract

def test_roles_and_facts(report):
 s=json.loads((report/'literature/SOURCES.json').read_text())['sources'];d={x['id']:x for x in s}
 assert {k:d[k]['role'] for k in 'ABCDE'}==contract()['literature_roles']
 assert d['A']['facts']['outer_cutoff_over_radius']==.05 and d['A']['facts']['lower_gap_order_over_radius']==1e-3
 assert all(x['doi'] and x['url'] and x['locator'] and x['scope_limit'] for x in s)
 text=(report/'literature/NEAR_FIELD_REGULARIZATION_LITERATURE.md').read_text()
 assert 'PROJECT MODEL CHOICE' in text and '2 nm' in text and '4 nm' in text
