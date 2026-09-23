from stage016_helpers import *
import pytest
from fem3d.planar_port import cap_gates
from fem3d.cap_remesh import triangle_quality

@pytest.mark.parametrize('name',CANDIDATES)
def test_actual_counts_recomputed_and_rejected(name):
 d=candidate(name);counts={n:int(np.count_nonzero(d['facet_tags']==p['entity_id'])) for n,p in CONTRACT['ports'].items()}
 q=triangle_quality(d['points_m'],d['triangles'][d['facet_tags']!=1])[0]
 gate=cap_gates(counts,q,POLICY)
 assert gate==read(OUT/name/'qc/surface_invariants.json')['cap_gate']
 assert gate['status']=='FAIL' and not gate['checks']['per_port_density'] and not gate['checks']['total_density']

def test_malicious_per_port_overbudget_despite_perfect_triangles():
 counts={'inlet':225,'outlet_01':44,'outlet_02':42,'outlet_03':49}
 gate=cap_gates(counts,np.ones(360),POLICY)
 assert gate['status']=='FAIL' and gate['checks']['P5'] and not gate['checks']['per_port_density']

def test_total_budget_cannot_be_hidden_by_good_quality():
 counts={k:201 for k in CONTRACT['ports']}
 assert not cap_gates(counts,np.ones(804),POLICY)['checks']['total_density']
