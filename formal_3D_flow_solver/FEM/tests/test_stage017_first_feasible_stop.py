from stage017_helpers import *
import pytest
from fem3d.adaptive_port import volume_acceptance,feedback

def test_quality_pass_never_triggers_more_refinement():
 m=measured_good();d=feedback(0,m,volume_acceptance(m,BASE,POLICY,True),POLICY)
 assert d['decision']=='ACCEPT' and d['stop'] and d['changed_port'] is None

def test_actual_first_feasible_was_final_production_iteration():
 r=result();assert r['termination_reason']=='FIRST_FEASIBLE_ACCEPTED'
 assert r['iterations'][-1]['iteration']==r['selected_iteration']
 assert sum(row['acceptance']['status']=='PASS' for row in r['iterations'])==1
 ledger=read(OUT/'volume_mesh_ledger.json')['entries']
 assert len([v for v in ledger if v['run']=='production'])==len(r['iterations'])
