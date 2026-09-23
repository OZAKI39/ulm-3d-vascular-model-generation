from stage016_helpers import *
import pytest
from fem3d.planar_port import select_candidate,volume_gates

def record(name,dofs,tetra=150000):
 v,q,p=good_volume();p['N_P2_velocity_proxy']=dofs;p['N_tetra']=tetra;g=volume_gates(v,q,p,POLICY)
 return {'candidate':name,'surface_status':'PASS','volume_status':g['status'],'volume_gate_checks':g['checks'],'proxy':p,'quality':q}

def test_malicious_both_quality_pass_requires_lowest_DOF():
 dense=record('denser',1000000,150000);sparse=record('sparser',900000,160000)
 assert select_candidate([dense,sparse])['selected_candidate']=='sparser'

def test_cost_rejected_candidate_cannot_win():
 cheap=record('bad_tetra',800000,200001);good=record('good',900000)
 assert select_candidate([cheap,good])['selected_candidate']=='good'

def test_tie_break_order_tetra_then_low_count_then_P1():
 a=record('a',900000,155000);b=record('b',900000,150000)
 assert select_candidate([a,b])['selected_candidate']=='b'
 a['proxy']['N_tetra']=150000;a['quality']['cap_adjacent_below_0_1']=9
 assert select_candidate([a,b])['selected_candidate']=='a'
 b['quality']['cap_adjacent_below_0_1']=9;b['quality']['total_below_0_1']=19
 assert select_candidate([a,b])['selected_candidate']=='b'
 a['quality']['total_below_0_1']=19;a['quality']['min_sicn']['P1']=.6
 assert select_candidate([b,a])['selected_candidate']=='a'

def test_actual_decision_has_no_survivor():
 d=read(REPORT/'candidate_selection.json');r=select_candidate(d['records'])
 assert r['selected_candidate'] is None and r['status']=='FAIL' and not r['survivors']
