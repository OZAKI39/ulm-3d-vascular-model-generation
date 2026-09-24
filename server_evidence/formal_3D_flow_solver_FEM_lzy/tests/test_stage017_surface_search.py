from stage017_helpers import *
import pytest
from fem3d.adaptive_port import surface_search,choose_surface

def evaluator(rule,counts=None):
 def run(H,index,reason):
  return {'trial':f'trial_{index:02d}','geometry_status':'PASS','quality':{'status':'PASS' if rule(H) else 'FAIL'},'triangle_count':counts(H) if counts else int(100/H)}
 return run

def test_pass_coarsens_then_log_bisects_within_both_limits():
 r=surface_search(evaluator(lambda H:H<=1.6),1.,.1,10.,POLICY)
 h=[t['H'] for t in r['trials']]
 assert h[:3]==[1.,1.25,1.5625] and h[3]>1.6
 assert r['refinement_evaluations']<=4 and len(r['trials'])<=8
 assert h[4]==pytest.approx(np.sqrt(h[2]*h[3]))

def test_initial_poor_quality_refines_then_brackets():
 r=surface_search(evaluator(lambda H:H<=.7),1.,.1,10.,POLICY)
 assert [t['H'] for t in r['trials'][:3]]==pytest.approx([1.,.8,.64])
 assert r['status']=='PASS' and r['selected_H']<=.7

def test_no_passing_surface_at_rim_scale_fails():
 r=surface_search(evaluator(lambda H:False),1.,.8,10.,POLICY)
 assert r['status']=='FAIL' and r['termination_reason']=='NO_PASS_AT_RIM_SCALE'

def test_actual_min_triangle_choice_not_monotone_H_assumption():
 records=[{'trial':str(i),'H':H,'triangle_count':n,'geometry_status':'PASS','quality':{'status':'PASS'}} for i,(H,n) in enumerate([(1.,20),(2.,16),(3.,19),(1.5,16)])]
 assert choose_surface(records)['H']==2.

def test_geometry_failure_stops_before_quality_decision():
 r=surface_search(lambda H,i,r:{'trial':'bad','triangle_count':10,'geometry_status':'FAIL','quality':None},1.,.1,10.,POLICY)
 assert len(r['trials'])==1 and r['status']=='FAIL'
