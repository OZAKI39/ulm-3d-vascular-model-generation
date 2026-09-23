from stage017_helpers import *
import pytest
from fem3d.adaptive_port import volume_acceptance,feedback,surface_search

@pytest.mark.parametrize('iteration',[4,5,6])
def test_failure_at_or_beyond_last_iteration_stops(iteration):
 m=measured_good();m['quality']['total_below_0_1']=200;m['quality']['cap_adjacent_below_0_1']=100;m['quality']['low_quality_nearest_boundary_counts']['INLET']=100
 d=feedback(iteration,m,volume_acceptance(m,BASE,POLICY,True),POLICY)
 assert d['stop'] and d['termination_reason']=='MAX_ADAPTIVE_ITERATIONS_REACHED'

def test_surface_search_never_exceeds_eight_actual_evaluations():
 calls=[]
 def evaluate(H,i,reason):
  calls.append(H);return {'trial':str(i),'geometry_status':'PASS','quality':{'status':'PASS'},'triangle_count':100}
 surface_search(evaluate,1.,.1,100.,POLICY)
 assert len(calls)==POLICY['limits']['maximum_surface_trials_per_port']
 assert len(read(OUT/'volume_mesh_ledger.json')['entries'])<=POLICY['limits']['maximum_total_volume_meshes']
