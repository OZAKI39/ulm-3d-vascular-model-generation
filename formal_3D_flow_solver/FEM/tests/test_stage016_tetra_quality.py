from stage016_helpers import *
import pytest
from fem3d.planar_port import volume_gates

def test_same_baseline_metric_and_nearest_boundary_classification():
 q=baseline_volume_audit()['quality'];old=read(ROOT/'outputs/stage01/medium/qc/geometry_qc.json')['quality']
 assert q['min_sicn']==old['gmsh_min_sicn']
 assert q['low_quality_nearest_boundary_counts']==old['low_quality_nearest_boundary_counts']
 assert q['cap_adjacent_below_0_1']==129 and q['total_below_0_1']==153

@pytest.mark.parametrize('field',['P1','P5','median','cap_adjacent','total','inverted'])
def test_each_tetra_quality_gate_is_hard(field):
 validity,q,proxy=good_volume()
 if field in ('P1','P5','median'):q['min_sicn'][field]=POLICY['tetra_quality'][field+'_min']-1e-5
 elif field=='cap_adjacent':q['cap_adjacent_below_0_1']=39
 elif field=='total':q['total_below_0_1']=77
 else:validity['negative_volume']=1
 assert volume_gates(validity,q,proxy,POLICY)['status']=='FAIL'
