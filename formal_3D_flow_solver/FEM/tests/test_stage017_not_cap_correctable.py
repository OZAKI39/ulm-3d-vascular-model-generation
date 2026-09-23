from stage017_helpers import *
import pytest
from fem3d.adaptive_port import volume_acceptance,feedback

def test_malicious_wall_dominated_failure_cannot_refine_caps():
 m=measured_good();m['quality']['total_below_0_1']=100;m['quality']['low_quality_nearest_boundary_counts']['WALL']=100
 m['quality']['worst_elements']=[{'nearest_boundary_patch':'WALL'}]*20
 d=feedback(0,m,volume_acceptance(m,BASE,POLICY,True),POLICY)
 assert d['stop'] and d['termination_reason']=='NOT_CAP_CORRECTABLE' and d['changed_port'] is None
