from stage017_helpers import *
import pytest
from fem3d.adaptive_port import volume_acceptance,feedback

def test_largest_positive_port_contribution_is_refined():
 m=measured_good();m['quality']['cap_adjacent_below_0_1']=100;m['quality']['total_below_0_1']=110
 m['quality']['low_quality_nearest_boundary_counts'].update(INLET=40,OUTLET_01=60,WALL=10)
 a=volume_acceptance(m,BASE,POLICY,True);d=feedback(0,m,a,POLICY)
 assert d['decision']=='REFINE_OUTLET_01' and d['H_divisor']==1.2 and not d['stop']

def test_zero_baseline_low_counts_must_remain_zero():
 b=deepcopy(BASE);b['quality']['cap_adjacent_below_0_1']=b['quality']['total_below_0_1']=0
 m=measured_good();m['quality']['cap_adjacent_below_0_1']=m['quality']['total_below_0_1']=1
 a=volume_acceptance(m,b,POLICY,True)
 assert not a['quality_checks']['cap_low'] and not a['quality_checks']['total_low']
