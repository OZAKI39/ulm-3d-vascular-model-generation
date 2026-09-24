from stage017_helpers import *
import pytest
from fem3d.adaptive_port import volume_acceptance,feedback

@pytest.mark.parametrize('key',['N_P2_velocity_proxy','N_tetra'])
def test_overbudget_even_with_excellent_quality_stops(key):
 m=measured_good();m['proxy'][key]=int(BASE['proxy'][key]*1.35)+1
 d=feedback(0,m,volume_acceptance(m,BASE,POLICY,True),POLICY)
 assert d['stop'] and d['termination_reason']=='FEM_COST_BUDGET_EXCEEDED'
