import pytest

@pytest.mark.parametrize('a',[.1e-6,.588015722765549e-6,1.2658269815352818e-6,2e-6,10e-6])
def test_definition(policy,a):
 d=policy.lower_handoff_gap(a);assert d['h_lower_m']==max(2e-9,.001*a)
 assert d['h_molecular_component_m']==2e-9 and d['h_suspension_component_m']==.001*a
