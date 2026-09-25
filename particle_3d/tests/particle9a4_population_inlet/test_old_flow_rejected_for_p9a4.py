import pytest
from particle_3d.continuous_infusion import ContinuousInfusionSource

def test_old_flow_refused(simple_source):
 with pytest.raises(ValueError,match='AUTHORITATIVE_NEW'):
  ContinuousInfusionSource(sampler=simple_source.sampler,distribution=simple_source.distribution,
   checker=simple_source.checker,concentration_m3=1.,master_seed=1,source_contract_sha256='x',
   flow_sha256='129ebb77396550b300b20ce29cde7b03919e6037b7a81de60fe75d6c6efb616d')
