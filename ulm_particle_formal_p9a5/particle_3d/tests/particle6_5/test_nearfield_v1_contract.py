import json,pytest
from particle_3d.nearfield_regularization import contract,NearFieldRegularizationV1

def test_exact_values(repo):
 c=json.loads((repo/'particle_3d/contracts/NEAR_FIELD_REGULARIZATION_V1.json').read_text());assert c==contract()
 assert (c['chi_full'],c['chi_off'],c['h_molecular_floor_m'],c['h_suspension_ratio'])==(.01,.05,2e-9,1e-3)
 assert c['status']=='USER_APPROVED_FOR_PARTICLE6_5'
 assert not any(c[k] for k in ['production_neighbor_cutoff_frozen','production_neighbor_skin_frozen','production_particle_timestep_frozen'])
@pytest.mark.parametrize('value',[1.5e-9,3e-9,0.,1e-9,float('nan')])
def test_no_silent_floor_override(value):
 with pytest.raises(ValueError):NearFieldRegularizationV1(value)
