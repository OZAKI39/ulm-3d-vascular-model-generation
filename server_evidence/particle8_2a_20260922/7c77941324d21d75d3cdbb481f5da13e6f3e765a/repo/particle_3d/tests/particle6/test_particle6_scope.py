import json
import pytest
from particle_3d.particle6_checkpoint import MODEL_FLAGS
from particle_3d.lammps_neighbors import ValidationNeighborPolicy
from particle_3d.lammps_state import schema_contract
from particle_3d.particle6_stepper import bind_query_dependency
from particle_3d import particle5_motion as p5

def test_particle6_scope(repo):
    assert MODEL_FLAGS['lammps_role']=='STATE_NEIGHBOR_CHECKPOINT_ONLY'
    for k in ['lammps_force_integration','lammps_time_integration','production_neighbor_cutoff_frozen','production_neighbor_skin_frozen','production_lubrication_cutoff_frozen','production_particle_timestep_frozen','particle7_started']:assert MODEL_FLAGS[k] is False
    assert MODEL_FLAGS['no_cfd_executed']
    assert json.loads((repo/'particle_3d/contracts/particle6_state_v1.json').read_text())==schema_contract()

def test_upstream_binding_is_private():
    original=p5.resistance_trial.__globals__['assemble_scene'];replacement=object()
    bound=bind_query_dependency(p5.resistance_trial,{'assemble_scene':replacement})
    assert bound.__code__ is p5.resistance_trial.__code__
    assert bound.__globals__['assemble_scene'] is replacement
    assert p5.resistance_trial.__globals__['assemble_scene'] is original

@pytest.mark.parametrize('cutoff,skin,role',[(0,0,'VALIDATION_ONLY'),(1,-1,'VALIDATION_ONLY'),(1,0,'PRODUCTION'),(float('nan'),0,'VALIDATION_ONLY')])
def test_policy_rejects_unfrozen_production(cutoff,skin,role):
    with pytest.raises(ValueError):ValidationNeighborPolicy(cutoff,skin,'test',role)
