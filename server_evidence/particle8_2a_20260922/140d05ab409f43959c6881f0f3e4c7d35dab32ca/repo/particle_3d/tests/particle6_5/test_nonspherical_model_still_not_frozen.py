import pytest
from particle_3d.particle6_cases import mixed_particles
from particle_3d.nearfield_handoff import require_admissible_initial
from particle_3d.hydrodynamic_resistance import NonsphericalLubricationNotFrozen

def test_no_equivalent_sphere(policy):
 records,_=mixed_particles()
 for p in records:
  if p.mode_code==1:continue
  with pytest.raises(NonsphericalLubricationNotFrozen):policy.reference_length(p.shape())
  with pytest.raises(NonsphericalLubricationNotFrozen):require_admissible_initial({p.particle_id:p.shape()},None,policy,.00345312)
