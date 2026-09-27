from dataclasses import replace
import pytest
from particle_3d.lammps_bridge import LammpsParticleBridge
from particle_3d.particle6_cases import mixed_particles

@pytest.mark.parametrize('field,value',[('particle_id',0),('particle_id',-1),('particle_id',2**31),('mode_code',4),('type_code',3)])
def test_invalid_identity_rejected(field,value):
    p,_=mixed_particles()
    with pytest.raises(ValueError):replace(p[0],**{field:value})

def test_duplicate_and_disagreeing_types():
    p,policy=mixed_particles()
    with pytest.raises(ValueError):LammpsParticleBridge([p[0],p[0]],policy)
    with LammpsParticleBridge(p[::-1],policy) as b:
        assert [x.particle_id for x in b.read()]==sorted(x.particle_id for x in p)
        b._lmp.numpy.extract_atom('type')[0]=1
        with pytest.raises(ValueError):b.read()
