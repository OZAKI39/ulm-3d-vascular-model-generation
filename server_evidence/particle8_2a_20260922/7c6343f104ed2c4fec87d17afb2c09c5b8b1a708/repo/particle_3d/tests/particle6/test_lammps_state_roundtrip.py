from dataclasses import replace
import pytest
from particle_3d.lammps_bridge import LammpsParticleBridge
from particle_3d.particle6_cases import mixed_particles,state_errors

def test_all_state_roundtrip():
    p,policy=mixed_particles()
    with LammpsParticleBridge(p[::-1],policy) as b:
        assert state_errors(p,b.read())['exact_equal']
        updated=[replace(x,position=x.position+[0,2e-6,0],velocity=x.velocity*2,omega=x.omega*3) for x in p]
        b.write(updated);b.rebuild()
        assert state_errors(updated,b.read())['exact_equal']
        with pytest.raises(ValueError):b.write(updated[:-1])
