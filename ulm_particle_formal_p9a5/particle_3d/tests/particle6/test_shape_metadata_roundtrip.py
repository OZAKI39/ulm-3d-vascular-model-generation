from dataclasses import replace
import numpy as np
import pytest
from particle_3d.lammps_bridge import LammpsParticleBridge
from particle_3d.particle6_cases import mixed_particles

@pytest.mark.parametrize('mode',[1,2,3])
def test_shape_metadata_roundtrip(mode):
    ps,policy=mixed_particles();p=next(p for p in ps if p.mode_code==mode)
    with LammpsParticleBridge([p],policy) as b:
        restored=b.read()[0]
        for direction in np.eye(3):np.testing.assert_array_equal(p.shape().support(direction),restored.shape().support(direction))
        if mode==3:
            q=replace(restored,q=np.array([1.,0,0,0]))
            np.testing.assert_array_equal(q.shape().axis_world,restored.capsule_axis)
            assert q.capsule_radius==restored.capsule_radius and q.capsule_length==restored.capsule_length
