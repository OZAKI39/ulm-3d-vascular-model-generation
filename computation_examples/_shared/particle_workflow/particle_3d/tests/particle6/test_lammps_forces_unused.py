import ast,inspect
import pytest
from particle_3d.lammps_bridge import LammpsParticleBridge
from particle_3d.particle6_cases import mixed_particles
from particle_3d import particle6_stepper

def test_force_zero_before_after_every_build(parity_cases):
    for case in parity_cases:
        assert len(case['force_audits'])>=400
        assert all(row['max_lammps_force']==0 and row['lammps_step']==0 and not row['used_as_physics_input'] for row in case['force_audits'])

def test_poison_is_detected_before_run_can_clear_force():
    ps,policy=mixed_particles()
    with LammpsParticleBridge(ps,policy) as b:
        b._lmp.numpy.extract_atom('f')[0,0]=1
        with pytest.raises(ValueError,match='NONZERO_UNUSED'):b.rebuild()

def test_physics_adapter_never_extracts_force_mass_torque():
    source=inspect.getsource(particle6_stepper)
    assert 'extract_atom' not in source and 'extract_global' not in source and 'mass' not in source
