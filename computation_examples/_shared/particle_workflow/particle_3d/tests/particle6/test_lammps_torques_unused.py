from particle_3d.lammps_bridge import LammpsParticleBridge
from particle_3d.particle6_cases import mixed_particles

def test_atomic_torque_absent_and_not_a_physics_input():
    p,policy=mixed_particles()
    with LammpsParticleBridge(p,policy) as b:
        assert b._lmp.numpy.extract_atom('torque') is None
        b.rebuild()
        assert all(x['torque_storage']=='UNDEFINED_FOR_ATOM_STYLE_ATOMIC' and x['max_lammps_torque']==0 and not x['used_as_physics_input'] for x in b.force_audits)
