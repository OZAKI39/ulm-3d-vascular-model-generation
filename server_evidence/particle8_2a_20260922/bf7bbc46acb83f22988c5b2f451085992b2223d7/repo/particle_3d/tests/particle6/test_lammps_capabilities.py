from particle_3d.lammps_bridge import LammpsParticleBridge
from particle_3d.particle6_cases import mixed_particles
from particle_3d.particle6_validation import environment

def test_actual_lammps_capabilities():
    p,policy=mixed_particles()
    with LammpsParticleBridge(p,policy) as b:
        e=environment(b)
        assert e['lammps_version']==20250722 and e['mpi_world_size']==1
        for k in ['python_module','pair_zero','property_atom','python_neighbor_access']:assert e['capabilities'][k]
        assert len(e['installed_packages'])>0
