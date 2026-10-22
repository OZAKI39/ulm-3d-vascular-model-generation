import numpy as np
from particle_3d.particle6_cases import mixed_particles
from particle_3d.lammps_bridge import LammpsParticleBridge
from particle_3d.lammps_neighbors import canonical_pairs

def test_half_full_duplicate_and_stable_ids():
    assert canonical_pairs([(8,2),(2,8),(8,2),(2,2)])==[(2,8)]
    p,policy=mixed_particles();expected=None
    for order in [p,p[::-1],list(np.random.default_rng(67).permutation(p))]:
        with LammpsParticleBridge(order,policy) as b:
            if expected is None:expected=b.raw_pairs
            assert b.raw_pairs==expected
            assert set(i for pair in b.raw_pairs for i in pair).issubset(x.particle_id for x in p)
