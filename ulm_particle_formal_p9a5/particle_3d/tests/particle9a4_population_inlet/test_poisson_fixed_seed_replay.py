from particle_3d.continuous_infusion import canonical_bytes
from particle_3d.population_inlet_p9a4 import generate

def test_exact_replay(simple_source):
 a=generate(simple_source,256);b=generate(simple_source,256)
 assert canonical_bytes(a.rows)==canonical_bytes(b.rows)
