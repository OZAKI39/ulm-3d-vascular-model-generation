from stage016_helpers import *
import pytest
from fem3d.mesh_input import surface_topology

@pytest.mark.parametrize('name',CANDIDATES)
def test_closed_single_surface_recomputed(name):
 d=candidate(name);r=surface_topology(d['points_m'],d['triangles'])
 assert r==read(OUT/name/'qc/surface_invariants.json')['topology']

def test_removed_surface_triangle_is_detected():
 d=candidate('sparse_C')
 with pytest.raises(ValueError):surface_topology(d['points_m'],d['triangles'][1:])
