import numpy as np
from particle_3d.interior_section import evaluate_candidate


def test_interior_perimeter_is_only_wall(tube):
    grid, candidate, topo, wall, caps = tube
    result, _ = evaluate_candidate(grid, candidate, topo, wall, caps, 2e-15, 1e-6)
    assert result['accepted']
    assert result['boundary_entirely_wall']
    assert not result['open_boundary_intersection']
    assert np.isclose(result['minimum_cap_distance_m'], .5e-6)


def test_oblique_section_touching_cap_rejected(tube):
    grid, candidate, topo, wall, caps = tube
    candidate.update(center_m=np.array([0, 0, .1e-6]),
                     normal=np.array([1., 0, 1.])/np.sqrt(2), arclength_m=.1e-6)
    result, _ = evaluate_candidate(grid, candidate, topo, wall, caps, 2e-15, 1e-6)
    assert result['open_boundary_intersection']
    assert not result['accepted']
