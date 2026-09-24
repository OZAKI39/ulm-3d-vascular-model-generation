import numpy as np
from particle_3d.interior_section import candidate_sequence, evaluate_candidate, root_topology


def test_sequence_never_passes_first_junction(tube):
    _, _, topo, _, _ = tube
    for refinement in (1, 2):
        rows = list(candidate_sequence(topo, .1e-6, refinement))
        assert all(0 < c['arclength_m'] < topo['arclength_m'][-1] for c in rows)
        assert np.allclose(rows[0]['center_m'], [0, 0, .05e-6/refinement], rtol=0, atol=1e-20)
    assert topo['junction_node_id'] == 2


def test_at_or_after_junction_rejected(tube):
    grid, candidate, topo, wall, caps = tube
    candidate.update(center_m=np.array([0, 0, .9e-6]), arclength_m=.9e-6)
    result, _ = evaluate_candidate(grid, candidate, topo, wall, caps, 2e-15, 1e-6)
    assert not result['accepted']
    assert 'NOT_BEFORE_FIRST_JUNCTION' in result['rejected_reason']
