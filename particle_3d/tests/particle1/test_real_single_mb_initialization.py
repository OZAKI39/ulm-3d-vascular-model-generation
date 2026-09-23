import json
from itertools import combinations
import numpy as np
from particle_3d.particle1_cases import select_inlet_centroid


def test_selection_uses_all_inlet_adjacent_cells_without_trial_tuning(p1_field,p1_audited,repo):
    inlet=p1_audited[3]['INLET']
    result,candidates=select_inlet_centroid(p1_field,inlet)
    ids=np.unique(np.asarray(inlet.cell_data['GlobalElementID'])-1)
    assert result['candidate_count']==len(ids)==len(candidates)
    inward=np.array(result['inward_normal'])
    scores=np.array([p1_field.velocity_nodes_m_s[p1_field.tetra[i]].mean(0)@inward for i in ids])
    winner=int(ids[np.lexsort((ids,-scores))[0]])
    assert result['initial_tetra_id']==winner
    np.testing.assert_array_equal(result['initial_position_m'],p1_field.points_m[p1_field.tetra[winner]].mean(0))
    stored=json.loads((repo/'particle_3d/reports/particle1/data/05_real_initialization.json').read_text())
    assert result==stored
    vertices=p1_field.points_m[p1_field.tetra[winner]]
    shortest=min(np.linalg.norm(vertices[a]-vertices[b]) for a,b in combinations(range(4),2))
    base=.25*shortest/np.linalg.norm(result['initial_velocity_m_s'])
    np.testing.assert_array_equal(result['validation_timesteps_s'],[base,base/2,base/4])
    assert not result['production_particle_timestep_frozen']
