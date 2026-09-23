from particle_3d.particle6_validation import rebuild_stress

def test_neighbor_rebuild_after_position_update():
    d=rebuild_stress();assert sum(x['previous_stale'] for x in d['rows'])>=4
    assert all(x['standalone']==x['lammps'] and x['mismatch_count']==0 for x in d['rows'])
