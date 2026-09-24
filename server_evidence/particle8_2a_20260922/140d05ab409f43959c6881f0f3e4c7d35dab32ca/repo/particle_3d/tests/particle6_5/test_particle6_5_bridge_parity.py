def test_actual_lammps_storage_neighbors_and_solver(bridge_result):
 d=bridge_result;assert d['status']=='PASS' and all(r['exact_equal'] for r in d['errors'])
 for r in d['matrix_rows']:
  assert all(r[k]==0 for k in ['R','b','U','J']) and r['constraints_equal']
  assert r['standalone_eligible_pairs']==r['filtered_nearfield_pairs']
  assert len(r['lammps_candidates'])>len(r['filtered_nearfield_pairs'])
 assert all(r['max_lammps_force']==0 for r in d['force_audits'])

def test_actual_binary_restart_continues(bridge_result):
 assert bridge_result['restart_error']['exact_equal']
 assert bridge_result['restart_at_step']==3 and bridge_result['total_steps']==6
