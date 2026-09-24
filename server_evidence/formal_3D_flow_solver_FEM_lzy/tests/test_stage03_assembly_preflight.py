from stage03_helpers import *

def test_measured_assembly_precedes_the_only_factorization():
    a=preflight()
    assert a['status']=='PASS' and a['matrix_shape']==[a['total_dofs']]*2
    assert a['total_dofs']==a['velocity_dofs']+a['pressure_dofs']+a['real_dofs']
    assert a['nnz']>a['total_dofs'] and 'Explicit LOCAL' in a['nnz_reduction']
    assert a['assembly_wall_time_s']>0 and a['sum_rank_peak_rss_kib']>0
    assert not a['matrix_factorized'] and not a['ksp_solve_called']
    assert a['timestamp']<read('outputs/stage03/reference/metadata/resource_gate.json')['fresh_memory']['timestamp']

def test_successful_formal_solution_uses_unchanged_direct_solver():
    s=solve()
    assert s['reference_pde_solve_count']==1 and s['factorization_status']=='SUCCESS'
    assert s['solver']['mpi_ranks']==4 and s['solver']['factor_solver']=='mumps'
    assert s['solver']['ksp_type']=='preonly' and s['solver']['pc_type']=='lu'


def test_actual_singular_matrix_has_verified_pressure_null_vectors():
    diag=read('outputs/stage03/reference/qc/singularity_diagnosis.json')
    failure=read('outputs/stage03/reference/metadata/failure.json')
    assert diag['diagnostic_status']=='CONFIRMED' and failure['status']=='FAIL'
    assert diag['zero_row_count']==diag['independent_exact_null_modes']==2
    assert all(row['actual_nonzero_entries']==0 and row['unit_pressure_mode_matrix_product_norm']==0 for row in diag['zero_rows'])
    assert not diag['matrix_factorized'] and not diag['ksp_solve_called'] and not diag['matrix_altered_for_repair']
    assert failure['formal_factorization_attempts']==1 and not failure['resource_blocked']
    assert not failure['pressure_pins_added'] and not failure['geometry_changed'] and not failure['formal_retry_performed']
