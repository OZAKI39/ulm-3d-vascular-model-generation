from particle_3d.particle2_cases import NORM_ATOL


def test_static_orientation_for_all_five_distribution_shapes(static_case_result):
    rows,summary=static_case_result
    assert summary["passed"] and summary["geometry_count"]==5
    assert all(row["rotation_matrix_error"]<=NORM_ATOL and row["omega_norm_s_inv"]==0 for row in rows)
