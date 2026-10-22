from particle_3d.particle2_cases import NORM_ATOL


def test_rigid_rotation_follows_exact_fluid_rotation_for_multiple_geometries(rotation_case_result):
    rows,summary=rotation_case_result
    assert summary["passed"] and summary["geometry_count"]==5
    assert summary["max_omega_error_s_inv"]==0
    assert max(r["axis_vector_error"] for r in rows)<=summary["axis_error_bound"]
    assert max(r["quaternion_norm_error"] for r in rows)<=NORM_ATOL
