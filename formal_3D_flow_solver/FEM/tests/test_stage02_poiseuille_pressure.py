from stage02_helpers import read_case,read_report


def test_pressure_gauge_and_multiplier_match_compatible_exact_solution():
    policy=read_report("acceptance_policy")
    r=read_case("pipe_fine_reference","analytic_comparison")
    assert r["lambda_relative_error"]<=policy["analytic_lambda_relative"]
    assert r["pressure_profile_error"]<=policy["analytic_pressure_relative"]
    assert r["pressure_global_L2_relative_error"]<=policy["analytic_pressure_relative"]
    assert len(r["pressure_sections"])>=5
    assert r["section_quadrature"]["radial_gauss_points"]*r["section_quadrature"]["azimuthal_points"]>100
    formal=read_case("pipe_fine_natural","analytic_comparison")
    assert not formal["exact_solution_for_this_boundary_model"]
