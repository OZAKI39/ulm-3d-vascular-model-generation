from stage02_helpers import read_case,read_report


def test_full_profiles_at_three_sections_match_compatible_exact_solution():
    policy=read_report("acceptance_policy")
    result=read_case("pipe_fine_reference","analytic_comparison")
    assert result["exact_solution_for_this_boundary_model"]
    assert {s["s_over_L"] for s in result["velocity_sections"]}=={.25,.5,.75}
    assert result["velocity_L2_relative_error"]<=policy["analytic_velocity_L2_relative"]
    assert result["velocity_global_L2_relative_error"]<=policy["analytic_velocity_L2_relative"]
    assert max(s["maximum_transverse_velocity_over_Umean"] for s in result["velocity_sections"])<policy["transverse_over_Umean"]
    assert len(result["velocity_profiles"][0]["r_over_R"])>50
