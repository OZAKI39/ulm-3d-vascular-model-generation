import inspect
from stage02_helpers import coefficient_scaling_error,read_case


def test_density_only_changes_reynolds_diagnostic():
    from fem3d.solver import solve_stokes
    assert "rho" not in inspect.signature(solve_stokes).parameters
    assert max(coefficient_scaling_error("rho_double",(1.,1.,1.)).values())==0
    a=read_case("rho_double","analytic_comparison")
    b=read_case("pipe_medium_natural","analytic_comparison")
    assert a["reynolds_number"]==2*b["reynolds_number"]
