import ast
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def test_no_pressure_boundary_dofs_or_constant_nullspace():
    for profile in ("pipe_coarse","pipe_medium","pipe_fine"):
        result=json.loads((ROOT/f"outputs/stage02/cases/{profile}_natural/qc/solver.json").read_text())
        bc=result["boundary_conditions"]
        assert bc["pressure_dirichlet_dofs"]==0 and bc["pressure_nullspace_registered"] is False
        assert bc["inlet_velocity_profile"] is None and bc["outlet_velocity_profile"] is None
        assert result["converged_reason"]>0
    source=(ROOT/"src/fem3d/solver.py").read_text()+(ROOT/"src/fem3d/boundary.py").read_text()
    tree=ast.parse(source)
    assert not any(isinstance(n,ast.Attribute) and n.attr in ("setNullSpace","setNearNullSpace") for n in ast.walk(tree))


def test_formal_api_rejects_nonpositive_Q_before_assembly():
    import pytest
    from fem3d.solver import solve_stokes
    for q in (0,-1e-14,float("nan")):
        with pytest.raises(ValueError,match="positive inflow"):
            solve_stokes(None,None,.003,q,5e-6,None)
