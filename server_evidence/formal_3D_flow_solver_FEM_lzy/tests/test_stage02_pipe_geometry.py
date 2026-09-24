import json
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("profile",["pipe_coarse","pipe_medium","pipe_fine"])
def test_synthetic_pipe_is_three_dimensional_si_and_correctly_tagged(profile):
    m=json.loads((ROOT/f"outputs/stage02/meshes/{profile}/metadata.json").read_text())
    assert m["status"]=="PASS" and not m["experimental_condition"] and not m["pde_solved"]
    assert m["units"]=="SI" and m["geometry"]["radius_m"]==5e-6 and m["geometry"]["length_m"]==50e-6
    assert m["boundary_tags"]=={"WALL":1,"OUTLET":2,"INLET":4} and m["cell_tags"]=={"FLUID":100}
    assert m["minimum_jacobian_determinant_m3"]>0 and m["tetrahedron_count"]>0
    assert max(m["area_relative_errors"].values())<1e-3 and m["volume_relative_error"]<1e-3
    assert m["normal_dot_axis"]["INLET"]<-.999999 and m["normal_dot_axis"]["OUTLET"]>.999999
