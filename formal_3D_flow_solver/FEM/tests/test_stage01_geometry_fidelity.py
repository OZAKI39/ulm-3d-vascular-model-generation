import numpy as np
import pytest
from fem3d.mesh_qc import fidelity


def test_boundary_geometry_and_ports_are_frozen(stage01_mesh):
    q=stage01_mesh["qc"]
    f=q["boundary_fidelity"]
    assert f["matched_tagged_source_triangles"]==67262
    assert f["coordinates_exactly_equal"] and f["tagged_triangle_connectivity_equal"]
    assert f["maximum_boundary_displacement_m"]==0
    # Reduction order/NumPy versions can differ by one ulp with identical triangles.
    assert f["source_surface_area_m2"]==pytest.approx(f["final_surface_area_m2"],rel=1e-14,abs=0)
    assert f["source_bounds_m"]==f["final_boundary_bounds_m"]
    for p in q["ports"].values():
        assert p["relative_error"]<1e-11
        assert p["normal_dot_product"]>.999
        assert p["centroid_displacement_m"]<f["roundoff_tolerance_m"]
        assert p["max_plane_deviation_m"]<=p["source_max_plane_deviation_m"]+f["roundoff_tolerance_m"]


@pytest.mark.parametrize("fault",["displacement","tag_swap","triangle_loss"])
def test_changed_geometry_or_semantics_fails_even_with_similar_bounds(stage01_mesh,fault):
    src=stage01_mesh["source"]
    points=src["points_m"].copy()
    triangles=src["triangles"].copy()
    tags=src["facet_tags"].copy()
    if fault=="displacement": points[len(points)//2,0]+=1e-10
    elif fault=="tag_swap": tags[tags==4]=3
    else: triangles,tags=triangles[:-1],tags[:-1]
    with pytest.raises(ValueError,match="Frozen"):
        fidelity(points,triangles,tags,src["points_m"],src["triangles"],src["facet_tags"],1e-17)
