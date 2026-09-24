from itertools import permutations
from dataclasses import fields
import numpy as np
import pyvista as pv
import pytest
from particle_3d.field import FrozenFEMField, FlowSample
from particle_3d.audit import FrozenIntegrityError


def grid(points, tetra):
    cells = np.column_stack((np.full(len(tetra), 4), tetra))
    return pv.UnstructuredGrid(cells.ravel(), np.full(len(tetra), 10, np.uint8), points)


def test_all_24_flow_vertex_permutations_leave_canonical_result(affine):
    f = affine["field"]
    mesh = grid(f.points_m, f.tetra)
    for permutation in permutations(range(4)):
        flow = grid(f.points_m, f.tetra[:, permutation])
        flow.point_data["Velocity"] = f.velocity_nodes_m_s
        flow.point_data["Pressure"] = f.pressure_nodes_pa
        actual = FrozenFEMField.from_grids(mesh, flow).sample_many(affine["positions"])
        for attribute in fields(FlowSample):
            np.testing.assert_array_equal(getattr(actual, attribute.name), getattr(affine["batch"], attribute.name))


def test_real_flow_connectivity_cannot_determine_interpolation(audited, real_field, faces):
    _, mesh, flow, _ = audited
    changed_flow = grid(flow.points, real_field.tetra[:, [3, 2, 0, 1]])
    changed_flow.point_data["Velocity"] = flow.point_data["Velocity"]
    changed_flow.point_data["Pressure"] = flow.point_data["Pressure"]
    other = FrozenFEMField.from_grids(mesh, changed_flow)
    points = np.vstack(([f["center"] for f in faces], mesh.points[[0, 100, 1000]], mesh.points[real_field.tetra[:10]].mean(1)))
    for p in points:
        a, b = real_field.sample(p), other.sample(p)
        for attribute in fields(FlowSample):
            np.testing.assert_array_equal(getattr(a, attribute.name), getattr(b, attribute.name))


def test_point_order_and_tetra_row_mismatch_rejected(audited):
    _, mesh, flow, _ = audited
    reordered = flow.copy(deep=True)
    points = reordered.points.copy(); points[[0, 1]] = points[[1, 0]]; reordered.points = points
    with pytest.raises(FrozenIntegrityError, match="coordinates/order"):
        FrozenFEMField.from_grids(mesh, reordered)
    tetra = mesh.cells.reshape(-1, 5)[:, 1:].copy(); tetra[[0, 1]] = tetra[[1, 0]]
    with pytest.raises(FrozenIntegrityError, match="tetra rows"):
        FrozenFEMField.from_grids(mesh, grid(mesh.points, tetra))
