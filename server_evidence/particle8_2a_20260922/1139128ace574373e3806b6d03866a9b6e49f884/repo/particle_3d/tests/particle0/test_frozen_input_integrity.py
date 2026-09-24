import numpy as np
import pytest
from particle_3d.audit import check_hash, FrozenIntegrityError, sha256


def test_scientific_hashes_counts_units_order_and_boundary_ids(audited):
    summary, mesh, flow, boundaries = audited
    assert summary["status"] == "PASS"
    assert summary["node_count"] == 70363
    assert summary["tetra_count"] == 371402
    assert summary["local_vertex_permutation"] == [1, 0, 2, 3]
    assert summary["permuted_tetra_count"] == summary["tetra_count"]
    assert summary["units"] == dict(position="m", velocity="m/s", pressure="Pa", derivatives="1/s")
    assert set(boundaries) == {"WALL", "INLET", "OUTLET_01", "OUTLET_02", "OUTLET_03"}
    assert np.array_equal(mesh.points, flow.points)


def test_hash_mismatch_reports_file_expected_actual(tmp_path):
    path = tmp_path / "damaged.vtu"
    path.write_bytes(b"not a frozen file")
    with pytest.raises(FrozenIntegrityError) as error:
        check_hash(path, "0" * 64)
    message = str(error.value)
    assert str(path) in message and "0" * 64 in message and sha256(path) in message


def test_boundary_facets_are_complete_disjoint_and_outward(audited, frozen_root):
    import pyvista as pv
    summary, mesh, _, boundaries = audited
    tetra = mesh.cells.reshape(-1, 5)[:, 1:]
    keys = []
    for role, surface in boundaries.items():
        triangles = surface.faces.reshape(-1, 4)[:, 1:]
        gids = np.asarray(surface.point_data["GlobalNodeID"]) - 1
        global_triangles = gids[triangles]
        owners = np.asarray(surface.cell_data["GlobalElementID"]) - 1
        assert np.all(np.any(global_triangles[:, :, None] == tetra[owners, None, :], axis=2))
        xyz = surface.points[triangles]
        normal = np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0])
        assert np.all(np.einsum("ij,ij->i", normal, xyz.mean(1) - mesh.points[tetra[owners]].mean(1)) > 0)
        keys.extend(map(tuple, np.sort(global_triangles, axis=1).tolist()))
    exterior = pv.read(frozen_root / summary["boundary_manifest"]["exterior"]["path"])
    exterior_triangles = (np.asarray(exterior.point_data["GlobalNodeID"]) - 1)[exterior.faces.reshape(-1, 4)[:, 1:]]
    assert len(keys) == len(set(keys)) == exterior.n_cells
    assert set(keys) == set(map(tuple, np.sort(exterior_triangles, axis=1).tolist()))
    for role, face_id in summary["boundary_face_ids"].items():
        assert np.count_nonzero(exterior.cell_data["ModelFaceID"] == face_id) == boundaries[role].n_cells
