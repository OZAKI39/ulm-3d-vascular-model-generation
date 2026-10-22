"""Read-only provenance gate. Every scientific input is resolved from manifests."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pyvista as pv

FEM_BRANCH = "sync/fem-simvascular-stage-q-particle-handoff-20260920"
EXPECTED_FACE_IDS = {"WALL": 1, "OUTLET_03": 2, "OUTLET_01": 3, "INLET": 4, "OUTLET_02": 5}


class FrozenIntegrityError(ValueError):
    pass


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def require(condition, message):
    if not condition:
        raise FrozenIntegrityError(message)


def checked_path(root, relative):
    path = (root / relative).resolve()
    require(path.is_relative_to(root), f"path escapes frozen root: {relative}")
    return path


def check_hash(path, expected):
    actual = sha256(path) if path.is_file() else "MISSING"
    require(actual == expected, f"{path}: expected SHA256={expected}; actual SHA256={actual}; frozen input blocker")
    return actual


def tetra_connectivity(grid):
    require(grid.n_cells > 0 and np.all(grid.celltypes == 10), "requires only VTK_TETRA cells")
    raw = grid.cells.reshape(-1, 5)
    require(np.all(raw[:, 0] == 4), "requires four vertices per tetra")
    return np.array(raw[:, 1:], dtype=np.int64)


def validate_pair(mesh, flow):
    """Flow connectivity is checked for correspondence, NEVER used to interpolate."""
    require(np.array_equal(mesh.points, flow.points), "flow point coordinates/order differ from canonical volume")
    tetra = tetra_connectivity(mesh)
    flow_tetra = tetra_connectivity(flow)
    require(tetra.shape == flow_tetra.shape and np.array_equal(np.sort(tetra, axis=1), np.sort(flow_tetra, axis=1)),
            "flow tetra rows do not correspond to canonical volume")
    for name, shape in (("Velocity", (mesh.n_points, 3)), ("Pressure", (mesh.n_points,))):
        require(name in flow.point_data, f"missing POINT {name}")
        array = flow.point_data[name]
        require(array.shape == shape and array.dtype == np.float64 and np.isfinite(array).all(),
                f"invalid POINT {name}: shape/dtype/nonfinite")
    return tetra, flow_tetra


def read_frozen(fem_root):
    root = Path(fem_root).resolve()
    sums = root / "frozen_reference/SHA256SUMS.txt"
    for line in sums.read_text().splitlines():
        digest, relative = line.split("  ", 1)
        check_hash(checked_path(root, relative), digest)
    read = lambda name: json.loads((root / name).read_text())
    mm = read("frozen_reference/mesh_manifest.json")
    fm = read("frozen_reference/flow/flow_field_manifest.json")
    bm = read(mm["boundary_manifest"])
    require(mm["units"] == fm["coordinate_units"] == bm["units"] == "m", "coordinate units must be m")
    require(mm["face_ids"] == EXPECTED_FACE_IDS, "boundary face IDs changed")
    require(fm["mesh"] == mm["path"] and fm["mesh_sha256"] == mm["sha256"], "mesh/flow manifest disagreement")
    for key, name, unit, components in (("velocity", "Velocity", "m/s", 3), ("pressure", "Pressure", "Pa", 1)):
        spec = fm[key]
        require((spec["name"], spec["association"], spec["units"], spec["dtype"], spec["components"]) ==
                (name, "POINT", unit, "float64", components), f"{key} manifest contract changed")
    mesh_path, flow_path = (checked_path(root, m["path"]) for m in (mm, fm))
    check_hash(mesh_path, mm["sha256"])
    check_hash(flow_path, fm["sha256"])
    mesh, flow = pv.read(mesh_path), pv.read(flow_path)
    tetra, flow_tetra = validate_pair(mesh, flow)
    require(mesh.n_points == mm["nodes"] == fm["points"], "point count disagreement")
    require(mesh.n_cells == mm["tetra"] == fm["cells"], "tetra count disagreement")
    require(np.array_equal(mesh.point_data["GlobalNodeID"], np.arange(1, mesh.n_points + 1)), "canonical node IDs")
    require(np.array_equal(mesh.cell_data["GlobalElementID"], np.arange(1, mesh.n_cells + 1)), "canonical cell IDs")
    require(mm["flow_points_identical"] and mm["flow_cell_rows_identical_ignoring_local_vertex_order"], "ordering manifest")
    permutation = mm["flow_local_vertex_permutation_from_volume"]
    require(np.array_equal(flow_tetra, tetra[:, permutation]), "local vertex permutation differs from manifest")
    boundaries = {}
    for role, spec in bm["boundaries"].items():
        check_hash(checked_path(root, spec["path"]), spec["sha256"])
        surface = pv.read(root / spec["path"])
        require(spec["sv_face_id"] == EXPECTED_FACE_IDS[role], f"{role} face ID")
        gids = np.asarray(surface.point_data["GlobalNodeID"]) - 1
        require(np.array_equal(surface.points, mesh.points[gids]), f"{role} coordinates/order")
        require(surface.n_points == spec["nodes"] and surface.n_cells == spec["facets"], f"{role} counts")
        boundaries[role] = surface
    summary = dict(stage="Particle-0", status="PASS", fem_branch=fm.get('source_branch',FEM_BRANCH), fem_root=str(root),
                   case_role=fm.get('case_role','LEGACY_0P352841_MMPS'),
                   mesh_path=str(mesh_path), flow_path=str(flow_path), mesh_sha256=mm["sha256"], flow_sha256=fm["sha256"],
                   node_count=mesh.n_points, tetra_count=mesh.n_cells, coordinates_identical=True,
                   point_order_identical=True, tetra_rows_correspond=True,
                   local_vertex_permutation=permutation,
                   permuted_tetra_count=int(np.count_nonzero(np.any(tetra != flow_tetra, axis=1))),
                   velocity_array="Velocity", pressure_array="Pressure", units=dict(position="m", velocity="m/s", pressure="Pa", derivatives="1/s"),
                   boundary_face_ids=mm["face_ids"], boundary_manifest=bm,
                   sonovue=dict(path="/home/lzy/projects/sonovue_size_distribution_v0", used=False,
                       branch="codex/sonovue-size-distribution-20260920_111258"),
                   old_2d_reference="/home/lzy/projects/ulm_microbubble_traj_gen_2D", copied_2d_code=False)
    return summary, mesh, flow, boundaries
