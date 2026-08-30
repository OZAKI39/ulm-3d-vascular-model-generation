"""Run the read-only repaired-BASE topology and unit-scaling investigation."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pyvista as pv


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from utils.cfd_flow.io import sha256_file, write_json  # noqa: E402
from utils.cfd_flow.repaired_topology_forensics import (  # noqa: E402
    CONNECTIVITY_LABELS,
    PINNED_SDR_SHA,
    PINNED_SEEDER_SHA,
    PINNED_TREELM_SHA,
    PORT_LABELS,
    SCALE_FACTOR,
    cell_centers,
    classify_component_topology,
    compare_scaled_meshes,
    component_identity_records,
    coordinate_set_difference,
    difference_cluster_records,
    intersect_raytriangle_callgraph,
    load_boundary_ids_by_cell,
    nearest_cell_membership,
    nearest_centerline_segment,
    nearest_component_gap,
    parse_seed_point,
    parse_uniform_lattice,
    port_component_membership,
    scale_binary_stl,
    scale_seeder_lua_geometry,
    sparse_component_labels,
    timed_runtime,
    unit_scaling_oracle_decision,
)
from utils.cfd_flow.musubi_boundary_mass_referee import (  # noqa: E402
    MeshContract,
    load_mesh_contract,
)
from utils.cfd_flow.restart_decode import tree_ids_to_ijk  # noqa: E402


RESEARCH_RUN = (
    "healthy_mouse_capillary_port_grid_sensitivity_research_anchor003274_20260830"
)
OLD_RUN = "axis_aligned_ideal_plane_inlet_preflight_anchor003274_20260829_120444"
REPAIRED_RUN = (
    "healthy_mouse_capillary_qvalue_repaired_base_preflight_anchor003274_20260830"
)
SCALED_SUBDIR = "qvalue_repair/unpatched_micrometer_scaled_base"


def _paths(root: Path) -> dict[str, Path]:
    research = root / "outputs/cfd_flow" / RESEARCH_RUN
    return {
        "root": root,
        "research": research,
        "qc": research / "qc",
        "visualization": research / "qvalue_repair/topology_forensics",
        "old": root / "outputs/cfd_flow" / OLD_RUN,
        "repaired": root / "outputs/cfd_flow" / REPAIRED_RUN,
        "scaled": research / SCALED_SUBDIR,
        "surface": root
        / "outputs/cfd_flow/axis_aligned_inlet_geometry_anchor003274_20260829_111451"
        / "geometry/cfd_surface_axis_aligned_inlet_m.stl",
        "transform": root
        / "outputs/cfd_flow/axis_aligned_inlet_geometry_anchor003274_20260829_111451"
        / "transform/anatomical_to_cfd_transform.json",
        "nodes": root
        / "outputs/cfd_preprocess/global_to_roi_anchor003274_20260825_183628"
        / "global_1d/nodes.csv",
        "edges": root
        / "outputs/cfd_preprocess/global_to_roi_anchor003274_20260825_183628"
        / "global_1d/edges.csv",
        "seeder_source": Path(
            "//wsl.localhost/Ubuntu/home/lzy/apes-pinned/seeder_official"
        ),
    }


def _head(root: Path) -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
    ).stdout.strip()


def _load_mesh(path: Path) -> tuple[MeshContract, Any]:
    mesh = path / "seeder/mesh"
    contract = load_mesh_contract(
        mesh, allow_zero_normals=True, require_runtime_order=False
    )
    lattice = parse_uniform_lattice((mesh / "header.lua").read_text(encoding="utf-8"))
    return contract, lattice


def _connectivity_records(contract: MeshContract, lattice: Any) -> tuple[dict, dict]:
    labels: dict[int, np.ndarray] = {}
    records: dict[str, Any] = {}
    for mode in CONNECTIVITY_LABELS:
        mode_labels, sizes = sparse_component_labels(
            contract.cell_ijk,
            mode,
            cells_per_axis=2**lattice.level,
        )
        labels[mode] = mode_labels
        records[str(mode)] = {
            "definition": CONNECTIVITY_LABELS[mode],
            "component_count": int(len(sizes)),
            "component_sizes": sizes.tolist(),
        }
        if mode == 6:
            records[str(mode)]["components"] = component_identity_records(
                contract.cell_ijk, mode_labels, sizes, lattice
            )
    return labels, records


def _surface_membership(
    surface_path: Path, points: np.ndarray, *, maximum_points: int = 4096
) -> dict[str, Any]:
    surface = pv.read(surface_path)
    connected = surface.connectivity()
    region_ids = np.asarray(connected.cell_data["RegionId"], dtype=np.int64)
    component_count = int(len(np.unique(region_ids)))
    if len(points):
        sample_count = min(len(points), int(maximum_points))
        sample_indices = np.linspace(
            0, len(points) - 1, num=sample_count, dtype=np.int64
        )
        sample = np.asarray(points)[sample_indices]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            selected = pv.PolyData(sample).select_enclosed_points(
                surface, tolerance=0.0, check_surface=True
            )
        inside = np.asarray(selected.point_data["SelectedPoints"], dtype=bool)
    else:
        inside = np.empty(0, dtype=bool)
    return {
        "surface_path": str(surface_path.resolve()),
        "surface_sha256": sha256_file(surface_path),
        "is_watertight": bool(surface.n_open_edges == 0),
        "surface_component_count": component_count,
        "population_points": int(len(points)),
        "points_tested": int(len(inside)),
        "sampling": "DETERMINISTIC_EVEN_INDEX_SAMPLE_MAX_4096",
        "inside_count": int(np.count_nonzero(inside)),
        "outside_count": int(np.count_nonzero(~inside)),
        "inside_fraction": float(np.mean(inside)) if len(inside) else None,
        "continuous_geometry_connected": bool(
            surface.n_open_edges == 0 and component_count == 1
        ),
    }


def run_zero(project_root: Path) -> dict[str, Any]:
    started = time.perf_counter()
    paths = _paths(project_root.resolve())
    paths["qc"].mkdir(parents=True, exist_ok=True)
    paths["visualization"].mkdir(parents=True, exist_ok=True)
    old, old_lattice = _load_mesh(paths["old"])
    repaired, repaired_lattice = _load_mesh(paths["repaired"])
    if not np.array_equal(
        old_lattice.origin, repaired_lattice.origin
    ) or not np.isclose(old_lattice.dx, repaired_lattice.dx, rtol=0.0, atol=0.0):
        raise ValueError("OLD and repaired BASE do not share an exact lattice")

    labels, connectivity = _connectivity_records(repaired, repaired_lattice)
    centers = cell_centers(repaired.cell_ijk, repaired_lattice)
    seed_point = parse_seed_point(
        (paths["repaired"] / "seeder/seeder.lua").read_text(encoding="utf-8")
    )
    seed = {
        "configured_seed_point_m": seed_point.tolist(),
        **nearest_cell_membership(seed_point, centers, labels[6]),
    }
    boundary_cells = {
        name: repaired.boundaries[name].cell_indices for name in PORT_LABELS
    }
    ports = port_component_membership(
        labels[6], boundary_cells, connectivity["6"]["component_count"]
    )
    gap = nearest_component_gap(
        repaired.cell_ijk,
        labels[6],
        np.asarray(connectivity["6"]["component_sizes"]),
        repaired_lattice,
    )
    classification = classify_component_topology(
        connectivity["6"]["component_sizes"],
        connectivity["18"]["component_sizes"],
        ports,
    )
    component_result = {
        "status": "PASS_ZERO_RUN_COMPONENT_FORENSICS",
        "actual_head_at_execution": _head(paths["root"]),
        "production_pipeline_modified": False,
        "seeder_calls": 0,
        "musubi_calls": 0,
        "fluid_cell_count": int(len(repaired.tree_ids)),
        "dx_m": repaired_lattice.dx,
        "connectivity": connectivity,
        "seed": seed,
        "ports": ports,
        "nearest_gap_between_two_largest_face_components": gap,
        "component_classification": classification,
        "face_disconnected_but_d3q19_connected": bool(
            connectivity["6"]["component_count"] > 1
            and connectivity["18"]["component_count"] == 1
        ),
        "topology_underresolved": classification == "D3Q19_DIAGONAL_ONLY_CONNECTION",
        "runtime_s": timed_runtime(started),
    }
    write_json(paths["qc"] / "repaired_base_component_forensics.json", component_result)

    differences = coordinate_set_difference(old.tree_ids, repaired.tree_ids)
    old_only_ijk = tree_ids_to_ijk(differences["old_only"])
    new_only_ijk = tree_ids_to_ijk(differences["new_only"])
    port_centers = {
        name: centers[repaired.boundaries[name].cell_indices] for name in PORT_LABELS
    }
    gap_centers = (
        np.asarray(
            [gap["component_a_center_m"], gap["component_b_center_m"]],
            dtype=float,
        )
        if gap is not None
        else np.empty((0, 3), dtype=float)
    )
    old_only_labels, old_clusters = difference_cluster_records(
        old_only_ijk,
        old_lattice,
        port_centers=port_centers,
        gap_centers=gap_centers,
    )
    new_only_labels, new_clusters = difference_cluster_records(
        new_only_ijk,
        repaired_lattice,
        port_centers=port_centers,
        gap_centers=gap_centers,
    )
    old_only_xyz = cell_centers(old_only_ijk, old_lattice)
    new_only_xyz = cell_centers(new_only_ijk, repaired_lattice)
    old_surface = _surface_membership(paths["surface"], old_only_xyz)
    new_surface = _surface_membership(paths["surface"], new_only_xyz)
    old_false_flooding = (
        "YES" if old_surface["outside_count"] > old_surface["inside_count"] else "NO"
    )
    spatial_pattern = (
        "CONCENTRATED_IN_LARGEST_BRANCH_REGIONS"
        if old_clusters and old_clusters[0]["fraction_of_difference_cells"] >= 0.5
        else "DISTRIBUTED_ACROSS_MULTIPLE_REGIONS"
    )
    npz_path = paths["visualization"] / "topology_point_clouds.npz"
    component0 = centers[labels[6] == 0]
    component1 = (
        centers[labels[6] == 1]
        if connectivity["6"]["component_count"] > 1
        else np.empty((0, 3))
    )
    np.savez_compressed(
        npz_path,
        old_only_cells_m=old_only_xyz,
        new_only_cells_m=new_only_xyz,
        component_0_cells_m=component0,
        component_1_cells_m=component1,
        nearest_gap_cells_m=gap_centers,
        old_only_cluster_id=old_only_labels,
        new_only_cluster_id=new_only_labels,
    )
    spatial = {
        "status": "PASS_ZERO_RUN_SPATIAL_DIFFERENCE",
        "old_cell_count": int(len(old.tree_ids)),
        "repaired_cell_count": int(len(repaired.tree_ids)),
        "common_cell_count": int(len(differences["common"])),
        "old_only_cell_count": int(len(differences["old_only"])),
        "new_only_cell_count": int(len(differences["new_only"])),
        "net_cell_difference": int(len(repaired.tree_ids) - len(old.tree_ids)),
        "old_only_largest_clusters": old_clusters,
        "new_only_largest_clusters": new_clusters,
        "spatial_pattern": spatial_pattern,
        "old_only_continuous_lumen_membership": old_surface,
        "new_only_continuous_lumen_membership": new_surface,
        "old_mesh_false_flooding": old_false_flooding,
        "visualization_npz": {
            "path": str(npz_path.resolve()),
            "sha256": sha256_file(npz_path),
            "arrays": [
                "old_only_cells_m",
                "new_only_cells_m",
                "component_0_cells_m",
                "component_1_cells_m",
                "nearest_gap_cells_m",
            ],
        },
        "runtime_s": timed_runtime(started),
    }
    write_json(paths["qc"] / "old_vs_repaired_spatial_difference.json", spatial)

    callgraph = intersect_raytriangle_callgraph(paths["seeder_source"])
    callgraph["actual_head_at_execution"] = _head(paths["root"])
    callgraph["seeder_calls"] = 0
    callgraph["musubi_calls"] = 0
    write_json(paths["qc"] / "intersect_raytriangle_callgraph.json", callgraph)

    if gap is not None:
        gap_midpoint = 0.5 * (
            np.asarray(gap["component_a_center_m"])
            + np.asarray(gap["component_b_center_m"])
        )
        centerline = nearest_centerline_segment(
            gap_midpoint,
            transform_json=paths["transform"],
            nodes_csv=paths["nodes"],
            edges_csv=paths["edges"],
        )
    else:
        centerline = None
    zero_summary = {
        "status": "PASS_ZERO_RUN_FORENSICS_COMPLETE",
        "component_classification": classification,
        "component_evidence": str(
            (paths["qc"] / "repaired_base_component_forensics.json").resolve()
        ),
        "spatial_difference_evidence": str(
            (paths["qc"] / "old_vs_repaired_spatial_difference.json").resolve()
        ),
        "callgraph_evidence": str(
            (paths["qc"] / "intersect_raytriangle_callgraph.json").resolve()
        ),
        "continuous_surface": old_surface,
        "nearest_centerline_segment": centerline,
        "runtime_s": timed_runtime(started),
    }
    write_json(paths["visualization"] / "zero_run_summary.json", zero_summary)
    return zero_summary


def prepare_scaled(project_root: Path) -> dict[str, Any]:
    started = time.perf_counter()
    paths = _paths(project_root.resolve())
    source_geometry = paths["repaired"] / "geometry"
    destination_geometry = paths["scaled"] / "geometry"
    scaled_files = []
    for source in sorted(source_geometry.rglob("*.stl")):
        relative = source.relative_to(source_geometry)
        scaled_files.append(
            scale_binary_stl(
                source,
                destination_geometry / relative,
                factor=SCALE_FACTOR,
            )
        )
    source_lua = paths["repaired"] / "seeder/seeder.lua"
    scaled_lua = paths["scaled"] / "seeder/seeder.lua"
    scaled_lua.parent.mkdir(parents=True, exist_ok=True)
    scaled_lua.write_text(
        scale_seeder_lua_geometry(source_lua.read_text(encoding="utf-8")),
        encoding="utf-8",
    )
    result = {
        "status": "PASS_SCALED_GEOMETRY_PREPARED",
        "scale_factor": SCALE_FACTOR,
        "coordinate_units_before": "m",
        "coordinate_units_after": "micrometer_numeric_units",
        "minlevel_changed": False,
        "relative_geometry_changed": False,
        "source_seeder_lua": str(source_lua.resolve()),
        "scaled_seeder_lua": str(scaled_lua.resolve()),
        "scaled_seeder_lua_sha256": sha256_file(scaled_lua),
        "stl_files": scaled_files,
        "pinned_unpatched_seeder_sha": PINNED_SEEDER_SHA,
        "pinned_unpatched_treelm_sha": PINNED_TREELM_SHA,
        "pinned_unpatched_sdr_sha": PINNED_SDR_SHA,
        "seeder_calls": 0,
        "musubi_calls": 0,
        "runtime_s": timed_runtime(started),
    }
    write_json(paths["scaled"] / "scaled_geometry_transform.json", result)
    return result


def finalize(project_root: Path) -> dict[str, Any]:
    started = time.perf_counter()
    paths = _paths(project_root.resolve())
    repaired, repaired_lattice = _load_mesh(paths["repaired"])
    scaled, scaled_lattice = _load_mesh(paths["scaled"])
    repaired_labels, repaired_connectivity = _connectivity_records(
        repaired, repaired_lattice
    )
    _, scaled_connectivity = _connectivity_records(scaled, scaled_lattice)
    repaired_tree_ids, repaired_bnd = load_boundary_ids_by_cell(
        paths["repaired"] / "seeder/mesh"
    )
    scaled_tree_ids, scaled_bnd = load_boundary_ids_by_cell(
        paths["scaled"] / "seeder/mesh"
    )
    if not np.array_equal(repaired_tree_ids, repaired.tree_ids):
        raise ValueError("Repaired boundary/tree-ID ordering mismatch")
    if not np.array_equal(scaled_tree_ids, scaled.tree_ids):
        raise ValueError("Scaled boundary/tree-ID ordering mismatch")
    comparison = compare_scaled_meshes(repaired, scaled, repaired_bnd, scaled_bnd)
    repaired_counts = {
        mode: row["component_count"] for mode, row in repaired_connectivity.items()
    }
    scaled_counts = {
        mode: row["component_count"] for mode, row in scaled_connectivity.items()
    }
    decision = unit_scaling_oracle_decision(comparison, repaired_counts, scaled_counts)
    scaled_boundary_counts = {
        name: int(len(boundary.cell_indices))
        for name, boundary in scaled.boundaries.items()
    }
    patched_boundary_counts = {
        name: int(len(boundary.cell_indices))
        for name, boundary in repaired.boundaries.items()
    }
    oracle = {
        "status": decision,
        "interpretation": (
            "SCALE_INVARIANCE_ORACLE_PASS"
            if decision == "PASS"
            else "SCALE_INVARIANCE_ORACLE_FAIL"
        ),
        "patched_meter": {
            "cell_count": int(len(repaired.tree_ids)),
            "component_counts": repaired_counts,
            "boundary_counts": patched_boundary_counts,
        },
        "unpatched_micrometer_scaled": {
            "cell_count": int(len(scaled.tree_ids)),
            "component_counts": scaled_counts,
            "boundary_counts": scaled_boundary_counts,
            "pinned_seeder_sha": PINNED_SEEDER_SHA,
        },
        "comparison": comparison,
        "seeder_calls": 1,
        "musubi_calls": 0,
        "runtime_s": timed_runtime(started),
    }
    write_json(paths["qc"] / "unit_scaling_oracle.json", oracle)

    components = json.loads(
        (paths["qc"] / "repaired_base_component_forensics.json").read_text(
            encoding="utf-8"
        )
    )
    spatial = json.loads(
        (paths["qc"] / "old_vs_repaired_spatial_difference.json").read_text(
            encoding="utf-8"
        )
    )
    callgraph = json.loads(
        (paths["qc"] / "intersect_raytriangle_callgraph.json").read_text(
            encoding="utf-8"
        )
    )
    zero_summary = json.loads(
        (paths["visualization"] / "zero_run_summary.json").read_text(encoding="utf-8")
    )
    classification = components["component_classification"]
    continuous_connected = bool(
        zero_summary["continuous_surface"]["continuous_geometry_connected"]
    )
    if decision == "FAIL":
        final_status = "CFD_FLOW_QVALUE_PATCH_TOPOLOGY_NOT_PROVEN"
        root_cause = "UNIT_SCALED_UNPATCHED_MESH_DID_NOT_MATCH_PATCHED_METER_MESH"
        next_step = "INVESTIGATE OTHER DIMENSIONAL PREDICATES AND GEOMETRY TRANSFORMS"
    elif classification == "D3Q19_DIAGONAL_ONLY_CONNECTION":
        final_status = "CFD_FLOW_REPAIRED_BASE_DIAGONAL_ONLY_CONNECTION"
        root_cause = "BASE_FACE_CONNECTIVITY_LOST_AT_A_D3Q19_DIAGONAL_ONLY_NECK"
        next_step = "PREDICT THROAT CELLS BEFORE ONE REPAIRED FINE SEEDER"
    elif classification == "PORTLESS_ISOLATED_POCKET":
        final_status = "CFD_FLOW_REPAIRED_BASE_PORTLESS_POCKET"
        root_cause = "SCALE_CORRECT_REPAIRED_BASE_CONTAINS_A_PORTLESS_ISOLATED_POCKET"
        next_step = "TRACE WHY FLOOD_BOUNDARY_OUTPUTS_THE_PORTLESS_POCKET"
    elif not continuous_connected:
        final_status = "CFD_FLOW_CONTINUOUS_GEOMETRY_DISCONNECTED"
        root_cause = "CONTINUOUS_WATERTIGHT_SURFACE_IS_NOT_A_SINGLE_CONNECTED_COMPONENT"
        next_step = "RETURN TO UPSTREAM GEOMETRY GENERATION"
    elif classification in {
        "MAJOR_NETWORK_SPLIT",
        "PORT_BEARING_SECONDARY_COMPONENT",
    }:
        final_status = "CFD_FLOW_BASE_RESOLUTION_TOPOLOGY_INSUFFICIENT"
        root_cause = (
            "CONTINUOUS_GEOMETRY_CONNECTED_BUT_BASE_FACE_GRID_SPLITS_THE_NETWORK"
        )
        next_step = "PREDICT THROAT CELLS BEFORE ONE REPAIRED FINE SEEDER"
    else:
        final_status = "CFD_FLOW_REPAIRED_TOPOLOGY_ROOT_CAUSE_UNRESOLVED"
        root_cause = (
            "SCALE_ORACLE_PASSED_BUT_COMPONENT_CLASSIFICATION_REMAINS_UNRESOLVED"
        )
        next_step = "REFINE CONTINUOUS_GEOMETRY_NECK_FORENSICS WITHOUT CFD"
    final = {
        "status": final_status,
        "actual_head_at_execution": _head(paths["root"]),
        "production_pipeline_modified": False,
        "seeder_calls": 1,
        "musubi_calls": 0,
        "harvester_calls": 0,
        "old_base_cells": 221309,
        "repaired_base_cells": int(len(repaired.tree_ids)),
        "scaled_unpatched_base_cells": int(len(scaled.tree_ids)),
        "component_classification": classification,
        "scale_invariance_oracle": decision,
        "old_mesh_false_flooding": spatial["old_mesh_false_flooding"],
        "continuous_geometry_connected": continuous_connected,
        "patch_affects_flooding_or_topology_source_proven": bool(
            callgraph["patch_affects_flooding_or_topology"]
        ),
        "break_location": zero_summary["nearest_centerline_segment"],
        "root_cause_final": root_cause,
        "first_failure": (
            oracle["comparison"]
            if decision == "FAIL"
            else components["nearest_gap_between_two_largest_face_components"]
        ),
        "evidence": {
            "component_forensics": "repaired_base_component_forensics.json",
            "spatial_difference": "old_vs_repaired_spatial_difference.json",
            "callgraph": "intersect_raytriangle_callgraph.json",
            "unit_scaling_oracle": "unit_scaling_oracle.json",
        },
        "runtime_s": timed_runtime(started),
        "next": next_step,
    }
    write_json(paths["qc"] / "repaired_topology_root_cause.json", final)
    return final


def finalize_failed_launch(project_root: Path) -> dict[str, Any]:
    """Record the exhausted one-call oracle without inventing scaled results."""

    started = time.perf_counter()
    paths = _paths(project_root.resolve())
    components = json.loads(
        (paths["qc"] / "repaired_base_component_forensics.json").read_text(
            encoding="utf-8"
        )
    )
    spatial = json.loads(
        (paths["qc"] / "old_vs_repaired_spatial_difference.json").read_text(
            encoding="utf-8"
        )
    )
    callgraph = json.loads(
        (paths["qc"] / "intersect_raytriangle_callgraph.json").read_text(
            encoding="utf-8"
        )
    )
    zero_summary = json.loads(
        (paths["visualization"] / "zero_run_summary.json").read_text(encoding="utf-8")
    )
    stdout = paths["scaled"] / "seeder/seeder_stdout.log"
    stderr = paths["scaled"] / "seeder/seeder_stderr.log"
    transform = paths["scaled"] / "scaled_geometry_transform.json"
    first_failure = (
        "The only permitted Seeder invocation did not receive its seeder.lua "
        "working-directory argument through PowerShell Start-Process -> WSL "
        "bash -lc. Seeder printed 'Cannot load configuration file: cannot open "
        "seeder.lua' and returned 0 without generating mesh/."
    )
    attempt = {
        "status": "FAIL_CONFIGURATION_NOT_LOADED",
        "seeder_calls": 1,
        "hard_max": 1,
        "retry_permitted": False,
        "binary": "/home/lzy/apes-pinned/seeder_official/build/seeder",
        "binary_sha256": (
            "178d01f153d01df49cbc16e3f6be2f98ebcc19922bf92dc5afd43c49c8a5e511"
        ),
        "pinned_seeder_sha": PINNED_SEEDER_SHA,
        "wall_time_s": 3.0271732,
        "process_exit_code": 0,
        "semantic_success": False,
        "mesh_generated": False,
        "stdout": str(stdout.resolve()),
        "stdout_sha256": sha256_file(stdout),
        "stderr": str(stderr.resolve()),
        "stderr_sha256": sha256_file(stderr),
        "first_failure": first_failure,
    }
    write_json(paths["scaled"] / "seeder_attempt_summary.json", attempt)
    oracle = {
        "status": "FAIL",
        "interpretation": "SCALE_INVARIANCE_ORACLE_NOT_PROVEN_NO_SCALED_MESH",
        "patched_meter": {
            "cell_count": int(components["fluid_cell_count"]),
            "component_counts": {
                mode: row["component_count"]
                for mode, row in components["connectivity"].items()
            },
        },
        "unpatched_micrometer_scaled": {
            "cell_count": None,
            "component_counts": None,
            "boundary_counts": None,
            "pinned_seeder_sha": PINNED_SEEDER_SHA,
        },
        "comparison": {
            "tree_id_set_exact_match": None,
            "boundary_id_exact_match": None,
            "common_q_links": None,
            "q_patched_minus_scaled": {
                "rms": None,
                "median": None,
                "p95": None,
                "max": None,
            },
        },
        "scaled_geometry_transform": str(transform.resolve()),
        "scaled_geometry_transform_sha256": sha256_file(transform),
        "seeder_attempt": attempt,
        "seeder_calls": 1,
        "musubi_calls": 0,
        "first_failure": first_failure,
        "runtime_s": timed_runtime(started),
    }
    write_json(paths["qc"] / "unit_scaling_oracle.json", oracle)
    finalization_s = timed_runtime(started)
    runtime = {
        "successful_zero_run_s": zero_summary["runtime_s"],
        "scaled_geometry_prepare_s": json.loads(transform.read_text(encoding="utf-8"))[
            "runtime_s"
        ],
        "failed_scaled_seeder_call_s": attempt["wall_time_s"],
        "finalization_s": finalization_s,
    }
    runtime["accepted_total_s"] = sum(runtime.values())
    final = {
        "status": "CFD_FLOW_QVALUE_PATCH_TOPOLOGY_NOT_PROVEN",
        "actual_head_at_execution": _head(paths["root"]),
        "production_pipeline_modified": False,
        "seeder_calls": 1,
        "musubi_calls": 0,
        "harvester_calls": 0,
        "old_base_cells": 221309,
        "repaired_base_cells": int(components["fluid_cell_count"]),
        "scaled_unpatched_base_cells": None,
        "component_classification": components["component_classification"],
        "scale_invariance_oracle": "FAIL_NOT_PROVEN_NO_SCALED_MESH",
        "old_mesh_false_flooding": spatial["old_mesh_false_flooding"],
        "continuous_geometry_connected": bool(
            zero_summary["continuous_surface"]["continuous_geometry_connected"]
        ),
        "patch_affects_flooding_or_topology_source_proven": bool(
            callgraph["patch_affects_flooding_or_topology"]
        ),
        "break_location": zero_summary["nearest_centerline_segment"],
        "root_cause_final": (
            "ZERO_RUN_PROVES_A_ONE_CELL_D3Q19_DIAGONAL_ONLY_OUTLET_01 "
            "COMPONENT_AND_OLD_FALSE_FLOODING, BUT THE REQUIRED UNIT-SCALING "
            "ORACLE REMAINS UNPROVEN BECAUSE ITS ONLY ALLOWED SEEDER CALL DID "
            "NOT LOAD THE CONFIGURATION"
        ),
        "first_failure": first_failure,
        "evidence": {
            "component_forensics": "repaired_base_component_forensics.json",
            "spatial_difference": "old_vs_repaired_spatial_difference.json",
            "callgraph": "intersect_raytriangle_callgraph.json",
            "unit_scaling_oracle": "unit_scaling_oracle.json",
            "seeder_attempt": str(
                (paths["scaled"] / "seeder_attempt_summary.json").resolve()
            ),
        },
        "runtime": runtime,
        "next": (
            "IN A NEW AUTHORIZED ROUND, FIX ONLY THE WSL ARGUMENT TRANSPORT AND "
            "RUN ONE UNPATCHED MICROMETER-SCALED BASE SEEDER ORACLE"
        ),
    }
    write_json(paths["qc"] / "repaired_topology_root_cause.json", final)
    return final


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "phase",
        choices=("zero-run", "prepare-scaled", "finalize", "finalize-failed-launch"),
    )
    parser.add_argument("--project-root", type=Path, default=ROOT)
    args = parser.parse_args()
    if args.phase == "zero-run":
        result = run_zero(args.project_root)
    elif args.phase == "prepare-scaled":
        result = prepare_scaled(args.project_root)
    elif args.phase == "finalize":
        result = finalize(args.project_root)
    else:
        result = finalize_failed_launch(args.project_root)
    print(result["status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
