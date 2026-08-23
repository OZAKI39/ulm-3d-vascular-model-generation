"""Minimal controls that isolate tube scalar, port, and Y-junction behavior."""

from __future__ import annotations

from types import SimpleNamespace
from typing import Any

import numpy as np

from .config import CFDLumenConfig
from .geometry_diagnostics import explicit_union_for_diagnostics
from .lumen_builder import build_lumen_primitives, build_variable_radius_tube
from .mesh_defects import diagnose_mesh_defects
from .surface_qc import _section_polygon
from .types import BranchGeometry, PortGeometry


def _branch(
    branch_id: int,
    local_ids: tuple[int, ...],
    points: np.ndarray,
    radii: np.ndarray,
) -> BranchGeometry:
    step = np.linalg.norm(np.diff(points, axis=0), axis=1)
    return BranchGeometry(
        branch_id=branch_id,
        local_node_ids=local_ids,
        source_global_nodes=local_ids,
        source_global_edges=tuple(range(len(points) - 1)),
        raw_points_um=points.copy(),
        raw_radius_um=radii.copy(),
        points_um=points.copy(),
        radius_um=radii.copy(),
        arc_length_um=np.concatenate(([0.0], np.cumsum(step))),
    )


def _equivalent_radius(mesh: Any, center: np.ndarray, normal: np.ndarray) -> float | None:
    section = _section_polygon(mesh, center, normal)
    return float(np.sqrt(section[0] / np.pi)) if section else None


def _straight_port_control(
    config: CFDLumenConfig,
    *,
    tapered: bool,
) -> dict[str, Any]:
    x = np.linspace(0.0, 10.0, 21)
    points = np.column_stack((x, np.zeros_like(x), np.zeros_like(x)))
    radii = np.linspace(1.0, 2.0, len(x)) if tapered else np.full(len(x), 2.0)
    branch = _branch(0, (0, 1), points, radii)
    port = PortGeometry(
        port_id=0,
        roi_id="synthetic_taper" if tapered else "synthetic_straight",
        cut_port_id="cut_000",
        local_node_id=1,
        global_edge_id=0,
        original_position_um=np.asarray((10.0, 0.0, 0.0)),
        cap_center_um=np.asarray((20.0, 0.0, 0.0)),
        radius_um=2.0,
        outward_tangent=np.asarray((1.0, 0.0, 0.0)),
        extension_length_um=10.0,
        overlap_length_um=1.0,
        cylinder_start_um=np.asarray((9.0, 0.0, 0.0)),
        cylinder_end_um=np.asarray((20.0, 0.0, 0.0)),
        boundary_face="x_max",
    )
    roi = SimpleNamespace(
        local_node_ids=np.arange(2),
        local_edges=np.asarray(((0, 1),), dtype=np.int64),
        local_node_positions_um=np.asarray(((0.0, 0.0, 0.0), (10.0, 0.0, 0.0))),
        local_node_radius_um=np.asarray((float(radii[0]), float(radii[-1]))),
    )
    primitives = build_lumen_primitives([branch], roi, [port], config)
    mesh, runtime = explicit_union_for_diagnostics(primitives)
    epsilon = 0.02
    branch_radius = _equivalent_radius(
        primitives.branch_tubes[0], np.asarray((10.0 - epsilon, 0.0, 0.0)), np.asarray((1.0, 0.0, 0.0))
    )
    extension_radius = _equivalent_radius(
        primitives.port_extensions[0], np.asarray((10.0 + epsilon, 0.0, 0.0)), np.asarray((1.0, 0.0, 0.0))
    )
    step = (
        abs(branch_radius - extension_radius) / ((branch_radius + extension_radius) * 0.5)
        if branch_radius is not None and extension_radius is not None
        else None
    )
    return {
        "name": "tapered_tube_constant_extension" if tapered else "straight_tube_straight_extension",
        "status": "PASS" if step is not None and step < 0.01 and mesh.is_watertight else "FAIL",
        "target_endpoint_radius_um": 2.0,
        "branch_mesh_radius_um": branch_radius,
        "extension_mesh_radius_um": extension_radius,
        "step_radius_relative_error": step,
        "watertight": bool(mesh.is_watertight),
        "triangle_count": int(len(mesh.faces)),
        "runtime_s": runtime,
    }


def _tube_scalar_control(config: CFDLumenConfig) -> dict[str, Any]:
    points = np.asarray(((0.0, 0.0, 0.0), (5.0, 0.0, 0.0), (10.0, 0.0, 0.0)))
    branch = _branch(0, (0, 1), points, np.full(3, 2.0))
    tube = build_variable_radius_tube(branch, config.geometry.tube_sides)
    measured = _equivalent_radius(tube, np.asarray((5.0, 0.0, 0.0)), np.asarray((1.0, 0.0, 0.0)))
    relative_error = abs(measured - 2.0) / 2.0 if measured is not None else None
    return {
        "name": "vtk_absolute_radius_scalar_semantics",
        "status": "PASS" if relative_error is not None and relative_error < 0.02 else "FAIL",
        "vtk_tube_set_radius": 1.0,
        "vtk_vary_radius_mode": "VaryRadiusByAbsoluteScalar",
        "input_scalar_semantics": "radius_um",
        "target_radius_um": 2.0,
        "measured_mesh_radius_um": measured,
        "absolute_relative_error": relative_error,
    }


def _y_control(config: CFDLumenConfig) -> dict[str, Any]:
    positions = np.asarray(
        ((0.0, 0.0, 0.0), (-8.0, 0.0, 0.0), (7.0, 5.0, 0.0), (7.0, -5.0, 0.0))
    )
    radii = np.asarray((2.0, 2.0, 2.0, 2.0))
    branches = [
        _branch(index, (0, endpoint), positions[[0, endpoint]], radii[[0, endpoint]])
        for index, endpoint in enumerate((1, 2, 3))
    ]
    roi = SimpleNamespace(
        local_node_ids=np.arange(4),
        local_edges=np.asarray(((0, 1), (0, 2), (0, 3)), dtype=np.int64),
        local_node_positions_um=positions,
        local_node_radius_um=radii,
    )
    primitives = build_lumen_primitives(branches, roi, [], config)
    mesh, runtime = explicit_union_for_diagnostics(primitives)
    defects, _ = diagnose_mesh_defects(mesh, [(0, positions[0], radii[0])])
    passed = (
        defects["watertight"]
        and defects["boundary_edge_count"] == 0
        and defects["non_manifold_edge_count"] == 0
        and defects["self_intersection_count"] == 0
        and defects["surface_connected_component_count"] == 1
    )
    return {
        "name": "simple_symmetric_y_bifurcation",
        "status": "PASS" if passed else "FAIL",
        "watertight": defects["watertight"],
        "boundary_edge_count": defects["boundary_edge_count"],
        "non_manifold_edge_count": defects["non_manifold_edge_count"],
        "self_intersection_count": defects["self_intersection_count"],
        "suspected_internal_face_count": defects["suspected_internal_face_count"],
        "surface_connected_component_count": defects["surface_connected_component_count"],
        "triangle_count": int(len(mesh.faces)),
        "runtime_s": runtime,
    }


def run_synthetic_controls(config: CFDLumenConfig) -> dict[str, Any]:
    controls = [
        _tube_scalar_control(config),
        _straight_port_control(config, tapered=False),
        _straight_port_control(config, tapered=True),
        _y_control(config),
    ]
    return {
        "status": "PASS" if all(row["status"] == "PASS" for row in controls) else "FAIL",
        "controls": controls,
    }
