"""v3 synthetic controls and three-method real-ROI validation helpers."""

from __future__ import annotations

from types import SimpleNamespace
from typing import Any

import numpy as np

from .config import CFDLumenConfig
from .hybrid_merge import build_hybrid_lumen
from .hybrid_qc import collar_radius_rows, evaluate_hybrid_surface_qc
from .lumen_builder import build_variable_radius_tube
from .surface_continuity_qc import port_continuity_report
from .surface_transition import CFD_DERIVED_EXTENSION, extend_branches_to_cfd_ports
from .types import BranchGeometry, PortGeometry


def _synthetic_branch(
    branch_id: int,
    endpoint_id: int,
    endpoint: np.ndarray,
    junction_radius: float,
    endpoint_radius: float,
) -> BranchGeometry:
    fractions = np.linspace(0.0, 1.0, 61)
    points = fractions[:, None] * endpoint[None, :]
    radii = (1.0 - fractions) * junction_radius + fractions * endpoint_radius
    arc = np.concatenate(([0.0], np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))))
    return BranchGeometry(
        branch_id=branch_id,
        local_node_ids=(0, endpoint_id),
        source_global_nodes=(0, endpoint_id),
        source_global_edges=(branch_id,),
        raw_points_um=points.copy(),
        raw_radius_um=radii.copy(),
        points_um=points,
        radius_um=radii,
        arc_length_um=arc,
    )


def _synthetic_case(
    name: str,
    endpoints: list[tuple[tuple[float, float, float], float]],
    junction_radius: float,
    config: CFDLumenConfig,
) -> dict[str, Any]:
    branches = [
        _synthetic_branch(
            branch_id,
            branch_id + 1,
            np.asarray(endpoint, dtype=float),
            junction_radius,
            endpoint_radius,
        )
        for branch_id, (endpoint, endpoint_radius) in enumerate(endpoints)
    ]
    positions = np.vstack((np.zeros(3), *[np.asarray(item[0], dtype=float) for item in endpoints]))
    radii = np.asarray((junction_radius, *[item[1] for item in endpoints]), dtype=float)
    edges = np.asarray([(0, index + 1) for index in range(len(endpoints))], dtype=np.int64)
    roi = SimpleNamespace(
        roi_id=name,
        node_count=len(positions),
        local_node_ids=np.arange(len(positions)),
        local_edges=edges,
        local_node_positions_um=positions,
        local_node_radius_um=radii,
        cut_ports=(),
    )
    mesh, details = build_hybrid_lumen(branches, roi, [], config)
    qc, _ = evaluate_hybrid_surface_qc(mesh, details, roi, config)
    collar_rows = collar_radius_rows(mesh, details, branches)
    max_collar_error = max(
        (
            row["absolute_radius_relative_error"]
            for row in collar_rows
            if row["absolute_radius_relative_error"] is not None
        ),
        default=None,
    )
    return {
        "name": name,
        "status": qc["status"],
        "watertight": qc["checks"]["watertight"],
        "connected": qc["checks"]["single_component"],
        "boundary_edges": qc["boundary_edge_count"],
        "nonmanifold_edges": qc["nonmanifold_edge_count"],
        "self_intersections": qc["self_intersection_pairs"],
        "internal_faces": qc["internal_face_count"],
        "internal_caps": qc["internal_cap_face_count"],
        "degenerate_triangles": qc["degenerate_triangle_count"],
        "triangle_count": qc["triangle_count"],
        "max_collar_radius_error": max_collar_error,
        "runtime_s": details.runtime_s["reconstruction_total"],
    }


def run_hybrid_synthetic_controls(config: CFDLumenConfig) -> dict[str, Any]:
    cases = [
        ("symmetric_y", [((-24.0, 0.0, 0.0), 2.0), ((22.0, 14.0, 0.0), 2.0), ((22.0, -14.0, 0.0), 2.0)], 2.0),
        ("asymmetric_y", [((-24.0, 0.0, 0.0), 2.6), ((22.0, 14.0, 0.0), 1.6), ((22.0, -14.0, 0.0), 1.1)], 2.0),
        ("acute_angle_bifurcation", [((-28.0, 0.0, 0.0), 2.0), ((28.0, 5.0, 0.0), 1.8), ((28.0, -5.0, 0.0), 1.8)], 2.0),
        ("large_radius_ratio", [((-28.0, 0.0, 0.0), 3.5), ((25.0, 15.0, 0.0), 1.0), ((25.0, -15.0, 0.0), 0.8)], 2.2),
        ("three_way_junction", [((-24.0, 0.0, 0.0), 2.0), ((22.0, 14.0, 0.0), 1.8), ((22.0, -14.0, 0.0), 1.8), ((0.0, 0.0, 24.0), 1.4)], 2.0),
        ("hybrid_collar_transition", [((-30.0, 0.0, 0.0), 2.4), ((28.0, 16.0, 0.0), 1.7), ((28.0, -16.0, 0.0), 1.3)], 2.0),
    ]
    rows: list[dict[str, Any]] = []
    for name, endpoints, radius in cases:
        try:
            rows.append(_synthetic_case(name, endpoints, radius, config))
        except Exception as exc:
            rows.append(
                {
                    "name": name,
                    "status": "FAIL",
                    "failure_reason": f"{type(exc).__name__}: {exc}",
                }
            )
    return {
        "status": "PASS" if all(row["status"] == "PASS" for row in rows) else "FAIL",
        "cases": rows,
    }


def _continuous_curved_port_case(config: CFDLumenConfig) -> dict[str, Any]:
    x = np.linspace(0.0, 18.0, 73)
    points = np.column_stack(
        (x, 1.8 * np.sin(np.pi * x / 36.0), 0.6 * np.sin(np.pi * x / 18.0))
    )
    radii = np.linspace(1.8, 1.4, len(points))
    branch = BranchGeometry(
        branch_id=0,
        local_node_ids=(0, 1),
        source_global_nodes=(0, 1),
        source_global_edges=(0,),
        raw_points_um=points.copy(),
        raw_radius_um=radii.copy(),
        points_um=points,
        radius_um=radii,
        arc_length_um=np.concatenate(
            ([0.0], np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1)))
        ),
    )
    tangent = points[-1] - points[-2]
    tangent /= np.linalg.norm(tangent)
    extension_length = 8.0
    port = PortGeometry(
        port_id=0,
        roi_id="curved_tube_straight_extension",
        cut_port_id="cut_curved",
        local_node_id=1,
        global_edge_id=0,
        original_position_um=points[-1].copy(),
        cap_center_um=points[-1] + extension_length * tangent,
        radius_um=float(radii[-1]),
        outward_tangent=tangent,
        extension_length_um=extension_length,
        overlap_length_um=0.0,
        cylinder_start_um=points[-1].copy(),
        cylinder_end_um=points[-1] + extension_length * tangent,
        boundary_face="x_max",
        boundary_role="CUT_PORT_OUTLET",
        source_core_cut_port_id="cut_curved",
        original_core_cut_position_um=tuple(map(float, points[-1])),
    )
    constructed, provenance = extend_branches_to_cfd_ports(
        [branch], [port], config
    )
    mesh = build_variable_radius_tube(constructed[0], config.geometry.tube_sides)
    rows, continuity = port_continuity_report(
        mesh, [branch], [port], config, continuous_centerline=True
    )
    derived_count = int(
        np.count_nonzero(
            constructed[0].construction_point_type == CFD_DERIVED_EXTENSION
        )
    )
    passed = (
        mesh.is_watertight
        and continuity["status"] == "PASS"
        and continuity["separate_cylinder_primitive_count"] == 0
        and derived_count > 0
        and len(provenance) > 0
    )
    return {
        "name": "curved_tube_to_straight_continuous_extension",
        "status": "PASS" if passed else "FAIL",
        "watertight": bool(mesh.is_watertight),
        "separate_cylinder_primitive_count": continuity[
            "separate_cylinder_primitive_count"
        ],
        "derived_extension_point_count": derived_count,
        "profile_sample_count": len(rows),
        "maximum_cut_area_absolute_relative_error": continuity[
            "maximum_cut_area_absolute_relative_error"
        ],
        "port_normal_jump_p95_deg": continuity["port_normal_jump"][
            "normal_jump_p95_deg"
        ],
        "port_normal_jump_p99_deg": continuity["port_normal_jump"][
            "normal_jump_p99_deg"
        ],
    }


def _curved_daughter_case(config: CFDLumenConfig) -> dict[str, Any]:
    fractions = np.linspace(0.0, 1.0, 81)
    curves = (
        np.column_stack((-26.0 * fractions, np.zeros_like(fractions), np.zeros_like(fractions))),
        np.column_stack(
            (
                24.0 * fractions,
                13.0 * fractions + 1.5 * np.sin(np.pi * fractions),
                2.0 * np.sin(np.pi * fractions),
            )
        ),
        np.column_stack(
            (
                23.0 * fractions,
                -12.0 * fractions - 1.2 * np.sin(np.pi * fractions),
                -1.5 * np.sin(np.pi * fractions),
            )
        ),
    )
    endpoint_radii = (2.4, 1.5, 1.0)
    junction_radius = 2.0
    branches: list[BranchGeometry] = []
    for branch_id, (points, endpoint_radius) in enumerate(
        zip(curves, endpoint_radii)
    ):
        radii = np.linspace(junction_radius, endpoint_radius, len(points))
        branches.append(
            BranchGeometry(
                branch_id=branch_id,
                local_node_ids=(0, branch_id + 1),
                source_global_nodes=(0, branch_id + 1),
                source_global_edges=(branch_id,),
                raw_points_um=points.copy(),
                raw_radius_um=radii.copy(),
                points_um=points,
                radius_um=radii,
                arc_length_um=np.concatenate(
                    ([0.0], np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1)))
                ),
            )
        )
    positions = np.vstack((np.zeros(3), *[curve[-1] for curve in curves]))
    roi = SimpleNamespace(
        roi_id="curved_daughter_transition",
        node_count=4,
        local_node_ids=np.arange(4),
        local_edges=np.asarray(((0, 1), (0, 2), (0, 3)), dtype=np.int64),
        local_node_positions_um=positions,
        local_node_radius_um=np.asarray(
            (junction_radius, *endpoint_radii), dtype=float
        ),
        cut_ports=(),
    )
    mesh, details = build_hybrid_lumen(branches, roi, [], config)
    qc, _ = evaluate_hybrid_surface_qc(mesh, details, roi, config)
    passed = (
        qc["status"] == "PASS"
        and details.transition_backend == "loop_stitch"
        and details.transition_fallback_reason is None
    )
    return {
        "name": "curved_daughter_loop_stitch_transition",
        "status": "PASS" if passed else "FAIL",
        "transition_backend": details.transition_backend,
        "fallback_reason": details.transition_fallback_reason,
        "self_intersections": qc["self_intersection_pairs"],
        "internal_faces": qc["internal_face_count"],
        "internal_caps": qc["internal_cap_face_count"],
        "boundary_edges": qc["boundary_edge_count"],
        "nonmanifold_edges": qc["nonmanifold_edge_count"],
        "components": qc["surface_component_count"],
        "transition_triangle_count": sum(
            int(row["transition_triangle_count"])
            for row in details.transition_rows
        ),
    }


def run_v5_synthetic_controls(config: CFDLumenConfig) -> dict[str, Any]:
    """Formal v5 controls covering continuous ports and collar stitching."""

    hybrid = run_hybrid_synthetic_controls(config)
    extra: list[dict[str, Any]] = []
    for function in (_continuous_curved_port_case, _curved_daughter_case):
        try:
            extra.append(function(config))
        except Exception as exc:
            extra.append(
                {
                    "name": function.__name__,
                    "status": "FAIL",
                    "failure_reason": f"{type(exc).__name__}: {exc}",
                }
            )
    all_cases = [*hybrid["cases"], *extra]
    return {
        "protocol": "v5 Surface Continuity Refinement synthetic controls",
        "status": "PASS" if all(row["status"] == "PASS" for row in all_cases) else "FAIL",
        "case_count": len(all_cases),
        "cases": all_cases,
    }
