"""Synthetic acceptance controls required by v6."""

from __future__ import annotations

from types import SimpleNamespace
from typing import Any

import numpy as np

from .config import CFDLumenConfig
from .continuous_field_transition import (
    build_continuous_field_hybrid,
    conditioned_port_profile,
)
from .hybrid_qc import evaluate_hybrid_surface_qc
from .types import BranchGeometry, PortGeometry


def _branch(
    branch_id: int,
    endpoint_id: int,
    points: np.ndarray,
    radii: np.ndarray,
) -> BranchGeometry:
    arc = np.concatenate(
        ([0.0], np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1)))
    )
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


def _port_case(name: str, curved: bool, config: CFDLumenConfig) -> dict[str, Any]:
    parameter = np.linspace(0.0, 1.0, 61)
    if curved:
        points = np.column_stack(
            (
                20.0 * parameter,
                3.0 * parameter**2,
                1.2 * np.sin(0.5 * np.pi * parameter),
            )
        )
    else:
        points = np.column_stack(
            (20.0 * parameter, np.zeros_like(parameter), np.zeros_like(parameter))
        )
    radii = 2.0 - 0.5 * parameter
    branch = _branch(0, 1, points, radii)
    endpoint_tangent = points[-1] - points[-2]
    endpoint_tangent /= np.linalg.norm(endpoint_tangent)
    if curved:
        endpoint_tangent = np.asarray((1.0, 0.0, 0.0))
    extension_length = 10.0
    port = PortGeometry(
        port_id=0,
        roi_id=name,
        cut_port_id=f"{name}_cut",
        local_node_id=1,
        global_edge_id=0,
        original_position_um=points[-1].copy(),
        cap_center_um=points[-1] + extension_length * endpoint_tangent,
        radius_um=float(radii[-1]),
        outward_tangent=endpoint_tangent,
        extension_length_um=extension_length,
        overlap_length_um=0.0,
        cylinder_start_um=points[-1].copy(),
        cylinder_end_um=points[-1] + extension_length * endpoint_tangent,
        boundary_face="x_max",
    )
    source_points = branch.points_um.copy()
    source_radius = branch.radius_um.copy()
    _, conditioned_radius, _, report = conditioned_port_profile(
        branch, -1, port, config
    )
    passed = (
        report["v6_tangent_jump_deg"] < report["v5_tangent_jump_deg"] + 1.0e-9
        and report["v6_radius_slope_jump_um_per_um"]
        < report["v5_radius_slope_jump_um_per_um"]
        and np.array_equal(branch.points_um, source_points)
        and np.array_equal(branch.radius_um, source_radius)
        and np.isclose(conditioned_radius[-1], port.radius_um)
    )
    return {
        "name": name,
        "status": "PASS" if passed else "FAIL",
        "v5_tangent_jump_deg": report["v5_tangent_jump_deg"],
        "v6_tangent_jump_deg": report["v6_tangent_jump_deg"],
        "v5_radius_slope_jump": report["v5_radius_slope_jump_um_per_um"],
        "v6_radius_slope_jump": report["v6_radius_slope_jump_um_per_um"],
        "source_points_unchanged": bool(np.array_equal(branch.points_um, source_points)),
        "source_radius_unchanged": bool(np.array_equal(branch.radius_um, source_radius)),
        "final_port_radius_preserved": bool(
            np.isclose(conditioned_radius[-1], port.radius_um)
        ),
    }


def _junction_case(
    name: str,
    curves: list[np.ndarray],
    endpoint_radii: tuple[float, float, float],
    config: CFDLumenConfig,
) -> dict[str, Any]:
    junction_radius = 2.0
    branches = [
        _branch(
            index,
            index + 1,
            points,
            np.linspace(junction_radius, endpoint_radii[index], len(points)),
        )
        for index, points in enumerate(curves)
    ]
    positions = np.vstack((np.zeros(3), *[points[-1] for points in curves]))
    roi = SimpleNamespace(
        roi_id=name,
        node_count=4,
        local_node_ids=np.arange(4),
        local_edges=np.asarray(((0, 1), (0, 2), (0, 3)), dtype=np.int64),
        local_node_positions_um=positions,
        local_node_radius_um=np.asarray((junction_radius, *endpoint_radii)),
        cut_ports=(),
    )
    mesh, details, _ = build_continuous_field_hybrid(
        branches, roi, [], config
    )
    qc, _ = evaluate_hybrid_surface_qc(mesh, details, roi, config)
    passed = (
        qc["status"] == "PASS"
        and details.transition_backend == "continuous_implicit_field"
        and all(
            patch.metadata["one_marching_cubes_extraction"]
            and not patch.metadata["surface_loop_stitching"]
            for patch in details.patches.values()
        )
    )
    return {
        "name": name,
        "status": "PASS" if passed else "FAIL",
        "transition_backend": details.transition_backend,
        "merge_backend": details.merge_steps[0]["merge_backend"],
        "self_intersections": qc["self_intersection_pairs"],
        "internal_faces": qc["internal_face_count"],
        "internal_caps": qc["internal_cap_face_count"],
        "boundary_edges": qc["boundary_edge_count"],
        "nonmanifold_edges": qc["nonmanifold_edge_count"],
        "components": qc["surface_component_count"],
        "triangle_count": qc["triangle_count"],
        "runtime_s": details.runtime_s["reconstruction_total"],
    }


def run_v6_synthetic_controls(config: CFDLumenConfig) -> dict[str, Any]:
    fraction = np.linspace(0.0, 1.0, 61)

    def line(endpoint: tuple[float, float, float]) -> np.ndarray:
        return fraction[:, None] * np.asarray(endpoint, dtype=float)[None, :]

    straight_y = [line((-26.0, 0.0, 0.0)), line((24.0, 14.0, 0.0)), line((24.0, -14.0, 0.0))]
    acute_y = [line((-28.0, 0.0, 0.0)), line((28.0, 5.5, 0.0)), line((28.0, -5.5, 0.0))]
    unequal_y = [line((-28.0, 0.0, 0.0)), line((25.0, 15.0, 0.0)), line((25.0, -15.0, 0.0))]
    curved = [
        line((-26.0, 0.0, 0.0)),
        np.column_stack(
            (
                24.0 * fraction,
                13.0 * fraction + 1.5 * np.sin(np.pi * fraction),
                2.0 * np.sin(np.pi * fraction),
            )
        ),
        np.column_stack(
            (
                23.0 * fraction,
                -12.0 * fraction - 1.2 * np.sin(np.pi * fraction),
                -1.5 * np.sin(np.pi * fraction),
            )
        ),
    ]
    cases: list[dict[str, Any]] = []
    specifications = (
        ("straight_branch_with_source_taper", lambda: _port_case("straight_branch_with_source_taper", False, config)),
        ("curved_branch_to_cfd_extension", lambda: _port_case("curved_branch_to_cfd_extension", True, config)),
        ("y_junction_continuous_implicit_transition", lambda: _junction_case("y_junction_continuous_implicit_transition", straight_y, (2.0, 1.8, 1.6), config)),
        ("acute_y_junction", lambda: _junction_case("acute_y_junction", acute_y, (2.0, 1.7, 1.7), config)),
        ("unequal_radius_y_junction", lambda: _junction_case("unequal_radius_y_junction", unequal_y, (3.2, 1.2, 0.9), config)),
        ("curved_daughter_branches", lambda: _junction_case("curved_daughter_branches", curved, (2.4, 1.5, 1.0), config)),
    )
    for name, function in specifications:
        try:
            cases.append(function())
        except Exception as exc:
            cases.append(
                {
                    "name": name,
                    "status": "FAIL",
                    "failure_reason": f"{type(exc).__name__}: {exc}",
                }
            )
    return {
        "protocol": "v6 Continuous-Field Transition Refinement synthetic controls",
        "status": "PASS" if all(row["status"] == "PASS" for row in cases) else "FAIL",
        "case_count": len(cases),
        "cases": cases,
    }
