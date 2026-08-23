"""Per-junction manifold merge of local implicit patches and explicit tubes."""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Any

import manifold3d
import numpy as np
import trimesh

from utils.sampling.sampling_types import ROIRecord

from .config import CFDLumenConfig
from .local_implicit_junction import (
    build_local_junction_patch,
    define_junction_collars,
    trim_explicit_branches,
)
from .lumen_builder import (
    balanced_manifold_union,
    build_port_extension_mesh,
    build_variable_radius_tube,
)
from .mesh_defects import diagnose_mesh_defects
from .surface_transition import (
    build_open_explicit_tubes,
    build_transition_strip,
    clip_implicit_patch_to_transition_inner,
    combine_stitched_surfaces,
    extend_branches_to_cfd_ports,
    trim_branches_to_transition_outer,
)
from .types import (
    BranchGeometry,
    GeometryValidationError,
    HybridBuildDetails,
    PortGeometry,
)


@dataclass(slots=True)
class _ActiveComponent:
    mesh: trimesh.Trimesh
    branch_ids: set[int]
    junction_ids: set[int]


def _remove_degenerate_faces(
    mesh: trimesh.Trimesh,
    cleanup_tolerance_um: float,
) -> trimesh.Trimesh:
    if cleanup_tolerance_um > 0:
        # Boolean unions can retain microscopic coplanar cap slivers even though
        # manifold topology is valid.  Manifold's bounded simplifier removes
        # those slivers while guaranteeing a displacement below this explicit
        # physical-unit tolerance (1e-7 um by default).
        manifold = manifold3d.Manifold(
            mesh=manifold3d.Mesh(
                vert_properties=np.asarray(mesh.vertices, dtype=np.float32),
                tri_verts=np.asarray(mesh.faces, dtype=np.uint32),
            )
        ).simplify(cleanup_tolerance_um)
        output = manifold.to_mesh()
        mesh = trimesh.Trimesh(
            vertices=np.asarray(output.vert_properties[:, :3], dtype=float),
            faces=np.asarray(output.tri_verts, dtype=np.int64),
            process=False,
        )
    cleaned = mesh.copy()
    cleaned.update_faces(cleaned.nondegenerate_faces(height=1.0e-12))
    cleaned.update_faces(cleaned.unique_faces())
    cleaned.remove_unreferenced_vertices()
    trimesh.repair.fix_normals(cleaned, multibody=True)
    return cleaned


def _step_metrics(
    mesh: trimesh.Trimesh,
    roi: ROIRecord,
    junction_node_id: int,
) -> dict[str, Any]:
    center = np.asarray(roi.local_node_positions_um[junction_node_id], dtype=float)
    radius = float(roi.local_node_radius_um[junction_node_id])
    defects, _ = diagnose_mesh_defects(
        mesh,
        [(junction_node_id, center, radius)],
        ray_sample_limit=512,
    )
    return {
        "triangles": int(len(mesh.faces)),
        "volume_um3": float(abs(mesh.volume)) if mesh.is_watertight else None,
        "self_intersection_pairs": defects["self_intersection_count"],
        "internal_face_count": defects["suspected_internal_face_count"],
        "boundary_edge_count": defects["boundary_edge_count"],
        "nonmanifold_edge_count": defects["non_manifold_edge_count"],
    }


def _build_manifold_boolean_hybrid(
    branches: list[BranchGeometry],
    roi: ROIRecord,
    ports: list[PortGeometry],
    config: CFDLumenConfig,
    *,
    tube_sides: int | None = None,
    cells_across_min_diameter: int | None = None,
    record_step_qc: bool = False,
    junction_node_ids: set[int] | None = None,
    controlled_local_implicit: bool = False,
    continuous_port_extensions: bool = False,
) -> tuple[trimesh.Trimesh, HybridBuildDetails]:
    """Build the v3 default: explicit branches/ports plus local implicit junctions."""

    started_total = time.perf_counter()
    sides = int(tube_sides or config.geometry.tube_sides)
    stage = time.perf_counter()
    collars = define_junction_collars(
        roi, branches, config, junction_node_ids=junction_node_ids
    )
    construction_branches = branches
    port_extension_rows: list[dict[str, Any]] = []
    if continuous_port_extensions:
        construction_branches, port_extension_rows = extend_branches_to_cfd_ports(
            branches, ports, config
        )
    branch_by_id = {branch.branch_id: branch for branch in branches}
    patches = {
        node_id: build_local_junction_patch(
            roi,
            branch_by_id,
            node_id,
            node_collars,
            config,
            cells_across_min_diameter=cells_across_min_diameter,
            controlled=controlled_local_implicit,
        )
        for node_id, node_collars in collars.items()
    }
    local_implicit_runtime = time.perf_counter() - stage

    stage = time.perf_counter()
    trimmed_branches, cap_planes = trim_explicit_branches(construction_branches, collars)
    branch_tubes = {
        branch.branch_id: build_variable_radius_tube(branch, sides)
        for branch in trimmed_branches
    }
    port_meshes = (
        []
        if continuous_port_extensions
        else [build_port_extension_mesh(port, sides) for port in ports]
    )
    explicit_runtime = time.perf_counter() - stage

    active = [
        _ActiveComponent(mesh=mesh, branch_ids={branch_id}, junction_ids=set())
        for branch_id, mesh in branch_tubes.items()
    ]
    merge_steps: list[dict[str, Any]] = []
    merged_junction_meshes: dict[int, trimesh.Trimesh] = {}
    merge_runtime = 0.0
    step_qc_runtime = 0.0
    for node_id, patch in patches.items():
        incident_ids = {collar.branch_id for collar in patch.collars}
        selected = [component for component in active if component.branch_ids & incident_ids]
        selected_identity = {id(component) for component in selected}
        active = [component for component in active if id(component) not in selected_identity]
        before_meshes = [patch.clean_mesh, *[component.mesh for component in selected]]
        stage = time.perf_counter()
        merged = _remove_degenerate_faces(
            balanced_manifold_union(before_meshes),
            config.hybrid_merge.cleanup_tolerance_um,
        )
        merge_runtime += time.perf_counter() - stage
        row: dict[str, Any] = {
            "junction_node_id": node_id,
            "merge_backend": "manifold",
            "input_solid_count": len(before_meshes),
            "incident_branch_ids": sorted(incident_ids),
            "before_triangles": int(sum(len(mesh.faces) for mesh in before_meshes)),
            "before_volume_sum_um3": float(sum(abs(mesh.volume) for mesh in before_meshes)),
            "after_triangles": int(len(merged.faces)),
            "after_volume_um3": float(abs(merged.volume)),
        }
        if record_step_qc:
            qc_started = time.perf_counter()
            row.update({f"after_{key}": value for key, value in _step_metrics(merged, roi, node_id).items()})
            step_qc_runtime += time.perf_counter() - qc_started
        merge_steps.append(row)
        merged_junction_meshes[node_id] = merged
        active.append(
            _ActiveComponent(
                mesh=merged,
                branch_ids=set().union(*(component.branch_ids for component in selected)),
                junction_ids={node_id}.union(*(component.junction_ids for component in selected)),
            )
        )

    stage = time.perf_counter()
    final_mesh = _remove_degenerate_faces(
        balanced_manifold_union([*[component.mesh for component in active], *port_meshes]),
        config.hybrid_merge.cleanup_tolerance_um,
    )
    final_merge_runtime = time.perf_counter() - stage
    details = HybridBuildDetails(
        patches=patches,
        merged_junction_meshes=merged_junction_meshes,
        trimmed_branches=trimmed_branches,
        merge_steps=merge_steps,
        explicit_cap_planes=cap_planes,
        runtime_s={
            "local_implicit": local_implicit_runtime,
            "explicit_primitives": explicit_runtime,
            "per_junction_merge": merge_runtime,
            "final_simple_merge": final_merge_runtime,
            "step_qc_excluded_from_reconstruction": step_qc_runtime,
            "reconstruction_total": (
                local_implicit_runtime + explicit_runtime + merge_runtime + final_merge_runtime
            ),
            "wall_total_with_step_qc": time.perf_counter() - started_total,
        },
        transition_backend="manifold_boolean",
        constructed_branches=construction_branches,
        port_extension_rows=port_extension_rows,
    )
    return final_mesh, details


def _build_loop_stitch_hybrid(
    branches: list[BranchGeometry],
    roi: ROIRecord,
    ports: list[PortGeometry],
    config: CFDLumenConfig,
    *,
    tube_sides: int,
    cells_across_min_diameter: int | None,
    controlled_local_implicit: bool,
) -> tuple[trimesh.Trimesh, HybridBuildDetails]:
    """Build v5 open surfaces joined by explicit transition collar strips."""

    started_total = time.perf_counter()
    stage = time.perf_counter()
    collars = define_junction_collars(roi, branches, config)
    branch_by_id = {branch.branch_id: branch for branch in branches}
    patches = {
        node_id: build_local_junction_patch(
            roi,
            branch_by_id,
            node_id,
            node_collars,
            config,
            cells_across_min_diameter=cells_across_min_diameter,
            controlled=controlled_local_implicit,
        )
        for node_id, node_collars in collars.items()
    }
    local_implicit_runtime = time.perf_counter() - stage

    stage = time.perf_counter()
    construction_branches, port_extension_rows = extend_branches_to_cfd_ports(
        branches, ports, config
    )
    trimmed_branches, collars_by_branch = trim_branches_to_transition_outer(
        construction_branches, collars
    )
    explicit_meshes, explicit_loops = build_open_explicit_tubes(
        trimmed_branches, collars_by_branch, tube_sides
    )
    explicit_runtime = time.perf_counter() - stage

    stage = time.perf_counter()
    clipped_patches: dict[int, trimesh.Trimesh] = {}
    implicit_loops: dict[
        tuple[int, int], tuple[np.ndarray, np.ndarray, np.ndarray]
    ] = {}
    for node_id, patch in patches.items():
        tolerance = max(float(patch.metadata["grid_spacing_um"]) * 1.0e-3, 1.0e-6)
        clipped, node_loops = clip_implicit_patch_to_transition_inner(
            patch.clean_mesh,
            patch.collars,
            branch_by_id,
            tolerance,
            tube_sides,
        )
        clipped_patches[node_id] = clipped
        patch.clean_mesh = clipped
        patch.metadata["v5_transition_inner_clip"] = True
        for branch_id, payload in node_loops.items():
            implicit_loops[(node_id, branch_id)] = payload

    transition_meshes: list[trimesh.Trimesh] = []
    transition_rows: list[dict[str, Any]] = []
    transition_by_junction: dict[int, list[trimesh.Trimesh]] = {}
    for node_id, node_collars in collars.items():
        for collar in node_collars:
            key = (node_id, collar.branch_id)
            if key not in implicit_loops or key not in explicit_loops:
                raise GeometryValidationError(
                    f"COLLAR_LOOP_EXTRACTION_FAILED: missing paired loops for J{node_id}/B{collar.branch_id}"
                )
            implicit_loop, implicit_center, _ = implicit_loops[key]
            explicit_loop, explicit_center, _ = explicit_loops[key]
            strip, row = build_transition_strip(
                clipped_patches[node_id],
                implicit_loop,
                implicit_center,
                explicit_meshes[collar.branch_id],
                explicit_loop,
                explicit_center,
                branch_by_id[collar.branch_id],
                collar,
                tube_sides,
                config,
            )
            row.update(
                {
                    "junction_node_id": node_id,
                    "branch_id": collar.branch_id,
                    "inner_distance_from_junction_um": collar.explicit_cap_distance_um,
                    "outer_distance_from_junction_um": collar.implicit_extent_um,
                    "source_collar_radius_um": collar.collar_radius_um,
                }
            )
            transition_meshes.append(strip)
            transition_rows.append(row)
            transition_by_junction.setdefault(node_id, []).append(strip)

    final_mesh, face_region = combine_stitched_surfaces(
        list(clipped_patches.values()),
        transition_meshes,
        list(explicit_meshes.values()),
    )
    stitch_runtime = time.perf_counter() - stage
    if not final_mesh.is_watertight:
        boundary_count = int(
            np.count_nonzero(
                np.unique(np.sort(final_mesh.edges, axis=1), axis=0, return_counts=True)[1]
                == 1
            )
        )
        raise GeometryValidationError(
            f"LOOP_STITCH_TOPOLOGY_FAILED: final surface has {boundary_count} boundary edges"
        )
    merged_junction_meshes = {
        node_id: trimesh.util.concatenate(
            [clipped_patches[node_id], *transition_by_junction.get(node_id, [])]
        )
        for node_id in clipped_patches
    }
    details = HybridBuildDetails(
        patches=patches,
        merged_junction_meshes=merged_junction_meshes,
        trimmed_branches=trimmed_branches,
        merge_steps=[
            {
                "junction_node_id": node_id,
                "merge_backend": "boundary_loop_stitch",
                "incident_branch_ids": [
                    collar.branch_id for collar in collars[node_id]
                ],
                "transition_triangle_count": sum(
                    int(row["transition_triangle_count"])
                    for row in transition_rows
                    if row["junction_node_id"] == node_id
                ),
            }
            for node_id in sorted(collars)
        ],
        explicit_cap_planes=[],
        runtime_s={
            "local_implicit": local_implicit_runtime,
            "explicit_primitives": explicit_runtime,
            "transition_stitch": stitch_runtime,
            "per_junction_merge": 0.0,
            "final_simple_merge": 0.0,
            "step_qc_excluded_from_reconstruction": 0.0,
            "reconstruction_total": time.perf_counter() - started_total,
            "wall_total_with_step_qc": time.perf_counter() - started_total,
        },
        transition_backend="loop_stitch",
        face_region=face_region,
        transition_rows=transition_rows,
        constructed_branches=construction_branches,
        port_extension_rows=port_extension_rows,
    )
    return final_mesh, details


def build_hybrid_lumen(
    branches: list[BranchGeometry],
    roi: ROIRecord,
    ports: list[PortGeometry],
    config: CFDLumenConfig,
    *,
    tube_sides: int | None = None,
    cells_across_min_diameter: int | None = None,
    record_step_qc: bool = False,
    junction_node_ids: set[int] | None = None,
    controlled_local_implicit: bool = False,
    transition_backend: str | None = None,
    continuous_port_extensions: bool | None = None,
) -> tuple[trimesh.Trimesh, HybridBuildDetails]:
    """Build v5 loop-stitch geometry, retaining v4 Boolean as explicit fallback."""

    sides = int(tube_sides or config.geometry.tube_sides)
    backend = transition_backend or config.hybrid_transition.backend
    continuous_ports = (
        config.port_transition.backend == "continuous_centerline"
        if continuous_port_extensions is None
        else continuous_port_extensions
    )
    if backend == "manifold_boolean":
        return _build_manifold_boolean_hybrid(
            branches,
            roi,
            ports,
            config,
            tube_sides=sides,
            cells_across_min_diameter=cells_across_min_diameter,
            record_step_qc=record_step_qc,
            junction_node_ids=junction_node_ids,
            controlled_local_implicit=controlled_local_implicit,
            continuous_port_extensions=continuous_ports,
        )
    try:
        return _build_loop_stitch_hybrid(
            branches,
            roi,
            ports,
            config,
            tube_sides=sides,
            cells_across_min_diameter=cells_across_min_diameter,
            controlled_local_implicit=controlled_local_implicit,
        )
    except Exception as exc:
        if config.hybrid_transition.fallback_backend != "manifold_boolean":
            raise
        mesh, details = _build_manifold_boolean_hybrid(
            branches,
            roi,
            ports,
            config,
            tube_sides=sides,
            cells_across_min_diameter=cells_across_min_diameter,
            record_step_qc=record_step_qc,
            junction_node_ids=junction_node_ids,
            controlled_local_implicit=controlled_local_implicit,
            continuous_port_extensions=continuous_ports,
        )
        details.transition_fallback_reason = f"{type(exc).__name__}: {exc}"
        return mesh, details


def build_hybrid_junction_component(
    branches: list[BranchGeometry],
    roi: ROIRecord,
    config: CFDLumenConfig,
    junction_node_id: int,
    *,
    tube_sides: int | None = None,
    cells_across_min_diameter: int | None = None,
    controlled_local_implicit: bool = False,
) -> tuple[trimesh.Trimesh, HybridBuildDetails]:
    """Build one real junction and its complete incident explicit branches for smoke tests.

    This deliberately does not alter or truncate the ROI graph.  It is a validation
    view of one junction component, used for the v3 junction-49 grid study when a
    different junction in the same ROI has an independently reported collar/port
    conflict.
    """

    started = time.perf_counter()
    sides = int(tube_sides or config.geometry.tube_sides)
    collars = define_junction_collars(
        roi,
        branches,
        config,
        junction_node_ids={int(junction_node_id)},
    )
    if junction_node_id not in collars:
        raise ValueError(f"Node {junction_node_id} is not a bifurcation in this ROI")
    branch_by_id = {branch.branch_id: branch for branch in branches}
    patch = build_local_junction_patch(
        roi,
        branch_by_id,
        junction_node_id,
        collars[junction_node_id],
        config,
        cells_across_min_diameter=cells_across_min_diameter,
        controlled=controlled_local_implicit,
    )
    implicit_finished = time.perf_counter()
    trimmed, cap_planes = trim_explicit_branches(branches, collars)
    incident_ids = {collar.branch_id for collar in patch.collars}
    incident_tubes = [
        build_variable_radius_tube(branch, sides)
        for branch in trimmed
        if branch.branch_id in incident_ids
    ]
    explicit_finished = time.perf_counter()
    merged = _remove_degenerate_faces(
        balanced_manifold_union([patch.clean_mesh, *incident_tubes]),
        config.hybrid_merge.cleanup_tolerance_um,
    )
    merge_finished = time.perf_counter()
    details = HybridBuildDetails(
        patches={junction_node_id: patch},
        merged_junction_meshes={junction_node_id: merged},
        trimmed_branches=[branch for branch in trimmed if branch.branch_id in incident_ids],
        merge_steps=[
            {
                "junction_node_id": junction_node_id,
                "merge_backend": "manifold",
                "input_solid_count": 1 + len(incident_tubes),
                "incident_branch_ids": sorted(incident_ids),
                "before_triangles": int(
                    len(patch.clean_mesh.faces)
                    + sum(len(tube.faces) for tube in incident_tubes)
                ),
                "after_triangles": int(len(merged.faces)),
                "after_volume_um3": float(abs(merged.volume)),
                "validation_scope": "single_junction_with_complete_incident_branches",
            }
        ],
        explicit_cap_planes=[
            row for row in cap_planes if int(row["branch_id"]) in incident_ids
        ],
        runtime_s={
            "local_implicit": implicit_finished - started,
            "explicit_primitives": explicit_finished - implicit_finished,
            "per_junction_merge": merge_finished - explicit_finished,
            "final_simple_merge": 0.0,
            "step_qc_excluded_from_reconstruction": 0.0,
            "reconstruction_total": merge_finished - started,
            "wall_total_with_step_qc": merge_finished - started,
        },
    )
    return merged, details
