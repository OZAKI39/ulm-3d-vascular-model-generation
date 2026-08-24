"""Optional constrained, topology-preserving local fairing for v9."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import trimesh

from .config import CFDLumenConfig
from .types import PatchResult, PortGeometry
from .unified_polyball import JunctionBlendSpec


@dataclass(slots=True)
class JunctionFairingResult:
    mesh: trimesh.Trimesh
    moved_vertex_ids: np.ndarray
    fixed_vertex_ids: np.ndarray
    report: dict[str, Any]


def _vertex_neighbors(mesh: trimesh.Trimesh) -> list[np.ndarray]:
    output: list[list[int]] = [[] for _ in range(len(mesh.vertices))]
    for first, second in np.asarray(mesh.edges_unique, dtype=np.int64):
        output[int(first)].append(int(second))
        output[int(second)].append(int(first))
    return [np.asarray(sorted(set(row)), dtype=np.int64) for row in output]


def constrained_junction_fairing(
    mesh: trimesh.Trimesh,
    defect_edges: np.ndarray,
    specs: tuple[JunctionBlendSpec, ...],
    ports: list[PortGeometry],
    patch: PatchResult,
    config: CFDLumenConfig,
) -> JunctionFairingResult:
    """Fair only the detected defect band; ports and the outer support ring are fixed."""

    vertices = np.asarray(mesh.vertices, dtype=float).copy()
    original = vertices.copy()
    edges = np.asarray(defect_edges, dtype=np.int64).reshape((-1, 2))
    if not len(edges):
        return JunctionFairingResult(
            mesh=mesh.copy(),
            moved_vertex_ids=np.empty(0, dtype=np.int64),
            fixed_vertex_ids=np.arange(len(vertices), dtype=np.int64),
            report={
                "status": "NOT_NEEDED",
                "reason": "no residual segment-switch defect band",
                "newton_projection_after_fairing": False,
            },
        )
    neighbors = _vertex_neighbors(mesh)
    core = set(map(int, np.unique(edges)))
    support = set(core)
    outer_ring: set[int] = set()
    frontier = set(core)
    for _ in range(config.v9.fairing_neighborhood_rings):
        following = {
            neighbor
            for vertex in frontier
            for neighbor in neighbors[vertex].tolist()
            if neighbor not in support
        }
        support.update(following)
        outer_ring = following
        frontier = following

    inside_junction = np.zeros(len(vertices), dtype=bool)
    local_radius = np.full(len(vertices), np.inf, dtype=float)
    for spec in specs:
        distance = np.linalg.norm(
            vertices - np.asarray(spec.center_world_um)[None, :], axis=1
        )
        selected = distance < spec.blend_length_um
        inside_junction |= selected
        local_radius[selected] = np.minimum(local_radius[selected], spec.radius_um)
    local_radius[~np.isfinite(local_radius)] = min(spec.radius_um for spec in specs)

    port_vertices: set[int] = set()
    port_faces = np.flatnonzero(np.asarray(patch.patch_type) != 0)
    if len(port_faces):
        port_vertices.update(
            map(int, np.unique(np.asarray(mesh.faces, dtype=np.int64)[port_faces]))
        )
    # Also protect the exact port planes even if a caller supplies a wall-only patch.
    for port in ports:
        relative = vertices - np.asarray(port.cap_center_um)[None, :]
        axial = np.abs(relative @ np.asarray(port.outward_tangent))
        radial_vector = relative - (
            relative @ np.asarray(port.outward_tangent)
        )[:, None] * np.asarray(port.outward_tangent)[None, :]
        radial = np.linalg.norm(radial_vector, axis=1)
        port_vertices.update(
            map(
                int,
                np.flatnonzero(
                    (axial <= max(config.ports.plane_tolerance_um, 1.0e-8))
                    & (radial <= 1.1 * port.radius_um)
                ),
            )
        )
    movable = np.asarray(
        sorted(
            vertex
            for vertex in support
            if vertex not in outer_ring
            and vertex not in port_vertices
            and inside_junction[vertex]
        ),
        dtype=np.int64,
    )
    fixed = np.setdiff1d(np.arange(len(vertices), dtype=np.int64), movable)
    maximum_displacement = (
        config.v9.fairing_max_displacement_radius_fraction * local_radius
    )
    iteration_rows: list[dict[str, Any]] = []
    for iteration in range(config.v9.fairing_iterations):
        update = vertices.copy()
        displacements: list[float] = []
        for vertex in movable:
            adjacent = neighbors[int(vertex)]
            if not len(adjacent):
                continue
            target = np.mean(vertices[adjacent], axis=0)
            candidate = vertices[vertex] + config.v9.fairing_relaxation * (
                target - vertices[vertex]
            )
            delta = candidate - original[vertex]
            norm = float(np.linalg.norm(delta))
            limit = float(maximum_displacement[vertex])
            if norm > limit:
                candidate = original[vertex] + delta * (limit / norm)
                norm = limit
            update[vertex] = candidate
            displacements.append(norm)
        vertices = update
        iteration_rows.append(
            {
                "iteration": iteration + 1,
                "maximum_displacement_um": max(displacements, default=0.0),
                "mean_displacement_um": float(np.mean(displacements))
                if displacements
                else 0.0,
            }
        )
    output = trimesh.Trimesh(
        vertices=vertices,
        faces=np.asarray(mesh.faces, dtype=np.int64),
        process=False,
    )
    trimesh.repair.fix_normals(output, multibody=True)
    fixed_error = float(np.max(np.linalg.norm(vertices[fixed] - original[fixed], axis=1)))
    moved_distance = np.linalg.norm(vertices[movable] - original[movable], axis=1)
    return JunctionFairingResult(
        mesh=output,
        moved_vertex_ids=movable,
        fixed_vertex_ids=fixed,
        report={
            "status": "APPLIED",
            "scope": "detected defect band plus constrained ring neighborhood",
            "defect_band_vertex_count": len(core),
            "support_vertex_count": len(support),
            "outer_band_fixed_vertex_count": len(outer_ring),
            "port_fixed_vertex_count": len(port_vertices),
            "moved_vertex_count": int(len(movable)),
            "fixed_vertex_count": int(len(fixed)),
            "fixed_vertex_maximum_displacement_um": fixed_error,
            "moved_vertex_displacement_um": {
                "mean": float(np.mean(moved_distance)) if len(moved_distance) else 0.0,
                "p95": float(np.percentile(moved_distance, 95)) if len(moved_distance) else 0.0,
                "max": float(np.max(moved_distance)) if len(moved_distance) else 0.0,
            },
            "iterations": iteration_rows,
            "newton_projection_after_fairing": False,
            "hard_min_field_reprojection_forbidden": True,
        },
    )
