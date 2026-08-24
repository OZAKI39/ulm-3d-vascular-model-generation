"""Conservative, topology-fixed Mask-assisted junction refinement experiment."""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Any

import numpy as np
import trimesh

from .mask_surface import sample_sdf_and_gradient


@dataclass(frozen=True, slots=True)
class MaskInfluenceConfig:
    name: str
    lambda_mask: float
    lambda_radius: float
    lambda_smooth: float
    step_size: float
    iterations: int
    maximum_displacement_um: float


INFLUENCE_CONFIGS = (
    MaskInfluenceConfig("weak", 0.15, 1.0, 0.08, 0.08, 10, 0.15),
    MaskInfluenceConfig("medium", 0.35, 1.0, 0.10, 0.08, 12, 0.30),
    MaskInfluenceConfig("strong", 0.65, 1.0, 0.12, 0.08, 15, 0.50),
)


def spatial_influence_weights(
    points_xyz_um: np.ndarray,
    *,
    junction_center_xyz_um: np.ndarray,
    core_radius_um: float,
    outer_radius_um: float,
) -> np.ndarray:
    """One in the core, smoothstep decay through the collar, zero outside."""

    distance = np.linalg.norm(
        np.asarray(points_xyz_um, dtype=float)
        - np.asarray(junction_center_xyz_um, dtype=float)[None, :],
        axis=1,
    )
    denominator = max(float(outer_radius_um - core_radius_um), 1.0e-12)
    coordinate = np.clip((distance - core_radius_um) / denominator, 0.0, 1.0)
    smoothstep = coordinate * coordinate * (3.0 - 2.0 * coordinate)
    return np.where(distance <= core_radius_um, 1.0, 1.0 - smoothstep)


def _movable_vertices(mesh: trimesh.Trimesh, face_region: np.ndarray) -> np.ndarray:
    """Move only vertices whose incident faces are all core/collar (codes 1/2)."""

    labels = np.asarray(face_region, dtype=np.uint8)
    if len(labels) != len(mesh.faces):
        raise ValueError("face_region length does not match the surface")
    incident: list[list[int]] = [[] for _ in range(len(mesh.vertices))]
    for face_id, face in enumerate(np.asarray(mesh.faces, dtype=np.int64)):
        for vertex_id in face:
            incident[int(vertex_id)].append(face_id)
    movable = np.zeros(len(mesh.vertices), dtype=bool)
    for vertex_id, face_ids in enumerate(incident):
        if face_ids and np.all(np.isin(labels[np.asarray(face_ids)], (1, 2))):
            movable[vertex_id] = True
    return movable


def _energies(
    mesh: trimesh.Trimesh,
    original_vertices: np.ndarray,
    sdf: np.ndarray,
    *,
    origin_xyz_um: tuple[float, float, float],
    active: np.ndarray,
) -> dict[str, float]:
    vertices = np.asarray(mesh.vertices, dtype=float)
    values, _ = sample_sdf_and_gradient(sdf, vertices[active], origin_xyz_um=origin_xyz_um)
    displacement = vertices[active] - original_vertices[active]
    smooth_values: list[float] = []
    for vertex_id in np.flatnonzero(active):
        neighbors = np.asarray(mesh.vertex_neighbors[int(vertex_id)], dtype=np.int64)
        if len(neighbors):
            smooth_values.append(
                float(np.sum((vertices[vertex_id] - vertices[neighbors].mean(axis=0)) ** 2))
            )
    return {
        "E_mask_um2": float(np.mean(values**2)) if len(values) else 0.0,
        "E_radius_proxy_um2": float(np.mean(np.sum(displacement**2, axis=1))) if len(displacement) else 0.0,
        "E_smooth_um2": float(np.mean(smooth_values)) if smooth_values else 0.0,
    }


def refine_surface_with_mask(
    mesh: trimesh.Trimesh,
    face_region: np.ndarray,
    sdf: np.ndarray,
    *,
    origin_xyz_um: tuple[float, float, float],
    junction_center_xyz_um: np.ndarray,
    junction_radius_um: float,
    influence_outer_radius_um: float,
    config: MaskInfluenceConfig,
) -> tuple[trimesh.Trimesh, dict[str, Any]]:
    """Refine positions while preserving faces, connectivity, and explicit vertices."""

    started = time.perf_counter()
    refined = mesh.copy()
    original = np.asarray(mesh.vertices, dtype=float).copy()
    movable = _movable_vertices(mesh, face_region)
    weights = spatial_influence_weights(
        original,
        junction_center_xyz_um=np.asarray(junction_center_xyz_um, dtype=float),
        core_radius_um=max(2.0 * float(junction_radius_um), 2.0),
        outer_radius_um=float(influence_outer_radius_um),
    )
    active = movable & (weights > 0.0)
    fixed = ~active
    history: list[dict[str, float | int]] = []
    if not np.any(active):
        raise ValueError("Mask-assisted refinement has no active core/collar vertices")

    for iteration in range(config.iterations):
        vertices = np.asarray(refined.vertices, dtype=float)
        values, gradient = sample_sdf_and_gradient(
            sdf, vertices, origin_xyz_um=origin_xyz_um
        )
        mask_gradient = 2.0 * values[:, None] * gradient
        radius_gradient = 2.0 * (vertices - original)
        smooth_gradient = np.zeros_like(vertices)
        for vertex_id in np.flatnonzero(active):
            neighbors = np.asarray(refined.vertex_neighbors[int(vertex_id)], dtype=np.int64)
            if len(neighbors):
                smooth_gradient[vertex_id] = 2.0 * (
                    vertices[vertex_id] - vertices[neighbors].mean(axis=0)
                )
        total_gradient = (
            config.lambda_mask * mask_gradient
            + config.lambda_radius * radius_gradient
            + config.lambda_smooth * smooth_gradient
        )
        proposed = vertices.copy()
        proposed[active] -= (
            config.step_size * weights[active, None] * total_gradient[active]
        )
        displacement = proposed - original
        magnitude = np.linalg.norm(displacement, axis=1)
        over = active & (magnitude > config.maximum_displacement_um)
        proposed[over] = original[over] + displacement[over] * (
            config.maximum_displacement_um / magnitude[over]
        )[:, None]
        proposed[fixed] = original[fixed]
        refined.vertices = proposed
        energy = _energies(
            refined,
            original,
            sdf,
            origin_xyz_um=origin_xyz_um,
            active=active,
        )
        energy["iteration"] = iteration + 1
        energy["E_total"] = (
            config.lambda_mask * energy["E_mask_um2"]
            + config.lambda_radius * energy["E_radius_proxy_um2"]
            + config.lambda_smooth * energy["E_smooth_um2"]
        )
        history.append(energy)

    trimesh.repair.fix_normals(refined, multibody=True)
    displacement = np.asarray(refined.vertices) - original
    report: dict[str, Any] = {
        "experimental_only": True,
        "configuration": config.name,
        "lambda_mask": config.lambda_mask,
        "lambda_radius": config.lambda_radius,
        "lambda_smooth": config.lambda_smooth,
        "iterations": config.iterations,
        "step_size": config.step_size,
        "maximum_displacement_um": config.maximum_displacement_um,
        "active_vertex_count": int(np.count_nonzero(active)),
        "fixed_vertex_count": int(np.count_nonzero(fixed)),
        "explicit_branch_vertex_count_moved": 0,
        "maximum_actual_displacement_um": float(np.linalg.norm(displacement, axis=1).max()),
        "mean_active_displacement_um": float(
            np.linalg.norm(displacement[active], axis=1).mean()
        ),
        "source_vertex_count_unchanged": len(refined.vertices) == len(mesh.vertices),
        "source_face_count_unchanged": len(refined.faces) == len(mesh.faces),
        "source_faces_byte_identical": bool(np.array_equal(refined.faces, mesh.faces)),
        "source_branch_count_unchanged": True,
        "source_junction_connectivity_unchanged": True,
        "topology_definition": "corrected SWC; Mask is a local shape soft reference only",
        "energy_history": history,
        "runtime_s": time.perf_counter() - started,
    }
    return refined, report
