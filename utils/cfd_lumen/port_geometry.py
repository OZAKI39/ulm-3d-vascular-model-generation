"""Geometric outward directions and straight CUT_PORT extensions."""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from utils.sampling.sampling_types import ROIRecord

from .config import CFDLumenConfig
from .types import GeometryValidationError, PortGeometry


_FACE_NORMALS = {
    "x_min": np.asarray((-1.0, 0.0, 0.0)),
    "x_max": np.asarray((1.0, 0.0, 0.0)),
    "y_min": np.asarray((0.0, -1.0, 0.0)),
    "y_max": np.asarray((0.0, 1.0, 0.0)),
    "z_min": np.asarray((0.0, 0.0, -1.0)),
    "z_max": np.asarray((0.0, 0.0, 1.0)),
}


def _outside_score(point: np.ndarray, minimum: np.ndarray, maximum: np.ndarray) -> float:
    below = np.maximum(minimum - point, 0.0)
    above = np.maximum(point - maximum, 0.0)
    return float(np.linalg.norm(below + above))


def construct_port_geometry(roi: ROIRecord, config: CFDLumenConfig) -> list[PortGeometry]:
    edges = np.asarray(roi.local_edges, dtype=np.int64).reshape((-1, 2))
    positions = np.asarray(roi.local_node_positions_um, dtype=float)
    minimum = np.asarray(roi.bbox_min_um, dtype=float)
    maximum = np.asarray(roi.bbox_max_um, dtype=float)
    output: list[PortGeometry] = []
    for port_id, cut in enumerate(roi.cut_ports):
        local_id = int(cut.local_node_id)
        incident = edges[np.any(edges == local_id, axis=1)]
        if len(incident) != 1:
            raise GeometryValidationError(
                f"CUT_PORT {cut.cut_port_id} must be a local degree-one node; got degree {len(incident)}"
            )
        neighbor = int(incident[0, 1] if incident[0, 0] == local_id else incident[0, 0])
        cut_position = np.asarray(cut.intersection_position_um, dtype=float)
        tangent = cut_position - positions[neighbor]
        norm = float(np.linalg.norm(tangent))
        if norm <= config.geometry.minimum_edge_length_um:
            raise GeometryValidationError(f"CUT_PORT {cut.cut_port_id} has an undefined tangent")
        tangent /= norm
        epsilon = max(0.1 * float(cut.radius_at_cut_um), 1.0e-3)
        candidates = (tangent, -tangent)
        scores = [
            _outside_score(cut_position + epsilon * candidate, minimum, maximum)
            for candidate in candidates
        ]
        face_normal = _FACE_NORMALS.get(cut.boundary_face)
        if abs(scores[0] - scores[1]) <= 1.0e-12 and face_normal is not None:
            outward = max(candidates, key=lambda candidate: float(np.dot(candidate, face_normal)))
        else:
            outward = candidates[int(scores[1] > scores[0])]
        if face_normal is not None and float(np.dot(outward, face_normal)) < 0:
            outward = -outward
        radius = float(cut.radius_at_cut_um)
        diameter = 2.0 * radius
        extension_length = config.ports.extension_diameters * diameter
        overlap_length = config.ports.overlap_diameters * diameter
        start = cut_position - overlap_length * outward
        cap_center = cut_position + extension_length * outward
        output.append(
            PortGeometry(
                port_id=port_id,
                roi_id=roi.roi_id,
                cut_port_id=cut.cut_port_id,
                local_node_id=local_id,
                global_edge_id=int(cut.global_edge_id),
                original_position_um=cut_position,
                cap_center_um=cap_center,
                radius_um=radius,
                outward_tangent=np.asarray(outward, dtype=float),
                extension_length_um=extension_length,
                overlap_length_um=overlap_length,
                cylinder_start_um=start,
                cylinder_end_um=cap_center,
                boundary_face=cut.boundary_face,
                boundary_role=cut.boundary_role,
                source_core_cut_port_id=cut.cut_port_id,
                original_core_cut_position_um=tuple(map(float, cut_position)),
            )
        )
    return output


def attach_core_port_provenance(
    ports: list[PortGeometry],
    port_mappings: list[dict[str, object]],
) -> list[PortGeometry]:
    """Attach immutable CORE boundary provenance to active CFD boundary ports."""

    provenance: dict[str, tuple[str, tuple[float, float, float], float]] = {}
    for mapping in port_mappings:
        original = mapping["original_cut_port"]
        if not isinstance(original, dict):
            continue
        position = (
            float(original["intersection_x_um"]),
            float(original["intersection_y_um"]),
            float(original["intersection_z_um"]),
        )
        for replacement in mapping["new_cfd_ports"]:
            if isinstance(replacement, dict):
                provenance[str(replacement["cut_port_id"])] = (
                    str(original["cut_port_id"]),
                    position,
                    float(mapping.get("added_centerline_length_um", 0.0)),
                )
    output: list[PortGeometry] = []
    for port in ports:
        source = provenance.get(port.cut_port_id)
        output.append(
            replace(
                port,
                source_core_cut_port_id=source[0] if source else port.cut_port_id,
                original_core_cut_position_um=(
                    source[1]
                    if source
                    else tuple(map(float, port.original_position_um))
                ),
                source_core_to_cfd_cut_length_um=source[2] if source else 0.0,
            )
        )
    return output
