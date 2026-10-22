"""Solver-independent port-plane directions and extension specifications."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from utils.sampling.sampling_types import CutPort, ROIRecord


@dataclass(frozen=True, slots=True)
class PortGeometry:
    simulation_tangent: np.ndarray
    outward_normal: np.ndarray
    extension_length_um: float
    extension_end_um: np.ndarray


def build_boundary_geometry(
    roi: ROIRecord,
    *,
    local_node_id: int,
    center_um: np.ndarray,
    radius_um: float,
    role: str,
    inlet_length_diameters: float,
    outlet_length_diameters: float,
    topology_error: str,
) -> PortGeometry:
    """Build common plane directions for either a cut or terminal boundary."""

    incident = [
        (int(upstream), int(downstream))
        for upstream, downstream in roi.local_edges
        if local_node_id in (int(upstream), int(downstream))
    ]
    if len(incident) != 1:
        raise ValueError(topology_error)
    upstream, downstream = incident[0]
    neighbor = downstream if upstream == local_node_id else upstream
    center = np.asarray(center_um, dtype=float)
    neighbor_position = np.asarray(roi.local_node_positions_um[neighbor], dtype=float)
    if role == "ASSUMED_INLET":
        vector = neighbor_position - center
        outward_sign = -1.0
        length_diameters = inlet_length_diameters
    elif role == "ASSUMED_OUTLET":
        vector = center - neighbor_position
        outward_sign = 1.0
        length_diameters = outlet_length_diameters
    else:
        raise ValueError(f"Unsupported port role: {role}")
    norm = float(np.linalg.norm(vector))
    if not np.isfinite(norm) or norm <= 0:
        raise ValueError(f"{topology_error}: tangent is not finite and non-zero")
    tangent = vector / norm
    outward = outward_sign * tangent
    extension_length = length_diameters * 2.0 * radius_um
    return PortGeometry(
        simulation_tangent=tangent,
        outward_normal=outward,
        extension_length_um=float(extension_length),
        extension_end_um=center + outward * extension_length,
    )


def build_port_geometry(
    roi: ROIRecord,
    port: CutPort,
    *,
    role: str,
    inlet_length_diameters: float,
    outlet_length_diameters: float,
) -> PortGeometry:
    """Compatibility wrapper for existing CUT_PORT geometry."""

    return build_boundary_geometry(
        roi,
        local_node_id=port.local_node_id,
        center_um=np.asarray(port.intersection_position_um, dtype=float),
        radius_um=port.radius_at_cut_um,
        role=role,
        inlet_length_diameters=inlet_length_diameters,
        outlet_length_diameters=outlet_length_diameters,
        topology_error="AMBIGUOUS_CUT_PORT_DIRECTION",
    )
