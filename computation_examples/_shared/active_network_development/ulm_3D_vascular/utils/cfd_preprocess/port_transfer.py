"""Classify ROI cut ports and transfer global 1D pressure/flow exactly."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from utils.sampling.sampling_types import CutPort, GlobalVascularModel, ROIRecord

from .one_d_flow import GlobalFlowResult, edge_resistance
from .port_geometry import (
    PortGeometry,
    build_boundary_geometry,
    build_port_geometry,
)


class PortTransferError(RuntimeError):
    """An exact cut-port/global-edge mapping validation failed."""


@dataclass(frozen=True, slots=True)
class PortTransfer:
    port_id: str
    local_node_id: int
    role: str
    boundary_origin: str
    global_node_id: int | None
    global_edge_id: int
    center_um: np.ndarray
    source_radius_um: float
    alpha_on_global_edge: float
    position_error_um: float
    radius_relative_error: float
    pressure_pa: float
    signed_parent_to_child_flow_m3_s: float
    role_flow_m3_s: float
    geometry: PortGeometry

    @property
    def pressure_mmhg(self) -> float:
        return self.pressure_pa / 133.32236842105263

    @property
    def flow_pl_s(self) -> float:
        return self.role_flow_m3_s * 1.0e15


def classify_cut_port(roi: ROIRecord, port: CutPort) -> str:
    incident = [
        (int(upstream), int(downstream))
        for upstream, downstream in roi.local_edges
        if port.local_node_id in (int(upstream), int(downstream))
    ]
    if len(incident) != 1:
        raise PortTransferError("AMBIGUOUS_CUT_PORT_DIRECTION")
    upstream, downstream = incident[0]
    if port.local_node_id == upstream:
        return "ASSUMED_INLET"
    if port.local_node_id == downstream:
        return "ASSUMED_OUTLET"
    raise PortTransferError("AMBIGUOUS_CUT_PORT_DIRECTION")


def transfer_cut_port(
    roi: ROIRecord,
    port: CutPort,
    model: GlobalVascularModel,
    flow: GlobalFlowResult,
    *,
    mu_pa_s: float,
    maximum_position_error_um: float,
    maximum_radius_relative_error: float,
    inlet_length_diameters: float,
    outlet_length_diameters: float,
) -> PortTransfer:
    if port.global_edge_id < 0 or port.global_edge_id >= model.edge_count:
        raise PortTransferError(
            "GLOBAL_TO_ROI_TRANSFER_FAILED: global edge ID is absent"
        )
    edge = model.edges[port.global_edge_id]
    if edge.edge_id != port.global_edge_id:
        raise PortTransferError(
            "GLOBAL_TO_ROI_TRANSFER_FAILED: edge index identity is invalid"
        )
    x0 = np.asarray(edge.upstream_position_um, dtype=float)
    x1 = np.asarray(edge.downstream_position_um, dtype=float)
    x_cut = np.asarray(port.intersection_position_um, dtype=float)
    edge_vector = x1 - x0
    denominator = float(np.dot(edge_vector, edge_vector))
    if denominator <= 0:
        raise PortTransferError(
            "GLOBAL_TO_ROI_TRANSFER_FAILED: zero-length global edge"
        )
    alpha = float(np.dot(x_cut - x0, edge_vector) / denominator)
    if alpha < -1.0e-10 or alpha > 1.0 + 1.0e-10:
        raise PortTransferError(
            "CUT_PORT_GLOBAL_EDGE_POSITION_MISMATCH: alpha outside [0,1]"
        )
    alpha = float(np.clip(alpha, 0.0, 1.0))
    reconstructed = x0 + alpha * edge_vector
    position_error = float(np.linalg.norm(reconstructed - x_cut))
    if position_error > maximum_position_error_um:
        raise PortTransferError(
            f"CUT_PORT_GLOBAL_EDGE_POSITION_MISMATCH: error={position_error:.6g} um"
        )
    radius_global = edge.upstream_radius_um + alpha * (
        edge.downstream_radius_um - edge.upstream_radius_um
    )
    radius_error = abs(radius_global - port.radius_at_cut_um) / radius_global
    if radius_error > maximum_radius_relative_error:
        raise PortTransferError(
            f"CUT_PORT_RADIUS_MAPPING_MISMATCH: relative_error={radius_error:.6g}"
        )
    role = classify_cut_port(roi, port)
    full_length_m = float(np.linalg.norm(edge_vector)) * 1.0e-6
    partial_resistance = (
        edge_resistance(
            alpha * full_length_m,
            edge.upstream_radius_um * 1.0e-6,
            radius_global * 1.0e-6,
            mu_pa_s,
        )
        if alpha > 0
        else 0.0
    )
    parent_pressure = flow.pressure_by_node_id[edge.upstream_node_id]
    signed_flow = flow.flow_by_edge_id[edge.edge_id]
    pressure = parent_pressure - signed_flow * partial_resistance
    role_flow = signed_flow if role == "ASSUMED_INLET" else signed_flow
    if not np.isfinite(pressure) or not np.isfinite(role_flow) or role_flow <= 0:
        raise PortTransferError("GLOBAL_TO_ROI_TRANSFER_FAILED: non-positive role flow")
    geometry = build_port_geometry(
        roi,
        port,
        role=role,
        inlet_length_diameters=inlet_length_diameters,
        outlet_length_diameters=outlet_length_diameters,
    )
    return PortTransfer(
        port_id=port.cut_port_id,
        local_node_id=port.local_node_id,
        role=role,
        boundary_origin="CUT_PORT",
        global_node_id=None,
        global_edge_id=edge.edge_id,
        center_um=x_cut,
        source_radius_um=float(radius_global),
        alpha_on_global_edge=alpha,
        position_error_um=position_error,
        radius_relative_error=float(radius_error),
        pressure_pa=float(pressure),
        signed_parent_to_child_flow_m3_s=float(signed_flow),
        role_flow_m3_s=float(role_flow),
        geometry=geometry,
    )


def transfer_true_terminal(
    roi: ROIRecord,
    model: GlobalVascularModel,
    flow: GlobalFlowResult,
    *,
    terminal_index: int,
    local_node_id: int,
    global_node_id: int,
    maximum_position_error_um: float,
    maximum_radius_relative_error: float,
    inlet_length_diameters: float,
    outlet_length_diameters: float,
) -> PortTransfer:
    """Map one saved ROI true terminal to its exact global structural leaf."""

    if local_node_id < 0 or local_node_id >= roi.node_count:
        raise PortTransferError("TRUE_TERMINAL_LOCAL_TOPOLOGY_INVALID")
    global_matches = np.flatnonzero(model.node_ids == global_node_id)
    if len(global_matches) != 1:
        raise PortTransferError("TRUE_TERMINAL_GLOBAL_MAPPING_INVALID")
    global_index = int(global_matches[0])
    if global_node_id == flow.root_node_id:
        raise PortTransferError("TRUE_TERMINAL_GLOBAL_MAPPING_INVALID")
    outgoing = [edge for edge in model.edges if edge.upstream_node_id == global_node_id]
    incoming = [
        edge for edge in model.edges if edge.downstream_node_id == global_node_id
    ]
    if outgoing or len(incoming) != 1:
        raise PortTransferError("TRUE_TERMINAL_GLOBAL_MAPPING_INVALID")
    edge = incoming[0]

    local_position = np.asarray(roi.local_node_positions_um[local_node_id], dtype=float)
    global_position = np.asarray(model.node_positions_um[global_index], dtype=float)
    position_error = float(np.linalg.norm(local_position - global_position))
    if position_error > maximum_position_error_um:
        raise PortTransferError(
            f"TRUE_TERMINAL_POSITION_MISMATCH: error={position_error:.6g} um"
        )
    local_radius = float(roi.local_node_radius_um[local_node_id])
    global_radius = float(model.node_radius_um[global_index])
    radius_error = abs(local_radius - global_radius) / global_radius
    if radius_error > maximum_radius_relative_error:
        raise PortTransferError(
            f"TRUE_TERMINAL_RADIUS_MISMATCH: relative_error={radius_error:.6g}"
        )
    pressure = flow.pressure_by_node_id[global_node_id]
    signed_flow = flow.flow_by_edge_id[edge.edge_id]
    if not np.isfinite(signed_flow) or signed_flow <= 0:
        raise PortTransferError("TRUE_TERMINAL_FLOW_INVALID")
    if not np.isfinite(pressure):
        raise PortTransferError(
            "TRUE_TERMINAL_GLOBAL_MAPPING_INVALID: non-finite pressure"
        )
    try:
        geometry = build_boundary_geometry(
            roi,
            local_node_id=local_node_id,
            center_um=local_position,
            radius_um=global_radius,
            role="ASSUMED_OUTLET",
            inlet_length_diameters=inlet_length_diameters,
            outlet_length_diameters=outlet_length_diameters,
            topology_error="TRUE_TERMINAL_LOCAL_TOPOLOGY_INVALID",
        )
    except ValueError as error:
        raise PortTransferError(str(error)) from error
    return PortTransfer(
        port_id=f"{roi.roi_id}__terminal_{terminal_index:03d}",
        local_node_id=local_node_id,
        role="ASSUMED_OUTLET",
        boundary_origin="TRUE_TERMINAL",
        global_node_id=global_node_id,
        global_edge_id=edge.edge_id,
        center_um=local_position,
        source_radius_um=global_radius,
        alpha_on_global_edge=1.0,
        position_error_um=position_error,
        radius_relative_error=radius_error,
        pressure_pa=float(pressure),
        signed_parent_to_child_flow_m3_s=float(signed_flow),
        role_flow_m3_s=float(signed_flow),
        geometry=geometry,
    )


def transfer_all_boundaries(
    roi: ROIRecord,
    model: GlobalVascularModel,
    flow: GlobalFlowResult,
    **kwargs: float,
) -> tuple[list[PortTransfer], float]:
    """Transfer every CUT_PORT and TRUE_TERMINAL CFD boundary exactly once."""

    if len(roi.true_terminal_local_ids) != len(roi.true_terminal_global_ids):
        raise PortTransferError("TRUE_TERMINAL_GLOBAL_MAPPING_INVALID")
    transfers = [
        transfer_cut_port(roi, port, model, flow, **kwargs) for port in roi.cut_ports
    ]
    terminal_kwargs = {key: value for key, value in kwargs.items() if key != "mu_pa_s"}
    transfers.extend(
        transfer_true_terminal(
            roi,
            model,
            flow,
            terminal_index=index,
            local_node_id=int(local_node_id),
            global_node_id=int(global_node_id),
            **terminal_kwargs,
        )
        for index, (local_node_id, global_node_id) in enumerate(
            zip(roi.true_terminal_local_ids, roi.true_terminal_global_ids)
        )
    )
    inlet_flow = sum(
        item.role_flow_m3_s for item in transfers if item.role == "ASSUMED_INLET"
    )
    outlet_flow = sum(
        item.role_flow_m3_s for item in transfers if item.role == "ASSUMED_OUTLET"
    )
    error = (
        abs(inlet_flow - outlet_flow) / inlet_flow if inlet_flow > 0 else float("inf")
    )
    return transfers, float(error)


def transfer_all_ports(
    roi: ROIRecord,
    model: GlobalVascularModel,
    flow: GlobalFlowResult,
    **kwargs: float,
) -> tuple[list[PortTransfer], float]:
    """Deprecated compatibility wrapper; use transfer_all_boundaries()."""

    return transfer_all_boundaries(roi, model, flow, **kwargs)
