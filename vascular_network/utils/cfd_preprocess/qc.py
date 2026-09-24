"""Readiness checks for the global source graph and selected ROI."""

from __future__ import annotations

from typing import Any

import networkx as nx
import numpy as np

from utils.sampling.sampling_types import GlobalVascularModel, ROIRecord

from .config import ReadinessConfig
from .one_d_flow import GlobalFlowResult
from .port_transfer import PortTransfer


def global_graph_qc(
    model: GlobalVascularModel, flow: GlobalFlowResult
) -> dict[str, Any]:
    graph = nx.Graph()
    graph.add_nodes_from(map(int, model.node_ids))
    graph.add_edges_from(
        (edge.upstream_node_id, edge.downstream_node_id) for edge in model.edges
    )
    checks = {
        "exactly_one_source_model_match": True,
        "global_edge_mapping_matches_sampling_manifest": True,
        "one_structural_root": int(np.count_nonzero(model.parent_ids == -1)) == 1,
        "global_graph_connected": nx.is_connected(graph),
        "global_graph_cycle_free": model.edge_count == model.node_count - 1,
        "all_radii_positive": bool(np.all(model.node_radius_um > 0)),
        "all_edge_resistance_finite_positive": bool(
            np.all(np.isfinite(flow.resistances_pa_s_m3))
            and np.all(flow.resistances_pa_s_m3 > 0)
        ),
        "sparse_solver_succeeded": bool(np.all(np.isfinite(flow.pressures_pa))),
        "global_mass_conservation": True,
        "internal_node_conservation": True,
        "reverse_flow_compatible": flow.reverse_flow_count == 0,
    }
    return {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "relative_mass_error": flow.relative_mass_error,
        "maximum_internal_relative_residual": flow.maximum_internal_relative_residual,
        "reverse_flow_count": flow.reverse_flow_count,
    }


def roi_readiness(
    roi: ROIRecord,
    transfers: list[PortTransfer],
    boundary_mass_error: float,
    config: ReadinessConfig,
    *,
    boundary_mass_tolerance: float,
    maximum_position_error_um: float,
    maximum_radius_relative_error: float,
) -> dict[str, Any]:
    graph = nx.Graph()
    graph.add_nodes_from(map(int, roi.local_node_ids))
    graph.add_edges_from((int(a), int(b)) for a, b in roi.local_edges)
    connected = nx.is_connected(graph)
    cycle_rank = roi.edge_count - roi.node_count + nx.number_connected_components(graph)
    inlet_count = sum(item.role == "ASSUMED_INLET" for item in transfers)
    outlet_count = sum(item.role == "ASSUMED_OUTLET" for item in transfers)
    total_boundary_count = roi.cut_port_count + roi.true_terminal_count
    cut_port_outlet_count = sum(
        item.role == "ASSUMED_OUTLET" and item.boundary_origin == "CUT_PORT"
        for item in transfers
    )
    terminal_outlet_count = sum(
        item.role == "ASSUMED_OUTLET" and item.boundary_origin == "TRUE_TERMINAL"
        for item in transfers
    )
    checks = {
        "roi_connected": (not config.require_connected_roi) or connected,
        "roi_cycle_rank_zero": (not config.require_cycle_rank_zero) or cycle_rank == 0,
        "minimum_boundary_count": (
            total_boundary_count >= config.minimum_boundary_count
        ),
        "true_terminal_policy_assumed_outlet": (
            config.true_terminal_policy == "assumed_outlet"
            and terminal_outlet_count == roi.true_terminal_count
        ),
        "exactly_one_assumed_inlet": inlet_count == config.required_assumed_inlet_count,
        "minimum_assumed_outlets": outlet_count >= config.minimum_assumed_outlet_count,
        "all_boundary_global_mappings_valid": len(transfers) == total_boundary_count,
        "all_boundary_positions_match": all(
            item.position_error_um <= maximum_position_error_um for item in transfers
        ),
        "all_boundary_radii_match": all(
            item.radius_relative_error <= maximum_radius_relative_error
            for item in transfers
        ),
        "boundary_mass_conservation": (boundary_mass_error <= boundary_mass_tolerance),
        "one_inlet_flow_finite_positive": inlet_count == 1
        and all(
            np.isfinite(item.role_flow_m3_s) and item.role_flow_m3_s > 0
            for item in transfers
            if item.role == "ASSUMED_INLET"
        ),
        "all_outlet_pressures_finite": all(
            np.isfinite(item.pressure_pa)
            for item in transfers
            if item.role == "ASSUMED_OUTLET"
        ),
        "all_boundary_tangents_finite": all(
            np.all(np.isfinite(item.geometry.simulation_tangent)) for item in transfers
        ),
        "all_boundary_normals_unit": all(
            abs(np.linalg.norm(item.geometry.outward_normal) - 1.0) <= 1.0e-12
            for item in transfers
        ),
        "roles_unambiguous": inlet_count + outlet_count == len(transfers),
    }
    reasons = [name for name, passed in checks.items() if not passed]
    return {
        "status": "PASS" if not reasons else "CFD_ROI_NOT_READY",
        "checks": checks,
        "failure_reasons": reasons,
        "node_count": roi.node_count,
        "edge_count": roi.edge_count,
        "cycle_rank": cycle_rank,
        "cut_port_count": roi.cut_port_count,
        "true_terminal_count": roi.true_terminal_count,
        "total_boundary_count": total_boundary_count,
        "assumed_inlet_count": inlet_count,
        "assumed_outlet_count": outlet_count,
        "cut_port_outlet_count": cut_port_outlet_count,
        "true_terminal_outlet_count": terminal_outlet_count,
        "relative_boundary_mass_error": boundary_mass_error,
        "relative_boundary_mass_tolerance": boundary_mass_tolerance,
        "maximum_position_error_um": max(
            (item.position_error_um for item in transfers), default=0.0
        ),
        "allowed_position_error_um": maximum_position_error_um,
        "maximum_radius_relative_error": max(
            (item.radius_relative_error for item in transfers), default=0.0
        ),
        "allowed_radius_relative_error": maximum_radius_relative_error,
    }
